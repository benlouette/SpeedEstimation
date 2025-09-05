function [fr_hat, OUT] = estimate_speed_hybrid(x, Fs, P)
% Estimation "tachless" hybride de la vitesse de rotation:
%  1) Détection de bandes résonantes via kurtosis par bande
%  2) Dans chaque bande:
%     - Enveloppe (Hilbert) + FFT + ACF + Cepstre
%     - Sidebands (peigne autour de la porteuse)
%     - CSC (Cyclic Spectral Coherence)
%     => Fusion robuste intra-bande (médiane filtrée + tie-break CSC)
%  3) Fusion inter-bandes (médiane / pondération simple)
%
% Requiert les helpers déjà fournis : find_impulsive_bands, bandpass_iir,
% envelope_fft, autocorr_envelope, cepstrum_peak, sideband_grid_score,
% stft_hann, compute_csc.

if size(x,2)>1, x = x(:); end
x = detrend(x,'linear');

% --- 1) bandes impulsives
[bands, bandScores] = find_impulsive_bands(x, Fs, P.searchBand, P.nBands);
[~, idxSort] = sort(bandScores, 'descend');
keep = idxSort(1:min(P.topK, numel(idxSort)));
bandsKeep = bands(keep,:);

fr_band = nan(numel(keep),1);
Band = struct([]);

for i = 1:numel(keep)
    f1 = bandsKeep(i,1); f2 = bandsKeep(i,2);
    y  = bandpass_iir(x, Fs, f1, f2, 4);

    % --- Enveloppe
    env = abs(hilbert(y));
    [Fenv, Aenv, fr_env] = envelope_fft(env, Fs, [P.minFr P.maxFr]);
    [fr_acf, acLag, acCurve] = autocorr_envelope(env, Fs, [P.minFr P.maxFr]);
    [fr_cep, qvec, Cmag]     = cepstrum_peak(env, Fs, [P.minFr P.maxFr]);

    % SNR sur enveloppe (grossier)
    if ~isempty(Aenv)
        [~,iPk] = max(Aenv); pk = Aenv(iPk);
        nb = max(3, round(0.01*numel(Aenv)));
        neigh = Aenv(max(1,iPk-nb):min(numel(Aenv),iPk+nb));
        noise = median(neigh);
        snr_env = 20*log10((pk+eps)/(noise+eps));
    else
        snr_env = 0;
    end

    % --- Sidebands (sur brut filtré bande porteuse)
    [fr_sb, sbScore] = sideband_grid_score(y, Fs, mean([f1 f2]), P.sidebandWinHz, [P.minFr P.maxFr]);

    % --- CSC (sur la même bande)
    alphaVec = P.minFr:P.alphaStep:P.maxFr;
    [Z, fVec] = stft_hann(y, Fs, P.stft.winLen, P.stft.hop, P.stft.nfft);
    [alphaVecOut, fValid, CSCabs, alphaScore] = compute_csc(Z, fVec, Fs, alphaVec, [f1 f2]); %#ok<ASGLU>
    [alphaPeak, iAlphaPk] = max(alphaScore);
    if isempty(alphaPeak) || all(~isfinite(alphaScore))
        fr_csc = 0; alphaPeak = 0;
    else
        fr_csc = alphaVecOut(iAlphaPk);
    end

    % --- Fusion intra-bande (robuste, simple et efficace)
    cand = [fr_env, fr_acf, fr_cep, fr_sb, fr_csc];
    cand = cand(cand>0 & isfinite(cand));
    if isempty(cand)
        fr_b = 0;
    else
        m0 = median(cand);
        if m0>0
            % garde les candidats à ±25%
            keepC = abs(cand - m0)./m0 <= 0.25;
            if any(keepC), cand = cand(keepC); end
            m1 = median(cand);
            % tie-break avec CSC si cohérent
            if fr_csc>0 && abs(fr_csc - m1)/max(m1,eps) <= 0.25
                fr_b = 0.5*(m1 + fr_csc);
            else
                fr_b = m1;
            end
        else
            fr_b = m0;
        end
    end
    fr_band(i) = fr_b;

    % --- stock diag
    Band(i).f1 = f1; Band(i).f2 = f2;
    Band(i).fr_env = fr_env; Band(i).fr_acf = fr_acf; Band(i).fr_cep = fr_cep;
    Band(i).fr_sb = fr_sb; Band(i).sbScore = sbScore;
    Band(i).fr_csc = fr_csc; Band(i).alphaPeak = alphaPeak;
    Band(i).snr_env = snr_env; Band(i).Fenv = Fenv; Band(i).Aenv = Aenv;
    Band(i).acLag = acLag; Band(i).acCurve = acCurve;
    Band(i).qvec = qvec; Band(i).Cmag = Cmag;
end

% --- 3) Fusion inter-bandes
% médiane par défaut; on peut pondérer par bandScores et alphaPeak si besoin
valid = fr_band(fr_band>0 & isfinite(fr_band));
if isempty(valid)
    fr_hat = 0;
else
    fr_hat = median(valid);
end

% --- sorties
OUT = struct();
OUT.bands = bands; OUT.bandScores = bandScores;
OUT.bandsKeep = bandsKeep; OUT.fr_band = fr_band; OUT.Band = Band;

if isfield(P,'plotting') && P.plotting
    visualize_hybrid_results(x, Fs, P, OUT, fr_hat);
end
end
