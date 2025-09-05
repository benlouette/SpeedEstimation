% DEMO — tachless speed estimation from vibration
% clear; close all; clc;

% % --------- Synthèse d’un signal type roulement (pour valider la chaîne) ----------
% Fs    = 51200;           % Hz
% T     = 10;              % s
% t     = (0:1/Fs:T-1/Fs).';
% fr0   = 17.3;            % Hz (vitesse de rotation "vraie" pour le test)
% f0    = 6500;            % Hz (résonance structurelle porteuse)
% Q     = 25;              % Facteur de qualité approx (bande ~ f0/Q)
% bw    = f0/Q;
% 
% % Impacts AM à fr0 modulant une résonance à f0
% x     = generate_bearing_like_signal(t, Fs, fr0, f0, bw);

% --------- Estimation tachless ----------
params = struct();
params.searchBand = [200 4000];     % zone de recherche des résonances
params.nBands     = 24;              % nb sous-bandes pour kurtosis
params.topK       = 10;               % nb bandes retenues
params.maxFr      = 50;             % Hz : vitesse max à chercher
params.minFr      = 0.1;             % Hz : vitesse min
params.sidebandWinHz = 600;          % largeur d’analyse autour de la porteuse
params.plotting   = true;

Fs = 48000;
fr0 = 29.95;
test1 = load('Normal Baseline Data 12K/Normal_0_1797_rpm.mat'); % 29.95 Hz
x=test1.X097_FE_time;
% [fr_hat, diagOut] = estimate_speed_from_vibration(x, Fs, params);

% ---------- Paramètres ----------
P = struct();
P.searchBand     = [16508 20000];   % zone de recherche des résonances
P.nBands         = 24;            % nb sous-bandes pour la kurtosis
P.topK           = 10;             % nb de bandes retenues
P.minFr          = 0.2;           % Hz
P.maxFr          = 50;           % Hz
P.sidebandWinHz  = 500;           % fenêtre autour de la porteuse pour sidebands
P.plotting       = true;

% CSC (nécessaire au mode hybride)
P.alphaStep      = 0.2;           % pas en Hz pour la grille alpha
P.stft.winLen    = 4096;
P.stft.hop       = 1024;
P.stft.nfft      = 8192;

% ---------- Estimation ----------
[fr_hat,OUT] = estimate_speed_hybrid(x, Fs, P); 
fprintf('\nHybride: fr_hat = %.3f Hz (vrai = %.3f Hz)\n', fr_hat, fr0);


fprintf('\nEstimation finale fr_hat = %.3f Hz (vrai=%.3f Hz) \n', fr_hat, fr0);

function [fr_hat] = m_estimate_speed_hybrid(x, Fs, P,f1,f2)
% Estimation "tachless" hybride de la vitesse de rotation:
%  1) Détection de bandes résonantes via kurtosis par bande
%  2) Dans chaque bande:
%     - Enveloppe (Hilbert) + FFT + ACF + Cepstre
%     - Sidebands (peigne autour de la porteuse)
%     - CSC (Cyclic Spectral Coherence)
%     => Fusion robuste intra-bande (médiane filtrée + tie-break CSC)
%  3) Fusion inter-bandes (médiane / pondération simple)
%ZZ
% Requiert les helpers déjà fournis : find_impulsive_bands, bandpass_iir,
% envelope_fft, autocorr_envelope, cepstrum_peak, sideband_grid_score,
% stft_hann, compute_csc.

if size(x,2)>1, x = x(:); end
x = detrend(x,'linear');

% --- 1) bandes impulsives
% [bands, bandScores] = find_impulsive_bands(x, Fs, P.searchBand, P.nBands);
% [~, idxSort] = sort(bandScores, 'descend');
% keep = idxSort(1:min(P.topK, numel(idxSort)));
% bandsKeep = bands(keep,:);
% 
% fr_band = nan(numel(keep),1);
Band = struct();

% for i = 1:numel(keep)
%     f1 = bandsKeep(i,1); f2 = bandsKeep(i,2);
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
    fr_band = fr_b;

    % --- stock diag
    Band.f1 = f1; Band.f2 = f2;
    Band.fr_env = fr_env; Band.fr_acf = fr_acf; Band.fr_cep = fr_cep;
    Band.fr_sb = fr_sb; Band.sbScore = sbScore;
    Band.fr_csc = fr_csc; Band.alphaPeak = alphaPeak;
    Band.snr_env = snr_env; Band.Fenv = Fenv; Band.Aenv = Aenv;
    Band.acLag = acLag; Band.acCurve = acCurve;
    Band.qvec = qvec; Band.Cmag = Cmag;
% end

% --- 3) Fusion inter-bandes
% médiane par défaut; on peut pondérer par bandScores et alphaPeak si besoin
valid = fr_band(fr_band>0 & isfinite(fr_band));
if isempty(valid)
    fr_hat = 0;
else
    fr_hat = median(valid);
end

end
