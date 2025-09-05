function [fr_hat, D] = estimate_speed_from_vibration(x, Fs, P)
% Orchestrateur : trouve bandes impulsives, démodule, estime fr
% In:
%   x : signal brut (colonne)
%   Fs: échantillonnage
%   P : params (voir demo)
% Out:
%   fr_hat : estimation finale de la vitesse de rotation (Hz)
%   D      : diagnostics (struct)

if size(x,2)>1, x = x(:); end
x = detrend(x,'linear');

% 1) Trouver bandes résonantes par kurtosis
[bands, bandScores] = find_impulsive_bands(x, Fs, P.searchBand, P.nBands);

% Garder topK
[~, idxSort] = sort(bandScores, 'descend');
keep = idxSort(1:min(P.topK, numel(idxSort)));
bandsKeep = bands(keep,:);

% 2) Pour chaque bande, démodulation + estimations
estimates = [];  % [fr_env fr_acf fr_cep fr_sideband SNR]
bandSummaries = cell(numel(keep),1);

for i = 1:numel(keep)
    f1 = bandsKeep(i,1); f2 = bandsKeep(i,2);
    y  = bandpass_iir(x, Fs, f1, f2, 4);           % IIR Butter + filtfilt
    env = abs(hilbert(y));
    [Fenv, Aenv, fr_env] = envelope_fft(env, Fs, [P.minFr P.maxFr]);
    [fr_acf, acLag, acCurve] = autocorr_envelope(env, Fs, [P.minFr P.maxFr]);
    [fr_cep, qvec, Cmag]     = cepstrum_peak(env, Fs, [P.minFr P.maxFr]);
    [fr_sb, sbScore]         = sideband_grid_score(y, Fs, mean([f1 f2]), P.sidebandWinHz, [P.minFr P.maxFr]);

    % SNR grossier de l’enveloppe à fr_env
    [~,iPk] = max(Aenv); pk = Aenv(iPk);
    noise = median(Aenv(max(iPk-10,1):min(iPk+10,numel(Aenv))));
    snr_env = 20*log10((pk+eps)/(noise+eps));

    estimates = [estimates; fr_env, fr_acf, fr_cep, fr_sb, snr_env]; %#ok<AGROW>

    bandSummaries{i} = struct('f1',f1,'f2',f2,'fr_env',fr_env,'fr_acf',fr_acf,...
        'fr_cep',fr_cep,'fr_sb',fr_sb,'snr_env',snr_env,'Fenv',Fenv,'Aenv',Aenv,...
        'acLag',acLag,'acCurve',acCurve,'qvec',qvec,'Cmag',Cmag);
end

% 3) Fusion robuste des estimations (médiane + cohérence)
%    On prend la médiane des valeurs valides (>0)
cand = estimates; cand(cand<=0) = NaN;
fr_candidates = cand(:);
fr_hat = median(fr_candidates,'omitnan');

% 4) Diagnostics + plots optionnels
D = struct();
D.bands = bands;
D.bandScores = bandScores;
D.bandsKeep = bandsKeep;
D.bandSummaries = bandSummaries;
D.estimates = estimates;

if isfield(P,'plotting') && P.plotting
    visualize_results(x, Fs, P, D, fr_hat);
end
end
