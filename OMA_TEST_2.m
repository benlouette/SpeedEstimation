close all;
clear;
Fs = 12000;
fmin        = 0.5;     % bande de recherche (Hz)
fmax        = 120;     % bande de recherche (Hz)
    
test1 = load('Normal Baseline Data 12K/Normal_0_1797_rpm.mat'); % 1796
freqRot1 = 1796/60.0; % tours/min → Hz
array1A=test1.X097_DE_time;
array1B=test1.X097_FE_time;

test2 = load('Normal Baseline Data 12K/Normal_1_1772_rpm.mat'); % 27.0880 Hz
freqRot2 = 1772/60.0; % tours/min → Hz
array2A=test2.X098_DE_time;
array2B=test2.X098_FE_time;

test3 = load('Normal Baseline Data 12K/Normal_2_1750_rpm.mat'); % 27.0880 Hz
freqRot3 = 1750/60.0; % tours/min → Hz
array3A=test3.X099_DE_time;
array3B=test3.X099_FE_time;

test4 = load('Normal Baseline Data 12K/Normal_3_1730_rpm.mat'); % 1725 Rpm
freqRot4 = 1725/60; % tours/min → Hz
array4A=test4.X100_DE_time;
array4B=test4.X100_FE_time;
% 
array   =   array1A;
freqRot = freqRot1 ;


sizeSlide=8192*2;
accOrigi=array(1:sizeSlide);
n=floor(numel(array)/sizeSlide);


[ c_q, q_axis] = coarse_cepstrum_f0(accOrigi, Fs, [0.5 120], 'Win','hann','Dither',1e-3);
c_q=c_q(1:min(numel(accOrigi),Fs/fmin));
q_axis=q_axis(1:min(numel(accOrigi),Fs/fmin));
f_axis = 1 ./ q_axis;
mask = (f_axis > fmax);
c_q(mask) = 0;
[~,idx]=max(c_q);
candf_Cepstrum = 1/q_axis(idx);


Opt = struct('rangeHz',[500 4000], 'welchWin',4096, 'welchOL',0.5, 'nfft',numel(accOrigi), ...
             'maxPeaks',20, 'minPromDB',6, 'minDistHz',10, ...
             'K',3, 'DeltaHz',[5 80], 'etaBox',0.15, 'thrSNR_dB',3, ...
             'useMedianBaseline',true, 'Delta0_small',10);

res = pick_best_carrier_comb(accOrigi, Fs, Opt);   % x: signal colonne
plot_pick_best_comb(res, Opt, 10);
M = 10;
% b) Version filtres FIR (zéro-phase)
[S_fir, info_fir] = split_top_candidates(accOrigi, Fs, res, ...
    'NumSignals',M, 'Method','fir', 'FIR_Order',512, ...
    'WidthFactor',3, 'TransWidthHz',0.15, 'K',3);
% a) Version rapide FFT/IFFT (offline)
[S_fft, info_fft] = split_top_candidates(accOrigi, Fs, res, ...
    'NumSignals',10, 'Method','fft', 'WidthFactor',3, 'K',3);

w = [0.25 0.25 0.20 0.15 0.15];   % [Qspec Qent Qsfm Qpksep Qprom]
TopK = 5;
metrics = repmat(struct( ...
    'Qspec',0,'Qent',0,'Qsfm',0,'Qpksep',0,'Qprom',0, ...
    'Score',0,'acf',[],'lags',[],'info',struct()), 1, M);

for m = 1:M
    accFir = S_fft(:,m);
    lags = 1:numel(accFir);
    signalEnvelope = getEnv(accFir, Fs, info_fft.fc_Hz(m));% res.best.fc BEST.band(1)-50);  %   extracts the envelope
    [r] = nsdf_mcleod_fast(signalEnvelope, Fs, 0, 160, numel(signalEnvelope));

    % --- métriques ACF (robustes) ---
    Mx = acf_quality_metrics2(r);         % fonction B ci-dessous
    score = acf_quality_score(Mx, w);     % fonction C ci-dessous
    metrics(m).Qspec = Mx.Qspec;
    metrics(m).Qent  = Mx.Qent;
    metrics(m).Qsfm  = Mx.Qsfm;
    metrics(m).Qpksep= Mx.Qpksep;
    metrics(m).Qprom = Mx.Qprom;
    metrics(m).Score = score;
    metrics(m).acf   = r;
    metrics(m).lags  = lags/Fs;           % en secondes
    metrics(m).info  = Mx.debug;
end

% --- tableau + tri ---
Qspec = [metrics.Qspec]';
Qent  = [metrics.Qent]';
Qsfm  = [metrics.Qsfm]';
Qpksep= [metrics.Qpksep]';
Qprom = [metrics.Qprom]';
Score = [metrics.Score]';

idx = (1:M).';
ranked = table(idx, Score, Qspec, Qent, Qsfm, Qpksep, Qprom, ...
    'VariableNames', {'Candidate','Score','Qspec','Qent','Qsfm','Qpksep','Qprom'});
ranked = sortrows(ranked, 'Score', 'descend');

% --- plots (optionnels) ---
if 1
    % 1) Barres des scores tous candidats
    figure('Color','w','Name','ACF quality – Scores');
    bar(ranked.Candidate, ranked.Score); grid on
    xlabel('Candidat'); ylabel('Score ACF (0..1)');
    title('Classement par qualité ACF (plus haut = mieux)');

    % 2) Comparaison des métriques pour TopK
    K = min(TopK, M);
    topIdx = ranked.Candidate(1:K);
    Mmat = [ [metrics(topIdx).Qspec]' [metrics(topIdx).Qent]' [metrics(topIdx).Qsfm]' ...
             [metrics(topIdx).Qpksep]' [metrics(topIdx).Qprom]' ];
    figure('Color','w','Name','ACF quality – Metrics (TopK)');
    bar(Mmat, 'grouped'); grid on
    xticklabels(string(topIdx));
    legend({'Qspec','Qent','Qsfm','Qpksep','Qprom'},'Location','northoutside','Orientation','horizontal');
    xlabel('Candidat (index colonne S)'); ylabel('Score métrique (0..1)');
    title('Comparaison des métriques ACF – TopK');

    % 3) Détail du #1 : ACF + spectre de l’ACF (débaseliné)
    c1 = ranked.Candidate(1);
    r  = metrics(c1).acf;
    tlag = metrics(c1).lags;
    figure('Color','w','Name','ACF quality – Detail Top1');
    tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

    % ACF
    nexttile; plot(tlag, r, 'LineWidth',1.2); grid on
    xlabel('Lag (s)'); ylabel('ACF (norm.)'); title(sprintf('ACF – candidat #%d', c1));

    % Spectre ACF (avec baseline médiane)
    nexttile;
    Rspec = abs(fft(r)).^2; Rspec = Rspec(1:floor(numel(r)/2));
    Lb = max(5, round(0.01*numel(Rspec))); base = movmedian(Rspec,Lb);
    Rs = max(Rspec - base, 0);
    fbin = (0:numel(Rs)-1)';  % en "bins" (pas de Fs ici: ACF de x)
    plot(fbin, 10*log10(Rs+eps), 'LineWidth',1.0); grid on
    xlabel('Bin'); ylabel('ACF-spectrum (dB, baseline-removed)');
    title('Spectre de l’ACF (débaseliné)');
end

disp(ranked(1:10,:));   % Top 10 classés par Score

%%
    %% cepstrum

    
for slide_idx = 1:M

     accFir = S_fft(:,slide_idx);

    signalEnvelope = getEnv(accFir, Fs, info_fft.fc_Hz(slide_idx));% res.best.fc BEST.band(1)-50);  %   extracts the envelope
    [AcfIn] = nsdf_mcleod_fast(signalEnvelope, Fs, 0, 160, numel(signalEnvelope));

    
    MaxLagSec =   1.0; % analyser l'ACF jusqu'à 1 s
    TopPeaks =   6; % # pics max pour la médiane des intervalles
    Prominence = 0.05;  % proéminence min (fraction du max hors lag 0)
    SmoothMs =   4.0;   % lissage avant détection de pics (ms)
    PsdFmax  =   300;

    Kharm       = 6;       % # d'harmoniques pour le peigne odd/even
    lambdaEven  = 0.95;     % pénalisation de l'énergie paire
    gammaHalf   = 1.05;    % si R(T/2) > gamma*R(T) ET even>odd => demi-tour probable
    Threshold   = 0.15;    % seuil YIN (0.1–0.2 typique)
    WeightMode  = '1/k';   % 'equal' ou '1/k'
    Eps         = 1e-3;    % epsilon dans log( eps + R+ )
    HPS_K      = 4;        % # harmoniques HPS
    HPS_dfHz   = 0.1;      % pas Hz pour HPS
    huberK = 0.8; %≈0.4–0.8
    Kharm_KSRH= 8;                               % # harmoniques
    
    Nlag   = numel(AcfIn);
    % ---------- bornes en lag ----------
    tauMinSamp = max(2, floor(Fs / fmax)); % >= 2 pour parabolique on ne cherche pas de période plus courte que 1/fmax, et jamais en dessous de 2 échantillons.
    tauMaxSamp = min(Nlag-2, ceil(Fs / fmin));% on ne cherche pas de période plus grande que 1/fmin, ni au-delà de ce que l’ACF permet.
    maxLagSamp = min(Nlag-1, round(MaxLagSec*Fs));% borne maximale imposée par un paramètre utilisateur (MaxLagSec).
    tauMaxSamp = min(tauMaxSamp, maxLagSamp);% la borne maximale effective est la plus restrictive des deux.
    
    if tauMinSamp+2 >= tauMaxSamp
        error('Plage de lag invalide: ajuster fmin/fmax/MaxLagSec.');
    end
    
    % seuil de proéminence (relatif au max hors lag 0)
    base = max(AcfIn(2:min(end,maxLagSamp+1)));
    prom = Prominence * max(base, eps);
    % impose une distance min entre pics pour éviter les doublons (≈ T/3 à fmax)
    minPkDist = max(2, round(0.33 * Fs / fmax));
    
    lags = (0:Nlag-1).';% [0 ,1 ,2 ... ,5999]
    tau  = lags / Fs;% [0 ,0.0002 ,0.0003,... ,0.999]

    AcfIn(1:floor(Fs/fmax)) = 0;
    AcfIn = normalize_0_1(AcfIn);
    


    %% Estimation de période par médiane des différences successives
%     [pkVal, pkLoc] = findpeaks(AcfIn, 'MinPeakProminence', prom, 'MinPeakDistance', minPkDist);   % pkLoc = lags, pkVal = amplitudes
%     medianVal = median(pkVal);
%     mask = (pkVal > medianVal);
%     pkVal=pkVal(mask);
%     pkLoc=pkLoc(mask);
%     L_best = median(diff(pkLoc));
%     [~, iNear] = min(abs(double(pkLoc) - L_best));
%     L = pkLoc(iNear);
%     [T_hat,L_ref] =interp_parabolic_acf(AcfIn, L, Fs, tauMaxSamp);
% 
%     candf_mediane=Fs/L_ref;

    Opt = struct('SmoothMs',4.0,'PromRel',0.2,'Parabolic',true,'AntiHalf',false);
[OUT] = period_from_acf_median_robust(AcfIn, Fs, tauMinSamp, tauMaxSamp, Opt);
candf_mediane = OUT.f_hat;

    %% Estimateur YIN/CMNDF
    [L, idxRange, dprime] = yin_pick_L_from_acf_cmndf(AcfIn, Fs, tauMinSamp, tauMaxSamp, Threshold, true, 0.10);
    % On choisit le pic le plus proche de T_samp et on récupère sa position entière L.
    [T_hat,Lref] =interp_parabolic_acf(dprime, L, Fs, tauMaxSamp);
    candf_yin=Fs/Lref;
    
    %% Estimateur par peigne-produit (log)
    opts = struct();
    opts.Kharm = Kharm;
    opts.WeightMode = WeightMode;
    opts.Eps = Eps;
    opts.gridStepSamp = 2;
    opts.needAtLeastHarm = 2;
    [Lbest, Sbest] = scorePeigneGridRefine(AcfIn, tauMinSamp, tauMaxSamp, opts);
    [T_hat,Lref] =interp_parabolic_acf(AcfIn, Lbest, Fs, tauMaxSamp);
    candf_comb = Fs/Lref;
    
    %% Estimateur HPS
    [candf_hps, Lref] = cand_hps_from_psd(AcfIn, Fs, fmin, fmax, HPS_dfHz, HPS_K);

    candf_hps = Fs/Lref;
    %% Decide
    fcandstring = ["Cepst  ","yin    ", "comb   ", "hps    ","mediane","delta  "];
    Delta = freqRot;
    fcands = [candf_Cepstrum,candf_yin, candf_comb, candf_hps,candf_mediane,Delta];
    


[f_twm, costs, grid, info] = twm_select_with_octave_grid( ...
    AcfIn, Fs, fcands, ...
    'Octaves',[0.5 1 2], 'RelPerturb',-0.04:0.01:0.04, ...
    'fminmax',[0.5 200], 'MaxGrid',250, ...
    'MaxPeaks',40, 'Kmax',18, 'rho',0.33, 'q',1.4, 'wMode','1/k');

fprintf('TWM (grid) picked: %.3f Hz (min cost=%.4f) over %d candidates\n', f_twm, min(costs), numel(grid));
%     [cSorted, idx] = sort(costs, 'ascend');
% disp(table(grid(idx(1:10)).', cSorted(1:10), 'VariableNames', {'Hz','Cost'}));

    
%     [f_cons, diag] = fuse_f0_candidates_fast(fcands, AcfIn, Fs, ...
%         'Weights', [], 'OctaveTol', 0.20, 'AntiHalf', true, 'GammaHalf', 1.05);
    % affichage des candidats
    fprintf('--- Candidats F0 ---\n');
    for i = 1:numel(fcands)
        if isfinite(fcands(i)) && fcands(i)>0
            fprintf('  cand[%s] = %8.3f Hz   (≈ %7.1f rpm)\n', ...
                fcandstring(i), fcands(i), 60*fcands(i));
        else
            fprintf('  cand[%d] =   NaN/invalid\n', i);
        end
    end
    
    %% Plot
    figureHandle =figure('Name','ACF + Candidats','Color','w');
%         subplot(4,1,1);
%     plot(diff(diff(complexityArray)));
%     xline(4);
    
    subplot(3,1,1);
    Nfft = length(c_q);
    q_axis = (0:Nfft-1)'/Fs;        % axe quefrency en s
    
    % borne max cohérente (par ex. 1/fmin ou max(tau))
    tmax = min([max(tau), max(q_axis)]);
    tau_lim = tau(tau <= tmax);
    acf_lim = AcfIn(1:numel(tau_lim));
    q_axis_lim = q_axis(q_axis <= tmax);
    c_q_lim = c_q(1:numel(q_axis_lim));
    
    plot(q_axis_lim, c_q_lim, 'r'); grid on;
%     xline(1/freqRot,'g--');
%     xline(1/candf_Cepstrum,'m:', ...
%         'Label',sprintf('Cepst %.2f Hz (%.0f rpm)', candf_Cepstrum,60*candf_Cepstrum), ...
%         'LabelOrientation','horizontal','LabelVerticalAlignment','bottom');
    xlabel('Quefrency (s)');
    ylabel('Cepstrum');
    title('Cepstrum (queferency domain)');

    
    % --- subplot 1 : ACF + lignes verticales ---
    subplot(3,1,2);
    plot(tau, AcfIn, 'b-','LineWidth',1.2); grid on; hold on;
    xline(1/freqRot,'g--','LineWidth',1.2, ...
        'Label',sprintf('Théorique %.2f Hz', freqRot), ...
        'LabelOrientation','horizontal','LabelVerticalAlignment','bottom');
    
    xlabel('\tau (s)');
    ylabel('r(\tau) (lissée)');
    title(sprintf('ACF %d + consensus final ',slide_idx));
    
    % --- subplot 2 : barres candidats ---
    subplot(3,1,3);
    vals = [freqRot,candf_Cepstrum,candf_yin, candf_comb, candf_hps, f_twm, candf_mediane, Delta];
    names= {'real','Cpst','YIN','Comb','HPS','TWM','MEDIANE','DELTA'};
    bh = bar(vals,'FaceColor',[0.2 0.6 0.8]); grid on;
    set(gca,'XTick',1:numel(vals),'XTickLabel',names);
    ylabel('Fréquence (Hz)');
    % ajouter texte au-dessus des barres
    xt = 1:numel(vals);
    for i=1:numel(vals)
        if isfinite(vals(i)) 
            
            
%             if(iscorrectionCepstrum && strcmp( names(i), 'Cpst') )
%                 text(xt(i), vals(i)+0.02*max(vals), ...
%                     sprintf('%.2f Hz\n%s',vals(i),"CorSRH"), ...
%                     'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
%             else
                text(xt(i), vals(i)+0.02*max(vals), ...
                    sprintf('%.2f Hz\n%.0f rpm', vals(i),60*vals(i)), ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
%             end
            
        else
            text(xt(i), 0, 'NaN','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
        end
    end
    set(figureHandle, 'Position', [100 100 1400 500]);
    
end

%%
function [lags, acf, info] = acf_fft_band(x, Fs, fmin, fmax, opts)
% ACF_FFT_BAND  Autocorrelation (FFT-based, linear) limited to [fmin..fmax]
%
% [lags, acf, info] = acf_fft_band(x, Fs, fmin, fmax, opts)
%
% Inputs:
%   x     : signal (vector)
%   Fs    : sampling rate [Hz]
%   fmin  : lowest fundamental of interest [Hz]   (sets tauMaxSamp)
%   fmax  : highest fundamental of interest [Hz]  (sets tauMinSamp)
%   opts  : struct with optional fields:
%           .unbiased  (logical, default=false)  % divide each lag by (L-k)
%           .normalize (logical, default=true)   % acf(0)=1
%           .window    ('none'|'hann', default='none') % window on x
%           .maxLagSec (scalar, optional)        % extra upper cap in seconds
%
% Outputs:
%   lags  : vector of lags [samples], starting at tauMinSamp..tauMaxSamp
%   acf   : ACF(tau) aligned with lags
%   info  : struct with fields:
%           .tauMinSamp, .tauMaxSamp, .L, .nfft
%
% Notes:
% - Removes DC (mean) before FFT.
% - Uses zero-padding >= 2L-1 for LINEAR ACF (not circular).
% - Use unbiased=true if you want comparability across window sizes.

%     arguments
%         x (:,1) double
%         Fs (1,1) double {mustBePositive}
%         fmin (1,1) double {mustBePositive}
%         fmax (1,1) double {mustBePositive}
%         opts.unbiased  logical = false
%         opts.normalize logical = true
%         opts.window    char    = 'none'   % or "hann"
%         opts.maxLagSec double  = Inf
%     end

% --- pre
x = x(:);
L = numel(x);
if fmax >= Fs/2
    warning('fmax >= Fs/2; setting fmax = 0.45*Fs');
    fmax = 0.45*Fs;
end
if fmin <= 0
    error('fmin must be > 0');
end
if fmax >= fmin
    warning('fmax >= fmin; swapping to enforce fmax < fmin.');
    tmp = fmax; fmax = fmin; fmin = tmp;
end

% --- tau bounds from band
tauMinSamp = max(1, floor(Fs / fmax));             % smallest lag we trust
tauMaxBand = floor(Fs / fmin);                      % largest lag from fmin
tauMaxUser = min(L-1, floor(opts.maxLagSec * Fs));  % user/size cap
tauMaxSamp = min([tauMaxBand, tauMaxUser, L-1]);

if tauMinSamp > tauMaxSamp
    % Degenerate band -> fall back to [1 .. min(ceil(L/2), L-1)]
    tauMinSamp = 1;
    tauMaxSamp = min(ceil(L/2), L-1);
    %         warning('Band limits collapsed; using fallback tau range [%d..%d].', ...
    %                  tauMinSamp, tauMaxSamp);
end

% --- optional window on x (avoid strong edge effects if needed)
switch lower(opts.window)
    case 'hann'
        w = hann(L);
        xw = (x - mean(x)) .* w;
    otherwise
        xw = x - mean(x);   % remove DC is crucial
end

% --- FFT-based linear ACF
nfft = 2^nextpow2(2*L - 1);
X = fft(xw, nfft);
S = X .* conj(X);
r = ifft(S, 'symmetric');   % linear ACF thanks to zero-padding
r = r(1:L);                 % lags >= 0

% --- biased / unbiased
if opts.unbiased
    k = (0:L-1).';
    denom = max(L - k, 1);
    r = r ./ denom;
end

% --- normalization (so acf(0)=1)
if opts.normalize
    r0 = max(r(1), eps);
    r = r / r0;
end

% --- crop to banded lags [tauMinSamp .. tauMaxSamp]
lags = (tauMinSamp:tauMaxSamp).';
acf  = r(lags + 1);   % +1 because r(1) is lag 0

% --- info
info = struct('tauMinSamp',tauMinSamp, 'tauMaxSamp',tauMaxSamp, ...
    'L',L, 'nfft',nfft, 'fmin',fmin, 'fmax',fmax);
end



function acf = getAcf(SigEnv)
% Compute the Autocorrelation Function (ACF) of the Signal Envelope

% Perform the Fourier Transform of the Signal Envelope
fftSignal = fft(SigEnv);

% Compute the Power Spectrum
% This is the modulus squared of the FFT, which gives the power at each frequency
powerSpectrum = abs(fftSignal) .^ 2;

% Compute the Inverse Fourier Transform of the Power Spectrum
% This yields the autocorrelation function of the signal
acfFull = ifft(powerSpectrum);

% Extract the first half of the ACF
% The ACF is symmetric, so the first half contains all the unique information
acf = acfFull(1:round(numel(acfFull)/2));
end

function xEnv = getEnv(signal, Fs, cutOff)
% Extract the Envelope of a Signal Using Complex Demodulation
% Parameters:
%   signal: The input signal from which to extract the envelope.
%   Fs: Sampling frequency of the signal.
%   cutOff: Cut-off frequency for the low-pass filter.

% Define the filter order
FilterOrder = 50;

% Create a time vector corresponding to the signal
t = (0:length(signal) - 1)' / Fs;

% Perform complex demodulation
% This shifts the signal frequency content around the cut-off frequency to baseband
x0 = signal .* exp(-1i * 2 * pi * cutOff * t);

% Design a low-pass FIR filter
% The effective cut-off frequency is halved as the signal is shifted to baseband
b = fir1(FilterOrder,  cutOff/2 / Fs/2);

% Apply the low-pass filter
% The filter is applied to extract the low-frequency components (envelope)
xEnv = conv2(x0, b(:), 'same') * 2;

% Compute the magnitude of the complex signal
% This step extracts the envelope of the original signal
xEnv = abs(xEnv);

% Remove the DC component from the envelope
% This centers the envelope around zero, removing any constant ofFset
xEnv = xEnv - mean(xEnv);
end


function normalizedData = normalize_0_1(data)
% Find the minimum and maximum values of the input data
dataMin = min(data);
dataMax = max(data);

% Check if dataMax and dataMin are the same (to avoid division by zero)
if dataMax == dataMin
    % If all data points are equal, return a zero array of the same size
    normalizedData = zeros(size(data));
else
    % Normalize the data to the range [0, 1]
    normalizedData = (data - dataMin) / (dataMax - dataMin);
end
end



function acfDetrended = detrendHodrickPrescott(acf)
% Hodrick-Prescott Filter for Detrending ACF
% -------------------------------------------

% Initialize the autocorrelation function
Y = acf(:);

% Smoothing Parameter
% Note: The smoothing value should be in accordance with the sampling rate.
smoothing = 1e+6;

% Number of Observations
numObs = numel(acf);

% Constructing the Hodrick-Prescott Filter Coefficients
% Normally this would be a sparse matrix for efficiency, but here we use a full matrix.
e = repmat([smoothing, -4 * smoothing, (1 + 6 * smoothing), -4 * smoothing, smoothing], numObs, 1);
A = spdiags(e, -2:2, numObs, numObs);
% Adjusting the first and last rows of matrix A for boundary conditions
A(1,1) = 1 + smoothing;
A(1,2) = -2 * smoothing;
A(2,1) = -2 * smoothing;
A(2,2) = 1 + 5 * smoothing;
A(numObs-1, numObs-1) = 1 + 5 * smoothing;
A(numObs-1, numObs) = -2 * smoothing;
A(numObs, numObs-1) = -2 * smoothing;
A(numObs, numObs) = 1 + smoothing;

% LU Decomposition for Solving the Filter
% L for lower triangular matrix and U for upper triangular matrix
[L, U, P] = lu(A);

% Solving the System Using LU Decomposition
y = L \ (P * Y);  % Solving Ly = Pb using forward substitution
trend = U \ y;    % Solving Ux = y using backward substitution

% Detrending the ACF Signal
acfDetrended = Y - trend;  % Subtracting the trend component
end

function acfHuber = huber(acf)
% Median and Mean Absolute Deviation (MAD) Calculation
mu = median(acf);  % Median of the detrended ACF
s = mean(abs(acf - mu));  % Mean absolute deviation around the median

% Applying the Huber Function
% The Huber function is a robust method of scaling data.
x = (acf' - mu) / s;
result = sign(x) .* min(abs(x), s);  % Applying Huber function to scale the data

% Applying the Huber Function Results to the ACF
acfHuber = acf' .* result;
end

% --------- sous-fonctions locales ---------
function y = movavg(x,w)
% moyenne mobile (filtfilt si possible, sinon conv), colonne en sortie
x = double(x(:));
k = ones(w,1)/w;
try
    y = filtfilt(k,1,x);
catch
    y = conv(x,k,'same');
end
y = cast(y, 'like', x);
end



function [pks, locs] = findpeaks_safe(x, varargin)
% utilise findpeaks si dispo (Signal Toolbox), sinon fallback
try
    [pks, locs] = findpeaks(x, varargin{:});
catch
    p = inputParser;
    p.addParameter('MinPeakProminence', 0);
    p.addParameter('MinPeakDistance', 1);
    p.parse(varargin{:});
    prom = p.Results.MinPeakProminence;
    mind = p.Results.MinPeakDistance;
    locs = [];
    for i=2:numel(x)-1
        if x(i)>x(i-1) && x(i)>x(i+1) && x(i)>prom
            if isempty(locs) || (i - locs(end)) >= mind
                locs(end+1)=i; %#ok<AGROW>
            end
        end
    end
    pks = x(locs);
end
end



% ===== PSD à partir de l'ACF (Wiener–Khinchin) =====
function [F, Pxx] = psd_from_acf(acf_n, Fs)
% ACF symétrisée : r[-k] = r[k]
acf_sym = [flipud(acf_n(2:end)); acf_n(:)];
% fenêtre Hann (évite les fuites dues au troncage)
Nlag = numel(acf_sym);
n = (0:Nlag-1).';
w = 0.5 - 0.5*cos(2*pi*n/(Nlag-1)); % fenetre de Hann
R = acf_sym .* w; % appliques cette fenêtre à l’ACF avant de faire la FFT :
% PSD = FFT(R) réelle (procès réel)
S = fft(R);
Pxx = real(S(1:floor(Nlag/2)+1));
F   = (0:floor(Nlag/2)).' * (Fs/Nlag);
end


function [T_hat, Lref] = interp_parabolic_acf(acf_s, L, Fs, tauMaxSamp)
% INTERP_PARABOLIC_ACF
% Raffine la position d'un pic ACF par interpolation parabolique (3 points).
% - Fallback sûr si proche des bords ou si la courbure n'est pas concave.
%
% In:
%   acf_s       : ACF (déjà lissée de préférence), vecteur colonne/ligne
%   L           : indice entier du pic (1-based, tel que acf_s(L) est le pic)
%   Fs          : fréquence d'échantillonnage (Hz)
%   tauMaxSamp  : borne max utile en lag (en échantillons)
%
% Out:
%   T_hat       : période raffinée (secondes)  = Lref / Fs
%   Lref        : position raffinée du pic (échantillons, fractionnaire)

acf_s = acf_s(:);                 % assure colonne
Nlag     = numel(acf_s);
tauMaxSamp = min(tauMaxSamp, Nlag-1);

% --- Cas bords : pas assez de points (il faut L-1, L, L+1)
if ~(L >= 2 && L <= tauMaxSamp-1)
    Lref  = double(L);            % fallback : pas d'interp
    T_hat = Lref / Fs;
    return;
end

% --- Parabole 3 points
y1 = double(acf_s(L-1));
y2 = double(acf_s(L));
y3 = double(acf_s(L+1));

% Conditions de validité :
% 1) y2 au moins aussi grand que ses voisins (pic local)
% 2) courbure concave : (y1 - 2*y2 + y3) < 0
if (y2 < y1) && (y2 < y3)
    % clair creux -> pas d'interp
    Lref  = double(L);
    T_hat = Lref / Fs;
    return;
end
den = (y1 - 2*y2 + y3);
if ~(den < 0)
    % pas un maximum concave exploitable -> fallback
    Lref  = double(L);
    T_hat = Lref / Fs;
    return;
end

% Décalage sous-échantillon (théorie parabole)
delta = 0.5 * (y1 - y3) / den;

% Clamp pour éviter sauts aberrants quand les 3 points sont plats/bruités
delta = max(min(delta, 0.5), -0.5);

% Position raffinée
Lref  = double(L) + delta;
% Sécurité bornes (en cas d'arrondi limite)
Lref  = max(1, min(Lref, double(tauMaxSamp)));

% Période (s)
T_hat = Lref / Fs;
end



function [L, idxRange, dprime] = yin_pick_L_from_acf_cmndf(r_ns, Fs, tauMinSamp, tauMaxSamp, thresh, useAdaptive, winProp)
% YIN_PICK_L_FROM_ACF_CMNDF
%  Sélection du lag fondamental via YIN/CMNDF depuis une ACF déjà
%  normalisée & lissée (r_ns). Ne fait PAS de parabole ici.
%
% In:
%   r_ns        : ACF normalisée ET lissée (colonne/ligne), r_ns(1)=1
%   Fs          : Hz (utilisé seulement pour cohérence d'unités; ici pas critique)
%   tauMinSamp  : lag mini (samples)   ~ floor(Fs/fmax), forcé >= 2
%   tauMaxSamp  : lag maxi (samples)   ~ ceil(Fs/fmin)
%   thresh      : seuil YIN (p.ex. 0.10–0.20)
%   useAdaptive : true/false → abaisse le seuil si pas de crossing
%   winProp     : largeur de fenêtre locale autour du crossing (ex: 0.10)
%
% Out:
%   L        : lag entier choisi (samples)
%   idxRange : [i1 i2] fenêtre locale utilisée pour le min (utile pour ta parabole)
%   dprime   : CMNDF (utile pour debug/plots en aval)
%
% Hypothèses:
%  - r_ns est déjà pré-traitée (centrage/normalisation, lissage) ailleurs.
%  - L’interpolation parabolique sera appliquée ailleurs à partir de L.

% --- garde-fous sur la bande ---
r_ns = r_ns(:);
Nlag = numel(r_ns);
tauMinSamp = max(2, floor(tauMinSamp));
tauMaxSamp = min(max(tauMinSamp+2, floor(tauMaxSamp)), Nlag-1);

% --- CMNDF: d(τ)=2*(1-r(τ)), d'(τ)=d(τ)*τ / cumsum(d) ---
tauIdx = (1:tauMaxSamp).';
d      = 2.0 * (1 - r_ns(tauIdx+1));
cum    = cumsum(d);
dprime = d .* (tauIdx ./ max(cum, eps));

% --- masque bande (ignore < tauMinSamp) ---
dprime(1:tauMinSamp-1) = +Inf;

% --- seuil absolu + option adaptative ---
idx = find(dprime < thresh, 1, 'first');
if isempty(idx) && useAdaptive
    % seuil adaptatif simple vers le minimum observé dans la bande
    dband = dprime;
    dband(~isfinite(dband)) = +Inf;
    mval  = min(dband);
    th2   = max(0.05, min(0.95, 0.5*(thresh + mval)));
    idx   = find(dprime < th2, 1, 'first');
end

% --- sélection du minimum local autour du premier crossing (ou min global) ---
if ~isempty(idx)
    firstIdx = idx;                                   % déjà en indices absolus car dprime est tronqué seulement en <tauMin
    win      = max(3, round(winProp * firstIdx));     % ex: 10% de τ
    i1       = max(tauMinSamp, firstIdx - win);
    i2       = min(tauMaxSamp, firstIdx + win);
    [~,rel]  = min(dprime(i1:i2));
    L        = i1 + rel - 1;
    idxRange = [i1 i2];
else
    % pas de crossing → minimum global dans la bande
    [~,L] = min(dprime);
    if L<tauMinSamp || L>tauMaxSamp
        [~,L] = min(dprime(tauMinSamp:tauMaxSamp));
        L = L + tauMinSamp - 1;
    end
    idxRange = [max(tauMinSamp, L-3), min(tauMaxSamp, L+3)]; % petite fenêtre par défaut
end
end







function [Lbest, Sbest] = scorePeigneGridRefine(acf_s, tauMinSamp, tauMaxSamp, opts)
% SCOREPEIGNEGRIDREFINE
%  Peigne log sur ACF lissée/normalisée (positif), échantillonné sur une grille grossière,
%  puis raffinement local à pas 1 et interpolation parabolique du score.
%
% In:
%   acf_s        : ACF (déjà lissée de préférence), colonne/ligne
%   tauMinSamp   : lag min (samples)
%   tauMaxSamp   : lag max (samples)
%   opts struct  :
%       .Kharm        (default 6)
%       .WeightMode   (default '1/k')   % or 'equal'
%       .Eps          (default 1e-3)
%       .gridStepSamp (default 2)       % grille grossière (>=1)
%       .needAtLeastHarm (default 2)    % nb mini d'harmoniques
%
% Out:
%   Lbest  : lag entier (échantillons) après raffinement local
%   Sbest  : score peigne au meilleur L (sans parabolique)

if nargin<4, opts = struct(); end
if ~isfield(opts,'Kharm'),        opts.Kharm = 6; end
if ~isfield(opts,'WeightMode'),   opts.WeightMode = '1/k'; end
if ~isfield(opts,'Eps'),          opts.Eps = 1e-3; end
if ~isfield(opts,'gridStepSamp'), opts.gridStepSamp = 2; end
if ~isfield(opts,'needAtLeastHarm'), opts.needAtLeastHarm = 2; end

acf_s = acf_s(:);
Nlag  = numel(acf_s);
tauMaxSamp = min(tauMaxSamp, Nlag-1);

% ---- Poids ----
K  = opts.Kharm;
if strcmpi(opts.WeightMode,'equal')
    wk = ones(1,K);
else
    wk = 1./(1:K);
end
wk = wk / sum(wk);

% ---- Préparation ----
rpos   = max(acf_s, 0);
logeps = @(x) log(opts.Eps + x);

% ---- Grille grossière ----
gridStep = max(1, round(opts.gridStepSamp));
tauGrid  = (tauMinSamp:gridStep:tauMaxSamp).';
Sg = numel(tauGrid);
ScoreG = -inf(Sg,1);

for s=1:Sg
    Lg   = tauGrid(s);
    kmax = min(K, floor(tauMaxSamp / Lg));
    if kmax < opts.needAtLeastHarm, continue; end
    used = 1:kmax;
    wloc = wk(used); wloc = wloc / sum(wloc);
    idx  = used.*Lg + 1;        % +1 : r(1) = lag0
    vals = rpos(idx);
    ScoreG(s) = dot(wloc(:), logeps(vals(:)));
end

% Si rien d'exploitable sur la grille → fallback pic ACF
[Sm, imax] = max(ScoreG);
if ~isfinite(Sm)
    [~, Lbest] = max(acf_s(tauMinSamp:tauMaxSamp));
    Lbest = Lbest + tauMinSamp - 1;
    Lref  = double(Lbest);
    Sbest = -inf;
    return;
end
Lcand = tauGrid(imax);

% ---- Raffinement local à pas 1 (±gridStep) ----
L1 = max(tauMinSamp, Lcand - gridStep);
L2 = min(tauMaxSamp, Lcand + gridStep);
Lvec = (L1:L2);
ScoreFine = -inf(numel(Lvec),1);
for ii=1:numel(Lvec)
    L = Lvec(ii);
    kmax = min(K, floor(tauMaxSamp / L));
    if kmax < opts.needAtLeastHarm, continue; end
    used = 1:kmax;
    wloc = wk(used); wloc = wloc / sum(wloc);
    idx  = used.*L + 1;
    vals = rpos(idx);
    ScoreFine(ii) = dot(wloc(:), logeps(vals(:)));
end

[Sm2, iFine] = max(ScoreFine);
Lbest = Lvec(iFine);
Sbest = Sm2;
end

function V = compute_correntropy_signal(x, tauMax, sigma)
% V(τ) = mean( exp( - (x[n]-x[n-τ])^2 / (2σ^2) ) )
x = double(x(:));
N = numel(x);
V = zeros(tauMax+1,1);
V(1) = 1;  % à τ=0, la correntropy vaut 1 (kernel à distance nulle)
for t = 1:tauMax
    M = N - t;
    if M<=0, break; end
    d = x(1+t:N) - x(1:N-t);
    V(1+t) = mean( exp( - (d.^2) / (2*sigma^2) ) );
end
end

function nsdf = nsdf_mcleod_fast(x, Fs, fmin, fmax, MaxLagSec)
% NSDF_MCLEOD_FAST  —  NSDF (McLeod) rapide via ACF par FFT (linéaire).
%
% [nsdf, tau, OUT] = nsdf_mcleod_fast(x, Fs, fmin, fmax, 'Name',Value,...)
%
% In:
%   x     : signal (colonne/ligne)
%   Fs    : Hz
%   fmin  : Hz (borne basse)   -> tauMax = ceil(Fs/fmin)
%   fmax  : Hz (borne haute)   -> tauMin = max(2,floor(Fs/fmax))
%   MaxLagSec  (default 1.0)  : borne dure en secondes
%
% Out:
%   nsdf : vecteur NSDF(τ), τ=0..tauMax (nsdf(1)=0 par convention ici)

%
% Complexité: O(N log N) (ACF via FFT) + O(tauMax) (dénominateur).


% ---- prep ----
x = double(x(:));

x = x - mean(x);

N = numel(x);

% ---- ACF linéaire via FFT (Wiener–Khinchin + zero-pad >= 2N-1) ----
nfft = 2^nextpow2(2*N-1);
X = fft(x, nfft);
S = X .* conj(X);
rfull = ifft(S, 'symmetric');        % taille nfft, corrélation circulaire zero-paddée
R = rfull(1:N);                      % r(τ), τ=0..N-1 (ACF linéaire pour τ>=0)

% ---- bornes lag utiles ----
tauMax_from_fmin = ceil(Fs / max(fmin, eps));
tauMax_hard      = min(N-1, round(MaxLagSec * Fs));
tauMax = min(tauMax_hard, tauMax_from_fmin);

% ---- cumuls d'énergie pour dénominateur ----
x2 = x.^2;
E  = cumsum(x2);          % E(k) = sum_{n=1..k} x[n]^2
E0 = E(end);

% ---- NSDF = 2*R(τ) / (E0 + E(N-τ))  pour τ=1..tauMax, nsdf(0)=1 ----
nsdf = zeros(tauMax+1,1);
nsdf(1) = 1;  % τ=0
for T = 1:tauMax
    denom = E0 + E(N-T);               % somme d'énergies “alignées”
    nsdf(1+T) = (2*R(1+T)) / max(denom, eps);
end

end



function [fhat , L0]= cand_hps_from_psd(r, Fs, fmin, fmax, dfHz, K)
% PSD via Wiener–Khinchin (ACF -> FFT), puis Harmonic Sum sur une grille f0.
[F,Pxx] = psd_from_acf(r, Fs);
mask = (F >= fmax-1);          % indices au-dessus de fmax
Pxx(mask) = 0;              % mise à zéro
fgrid = (fmin:dfHz:fmax).';
Spec = @(f) interp1(F, Pxx, f, 'linear', 0);
S = zeros(size(fgrid));
for i=1:numel(fgrid)
    if(fgrid(i) > fmax)
        break;
    end
    f0 = fgrid(i); s=0; cnt=0;
    for k=1:K
        fk = k*f0; if fk > F(end), break; end
        s = s + Spec(fk); cnt = cnt + 1;
    end
    if cnt>=2, S(i)=s; else, S(i)=-Inf; end
end

[~,ix] = max(S);
% raffinement parabolique
if ix>1 && ix<numel(S)
    y1=S(ix-1); y2=S(ix); y3=S(ix+1);
    delta = 0.5*(y1-y3)/max(y1-2*y2+y3, eps); delta = max(min(delta,0.5),-0.5);
else
    delta = 0;
end
fhat = fgrid(ix) + delta*dfHz;
L0 = (Fs/fhat);
end


function [ c_q, q_axis] = coarse_cepstrum_f0(x, Fs, frange, varargin)
% COARSE_CEPSTRUM_F0  — estimation grossière de f0 par cepstre réel
% In:
%   x       : signal (colonne/ligne)
%   Fs      : Hz
%   frange  : [fmin fmax] (Hz) bande plausible de f0
% Options (name/value):
%   'Nfft'    ([])   : longueur FFT (def: nextpow2(length(x))*2)
%   'Win'     ('hann'): fenêtre ('hann','hamming', vecteur,...)
%   'Dither'  (0)    : bruit blanc additif RMS relatif (ex 1e-3) pour stabiliser log
%
% Out:
%   candf_Cepstrum, T_coarse
%   conf    : confiance (0..1) basée sur rapport pic/baseline
%   q_peak  : indice quefrency du pic
%   c_q     : cepstre (vector)
%   q_axis  : axe quefrency en secondes

p = inputParser;
addParameter(p,'Nfft',[]);
addParameter(p,'Win','hann');
addParameter(p,'Dither',0);
parse(p,varargin{:});
Nfft   = p.Results.Nfft;
winArg = p.Results.Win;
dith   = p.Results.Dither;

x = x(:);
N = numel(x);
if isempty(Nfft)
    Nfft = 2^nextpow2(N*2); % sur-échantillonnage fréquentiel léger
end

% Fenêtre
if ischar(winArg)
    switch lower(winArg)
        case 'hann',    w = hann(N,'periodic');
        case 'hamming', w = hamming(N,'periodic');
        otherwise,      w = hann(N,'periodic');
    end
elseif isnumeric(winArg)
    if numel(winArg)==N, w = winArg(:);
    else, w = hann(N,'periodic'); end
else
    w = hann(N,'periodic');
end

xw = x .* w;

% Dither optionnel pour stabiliser log(|X|)
if dith>0
    xw = xw + dith*std(xw)*randn(size(xw));
end

% FFT -> log magnitude -> cepstre réel
X   = fft(xw, Nfft);
mag = abs(X);
mag(mag==0) = eps;               % éviter log(0)
C   = real(ifft(log(mag)));      % cepstre réel
c_q = C(:);
q_axis = (0:Nfft-1).' / Fs;      % secondes

% % Fenêtre de recherche en quefrency
% fmin = max(frange(1), 1/N);      % garde-fou
% fmax = min(frange(2), Fs/2);
% qmin = 1/fmax; qmax = 1/fmin;
% 
% i1 = max(2, floor(qmin*Fs));     % éviter q=0
% i2 = min(floor(qmax*Fs), floor(Nfft/2)-1);
% 
% if i2 <= i1
%     % fallback trivial
%     [~, im] = max(c_q(2:floor(Nfft/2)));
%     q_peak  = im+1;
% else
%     % lissage léger pour robustesse
%     Cw = movmean(c_q, 3);
%     [~, im] = max(Cw(i1:i2));
%     q_peak  = i1 + im - 1;
% end
% 
% T_coarse = q_peak / Fs;
% candf_Cepstrum = 1 / max(T_coarse, eps);
% 
% % Confiance: contraste pic / baseline locale
% band = max(5, round(0.02*(i2-i1)));
% lo = max(i1, q_peak - band); hi = min(i2, q_peak + band);
% win = true(i2-i1+1,1); win((lo-i1+1):(hi-i1+1)) = false;
% seg = c_q(i1:i2);          % extrait la tranche utile
% baseline = median(seg(win));
% peakval  = c_q(q_peak);
% conf = max(0, min(1, (peakval - baseline) / (abs(peakval)+eps)));

end


function plotFFT(acc,Fs)
figure;
Nlag = length(acc);

% FFT
Y = fft(acc);
f = ((0:Nlag-1) * Fs / Nlag);
P = abs(Y) / Nlag;

% Affichage (jusqu'à Nyquist)
plot(f(1:Nlag/2), P(1:Nlag/2))
end

function  [mob , com] = m_hjorth_params( x)
%Calculate Hjorth mobility and complexity on given axis.
%     .. versionadded: 0.1.3
%     Parameters
%     ----------
%     x : list or np.array
%         1D or N-D data.
%     axis : int
%         The axis along which to perform the computation. Default is -1 (last).
%     Returns
%     -------
%     mobility, complexity : float
%         Mobility and complexity parameters.
%     Notes
%     -----
%     Hjorth Parameters are indicators of statistical properties used in signal
%     processing in the time domain introduced by Bo Hjorth in 1970. The
%     parameters are activity, mobility, and complexity. EntroPy only returns the
%     mobility and complexity parameters, since activity is simply the variance
%     of :math:`x`, which can be computed easily with :py:func:`numpy.var`.
%     The **mobility** parameter represents the mean frequency or the proportion
%     of standard deviation of the power spectrum. This is defined as the square
%     root of variance of the first derivative of :math:`x` divided by the
%     variance of :math:`x`.
%     The **complexity** gives an estimate of the bandwidth of the signal, which
%     indicates the similarity of the shape of the signal to a pure sine wave
%     (where the value converges to 1). Complexity is defined as the ratio of
%     the mobility of the first derivative of :math:`x` to the mobility of
%     :math:`x`.
%     References
%     ----------
%     - https://en.wikipedia.org/wiki/Hjorth_parameters
%     - https://doi.org/10.1016%2F0013-4694%2870%2990143-4
% Calculate derivatives
dx = diff(x);
ddx = diff(dx);
% Calculate variance
x_var = var(x);  % = activity
dx_var = var(dx);
ddx_var = var(ddx);
% Mobility and complexity
mob = sqrt(dx_var / x_var);
com = sqrt(ddx_var / dx_var) / mob;
end



function [f_hat, costs, diag] = twm_select_f0(x, Fs, fcands, varargin)
% TWM_SELECT_F0  — arbitre f0 via Two-Way Mismatch (Maher & Beauchamp)
% In:
%   x       : signal
%   Fs      : Hz
%   fcands  : vecteur de candidats f0 (Hz)
% Options:
%   'Nfft'      ([])      : taille FFT (auto: 2^nextpow2( len*2 ))
%   'MaxPeaks'  (40)      : nb max de pics spectre
%   'Kmax'      (20)      : nb max d'harmoniques simulés
%   'rho'       (0.33)    : pondération terme amplitude
%   'q'         (1.4)     : puissance pour l’écart fréquentiel
%   'wMode'     ('1/k')   : poids harmoniques ('1/k' ou '1')
%
% Out:
%   f_hat : f0 choisi (coût minimal)
%   costs : coût TWM pour chaque candidat (même ordre que fcands)
%   diag  : struct (pics spectraux, paramètres, etc.)

p = inputParser;
addParameter(p,'Nfft',[]);
addParameter(p,'MaxPeaks',40);
addParameter(p,'Kmax',20);
addParameter(p,'rho',0.33);
addParameter(p,'q',1.4);
addParameter(p,'wMode','1/k');
parse(p,varargin{:});
Nfft = p.Results.Nfft; MaxPeaks=p.Results.MaxPeaks;
Kmax = p.Results.Kmax; rho=p.Results.rho; q=p.Results.q; wMode=p.Results.wMode;

x = x(:);
x = x - mean(x); s = std(x); if s>0, x = x/s; end

% --- pics spectraux (fréquence/amplitude) ---
[F, A] = spectrum_peaks(x, Fs, Nfft, MaxPeaks);   % F: Hz (croissant), A: linéaire
if isempty(F), f_hat = NaN; costs = []; diag = struct(); return; end
A1 = max(A);  % ref amplitude

% --- coût TWM pour chaque candidat ---
fc = fcands(:);
fc = fc(isfinite(fc) & fc>0);
if isempty(fc), f_hat = NaN; costs = []; diag = struct(); return; end

costs = zeros(numel(fc),1);
for i=1:numel(fc)
    f0 = fc(i);
    costs(i) = twm_cost_one(f0, F, A, Fs, Kmax, rho, q, wMode, A1);
end

% meilleur
[~,ix] = min(costs);
f_hat = fc(ix);

diag = struct('F',F,'A',A,'params',p.Results,'fcands',fc,'costs',costs);
end

% -------------------------------------------------------------------------

function C = twm_cost_one(f0, F, A, Fs, Kmax, rho, q, wMode, Aref)
% Coût TWM pour un candidat f0
% (1) Mismatch des harmoniques prédits -> pics mesurés (PM)
% (2) Mismatch des pics mesurés -> harmoniques prédits (MP)
% coût = (PM + MP)/2
if f0<=0, C = inf; return; end
H = (1:Kmax)';                 % indices d’harmoniques
Hf = H * f0;                   % fréquences harmoniques (Hz)
Hf = Hf(Hf < Fs/2);            % on ne garde que sous Nyquist
if isempty(Hf), C = inf; return; end

% poids
switch lower(wMode)
    case '1/k', w = 1./(1:numel(Hf))';
    otherwise,  w = ones(numel(Hf),1);
end
w = w / sum(w);

% (1) harmoniques -> pics : pour chaque h, cherche pic le plus proche
PM = 0;
for k = 1:numel(Hf)
    [df, j] = min(abs(F - Hf(k)));      % écart freq
    % terme fréquence (normalisé) + terme amplitude
    Ef = (df / max(Hf(k), eps))^q;
    Ea = rho * abs( (A(j) - Aref) / max(Aref,eps) ); % simple écart relatif
    PM = PM + w(k)*(Ef + Ea);
end

% (2) pics -> harmoniques : pour chaque pic, proche harmonique
MP = 0;
for j = 1:numel(F)
    [df, k] = min(abs(Hf - F(j)));
    Ef = (df / max(F(j), eps))^q;
    Ea = rho * abs( (A(j) - Aref) / max(Aref,eps) );
    % poids côté harmonic index
    wk = (k<=numel(w)) * w(k) + (k>numel(w)) * w(end);
    MP = MP + wk*(Ef + Ea);
end

C = 0.5*(PM + MP);
end

% -------------------------------------------------------------------------

function [Fpeaks, Apeaks] = spectrum_peaks(x, Fs, Nfft, MaxPeaks)
% Spectre + détection simple de pics dominants (pour TWM)
if isempty(Nfft)
    Nfft = 2^nextpow2(numel(x)*2);
end
w = hann(numel(x),'periodic');
X = fft(w(:).*x, Nfft);
mag = abs(X(1:floor(Nfft/2)+1));
fax = (0:floor(Nfft/2))' * (Fs/Nfft);

% lissage léger pour robustesse
magS = movmean(mag, 3);

% pics (exclut DC)
[pv,pl] = findpeaks(magS, 'NPeaks', MaxPeaks, 'SortStr','descend', ...
                    'MinPeakDistance', max(3, round(0.005*Nfft/Fs*Fs))); %#ok<NASGU>
pl = sort(pl, 'ascend');
Fpeaks = fax(pl);
Apeaks = mag(pl);

% sécurité: enlever DC et > Nyquist (déjà géré)
m = Fpeaks>0 & isfinite(Fpeaks) & isfinite(Apeaks);
Fpeaks = Fpeaks(m); Apeaks = Apeaks(m);
end

function [f_best, costs, candGrid, diag] = twm_select_with_octave_grid(x, Fs, fcands, varargin)
% TWM_SELECT_WITH_OCTAVE_GRID
%  Étend fcands avec {×0.5, ×1, ×2} et petites perturbations relatives,
%  filtre par [fmin,fmax], déduplique, puis arbitre via TWM.
%
% In:
%   x, Fs     : signal, Hz
%   fcands    : vecteur des candidats initiaux (Hz)
% Options:
%   'Octaves'      ([0.5 1 2])     : multiplicateurs d’octave
%   'RelPerturb'   ([-0.05:0.01:0.05]) : perturbations relatives (±5%)
%   'fminmax'      ([0.2 500])     : bande valide en Hz
%   'MaxGrid'      (200)           : nb max de candidats après dédup
%   'Nfft'         ([])            : passe à twm_select_f0
%   'MaxPeaks'     (40)            : "
%   'Kmax'         (20)            : "
%   'rho'          (0.33)          : "
%   'q'            (1.4)           : "
%   'wMode'        ('1/k')         : "
%
% Out:
%   f_best   : f0 choisi (Hz)
%   costs    : coûts TWM alignés à candGrid
%   candGrid : grille finale de candidats (Hz)
%   diag     : struct (grille brute, masque, diag TWM)

p = inputParser;
addParameter(p,'Octaves',[0.5 1 2]);
addParameter(p,'RelPerturb',-0.05:0.01:0.05);
addParameter(p,'fminmax',[0.2 500]);
addParameter(p,'MaxGrid',200);
addParameter(p,'Nfft',[]);
addParameter(p,'MaxPeaks',40);
addParameter(p,'Kmax',20);
addParameter(p,'rho',0.33);
addParameter(p,'q',1.4);
addParameter(p,'wMode','1/k');
parse(p,varargin{:});
octs  = p.Results.Octaves;
drels = p.Results.RelPerturb;
fmm   = p.Results.fminmax;
MaxG  = p.Results.MaxGrid;

% 0) nettoyage
fc = fcands(:);
fc = fc(isfinite(fc) & fc>0);

% 1) grille (octaves × perturbations relatives, en log pour symétrie)
candList = [];
for m = octs(:).'
    fbase = fc * m;
    % perturbations relatives en log -> f * exp(log(1+dr))
    for dr = drels(:).'
        candList = [candList; fbase .* (1+dr)]; %#ok<AGROW>
    end
end

% 2) filtrage bande & dé-duplication (en log-f)
candList = candList(isfinite(candList) & candList>0);
candList = candList(candList>=fmm(1) & candList<=fmm(2));
if isempty(candList), f_best=NaN; costs=[]; candGrid=[]; diag=struct(); return; end
lf = log(candList);
[lfu, ia] = uniquetol(lf, 1e-3);     % ~0.1% en fréquence
candGrid = exp(lfu);                  % grille unique
% limiter la taille si nécessaire (garde les plus “centrés”)
if numel(candGrid) > MaxG
    mid = median(log(fc)); 
    [~,ord] = sort(abs(candGrid - exp(mid)),'ascend');
    candGrid = candGrid(ord(1:MaxG));
end

% 3) TWM sur la grille
[f_best, costs, d] = twm_select_f0(x, Fs, candGrid, ...
    'Nfft',p.Results.Nfft, 'MaxPeaks',p.Results.MaxPeaks, ...
    'Kmax',p.Results.Kmax, 'rho',p.Results.rho, ...
    'q',p.Results.q, 'wMode',p.Results.wMode);

% 4) diag
diag = struct();
diag.grid_raw      = candList;
diag.grid_unique   = candGrid;
diag.twm_diag      = d;
diag.params        = p.Results;
end


% 
% function out = acf_quality_metrics(r, fs, varargin)
% % out = acf_quality_metrics(r, fs, 'Tau', tau, 'SearchHz', [fmin fmax], ...
% %                           'Normalize', true, 'MaxHarm', 3)
% %
% % Calcule des métriques de "qualité" à partir de l'ACF uniquement.
% %
% % ENTRÉES
% %   r    : vecteur ACF. Peut être (i) symétrique (lags négatifs->positifs)
% %          ou (ii) limité aux lags >= 0, avec le premier élément à lag=0.
% %   fs   : fréquence d'échantillonnage (Hz). Si vide, fournir 'Tau' en secondes.
% %
% % Options (Name-Value)
% %   'Tau'        : vecteur des lags en secondes aligné avec la partie lags>=0 de r.
% %   'SearchHz'   : [fmin fmax] en Hz pour contraindre la recherche du pic fondamental.
% %                  (Ex: [200 2000]). Si vide, cherche sur toute la plage.
% %   'Normalize'  : true/false, normaliser r par r(0). Défaut: true.
% %   'MaxHarm'    : entier >= 2, nombre de multiples à tester pour "harmonicité". Défaut: 3.
% %
% % SORTIE
% %   out : struct avec champs:
% %     peak         : hauteur du pic principal (hors lag 0), ACF normalisée si Normalize=true
% %     lag_pk_samp  : lag du pic (échantillons)
% %     period_s     : période estimée (s)
% %     f0_Hz        : fréquence fondamentale estimée (Hz) si fs fourni
% %     FWHM_s       : largeur à mi-hauteur (s) (interp linéaire)
% %     FWHM_samp    : largeur à mi-hauteur (échantillons)
% %     PSR_mean     : Peak-to-Sidelobe Ratio (pic / mean lobe secondaire)
% %     PSR_max      : pic / max lobe secondaire
% %     sharpness    : pic / FWHM_s (grande valeur => pic étroit et haut)
% %     prominence   : pic - baseline locale (hors pic)
% %     jitter_CV    : coefficient de variation des espacements entre pics ACF successifs
% %     harmonicity  : moyenne des amplitudes aux multiples (2..MaxHarm) / pic
% %     n_peaks      : nb de pics détectés (sur 0..5*period)
% %     ACFScore     : score agrégé (0..1) basé sur peak, PSR, FWHM
% %     details      : sous-struct avec éléments intermédiaires
% %
% % Remarque: nécessite 'findpeaks' (Signal Processing Toolbox) pour certaines métriques.
% 
% % ------------ Parse options ------------
% p = inputParser;
% p.addParameter('Tau', [], @(v)isvector(v));
% p.addParameter('SearchHz', [], @(v)isvector(v) && numel(v)==2);
% p.addParameter('Normalize', true, @(v)islogical(v) || ismember(v,[0 1]));
% p.addParameter('MaxHarm', 3, @(v)isnumeric(v) && v>=2);
% p.parse(varargin{:});
% tau = p.Results.Tau;
% searchHz = p.Results.SearchHz;
% doNorm = p.Results.Normalize;
% Kmax = p.Results.MaxHarm;
% 
% r = r(:);
% 
% % ------------ Isoler lags >= 0 et normaliser ------------
% % Si r est symétrique, détecter l'indice du lag 0 comme argmax (souvent au centre).
% [~, i0] = max(r);
% if i0 > 1 && i0 < numel(r)
%     rpos = r(i0:end);
% else
%     rpos = r; % suppose déjà lags >= 0
%     i0 = 1;
% end
% 
% % Normalisation par r(0) si demandé
% if doNorm
%     r0 = rpos(1);
%     if r0 ~= 0
%         rpos = rpos / r0;
%     end
% end
% 
% Npos = numel(rpos);
% 
% % Définir l'axe des lags (en échantillons) et en secondes si possible
% lags_samp = (0:Npos-1).';
% if ~isempty(tau)
%     % tau = temps (s) pour lags >= 0
%     if numel(tau) ~= Npos
%         error('La longueur de ''Tau'' doit correspondre à la partie lags>=0 de r.');
%     end
%     tlag = tau(:);
%     dt = mean(diff(tlag));
% elseif ~isempty(fs)
%     dt = 1/fs;
%     tlag = lags_samp * dt;
% else
%     dt = NaN;
%     tlag = NaN(size(lags_samp));
% end
% 
% % ------------ Définir fenêtre de recherche du pic ------------
% % Exclure lag 0
% imin = 2; imax = Npos;
% if ~isempty(searchHz) && ~isempty(fs)
%     fmin = max(1e-6, searchHz(1)); fmax = min(fs/2*0.999, searchHz(2));
%     % Convertir en lags: période = 1/f -> lag = période/dt = fs/f
%     qmin = max(2, floor(fs / fmax));
%     qmax = min(Npos, ceil(fs / fmin));
%     if qmax > qmin
%         imin = qmin; imax = qmax;
%     end
% end
% 
% if imax <= imin
%     % Pas de fenêtre valide -> retourner neutre
%     out = default_out();
%     return;
% end
% 
% seg = rpos(imin:imax);
% [peak, krel] = max(seg);
% lag_pk = (imin - 1) + (krel - 1);   % en échantillons
% period_s = lag_pk * dt;
% if ~isempty(fs)
%     f0 = 1/max(period_s, eps);
% else
%     f0 = NaN;
% end
% 
% % ------------ FWHM (largeur à mi-hauteur) ------------
% half = peak/2;
% % Chercher passages demi-hauteur autour du pic
% iL = lag_pk;
% while iL > 2 && rpos(iL) > half
%     iL = iL - 1;
% end
% iR = lag_pk;
% while iR < Npos && rpos(iR) > half
%     iR = iR + 1;
% end
% % Interpolation linéaire pour précision sub-échantillon
% if iL < lag_pk && iR > lag_pk
%     % côté gauche
%     tL = interp1([rpos(iL) rpos(iL+1)], [iL iL+1], half, 'linear', 'extrap');
%     % côté droit
%     tR = interp1([rpos(iR-1) rpos(iR)], [iR-1 iR], half, 'linear', 'extrap');
%     FWHM_samp = max(tR - tL, 1);
% else
%     FWHM_samp = 1;
% end
% FWHM_s = FWHM_samp * dt;
% 
% % ------------ PSR (peak-to-sidelobe ratio) ------------
% w = max(3, round(0.05*lag_pk)); % fenêtre exclue autour du pic
% idxSL = true(size(rpos));
% idxSL(1:max(1,lag_pk-w):min(Npos,lag_pk+w)) = false; % garder lobes hors ±w
% sidelobes = rpos(idxSL & ( (1:Npos)' >= imin & (1:Npos)' <= imax ));
% if isempty(sidelobes)
%     PSR_mean = 1; PSR_max = 1;
% else
%     PSR_mean = peak / max(mean(max(sidelobes, 1e-6)), 1e-6);
%     PSR_max  = peak / max(max(sidelobes), 1e-6);
% end
% 
% % ------------ Prominence locale (baseline) ------------
% ringL = rpos(max(1,lag_pk-2*w):max(1,lag_pk-w));
% ringR = rpos(min(Npos,lag_pk+w):min(Npos,lag_pk+2*w));
% ring  = [ringL; ringR];
% if isempty(ring)
%     baseline = 0;
% else
%     baseline = median(ring);
% end
% prominence = max(peak - baseline, 0);
% 
% % ------------ Pics successifs & jitter des périodes ------------
% % On recherche les pics sur [1 .. min(Npos, ceil(5*lag_pk))]
% imaxSearch = min(Npos, max(lag_pk*5, lag_pk + 5*w));
% region = rpos(2:imaxSearch); % ignorer lag 0
% % 'MinPeakDistance' ~ 60% du lag fondamental pour éviter les faux doublons
% minDist = max(2, round(0.6*lag_pk));
% try
%     [pks, locs] = findpeaks(region, 'MinPeakDistance', minDist); %#ok<*ASGLU>
%     locs = locs + 1; % décalage car region démarre à index 2
% catch
%     % fallback: détection manuelle très simple
%     locs = simple_findpeaks(region) + 1;
%     pks = rpos(locs);
% end
% 
% if numel(locs) >= 3
%     periods_samp = diff(locs);
%     jitter_CV = std(periods_samp) / max(mean(periods_samp), eps);
% else
%     jitter_CV = NaN;
% end
% n_peaks = numel(locs);
% 
% % ------------ Harmonicité (amplitude aux multiples du lag_pk) ------------
% harm_vals = [];
% for k = 2:Kmax
%     mu = k*lag_pk;
%     wh = max(2, round(0.08*lag_pk));
%     i1 = max(2, mu - wh); i2 = min(Npos, mu + wh);
%     if i2 > i1
%         [hk, ~] = max(rpos(i1:i2));
%         harm_vals(end+1) = hk; %#ok<AGROW>
%     end
% end
% if ~isempty(harm_vals)
%     harmonicity = mean(harm_vals) / max(peak, eps);
% else
%     harmonicity = NaN;
% end
% 
% % ------------ Score agrégé (0..1) basé sur ACF uniquement ------------
% % squashers vers [0,1]
% peak01 = min(max(peak,0),1);
% psr01  = 1 - exp(-PSR_mean/5);        % PSR ~ [5..20] → 0.6..0.98
% width01 = 1 - min(FWHM_s / max(0.2, 5*dt), 1); % plus étroit = mieux
% ACFScore = 0.60*peak01 + 0.25*psr01 + 0.15*width01;
% 
% % ------------ Pack sortie ------------
% out = struct();
% out.peak        = peak;
% out.lag_pk_samp = lag_pk;
% out.period_s    = period_s;
% out.f0_Hz       = f0;
% out.FWHM_s      = FWHM_s;
% out.FWHM_samp   = FWHM_samp;
% out.PSR_mean    = PSR_mean;
% out.PSR_max     = PSR_max;
% out.sharpness   = peak / max(FWHM_s, eps);
% out.prominence  = prominence;
% out.jitter_CV   = jitter_CV;
% out.harmonicity = harmonicity;
% out.n_peaks     = n_peaks;
% out.ACFScore    = ACFScore;
% 
% % Détails utiles pour debug/plots
% out.details = struct('rpos', rpos, 'tlag', tlag, 'lags_samp', lags_samp, ...
%                      'imin', imin, 'imax', imax, 'w_excl', w);
% 
% end

% ---------- Helpers ----------

function out = default_out()
out = struct('peak',0,'lag_pk_samp',NaN,'period_s',NaN,'f0_Hz',NaN, ...
             'FWHM_s',NaN,'FWHM_samp',NaN,'PSR_mean',NaN,'PSR_max',NaN, ...
             'sharpness',NaN,'prominence',NaN,'jitter_CV',NaN, ...
             'harmonicity',NaN,'n_peaks',0,'ACFScore',0,'details',struct());
end

function locs = simple_findpeaks(x)
% Détection de pics très simple si findpeaks indisponible
locs = find(x(2:end-1) > x(1:end-2) & x(2:end-1) >= x(3:end)) + 1;
end

function [pks, locs] = simple_findpeaks_robust(x, varargin)
% Détection de pics sans Signal Processing Toolbox, avec options:
%   'MinPeakDistance' (entier >=1, défaut 1)
%   'MinPeakHeight'   (défaut -Inf)
%   'MinProminence'   (défaut 0)
%
% Renvoie:
%   pks  : amplitudes des pics retenus
%   locs : indices (1-based) des pics retenus

p = inputParser;
p.addParameter('MinPeakDistance', 1, @(v)isnumeric(v) && v>=1);
p.addParameter('MinPeakHeight', -Inf, @(v)isnumeric(v) && isscalar(v));
p.addParameter('MinProminence', 0, @(v)isnumeric(v) && isscalar(v) && v>=0);
p.parse(varargin{:});
minDist = round(p.Results.MinPeakDistance);
minH   = p.Results.MinPeakHeight;
minProm= p.Results.MinProminence;

x = x(:);
N = numel(x);
pks = []; locs = [];
if N < 3, return; end

% 1) Candidats = maxima locaux
cand = find(x(2:end-1) > x(1:end-2) & x(2:end-1) >= x(3:end)) + 1;
if isempty(cand), return; end

% 2) Filtre par hauteur minimale
cand = cand(x(cand) >= minH);
if isempty(cand), return; end

% 3) Prominence (approx) : pic - max(min gauche, min droite)
if minProm > 0
    keep = true(size(cand));
    for ii = 1:numel(cand)
        c = cand(ii);

        % Chercher le "bassin" gauche
        i = c; leftMin = x(c);
        while i > 1 && x(i-1) <= x(i)
            i = i - 1;
            leftMin = min(leftMin, x(i));
            % Stop si rencontre un point plus haut que le pic courant
            if x(i) > x(c), break; end
        end

        % Chercher le "bassin" droite
        j = c; rightMin = x(c);
        while j < N && x(j+1) <= x(j)
            j = j + 1;
            rightMin = min(rightMin, x(j));
            if x(j) > x(c), break; end
        end

        prom = x(c) - max(leftMin, rightMin);
        if prom < minProm
            keep(ii) = false;
        end
    end
    cand = cand(keep);
    if isempty(cand), return; end
end

% 4) Distance minimale: garder les plus hauts, éliminer voisins trop proches
[~, order] = sort(x(cand), 'descend');
cand_sorted = cand(order);

sel = false(size(cand_sorted));
taken = [];
for k = 1:numel(cand_sorted)
    c = cand_sorted(k);
    if isempty(taken) || all(abs(c - taken) >= minDist)
        sel(k) = true;
        taken(end+1) = c; %#ok<AGROW>
    end
end
locs = sort(cand_sorted(sel));
pks  = x(locs);
end

function BEST = pick_best_carrier_fast(x, Fs, Opt)
% PICK_BEST_CARRIER_FAST  (une FFT, cepstre local)
% Trouve la porteuse la plus exploitable (résonance + AM sidebands).
%
% In:
%   x   : signal
%   Fs  : Hz
%   Opt : struct optionnel
%     .rangeHz   = [200 5000]   % bande de recherche des porteuses
%     .nfft      = []           % auto: 2^nextpow2(N)
%     .maxPeaks  = 15           % nbre max de pics à évaluer
%     .minPromDB = 6            % proéminence min (dB)
%     .K         = 3            % nb d’ordres de sidebands (±kΔ)
%     .DeltaHz   = [5 150]      % fenêtre de Δ à considérer (Hz)
%     .cepBand   = 4            % largeur (en multiples de Δ estimée grossière) pour le cepstre local
%     .wProm     = 0.35         % pondération proéminence
%     .wQ        = 0.25         % pondération Q
%     .wSB       = 0.40         % pondération sidebands
%
% Out:
%   BEST.fc, BEST.Delta, BEST.band = [fc-3Δ fc+3Δ], BEST.score
%   BEST.F, BEST.Pdb (pour plots), BEST.candidates (diag)

if nargin<3, Opt = struct; end
DEF = struct('rangeHz',[200 5000],'nfft',[], 'maxPeaks',15, ...
             'minPromDB',6,'K',3,'DeltaHz',[5 150], 'cepBand',4, ...
             'wProm',0.35,'wQ',0.25,'wSB',0.40);
fn = fieldnames(DEF); for k=1:numel(fn), if ~isfield(Opt,fn{k}), Opt.(fn{k})=DEF.(fn{k}); end, end

x = x(:) - mean(x(:));
N = numel(x);
if isempty(Opt.nfft), nfft = 2^nextpow2(N); else, nfft = Opt.nfft; end

% ===== 1) Periodogram (une FFT) =====
X   = fft(x, nfft);
F   = (0:nfft-1).' * (Fs/nfft);
P   = abs(X).^2;             % puissance
P   = P(1:floor(nfft/2)+1);
F   = F(1:numel(P));

% bande utile
m    = (F>=Opt.rangeHz(1) & F<=Opt.rangeHz(2));
Fsr  = F(m);
Psr  = P(m);
Pdb  = 10*log10(Psr + eps);

% ===== 2) Pic candidats via proéminence (dB) =====
[pks,locs,w,prom] = findpeaks(Pdb, Fsr, 'MinPeakProminence', Opt.minPromDB);
if isempty(locs)
    BEST = struct('fc',NaN,'Delta',NaN,'band',[NaN NaN],'score',0,'F',F,'Pdb',10*log10(P+eps),'candidates',[]);
    return;
end
% garder les plus proéminents
if numel(locs) > Opt.maxPeaks
    [~,ord] = maxk(prom, Opt.maxPeaks);
    locs = locs(ord); pks=pks(ord); w=w(ord); prom=prom(ord);
end

% ===== 3) Pour chaque candidat: Q + Δ (cepstre local) + score sidebands =====
K = Opt.K; 
DeltaMin = Opt.DeltaHz(1); DeltaMax = Opt.DeltaHz(2);
nC = numel(locs);
Qval = zeros(nC,1);
Delta = zeros(nC,1);
SB    = zeros(nC,1);

% Loganalyse: indice d’une fréquence
F2idx = @(f) max(1, min(numel(Fsr), round( 1 + (f - Fsr(1))*(numel(Fsr)-1)/(Fsr(end)-Fsr(1)) )));

for i=1:nC
    fc = locs(i);

    % --- Q local par -3 dB ---
    [~,i0] = min(abs(Fsr - fc));
    pk  = Pdb(i0);
    th  = pk - 3;
    iL=i0; while iL>1             && Pdb(iL) > th, iL=iL-1; end
    iR=i0; while iR<numel(Pdb)    && Pdb(iR) > th, iR=iR+1; end
    BW = max(1e-6, Fsr(iR) - Fsr(iL));
    Qval(i) = fc / BW;

    % --- Cepstre local (estime Δ) ---
    % on prend une fenêtre spectrale étroite autour de fc : [fc - W, fc + W]
    % W est initialisé grossièrement pour contenir quelques Δ possibles
    W = max(3*DeltaMax, 30);         % largeur en Hz (suffisante)
    f1 = max(Fsr(1), fc - W); f2 = min(Fsr(end), fc + W);
    i1 = F2idx(f1); i2 = F2idx(f2);
    Sl = Pdb(i1:i2);
    % enlève la pente/baseline (detrend simple)
    Sl = Sl - movmean(Sl, max(3, round(0.05*numel(Sl))));
    % cepstre (IFFT du log-magnitude — ici déjà en dB, approximation OK)
    C = real(ifft( max(Sl, -200) ));      % sécurité sur dB
    % l’index cepstral ~ période en "bins", qu’on convertit en Δ (Hz)
    % calage: quefrency q bins -> fréquence spacing ≈ Fs_band / (q * Nloc)
    % plus simple pratique: balaye un petit set Δ & prend le meilleur en corrélant ±Δ
    dGrid = linspace(DeltaMin, DeltaMax, 60);
    SBscore = -Inf(numel(dGrid),1);
    for j=1:numel(dGrid)
        d = dGrid(j);
        s = 0; used=0;
        for k2=1:K
            fL = fc - k2*d; fH = fc + k2*d;
            if fL<Fsr(1) || fH>Fsr(end), break; end
            vL = Pdb(F2idx(fL));
            vR = Pdb(F2idx(fH));
            s  = s + min(vL, vR); used=used+1;
        end
        if used>=2, SBscore(j)=s; end
    end
    [SB(i), jbest] = max(SBscore);
    Delta(i) = dGrid(jbest);
end

% normalisations robustes
PromN = 1 - 10.^(-prom/20);                % 0..1 depuis dB
Qn    = robust_unit(Qval);
SBn   = robust_unit(SB);

Score = Opt.wProm*PromN + Opt.wQ*Qn + Opt.wSB*SBn;
[~, ibest] = max(Score);

fc  = locs(ibest);
dlt = Delta(ibest);
band = [max(Fsr(1), fc-3*dlt), min(Fsr(end), fc+3*dlt)];

BEST = struct();
BEST.fc   = fc;
BEST.Delta= dlt;
BEST.band = band;
BEST.score= Score(ibest);
BEST.F    = F;
BEST.Pdb  = 10*log10(P+eps);
BEST.candidates = table(locs(:), prom(:), Qval(:), Delta(:), SB(:), Score(:), ...
    'VariableNames', {'fc_Hz','prom_dB','Q','Delta_Hz','SBscore','Score'});
end

function y = robust_unit(x)
x = x(:);
ix = isfinite(x);
if ~any(ix), y = zeros(size(x)); return; end
lo = prctile(x(ix), 5);
hi = prctile(x(ix),95);
if hi<=lo, y = 0*x; return; end
y = (x - lo) / (hi - lo);
y = max(0, min(1, y));
end

function plot_metrics_by_fc_and_hw(T, metrics)
    % Détecter les demi-largeurs disponibles (ou valeur unique NaN -> on ignore)
    hasHW = ismember('hw', T.Properties.VariableNames) && any(~isnan(T.hw));
    if hasHW
        UHW = unique(T.hw(~isnan(T.hw)));
    else
        UHW = NaN;
    end
    nM = numel(metrics);
    nRows = ceil(nM/2);
    figure('Name','Courbes des métriques vs fc','Color','w');
    tl = tiledlayout(nRows, 2, 'Padding','compact', 'TileSpacing','compact');

    cmap = lines(max(1, numel(UHW)));
    for m = 1:nM
        met = metrics{m};
        nexttile; hold on
        if hasHW
            for i = 1:numel(UHW)
                sel = T.hw == UHW(i);
                % Il peut y avoir plusieurs points par fc : on trace tel quel
                plot(T.fc(sel), T.(met)(sel), '-o', 'Color', cmap(i,:), ...
                     'DisplayName', sprintf('hw=%g', UHW(i)));
            end
            legend('Location','bestoutside'); 
        else
            plot(T.fc, T.(met), '-o', 'DisplayName', met);
        end
        grid on
        xlabel('fc (Hz)');
        ylabel(met);
        title(met, 'Interpreter','none');
    end
    title(tl, 'Métriques par centre de bande (fc), couleurs = demi-largeur (hw)');
end

function plot_metrics_heatmaps_fc_by_hw(T, metrics)
    if ~(ismember('fc',T.Properties.VariableNames) && ismember('hw',T.Properties.VariableNames))
        warning('Pas de fc/hw -> heatmap ignorée'); return;
    end
    nM = numel(metrics);
    nRows = ceil(nM/2);
    figure('Name','Heatmaps des métriques (fc × hw)','Color','w');
    tiledlayout(nRows, 2, 'Padding','compact', 'TileSpacing','compact');

    for m = 1:nM
        met = metrics{m};
        % Agréger (au cas où plusieurs lignes par (fc,hw)) :
        G = varfun(@mean, T, 'InputVariables', met, ...
                   'GroupingVariables', {'fc','hw'});
        % Transformer pour heatmap (fc en lignes, hw en colonnes)
        wide = unstack(G, ['mean_', met], 'hw');
        nexttile;
        % Heatmap attend des tableaux, on passe par table "wide"
        h = heatmap(wide, 'hw', 'fc', 'ColorVariable', ['mean_', met]);
        h.Colormap = turbo; h.ColorbarVisible = 'on';
        h.Title = met; h.XLabel = 'hw (Hz)'; h.YLabel = 'fc (Hz)';
    end
end

function R2 = cosine_fit_r2_from_acf(tlag, rpos, tau0, K)
% Fit r(t) ~= a cos(2π t/tau0) + b sin(2π t/tau0) et calcule R² sur K périodes.
if nargin < 4, K = 5; end
if isnan(tau0) || tau0 <= 0 || any(isnan(rpos)), R2 = 0; return; end
mask = (tlag >= 0) & (tlag <= K*tau0);
t = tlag(mask); y = rpos(mask);
if numel(t) < 10, R2 = 0; return; end
w = 2*pi/tau0;
X = [cos(w*t) sin(w*t)];
beta = X\y;
yhat = X*beta;
ss_res = sum((y - yhat).^2);
ss_tot = sum((y - mean(y)).^2) + eps;
R2 = max(0, 1 - ss_res/ss_tot);
end

function cpp = cpp_from_acf(rpos, fs, fmin, fmax)
% ACF -> PSD (>=0) -> log -> cepstre -> pic dans [1/fmax, 1/fmin]
if nargin < 3, fmin = 50; end
if nargin < 4, fmax = fs/2*0.9; end
N = 2^nextpow2(numel(rpos));
S = real(fft(rpos, N));      % ~ PSD (Wiener-Khinchin, à une constante près)
S = max(S, eps);             % éviter log(0)
logS = log(S);
cep = real(ifft(logS));
q = (0:N-1)'/fs;
mask = q >= 1/fmax & q <= 1/fmin;
if ~any(mask), cpp = 0; return; end
cseg = cep(mask);
t = linspace(0,1,numel(cseg))';
A = [t, ones(size(t))]; ab = A\cseg; baseline = A*ab;
cpp = max(cseg - baseline); cpp = max(0, cpp);
end

function per = periodic_energy_ratio(tlag, rpos, tau0, K, bw)
% Energie autour des multiples k*tau0 (fenêtre ±bw), rapportée à l'énergie totale.
if nargin < 4, K = 5; end
if nargin < 5, bw = 0.1; end % 10% de tau0
if isnan(tau0) || tau0 <= 0, per = 0; return; end
mask = (tlag >= 0) & (tlag <= K*tau0);
t = tlag(mask); y = rpos(mask);
E = sum(y.^2) + eps;
Eper = 0;
for k = 1:K
    tk = k*tau0;
    sel = (t >= tk*(1-bw)) & (t <= tk*(1+bw));
    if any(sel), Eper = Eper + sum(y(sel).^2); end
end
per = min(1, Eper / E);
end

function R2 = damped_cos_fit_r2(tlag, rpos, tau0, K)
if nargin < 4, K = 5; end
if isnan(tau0) || tau0 <= 0, R2 = 0; return; end
mask = (tlag >= 0) & (tlag <= K*tau0);
t = tlag(mask); y = rpos(mask);
if numel(t) < 10, R2 = 0; return; end
w0 = 2*pi/tau0;
% Modèle: A*exp(-alpha*t).*cos(w*t); w libre autour de w0
model = @(p,t) p(1)*exp(-abs(p(2))*t).*cos((w0+p(3))*t);
cost  = @(p) sum((y - model(p,t)).^2);
p0 = [max(y), 0, 0]; % [A, alpha, delta_w]
opts = optimset('Display','off');
p = fminsearch(cost, p0, opts);
yhat = model(p, t);
ss_res = sum((y - yhat).^2);
ss_tot = sum((y - mean(y)).^2) + eps;
R2 = max(0, 1 - ss_res/ss_tot);
end


function RES = evaluate_carrier_sweep(x, Fs, fc_list, Opt)
% Evalue une liste de porteuses fc : enveloppe -> ACF -> métriques & score.
% Affiche un plot du score total et des trois métriques.

% --------- defaults ---------
D = struct('envLP_Hz',200, 'acfMaxSec',1.0, 'scoreW',[0.4 0.4 0.2], ...
           'smoothMs',1.0, 'nfftACF',[]);
fn=fieldnames(D); for k=1:numel(fn), if ~isfield(Opt,fn{k}), Opt.(fn{k})=D.(fn{k}); end, end

x = x(:);
x = x - mean(x);

% sorties
Qacf  = zeros(numel(fc_list),1);
Qspec = zeros(numel(fc_list),1);
Qent  = zeros(numel(fc_list),1);
fractal   = zeros(numel(fc_list),1);
Score = zeros(numel(fc_list),1);

% pré-calcul
maxLag = min(round(Opt.acfMaxSec*Fs), numel(x)-2);
wSamp  = max(1, round(Opt.smoothMs*1e-3*Fs));   % lissage léger

for i = 1:numel(fc_list)
    fc = fc_list(i);

    % 1) Enveloppe par démodulation complexe et LPF
%      env = envelope_demod(x, Fs, fc, Opt.envLP_Hz);
    env = getEnv(x, Fs, fc);  %   extracts the envelope
    % 2) ACF normalisée (FFT) jusqu’à maxLag
    r = acf_fft_norm(env, maxLag);
%     fractal(i) = sum(abs(diff(r)));
%     dx = gradient(r);
% ddx = gradient(dx);
% fractal(i)=max(abs(ddx));
%  MCosine = characterizeSmoothCosine(r, Fs);

%     fractal(i) = absFeatures_FractalDatas(numel(r),r);
    % 3) Lissage léger (stabilise les pics)
    if wSamp>1, r = movavg(r, wSamp); end

    % 4) Métriques de qualité
    M = acf_quality_metrics(r);
    Qacf(i)  = M.Qacf;
    Qspec(i) = M.Qspec;
    Qent(i)  = M.Qent;
%     figure;
%     plot(r);
%     title(sprintf('%d => Qacf %d Qspec %d Qent %d THD:%d',fc_list(i), Qacf(i),Qspec(i),Qent(i),MCosine.THD));

    % 5) Score pondéré
    w = Opt.scoreW(:).';
    Score(i) = max(0, w(1)*Qacf(i) + w(2)*Qspec(i) + w(3)*Qent(i));
end

% normalisation optionnelle du score (0..1) pour lecture facile
ScoreN = (Score - min(Score)) / max(eps, (max(Score)-min(Score)));

% meilleur fc
[~, imax] = max(Score);
fc_best = fc_list(imax);

% --------- plots ---------
figure('Name','ACF quality sweep');
tiledlayout(2,1);

% (a) Score global
nexttile;
plot(fc_list, ScoreN, 'k-', 'LineWidth',1.8); hold on; grid on;
xline(fc_best, 'r--', sprintf('best f_c = %.1f Hz', fc_best));
ylabel('Score (normalisé)');
title('Score global par porteuse');

% (b) Détail métriques
nexttile; 
plot(fc_list, Qacf,  '-', 'LineWidth',1.5); hold on;
plot(fc_list, Qspec, '-', 'LineWidth',1.5);
plot(fc_list, Qent,  '-', 'LineWidth',1.5);
grid on; xlabel('f_c (Hz)');
ylabel('Métriques');
legend('Q_{acf}','Q_{spec}','Q_{ent}','Location','best');
title('Détail des métriques');

% --------- retour ---------
RES = struct();
RES.fc_list = fc_list(:);
RES.Qacf    = Qacf;
RES.Qspec   = Qspec;
RES.Qent    = Qent;
RES.Score   = Score;
RES.ScoreN  = ScoreN;
RES.fc_best = fc_best;

end

function env = envelope_demod(x, Fs, fc, envLP_Hz)
% Démodulation complexe autour de fc puis enveloppe lissée
t = (0:numel(x)-1).' / Fs;
z = x .* exp(-1j*2*pi*fc*t);               % translate fc -> 0
% LPF passe-bas (Butter) sur I/Q
if envLP_Hz <= 0, env = abs(z); return; end
Wc = min(0.99, envLP_Hz/(Fs/2));
[b,a] = butter(4, Wc, 'low');
zi = filtfilt(b,a, real(z));
zq = filtfilt(b,a, imag(z));
env = sqrt(zi.^2 + zq.^2);                 % magnitude = enveloppe
% centrer et limiter les outliers
env = env - median(env);
end

function r = acf_fft_norm(x, maxLag)
% ACF (rapide) normalisée : r(0)=1, lags >= 0
x = x(:) - mean(x);
L = numel(x);
nfft = 2^nextpow2(2*L-1);
X = fft(x, nfft);
S = X .* conj(X);
rfull = ifft(S, 'symmetric');
r = rfull(1:maxLag+1);
r = r / max(r(1), eps);
end

function M = acf_quality_metrics(r)
% Trois métriques : Qacf (peak-to-background), Qspec (concentration),
% Qent (entropie "inverse").
r = r(:);
N = numel(r);

% 1) Spectral concentration de l’ACF
Rspec = abs(fft(r)).^2;
Rspec = Rspec(1:floor(N/2));
if sum(Rspec) <= 0, Qspec = 0; else, Qspec = max(Rspec)/sum(Rspec); end

% 2) Entropie du spectre de l’ACF (faible si pic dominant)
ps = Rspec / max(sum(Rspec), eps);
ps = ps(ps>0);
if isempty(ps)
    Qent = 0;
else
    H = -sum(ps .* log(ps));
    Hmax = log(numel(ps));
    Qent = 1 - H/max(Hmax, eps);
end

% 3) Peak-to-background sur l’ACF (hors lag 0)
if N >= 6
    [~,L1] = max(r(2:end)); L1=L1+1;
    mainPk = max(0, r(L1));                   % 1er pic > lag0
    bg = mean(abs(r(round(N*0.6):end)));      % fond sur fin de fenêtre
    Qacf = mainPk / max(bg, 1e-6);
    % comprimer pour rester ~[0..1]
    Qacf = Qacf / (1 + Qacf);
else
    Qacf = 0;
end

M = struct('Qacf', Qacf, 'Qspec', Qspec, 'Qent', Qent);
end

function BEST = robust_best_carrier(x, Fs, fc_list, Opt)
% Etage 1: balayage métriques
R = evaluate_carrier_sweep(x, Fs, fc_list, Opt.sweep);  % ta fonction existante
ScoreN = R.ScoreN;

% Fenêtre adaptative: garder les fc avec ScoreN >= theta
theta = 0.80;                              % seuil de confiance
idxGood = find(ScoreN >= theta);
if isempty(idxGood)
    % fallback: prends top-K (ex. 3) et ouvre une petite fenêtre autour
    [~,ord] = maxk(ScoreN, min(3,numel(ScoreN)));
    fc_win = sort(R.fc_list(ord));
    f_low  = max( min(fc_win)-200, 500);
    f_high = min( max(fc_win)+200, 2000);
else
    f_low  = max( min(R.fc_list(idxGood)) - 200, 500);
    f_high = min( max(R.fc_list(idxGood)) + 200, 2000);
end

% Petit handicap optionnel au-dessus de 1.5 kHz (si utile chez toi)
if isfield(Opt,'preferLow') && Opt.preferLow
    Opt.pick.penalizeAboveHz = 1500;  % simple indicateur pour tie-break
else
    Opt.pick.penalizeAboveHz = inf;
end

% Etage 2: affinage
Opt.pick.rangeHz = [f_low f_high];
Best1 = pick_best_carrier_fast(x, Fs, Opt.pick);

% Tie-break si Best est au-dessus de penalizeAboveHz et un concurrent plus bas est très proche
% -> on rescanne la sous-bande [500, penalizeAboveHz]
if Best1.fc > Opt.pick.penalizeAboveHz
    Opt2 = Opt.pick; Opt2.rangeHz = [500 Opt.pick.penalizeAboveHz];
    BestLow = pick_best_carrier_fast(x, Fs, Opt2);
    if isfinite(BestLow.score) && (BestLow.score >= 0.97*Best1.score)
        BEST = BestLow;  % préférer la basse fréquence si scores ~égaux
        BEST.note = 'tie-break->lower-fc';
        return;
    end
end

BEST = Best1;
BEST.note = 'adaptive-window';
end

% function y = movavg(x,w)
% if w<=1, y=x; return; end
% k = ones(w,1)/w;
% try, y = filtfilt(k,1,double(x)); catch, y = conv(double(x),k,'same'); end
% end

function [OUT ] = period_from_acf_median_robust(acf, Fs, tauMinSamp, tauMaxSamp, Opt)
% PERIOD_FROM_ACF_MEDIAN_ROBUST
%  Estime le fondamental par "médiane des intervalles" + durcissements :
%   - lissage léger
%   - détection de pics avec proéminence relative dans [tauMin..tauMax]
%   - rejet d'outliers (MAD) sur les intervalles
%   - correction de multiples (ramène 2T,3T->T via score peigne)
%   - interpolation parabolique locale
%   - anti-demi-tour (R(T/2) vs R(T) + odd/even)
%
% In:
%   acf : ACF (lags>=0), pas nécessairement normalisée
%   Fs  : Hz
%   tauMinSamp, tauMaxSamp : bornes de recherche en échantillons
%   Opt (struct) champs optionnels:
%     .SmoothMs    (default 1.0)
%     .PromRel     (default 0.05)  % proéminence relative vs max hors lag0
%     .MinPkDist   (default round(0.33*Fs/(Fs/tauMaxSamp))) % ~T/3 @ fmax
%     .UseMAD      (default true)
%     .MAD_k       (default 3.5)   % seuil outliers
%     .CheckMults  (default true)  % essayer T/2, T, 2T, 3T...
%     .MaxMult     (default 3)
%     .Parabolic   (default true)
%     .AntiHalf    (default true)
%     .Kharm       (default 5)
%     .lambdaEven  (default 0.8)
%     .gammaHalf   (default 1.12)
%
% Out (struct OUT):
%   .f_hat, .T_hat, .L_hat           : freq, période(s), lag (samples)
%   .L_med, .L_candidates, .L_choice : médiane, candidats multiples, choisi
%   .pkLoc, .pkVal                   : pics retenus (indices & amplitudes)
%   .antiHalfTriggered (bool), .combContrast
%   .status ('ok'|'no_peaks'|'bad_range')

% -------- Defaults --------
D = struct('SmoothMs',1.0,'PromRel',0.05,'MinPkDist',[], ...
           'UseMAD',true,'MAD_k',3.5,'CheckMults',true,'MaxMult',1, ...
           'Parabolic',true,'AntiHalf',true,'Kharm',5,'lambdaEven',0.8,'gammaHalf',1.12);
if nargin<5, Opt = struct; end
fn = fieldnames(D); for k=1:numel(fn), if ~isfield(Opt,fn{k}), Opt.(fn{k})=D.(fn{k}); end, end

r = acf(:);

Nlag = numel(r);
tauMinSamp = max(2, floor(tauMinSamp));
tauMaxSamp = min(max(tauMinSamp+2, floor(tauMaxSamp)), Nlag-1);

if tauMinSamp+2 >= tauMaxSamp
    OUT = struct('status','bad_range','f_hat',0,'T_hat',NaN,'L_hat',NaN);
    return;
end

% -------- Lissage léger --------
w = max(1, round(Opt.SmoothMs*1e-3*Fs));
if w>1
    r_s = movavg(r, w);
else
    r_s = r;
end

% -------- Détection de pics dans [tauMin..tauMax] --------
base = max(r_s(2:tauMaxSamp+1));
prom = Opt.PromRel * max(base, eps);
if isempty(Opt.MinPkDist)
    % approx: T @ fmax ~ Fs/tauMin -> distance min ~ T/3
    Tmin = Fs / tauMaxSamp;                     % ~ 1/fmax
    minPkDist = max(2, round(0.33 * (Fs / (Fs/tauMaxSamp)) )); %#ok<NASGU>
    % plus robuste : en samples, ~ 0.33*tauMin
    minPkDist = max(2, round(0.33 * tauMinSamp));
else
    minPkDist = Opt.MinPkDist;
end

[pkVal0, pkLoc0] = findpeaks(r_s, 'MinPeakProminence', prom, 'MinPeakDistance', minPkDist);
% garde dans bande & enlève lag0
mask = (pkLoc0>= (tauMinSamp+1)) & (pkLoc0 <= (tauMaxSamp+1));
pkLoc = pkLoc0(mask);
pkVal = pkVal0(mask);

if isempty(pkLoc)
    OUT = struct('status','no_peaks','f_hat',0,'T_hat',NaN,'L_hat',NaN, ...
                 'pkLoc',[],'pkVal',[]); 
    return;
end

% -------- Intervalles & rejet d’outliers (MAD) --------
d = diff(pkLoc);                 % en samples
d = double(d(:));
if numel(d) == 0
    L_med = double(pkLoc(1)-1);  % seul pic -> fallback
else
    if Opt.UseMAD && numel(d)>=4
        medD = median(d);
        madD = median(abs(d - medD)) + eps;
        keep = abs(d - medD) <= Opt.MAD_k * madD;
        d = d(keep);
    end
    if isempty(d)
        d = diff(pkLoc); d = double(d(:));
    end
    L_med = median(d);
end

% -------- Correction de multiples (ramener à T) --------
L_candidates = L_med ./ (1:Opt.MaxMult);   % [L, L/2, L/3, ...] -> candidats T
L_scores     = -Inf(size(L_candidates));
for i=1:numel(L_candidates)
    Lc = L_candidates(i);
    tol = max(1.0, 0.08 * Lc);             % tolérance ±8% (amplitude à ajuster)
    % score peigne: somme des r_s(k*Lc) (k=1..K) bornés au domaine
    S = 0; used=0;
    for k2=1:ceil(tauMaxSamp / Lc)
        q = k2*Lc;
        if q>tauMaxSamp, break; end
        % interp linéaire sur r_s (indices fractionnaires)
        val = safe_interp((0:Nlag-1).', r_s, q);
        S = S + max(val,0);
        used = used + 1;
    end
    if used>=2, L_scores(i) = S/used; end
end
[~, ibest] = max(L_scores);
L0 = L_candidates( max(1, ibest) );

% --- Affinage autour d’un pic près du FONDAMENTAL (L0) ---
lags_pk = pkLoc - 1;             % passer en lags 0-based
tol = max(2, round(0.20*L0));    % tolérance ±20% (ajuste selon tes données)

% 1) Cherche un pic autour de L0
[delta1, i1] = min(abs(lags_pk - L0));

if delta1 <= tol
    L_int = lags_pk(i1);         % pic trouvé ~L0 (fondamental)
    L_used = L0;
else
    % 2) Sinon, essaie 2*L0 (harmonique), puis 0.5*L0 (sous-harmonique)
    [delta2, i2] = min(abs(lags_pk - 2*L0));
    [deltah, ih] = min(abs(lags_pk - 0.5*L0));
    [dmin, which] = min([delta2 deltah]);

    if dmin <= tol
        if which==1
            % accroché à 2T -> revient au fondamental
            L_int = lags_pk(i2);
            L_used = 2*L0;
            % après interpolation, on divisera par 2
        else
            % accroché à T/2 -> remonte au fondamental
            L_int = lags_pk(ih);
            L_used = 0.5*L0;
            % après interpolation, on multipliera par 2
        end
    else
        % 3) fallback: prends le meilleur autour de L0, même si hors tol
        L_int = lags_pk(i1);
        L_used = L0;
    end
end

% --- Interpolation parabolique autour de L_int ---
if L_int>=2 && L_int<=tauMaxSamp-1
    y1 = r_s(L_int-1+1); y2 = r_s(L_int+1); y3 = r_s(L_int+1+1); % +1 car r_s(1)=lag0
    denom = (y1 - 2*y2 + y3);
    if denom < 0
        delta = 0.5*(y1 - y3)/denom;
        delta = max(min(delta, 0.5), -0.5);
    else
        delta = 0;
    end
else
    delta = 0;
end

L_ref = L_int + delta;  % lag (samples) près de L_used (peut être 0.5, 1, 2 * L0)

% % --- Si on a locké 2T ou T/2, re-projette vers le fondamental ---
% if abs(L_used - 2*L0) < 1e-9
%     L_ref = L_ref_local / 2;   % revenir au T fondamental
% elseif abs(L_used - 0.5*L0) < 1e-9
%     L_ref = L_ref_local * 2;   % revenir au T fondamental
% else
%     L_ref = L_ref_local;       % déjà sur T
% end

T_hat = L_ref / Fs;

% % -------- Anti-demi-tour (ACF) --------
% antiHalfTriggered = false; combContrast = NaN;
% if Opt.AntiHalf
%     interpR = @(q) safe_interp((0:Nlag-1).', r_s, q);
%     T     = T_hat * Fs;
%     R_T   = interpR(T);
%     R_T2  = interpR(T/2);
%     % peigne odd/even à T
%     S_odd=0; S_even=0;
%     for k2=1:Opt.Kharm
%         q = k2*T;
%         if q>tauMaxSamp, break; end
%         v = interpR(q);
%         if mod(k2,2)==1, S_odd=S_odd+max(v,0); else, S_even=S_even+max(v,0); end
%     end
%     combContrast = (S_odd - Opt.lambdaEven*S_even)/max(S_odd + Opt.lambdaEven*S_even, eps);
%     isHalf = (R_T2 > Opt.gammaHalf*R_T) && (S_even > S_odd);
%     if isHalf
%         T_hat = T_hat/2;
%         L_ref = L_ref/2;
%         antiHalfTriggered = true;
%     end
% end

OUT = struct();
OUT.status = 'ok';
OUT.L_med  = L_med;
OUT.L_candidates = L_candidates;
OUT.L_scores     = L_scores;
OUT.L_choice     = L0;
OUT.L_hat  = L_ref;
OUT.T_hat  = T_hat;
OUT.f_hat  = 1/max(T_hat,eps);
OUT.pkLoc  = pkLoc-1;  % en samples (0-based lag)
OUT.pkVal  = pkVal;
% OUT.antiHalfTriggered = antiHalfTriggered;
% OUT.combContrast = combContrast;
end


function v = safe_interp(x, y, q)
q = q(:); v = nan(size(q)); N = numel(x);
for i=1:numel(q)
    if q(i) < 1 || q(i) > N-2
        v(i)=NaN;
    else
        i0=floor(q(i)); a=q(i)-i0; v(i)=(1-a)*y(i0+1)+a*y(i0+2);
    end
end
if isscalar(v), v=v(1); end
end

function fractal = absFeatures_FractalDatas(SIZESIGNAL,vibrationSamples)
POW16=65536;
Area = 0;
dilations = 0;
erosions = 0;
dim = 16;
isM = 1;

for  i = 0:SIZESIGNAL - dim -isM
    
    dilations = vibrationSamples(i+isM); % Max
    erosions = vibrationSamples(i+isM); % Min
    for j = 0: dim-isM
        if (dilations < vibrationSamples(i + j + 1 + isM))
            dilations = vibrationSamples(i + j + 1 + isM);
        end
        if (erosions > vibrationSamples(i + j + 1 + isM))
            erosions = vibrationSamples(i + j + 1 + isM);
        end
    end
    Area = Area + (dilations - erosions);
end
fractal= fix(Area * (6 * POW16) / (Area + 4 * fix(absFeatures_Msqrt(double(Area) * 2))));
end

function out = absFeatures_Msqrt( in)
temp = in;
isNeg = 0;
if (temp < 0)
    isNeg = 1;
    temp = temp * (-1);
end

out = sqrt(temp);
if (isNeg == 1)
    out = out * (-1);
end
end


function M = characterizeSmoothCosine(x, Fs)
% CHARACTERIZESMOOTHCOSINE  Smoothness & cosine-likeness metrics for a 1D signal
% x : signal (vector)
% Fs: sampling rate (Hz)
%
% Returns struct M avec :
%   .TV                  variation totale normalisée
%   .roughnessL2         énergie du gradient (diff^2)
%   .bendingEnergy       énergie de la courbure (diff2^2)
%   .spectralCentroid    centroïde spectral (Hz)
%   .spectralBandwidth   écart-type spectral (Hz)
%   .f0_est              fondamentale estimée (Hz)
%   .lowFreqFrac         fraction d’énergie sous 2.5*f0
%   .cosineProjFrac      fraction d’énergie projetée sur {cos,sin}@f0 (0..1)
%   .R2_cosfit           R² du fit [cos,sin,DC] à f0
%   .THD                 Total Harmonic Distortion (≈ énergie harmoniques / fondamentale)

    % --- prétraitements ---
    x = x(:);
    x = detrend(x, 'linear');           % enlève la tendance
    N = numel(x);
    t = (0:N-1)'/Fs;
    Ex = sum(x.^2) + eps;

    % --- métriques “lissité” dans le domaine temporel ---
    dx  = diff(x);
    ddx = diff(x,2);
    M.TV            = sum(abs(dx)) / max(1,N-1);           % variation totale (plus petit = plus lisse)
    M.roughnessL2   = mean(dx.^2);                         % énergie du gradient
    M.bendingEnergy = mean(ddx.^2);                        % énergie de la courbure (punit les changements de pente)

    % --- spectre (énergie haute fréquence faible = signal lisse) ---
    [Pxx, f] = periodogram(x, [], max(1024,2^nextpow2(N)), Fs, 'onesided');
    % centroïde & largeur spectrale
    S = sum(Pxx) + eps;
    M.spectralCentroid  = sum(f .* Pxx) / S;
    M.spectralBandwidth = sqrt( sum(((f - M.spectralCentroid).^2).*Pxx) / S );

    % --- estimation f0 via pic principal (hors DC) ---
    if numel(Pxx) >= 3
        Pxx_noDC = Pxx;    Pxx_noDC(1) = 0;
        [~, imax] = max(Pxx_noDC);
        M.f0_est = f(imax);
    else
        M.f0_est = 0;
    end

    % fraction d’énergie sous 2.5*f0 (mesure de “concentration basse fréquence”)
    if M.f0_est > 0
        bandMax = min(2.5*M.f0_est, f(end));
    else
        bandMax = min(Fs*0.1, f(end)); % fallback : 10% de Fs si pas de f0
    end
    M.lowFreqFrac = sum(Pxx(f <= bandMax)) / S;

    % --- ressemblance à un cosinus : projection & R² ---
    if M.f0_est > 0
        C = cos(2*pi*M.f0_est*t);
        Sg = sin(2*pi*M.f0_est*t);
        Phi = [C Sg ones(N,1)];
        theta = Phi \ x;                 % LS fit
        xhat  = Phi * theta;             % reconstruction cos+sin+DC
        proj  = [C Sg] * theta(1:2);     % projection sur sous-espace {cos,sin} uniquement

        SSE = sum((x - xhat).^2);
        SST = sum( (x - mean(x)).^2 ) + eps;
        M.R2_cosfit      = max(0, 1 - SSE/SST);            % qualité du fit global (cos+sin+DC)
        M.cosineProjFrac = sum(proj.^2) / Ex;              % part d’énergie “expliquée” par cosinus pur (0..1)
    else
        M.R2_cosfit = 0;  M.cosineProjFrac = 0;
    end

    % --- THD approximatif (mesure d’harmonicité : plus petit ≈ plus “cosinus pur”) ---
    if M.f0_est > 0
        % énergie à la fondamentale (± un bin)
        E1 = localBinEnergy(f, Pxx, M.f0_est);
        % énergie des k = 2..5 harmoniques
        Ek = 0;
        for k = 2:5
            Ek = Ek + localBinEnergy(f, Pxx, k*M.f0_est);
        end
        M.THD = sqrt(Ek) / (sqrt(E1) + eps);
    else
        M.THD = NaN;
    end
end

function E = localBinEnergy(f, Pxx, ftarget)
    if ftarget <= 0 || ftarget > f(end)
        E = 0; return;
    end
    [~, i0] = min(abs(f - ftarget));
    i = unique(max(1, i0-1) : min(numel(f), i0+1)); % petite bande ±1 bin
    E = sum(Pxx(i));
end

function state = oma_init(nChannels, fs, varargin)
% OMA init – PCA avec facteur d'oubli (EW covariance), suivi modal par FFT
% Args:
%   nChannels : nb de capteurs (colonnes)
%   fs        : fréquence d'échantillonnage (Hz)
% Options (Name,Value):
%   'ForgettingFactor' (mu) : 0.9990 par défaut (0<mu<1, proche de 1)
%   'NumModes'        (M)   : 3    (nb de modes à suivre)
%   'FftLen'          (L)   : 4096 (longueur FFT pour estimer f_n)
%   'Decim'           (D)   : 1    (échantillons traités sur D; accélère)
%
% state contient:
%   mu, fs, M, L, D, t, m (moyenne), C (cov), W,Dl (VP), 
%   qBuf (L×M) buffers coord. modales, qIdx, win, 
%   f_est (M×1), shapes (nChannels×M), varExpl (M×1)

p = inputParser;
addParameter(p,'ForgettingFactor',0.9990);
addParameter(p,'NumModes',3);
addParameter(p,'FftLen',4096);
addParameter(p,'Decim',1);
parse(p,varargin{:});
mu  = p.Results.ForgettingFactor;
M   = p.Results.NumModes;
L   = p.Results.FftLen;
D   = p.Results.Decim;

state.mu = mu;
state.fs = fs;
state.M  = M;
state.L  = L;
state.D  = D;

state.t  = 0;                        % compteur d'échantillons traités
state.m  = zeros(nChannels,1);       % moyenne EW
state.C  = eye(nChannels);           % covariance EW (initialisation)
state.W  = eye(nChannels);           % vecteurs propres (colonnes)
state.Dl = ones(nChannels,1);        % valeurs propres

state.qBuf = zeros(L,M);             % buffers coordonnées modales
state.qIdx = 0;                      % index circulaire
state.win  = hann(L);                % fenêtre pour FFT
state.f_est   = nan(M,1);            % fréquences propres estimées
state.shapes  = state.W(:,1:M);      % formes modales courantes
state.varExpl = nan(M,1);            % variance expliquée (%) par mode
end

function [state,out] = oma_step(state, x)
% OMA step – met à jour moyenne/cov EW, recalcule PCA, pousse la coord. modale,
% estime les fréquences des M premiers modes via FFT de q_k (fenêtrée).
% Args:
%   x : vecteur colonne (nChannels×1) à l'instant courant
%
% out:
%   out.t        : échantillon absolu (entier)
%   out.freqs    : (M×1) dernières f_n estimées (Hz) [NaN si pas assez d'hist.]
%   out.shapes   : (nChannels×M) formes modales (colonnes normalisées)
%   out.varExpl  : (M×1) % variance expliquée cumulée par mode
%   out.eigs     : (nChannels×1) valeurs propres triées (énergie modale)

% décimation (facultatif)
if state.D > 1 && mod(state.t, state.D) ~= 0
    state.t = state.t + 1;
    out = make_out(state, false);
    return
end

mu = state.mu;
x  = x(:);

% --- MAJ moyenne EW ---
% m_{t+1} = mu*m_t + (1-mu)*x_t
m_old  = state.m;
state.m = mu*state.m + (1-mu)*x;

% --- MAJ covariance EW (autour de la moyenne courante) ---
% C_{t+1} = mu*C_t + (1-mu)*(x - m)(x - m)^T
xc = x - state.m;
state.C = mu*state.C + (1-mu)*(xc*xc.');

% --- Eigen-decomp petite matrice (nCh×nCh) ---
[W, Dl] = eig((state.C + state.C.')/2);      % symétrise par sûreté
[lam, idx] = sort(real(diag(Dl)), 'descend');
W = W(:, idx);

% Signage cohérent (optionnel): fixe la direction des vecteurs
for k = 1:size(W,2)
    [~,imax] = max(abs(W(:,k)));
    if W(imax,k) < 0, W(:,k) = -W(:,k); end
end

state.W  = W;
state.Dl = lam;

% --- Coordonnées modales instantanées (scores) ---
% q = W^T (x - m)
q = W.' * (x - state.m);

% --- pousse q_k dans buffers FFT pour les M premiers modes ---
state.qIdx = mod(state.qIdx, state.L) + 1;
for k = 1:state.M
    state.qBuf(state.qIdx, k) = q(k);
end

% --- estime f_n si buffer plein ---
if state.t >= state.L
    nfft = state.L;
    for k = 1:state.M
        y = state.qBuf(:,k) .* state.win;
        Y = fft(y, nfft);
        % On ne regarde que [0, fs/2)
        half = 1:floor(nfft/2);
        [~, imax] = max(abs(Y(half)));
        state.f_est(k) = (imax-1) * state.fs / nfft;
    end
end

% --- variance expliquée (%)
varTot = sum(state.Dl);
vexp   = state.Dl(1:state.M) / max(varTot, eps) * 100;

state.shapes  = state.W(:,1:state.M);
state.varExpl = vexp;

state.t = state.t + 1;
out = make_out(state, true);
end

function out = make_out(state, haveNew)
out.t       = state.t;
out.freqs   = state.f_est;
out.shapes  = state.shapes;
out.varExpl = state.varExpl;
out.eigs    = state.Dl;
out.updated = haveNew;
end

function res = oma_finalize(state)
% Récupère les derniers résultats consolidés
res = struct();
res.freqs   = state.f_est;
res.shapes  = state.shapes;
res.varExpl = state.varExpl;
res.eigs    = state.Dl;
end

function [st, out] = sc_step(st, x)
% x : scalaire (un échantillon)
% out.f_est (Hz) : estimation courante (NaN tant que pas assez d'hist)
% out.updated : true si nouvelle estimation
out = struct('f_est', st.f_est, 'updated', false);

% décimation (optionnel)
st.nproc = st.nproc + 1;
if st.D > 1 && mod(st.nproc-1, st.D) ~= 0
    return
end

% push dans buffer circulaire
st.idx = st.idx + 1;
if st.idx > st.N, st.idx = 1; end
st.buf(st.idx) = x;

% pas encore assez d'hist ?
persistent filled;
if isempty(filled), filled = false; end
if ~filled && st.nproc < st.N*st.D
    return
else
    filled = true;
end

% ré-indexer en segment temporel (ordre croissant)
if st.idx == st.N
    seg = st.buf;
else
    seg = [st.buf(st.idx+1:end); st.buf(1:st.idx)];
end

% fenêtre + FFT
segw = seg .* st.win;
Y = fft(segw);
H = Y(1:floor(st.N/2));              % 0..fs/2
P = (abs(H).^2) / sum(st.win.^2);    % périodogramme normalisé

% lissage exponentiel du spectre
mu = st.mu;
st.Pavg = mu*st.Pavg + (1-mu)*P;

% ignore DC (bin 1) pour la recherche de pic
[~, kmax] = max(st.Pavg(2:end));
kmax = kmax + 1;

% interpolation parabolique autour du pic (±1 bin) si possible
if kmax>2 && kmax<length(st.Pavg)
    a = st.Pavg(kmax-1); b = st.Pavg(kmax); c = st.Pavg(kmax+1);
    denom = (a - 2*b + c);
    if denom ~= 0
        delta = 0.5*(a - c) / denom;   % décalage sous-bin [-0.5,0.5]
    else
        delta = 0;
    end
else
    delta = 0;
end

kref = (kmax-1) + delta;                  % index réel (0-based)
st.f_est = kref * (st.fs / st.N);         % fréquence en Hz

out.f_est = st.f_est;
out.updated = true;
end

function st = sc_init(fs, varargin)
% Single-channel dominant frequency tracker (exponential-averaged periodogram)
% Options:
%   'FftLen' (default 4096)     – taille FFT
%   'Forgetting' (default 0.9)  – lissage expo. du spectre (0<mu<1)
%   'Decim' (default 1)         – traiter 1 échantillon sur D
p = inputParser;
addParameter(p,'FftLen',4096);
addParameter(p,'Forgetting',0.9);
addParameter(p,'Decim',1);
parse(p,varargin{:});
st.fs   = fs;
st.N    = p.Results.FftLen;
st.mu   = p.Results.Forgetting;
st.D    = p.Results.Decim;

st.buf  = zeros(st.N,1);
st.idx  = 0;
st.win  = hann(st.N);
st.Pavg = zeros(st.N/2,1);     % périodogramme lissé (0..fs/2)
st.f_est = NaN;
st.nproc = 0;
end

function st = sc_acf_init(fs, varargin)
% Autocorrelation-based tracker (fenêtre courte + peak picking)
p = inputParser;
addParameter(p,'WinLen',4096);   % taille fenêtre d'analyse
addParameter(p,'Forgetting',0.9);
addParameter(p,'Decim',1);
parse(p,varargin{:});

st.fs = fs; st.N = p.Results.WinLen; st.mu = p.Results.Forgetting; st.D = p.Results.Decim;
st.buf = zeros(st.N,1); st.idx=0; st.nproc=0;
st.Ravg = zeros(st.N,1);
st.f_est = NaN;
end

function [st,out] = sc_acf_step(st, x)
out = struct('f_est', st.f_est, 'updated', false);
st.nproc = st.nproc + 1;
if st.D > 1 && mod(st.nproc-1, st.D) ~= 0, return; end

st.idx = st.idx + 1; if st.idx>st.N, st.idx=1; end
st.buf(st.idx) = x;

persistent filled; if isempty(filled), filled=false; end
if ~filled && st.nproc < st.N*st.D, return; else, filled=true; end

% segment
if st.idx==st.N, seg=st.buf; else, seg=[st.buf(st.idx+1:end); st.buf(1:st.idx)]; end
seg = seg - mean(seg);

% autocorr via FFT (rapide)
M = 2^nextpow2(2*st.N-1);
S = fft(seg, M);
R = ifft(abs(S).^2);
R = real(R(1:st.N));
R = R / max(R+eps);       % normalisation

% lissage expo. de l'autocorr
st.Ravg = st.mu*st.Ravg + (1-st.mu)*R;

% chercher premier pic significatif après un lag min (élimine DC)
lagMin = round(st.fs/1000);  % ignore <1 ms (à adapter)
[~, locs] = findpeaks(st.Ravg,'MinPeakDistance',lagMin,'MinPeakProminence',0.05);
if ~isempty(locs)
    % premier pic non nul
    L = locs(1);
    st.f_est = st.fs / L;
    out.f_est = st.f_est;
    out.updated = true;
end
end

function BEST = pick_best_carrier_comb(x, Fs, Opt)
% Carrier picker robuste avec estimation de Δ par "comb-correlation".
% Sort:
%   BEST.best : struct du gagnant (fc, Delta, band, score, Q, etc.)
%   BEST.candidates : table triée de tous les candidats
%   BEST.psd : struct(F, Pxx_dB) pour debug/plots
%
% Requiert: Signal Processing Toolbox (pwelch)

% ---- défauts ----
DEF = struct( ...
  'rangeHz',[200 5000], ...
  'nfft',[], ...
  'welchWin',4096, ...      % fenêtre Welch
  'welchOL',0.5, ...        % overlap 0..1
  'maxPeaks',20, ...
  'minPromDB',6, ...
  'minDistHz',10, ...
  'K',3, ...                % ordres de sidebands ±kΔ
  'DeltaHz',[5 150], ...    % grille Δ (Hz)
  'etaBox',0.15, ...        % largeur boîte = max(2*df, eta*Δ)
  'thrSNR_dB',3, ...        % seuil SNR (dB) pour "couverture"
  'useMedianBaseline',true, ...
  'wProm',0.30,'wQ',0.20,'wSB',0.35,'wClean',0.15, ... % pondérations
  'Delta0_small',10 ...     % Δ0 pour pénalité douce des petits Δ
);
if nargin<3, Opt = struct; end
fn = fieldnames(DEF); for k=1:numel(fn), if ~isfield(Opt,fn{k}), Opt.(fn{k})=DEF.(fn{k}); end, end

% ---- pré-traitement ----
x = x(:); x = x - mean(x,'omitnan'); x(~isfinite(x)) = 0;
N = numel(x);
if isempty(Opt.nfft), nfft = 2^nextpow2(max(N, Opt.welchWin)); else, nfft = Opt.nfft; end

% ---- PSD Welch ----
welchWin = min(Opt.welchWin, N);
welchWin = 2^floor(log2(max(64, welchWin)));
noverlap = max(0, min(welchWin-1, round(Opt.welchOL * welchWin)));
win = hann(welchWin,'periodic');
[Pxx, F] = pwelch(x, win, noverlap, nfft, Fs, 'psd');   % V^2/Hz
Pdb = 10*log10(Pxx + eps);

% ---- bande utile ----
m = (F>=Opt.rangeHz(1) & F<=Opt.rangeHz(2));
Fsr = F(m); Psr = Pxx(m); Pdbsr = Pdb(m);
if numel(Fsr)<16, BEST = empty_out(F,Pdb); return; end
df = mean(diff(Fsr));

% ---- retrait baseline (option) ----
if Opt.useMedianBaseline
    Lb = max(5, round(0.01*numel(Pdbsr)));         % ~1% de la bande
    base = movmedian(Pdbsr, Lb, 'omitnan');
    Pdb_f = Pdbsr - base;                          % pour la détection
else
    Pdb_f = Pdbsr;
end

% ---- pics candidats ----
minDistBins = max(1, round(Opt.minDistHz/df));
[pks,locs,~,prom] = findpeaks(Pdb_f, 'MinPeakProminence',Opt.minPromDB, ...
                                        'MinPeakDistance',minDistBins);
if isempty(locs), BEST = empty_out(F,Pdb); return; end

fcand = Fsr(locs);
if numel(fcand) > Opt.maxPeaks
    [~,ord] = maxk(prom, Opt.maxPeaks);
    fcand = fcand(ord); prom = prom(ord);
end
nC = numel(fcand);

% ---- métriques par candidat ----
Qval   = zeros(nC,1);
Delta  = zeros(nC,1);
SBcomb = zeros(nC,1);   % score "comb"
Cov    = zeros(nC,1);   % couverture multi-ordres (0..1)
Asym   = zeros(nC,1);   % asymétrie moyenne (dB)
NoiseL = zeros(nC,1);   % bruit local médian (dB)

F2i = @(f) max(1, min(numel(Fsr), round( 1 + (f - Fsr(1))*(numel(Fsr)-1)/(Fsr(end)-Fsr(1)) )));

for i = 1:nC
    fc = fcand(i);

    % --- Q (-3 dB) sur PSD "brut" (Pdbsr) ---
    [~,i0] = min(abs(Fsr - fc));
    pk_db = Pdbsr(i0);
    th = pk_db - 3;
    iL = i0; while iL>1 && Pdbsr(iL)>th, iL=iL-1; end
    iR = i0; while iR<numel(Pdbsr) && Pdbsr(iR)>th, iR=iR+1; end
    BW = max(1e-6, Fsr(iR) - Fsr(iL));
    Qval(i) = fc / BW;

    % --- recherche Δ par corrélation peigne ---
    [Delta(i), SBcomb(i), Cov(i), Asym(i), NoiseL(i)] = ...
        best_delta_comb(Fsr, Pdbsr, fc, Opt.K, Opt.DeltaHz, df, Opt.etaBox, Opt.thrSNR_dB, Opt.Delta0_small, F2i);
end

% ---- normalisations robustes ----
PromN = soft01(prom(:));
Qn    = soft01(Qval(:));
CombN = soft01(SBcomb(:));
AsymN = soft01(max(Asym(:),0));
NoiseN= soft01(max(NoiseL(:), median(Pdbsr)));
CovN  = soft01(Cov(:));

% ---- score global "propre" ----
Score = Opt.wProm*PromN + Opt.wQ*Qn + Opt.wSB*CombN ...
      + Opt.wClean*(0.5*(1-AsymN) + 0.5*(1-NoiseN)) ...
      + 0.20*CovN;

% ---- tri & sorties ----
[Score, idx] = sort(Score, 'descend');
fcand = fcand(idx); prom = prom(idx); Qval = Qval(idx);
Delta = Delta(idx); SBcomb = SBcomb(idx); Cov = Cov(idx); Asym = Asym(idx); NoiseL = NoiseL(idx);

candTab = table(fcand(:), prom(:), Qval(:), Delta(:), SBcomb(:), Cov(:), Asym(:), NoiseL(:), Score(:), ...
  'VariableNames', {'fc_Hz','prom_dB','Q','Delta_Hz','CombScore','Coverage','Asym_dB','NoiseLocal_dB','Score'});

best.fc    = fcand(1);
best.Delta = Delta(1);
best.band  = [max(Fsr(1), best.fc-3*best.Delta), min(Fsr(end), best.fc+3*best.Delta)];
best.score = Score(1);
best.Q     = Qval(1);
best.prom_dB = prom(1);
best.CombScore = SBcomb(1);
best.Coverage  = Cov(1);
best.Asym_dB   = Asym(1);
best.Noise_dB  = NoiseL(1);

BEST.best = best;
BEST.candidates = candTab;
BEST.psd = struct('F',F,'Pxx_dB',Pdb);
end

% ---------- helpers locaux ----------
function [DeltaBest, CombBest, CovBest, AsymBest, noiseMed] = best_delta_comb(F, Pdb, fc, K, DeltaHz, df, eta, thrSNR, Delta0, F2i)
dGrid = linspace(DeltaHz(1), DeltaHz(2), max(25, ceil((DeltaHz(2)-DeltaHz(1))/2)));

CombBest = -inf; 
DeltaBest = NaN; 
CovBest = 0; 
AsymBest = Inf; 
noiseMed = median(Pdb, 'omitnan');

for d = dGrid
    boxHz = max(2*df, eta*d);              % fenêtre autour de chaque sideband
    w_k   = @(k) (1 + 0.15*(k-1));         % poids croissant avec l'ordre k

    % --- bruit local (médiane hors boîtes porteuse/sidebands) ---
    iN1 = F2i(fc - 3*d - 2*boxHz);  iN1 = max(1, min(iN1, numel(F)));
    iN2 = F2i(fc + 3*d + 2*boxHz);  iN2 = max(1, min(iN2, numel(F)));
    if iN2 < iN1, [iN1,iN2] = deal(iN2,iN1); end

    maskNoise = true(iN2 - iN1 + 1, 1);

    % masque porteuse
    iC1 = F2i(fc - boxHz);  iC2 = F2i(fc + boxHz);
    iC1 = max(iN1, min(iC1, numel(F))); 
    iC2 = max(1,   min(iC2, numel(F)));
    iC1 = max(iN1, iC1); iC2 = min(iN2, iC2);
    if iC2 >= iC1
        maskNoise((iC1-iN1+1):(iC2-iN1+1)) = false;
    end

    % masque sidebands
    for k = 1:K
        fL = fc - k*d; fR = fc + k*d;
        if fL<F(1) || fR>F(end), break; end
        iL1 = F2i(fL - boxHz); iL2 = F2i(fL + boxHz);
        iR1 = F2i(fR - boxHz); iR2 = F2i(fR + boxHz);
        iL1 = max(iN1, min(iL1, numel(F))); iL2 = max(1, min(iL2, numel(F)));
        iR1 = max(iN1, min(iR1, numel(F))); iR2 = max(1, min(iR2, numel(F)));
        iL1 = max(iN1, iL1); iL2 = min(iN2, iL2);
        iR1 = max(iN1, iR1); iR2 = min(iN2, iR2);
        if iL2 >= iL1, maskNoise((iL1-iN1+1):(iL2-iN1+1)) = false; end
        if iR2 >= iR1, maskNoise((iR1-iN1+1):(iR2-iN1+1)) = false; end
    end

    segNoise = Pdb(iN1:iN2);
    if any(maskNoise)
        noiseMed = median(segNoise(maskNoise), 'omitnan');   % <-- FIX MATLAB
    else
        noiseMed = median(Pdb, 'omitnan');
    end

    % --- accumulation "comb" ---
    Esum = 0; Asum = 0; used = 0; covered = 0;
    for k = 1:K
        fL = fc - k*d; fR = fc + k*d;
        if fL<F(1) || fR>F(end), break; end
        iL1 = F2i(fL - boxHz); iL2 = F2i(fL + boxHz);
        iR1 = F2i(fR - boxHz); iR2 = F2i(fR + boxHz);
        iL1 = max(1, min(iL1, numel(F))); iL2 = max(1, min(iL2, numel(F)));
        iR1 = max(1, min(iR1, numel(F))); iR2 = max(1, min(iR2, numel(F)));
        if iL2 < iL1 || iR2 < iR1, continue; end

        EL = mean(Pdb(iL1:iL2), 'omitnan');
        ER = mean(Pdb(iR1:iR2), 'omitnan');
        SNR_L = max(0, EL - noiseMed);
        SNR_R = max(0, ER - noiseMed);

        Esym = min(SNR_L, SNR_R);
        Asym1 = abs(SNR_L - SNR_R);

        used   = used + 1;
        Esum   = Esum + w_k(k)*Esym;
        Asum   = Asum + Asym1;
        if SNR_L > thrSNR && SNR_R > thrSNR
            covered = covered + 1;
        end
    end

    if used < 2
        continue;  % besoin d'au moins 2 paires crédibles
    end

    coverage = covered / used;             % 0..1
    Psmall   = 1 - exp(-d/Delta0);         % pénalité douce petits Δ

    combScore = (Esum/used) * (1 - 0.4*(Asum/used)) * (0.5 + 0.5*coverage) * Psmall;

    if combScore > CombBest
        CombBest = combScore;
        DeltaBest = d;
        CovBest   = coverage;
        AsymBest  = Asum/used;
    end
end
end

function y = soft01(x)
x = x(:);
q1 = quantile(x,0.10); q9 = quantile(x,0.90);
if q9>q1, y = (x - q1) / (q9 - q1);
else,     y = x - min(x); y = y / max(y + eps);
end
y = max(0, min(1, y));
end

function BEST = empty_out(F,Pdb)
BEST.best = struct('fc',NaN,'Delta',NaN,'band',[NaN NaN],'score',0,'Q',NaN, ...
                   'prom_dB',NaN,'CombScore',NaN,'Coverage',NaN,'Asym_dB',NaN,'Noise_dB',NaN);
BEST.candidates = table();
BEST.psd = struct('F',F,'Pxx_dB',Pdb);
end

function scan = scan_delta_curve(psdS, best, Opt, spanFactor)
% Calcule la courbe CombScore(Δ) autour du Δ optimal
% psdS.F (Hz), psdS.Pxx_dB (dB), best.fc, Opt.K, Opt.etaBox, Opt.thrSNR_dB, Opt.Delta0_small
% spanFactor : largeur de balayage relative (ex: 0.6 -> ±60%)
if nargin<4, spanFactor = 0.6; end

F = psdS.F(:);
Pdb = psdS.Pxx_dB(:);
fc = best.fc;

df = mean(diff(F));

% bande locale utile
Fmin = max(min(F), fc - 3*best.Delta - 80);
Fmax = min(max(F), fc + 3*best.Delta + 80);
m = (F>=Fmin & F<=Fmax);
Fl = F(m); Pl = Pdb(m);
if numel(Fl) < 16
    scan = struct('Delta', [], 'CombScore', [], 'Coverage', [], 'Asym', []);
    return;
end

% grille Δ autour du meilleur
d0   = best.Delta;
dLo  = max( (1-spanFactor)*d0, max(df, 1) );
dHi  = (1+spanFactor)*d0;
dGrid = linspace(dLo, dHi, 61);

% wrappers
F2i = @(f) max(1, min(numel(Fl), round( 1 + (f - Fl(1))*(numel(Fl)-1)/(Fl(end)-Fl(1)) )));
eta = Opt.etaBox; K = Opt.K; thr = Opt.thrSNR_dB; D0s = Opt.Delta0_small;

Comb = -inf(numel(dGrid),1);
Cov  = zeros(numel(dGrid),1);
Asym = zeros(numel(dGrid),1);

for j = 1:numel(dGrid)
    d = dGrid(j);
    boxHz = max(2*df, eta*d);

    % bruit local (médiane hors boîtes)
    iN1 = F2i(fc - 3*d - 2*boxHz);  iN1 = max(1, min(iN1, numel(Fl)));
    iN2 = F2i(fc + 3*d + 2*boxHz);  iN2 = max(1, min(iN2, numel(Fl)));
    if iN2 < iN1, [iN1,iN2] = deal(iN2,iN1); end

    maskNoise = true(iN2-iN1+1,1);

    iC1 = F2i(fc - boxHz); iC2 = F2i(fc + boxHz);
    iC1 = max(iN1, min(iC1, numel(Fl)));
    iC2 = max(1,   min(iC2, numel(Fl)));
    iC1 = max(iN1, iC1); iC2 = min(iN2, iC2);
    if iC2 >= iC1
        maskNoise((iC1-iN1+1):(iC2-iN1+1)) = false;
    end

    for k = 1:K
        fL = fc - k*d; fR = fc + k*d;
        if fL<Fl(1) || fR>Fl(end), continue; end
        iL1 = F2i(fL - boxHz); iL2 = F2i(fL + boxHz);
        iR1 = F2i(fR - boxHz); iR2 = F2i(fR + boxHz);
        iL1 = max(iN1, min(iL1, numel(Fl))); iL2 = max(1, min(iL2, numel(Fl)));
        iR1 = max(iN1, min(iR1, numel(Fl))); iR2 = max(1, min(iR2, numel(Fl)));
        iL1 = max(iN1, iL1); iL2 = min(iN2, iL2);
        iR1 = max(iN1, iR1); iR2 = min(iN2, iR2);
        if iL2 >= iL1, maskNoise((iL1-iN1+1):(iL2-iN1+1)) = false; end
        if iR2 >= iR1, maskNoise((iR1-iN1+1):(iR2-iN1+1)) = false; end
    end

    segNoise = Pl(iN1:iN2);
    if any(maskNoise)
        noiseMed = median(segNoise(maskNoise), 'omitnan');   % <- FIX
    else
        noiseMed = median(Pl, 'omitnan');
    end

    Esum=0; Asum=0; used=0; covered=0;
    for k = 1:K
        fL = fc - k*d; fR = fc + k*d;
        if fL<Fl(1) || fR>Fl(end), continue; end
        iL1 = F2i(fL - boxHz); iL2 = F2i(fL + boxHz);
        iR1 = F2i(fR - boxHz); iR2 = F2i(fR + boxHz);
        iL1 = max(1, min(iL1, numel(Fl))); iL2 = max(1, min(iL2, numel(Fl)));
        iR1 = max(1, min(iR1, numel(Fl))); iR2 = max(1, min(iR2, numel(Fl)));
        if iL2 < iL1 || iR2 < iR1, continue; end

        EL = mean(Pl(iL1:iL2), 'omitnan');
        ER = mean(Pl(iR1:iR2), 'omitnan');
        SNRL = max(0, EL - noiseMed);
        SNRR = max(0, ER - noiseMed);

        Esym  = min(SNRL, SNRR);
        Asym1 = abs(SNRL - SNRR);

        used   = used + 1;
        Esum   = Esum + (1 + 0.15*(k-1))*Esym;
        Asum   = Asum + Asym1;
        if SNRL>thr && SNRR>thr, covered = covered + 1; end
    end

    if used>=2
        coverage = covered/used;
        Psmall   = 1 - exp(-d/D0s);
        Comb(j)  = (Esum/used) * (1 - 0.4*(Asum/used)) * (0.5 + 0.5*coverage) * Psmall;
        Cov(j)   = coverage;
        Asym(j)  = Asum/used;
    end
end

scan = struct('Delta', dGrid(:), 'CombScore', Comb(:), 'Coverage', Cov(:), 'Asym', Asym(:));
end

function plot_pick_best_comb(result, Opt, Kshow)
% Visualisation 3 figures:
%  (1) Panorama PSD + Top-K candidats
%  (2) Zoom sur porteuse gagnante + sidebands ±kΔ
%  (3) Courbe CombScore(Δ) autour du meilleur Δ (avec couverture & asym)
if nargin<2, Opt = struct; end
if nargin<3, Kshow = 10; end

F = result.psd.F(:);
P = result.psd.Pxx_dB(:);
C = result.candidates;
B = result.best;

% === Figure 1: Panorama ===
figure('Name','Panorama – PSD + candidats','Color','w');
plot(F, P, 'LineWidth',1.0); grid on; hold on;
xlabel('Fréquence (Hz)'); ylabel('PSD (dB)');
title('Panorama : PSD (Welch) et candidats triés');

K = min(height(C), Kshow);
cmap = lines(K);
for i = 1:K
    fc = C.fc_Hz(i);
    y  = interp1(F, P, fc, 'linear', 'extrap');
    plot(fc, y, 'o', 'MarkerSize',6, 'LineWidth',1.2, 'Color', cmap(i,:), ...
        'DisplayName', sprintf('#%d  f=%.2f Hz  S=%.2f', i, fc, C.Score(i)));
end
legend('Location','best');

txt = sprintf('BEST: f_c=%.2f Hz, Δ=%.2f Hz, Score=%.2f, Q=%.1f', B.fc, B.Delta, B.score, B.Q);
yl = ylim; xlim([min(F) max(F)]);
text(min(F)+0.01*range(F), yl(2)-0.05*range(yl), txt, 'FontWeight','bold');

% === Figure 2: Zoom ===
figure('Name','Zoom – porteuse & sidebands','Color','w');
zoomHalf = max(7*B.Delta, 40);
fmin = max(min(F), B.fc - zoomHalf);
fmax = min(max(F), B.fc + zoomHalf);
mask = (F>=fmin & F<=fmax);
plot(F(mask), P(mask), 'LineWidth',1.2); grid on; hold on;
xlabel('Fréquence (Hz)'); ylabel('PSD (dB)');
title(sprintf('Zoom: f_c=%.2f Hz, Δ=%.2f Hz (K=±3)', B.fc, B.Delta));

% bande ±3Δ
draw_band([B.fc-3*B.Delta, B.fc+3*B.Delta], [0.92 0.96 1.00]);
yL = ylim;
plot([B.fc B.fc], yL, '--', 'LineWidth',1.2, 'Color',[0.2 0.2 0.9]);

Ksb = 3;
for k = 1:Ksb
    fL = B.fc - k*B.Delta; fR = B.fc + k*B.Delta;
    draw_sb(fL, yL, 2);
    draw_sb(fR, yL, 2);
end

% encart infos
txt2 = sprintf(['f_c=%.2f Hz\nΔ=%.2f Hz\nScore=%.2f\nQ=%.1f\n' ...
                'Comb=%.2f\nCov=%.2f\nAsym=%.2f dB\nNoise=%.1f dB'], ...
               B.fc, B.Delta, B.score, B.Q, B.CombScore, B.Coverage, B.Asym_dB, B.Noise_dB);
xBox = fmin + 0.02*(fmax - fmin);
yBox = yL(2) - 0.06*range(yL);
text(xBox, yBox, txt2, 'BackgroundColor',[1 1 1]*0.95, 'Margin',6, 'EdgeColor',[0.7 0.7 0.7]);

xlim([fmin fmax]);

% === Figure 3: Courbe Δ ===
figure('Name','Courbe Δ – comb score','Color','w');
scan = scan_delta_curve(result.psd, result.best, Opt, 0.6);
yyaxis left;
plot(scan.Delta, scan.CombScore, '-', 'LineWidth',1.8); grid on; hold on;
ylabel('CombScore (a.u.)');
yyaxis right;
plot(scan.Delta, scan.Coverage, '--', 'LineWidth',1.2);
plot(scan.Delta, 1 - normalize01(scan.Asym), ':', 'LineWidth',1.2);
ylabel('Coverage / (1 - Asym norm)');
xlabel('\Delta (Hz)');
title(sprintf('Δ-scan autour de %.2f Hz (meilleur Δ=%.2f Hz)', B.fc, B.Delta));

legend({'CombScore','Coverage','1 - Asym (norm)'}, 'Location','best');

end

% --- helpers graphiques ---
function draw_band(band, rgb)
yl = ylim;
patch([band(1) band(2) band(2) band(1)], [yl(1) yl(1) yl(2) yl(2)], rgb, ...
      'FaceAlpha',0.15, 'EdgeColor','none'); 
end

function draw_sb(f, yL, halfWidthHz)
xline(f, '-', sprintf('%.1f Hz', f), 'Alpha',0.6);
patch([f-halfWidthHz f+halfWidthHz f+halfWidthHz f-halfWidthHz], ...
      [yL(1) yL(1) yL(2) yL(2)], [1.0 0.95 0.9], 'FaceAlpha',0.15, 'EdgeColor','none');
end

function y = normalize01(x)
x = x(:);
a = quantile(x(~isnan(x)),0.05); b = quantile(x(~isnan(x)),0.95);
if b>a, y = (x-a)/(b-a); else, y = zeros(size(x)); end
y = max(0, min(1, y));
end

function [S, info] = split_top_candidates(x, Fs, res, varargin)
% Sépare jusqu'à 10 candidats en signaux individuels (colonnes de S)
% x    : signal (Nx1)
% Fs   : Hz
% res  : sortie de pick_best_carrier_comb (res.candidates utilisé)
% Options Name-Value:
%   'NumSignals'    (default 10)   : nb max de candidats à extraire
%   'K'             (default 3)    : nb d'ordres de sidebands à couvrir (±K)
%   'WidthFactor'   (default 3)    : bande = [fc - W*Δ, fc + W*Δ]
%   'Method'        (default 'fft'): 'fft' ou 'fir'
%   'FIR_Order'     (default 512)  : ordre FIR si Method='fir'
%   'TransWidthHz'  (default 0.15): largeur de transition en fraction de Δ (ex: 0.15 => 15% de Δ)
%
% Out:
%   S    : (N x M) signaux filtrés, M = min(NumSignals, height(res.candidates))
%   info : table avec fc, Delta, bandHz et méthode

p = inputParser;
addParameter(p,'NumSignals',10);
addParameter(p,'K',3);
addParameter(p,'WidthFactor',3);
addParameter(p,'Method','fft');
addParameter(p,'FIR_Order',512);
addParameter(p,'TransWidthHz',0.15);
parse(p,varargin{:});
Opt = p.Results;

x = x(:);
N = numel(x);
M = min(Opt.NumSignals, height(res.candidates));
if M==0
    S = zeros(N,0); info = table(); return;
end

fc   = res.candidates.fc_Hz(1:M);
dlt  = res.candidates.Delta_Hz(1:M);
Wfac = Opt.WidthFactor;

S = zeros(N, M);
bands = zeros(M,2);

switch lower(Opt.Method)
    case 'fft'
        % --- 1) Masquage en fréquence et IFFT (offline) ---
        X = fft(x);
        F = (0:N-1)'*(Fs/N);
        % garder uniquement la demi-bande [0, Fs/2]; miroir pour négatifs
        for i = 1:M
            b1 = max(0, fc(i) - Wfac*dlt(i));
            b2 = min(Fs/2, fc(i) + Wfac*dlt(i));
            bands(i,:) = [b1 b2];

            mask = (F>=b1 & F<=b2) | (F>=Fs-b2 & F<=Fs-b1); % symétrie conjuguée
            Xf = zeros(size(X));
            Xf(mask) = X(mask);
            s = real(ifft(Xf));

            % Optionnel : atténuer chevauchements (non nécessaire pour test)
            S(:,i) = s;
        end

    case 'fir'
        % --- 2) Banc de filtres passe-bande FIR + filtfilt (zéro-phase) ---
        for i = 1:M
            b1 = max(0, fc(i) - Wfac*dlt(i));
            b2 = min(Fs/2, fc(i) + Wfac*dlt(i));
            bands(i,:) = [b1 b2];

            % transition proportionnelle à Δ (évite bandes trop "dures")
            tw = max(1.0, Opt.TransWidthHz * dlt(i));   % en Hz, >= 1 Hz
            fpass = [max(0, b1), min(Fs/2, b2)];
            fstop = [max(0, b1 - tw), min(Fs/2, b2 + tw)];

            % normalisation
            wpass = fpass/(Fs/2);
            wstop = fstop/(Fs/2);
            wstop(1) = max(0, wstop(1)); wstop(2) = min(1, wstop(2));
            wpass(1) = max(0, wpass(1)); wpass(2) = min(1, wpass(2));

            % FIR par fenêtre (kaiser ou hann via fir1)
            ord = Opt.FIR_Order;
            b = fir1(ord, wpass, 'bandpass', hann(ord+1));

            % zéro-phase (pas de retard de groupe)
            S(:,i) = filtfilt(b, 1, x);
        end

    otherwise
        error('Method must be ''fft'' or ''fir''.');
end

info = table(fc, dlt, bands(:,1), bands(:,2), repmat(string(Opt.Method),M,1), ...
    'VariableNames', {'fc_Hz','Delta_Hz','BandLow_Hz','BandHigh_Hz','Method'});
end

function [ranked, metrics] = rank_candidates_by_acf(S, Fs, varargin)
% Classe les signaux (colonnes de S) selon la "qualité ACF".
% In:
%   S  : (N x M) signaux (ex: sortis de split_top_candidates)
%   Fs : Hz
% Options:
%   'MaxLagSec'   : durée max pour ACF (par défaut 0.5 s)
%   'ACFWin'      : nb d'échantillons max pour l'ACF (par défaut 8192)
%   'TopK'        : nb de candidats à comparer en plot métriques (par défaut 5)
%   'DoPlots'     : true/false (par défaut true)
% Out:
%   ranked  : table triée (Score + métriques) par candidat
%   metrics : struct array avec détails et ACF (pour debug)

p = inputParser;
addParameter(p,'MaxLagSec', 0.5);
addParameter(p,'ACFWin', 8192);
addParameter(p,'TopK', 5);
addParameter(p,'DoPlots', true);
parse(p,varargin{:});
opt = p.Results;

[N, M] = size(S);
if M==0, ranked = table(); metrics = struct([]); return; end

% borne lags/longueur
Lmax = min([N, opt.ACFWin]);
lagMax = min([floor(opt.MaxLagSec*Fs), Lmax-1]);

% poids (tu peux ajuster)
w = [0.25 0.25 0.20 0.15 0.15];   % [Qspec Qent Qsfm Qpksep Qprom]

metrics = repmat(struct( ...
    'Qspec',0,'Qent',0,'Qsfm',0,'Qpksep',0,'Qprom',0, ...
    'Score',0,'acf',[],'lags',[],'info',struct()), 1, M);

for m = 1:M
    x = S(:,m);
    % --- ACF rapide via FFT (fenêtrée) ---
%     [r,lags] = local_acf_fft(x, lagMax, Lmax);
    signalEnvelope = getEnv(x, Fs, res.best.band(1));% res.best.fc BEST.band(1)-50);  %   extracts the envelope
    [AcfIn] = nsdf_mcleod_fast(signalEnvelope, Fs, 0, 160, numel(signalEnvelope));

    % --- métriques ACF (robustes) ---
    Mx = acf_quality_metrics2(r);         % fonction B ci-dessous
    score = acf_quality_score(Mx, w);     % fonction C ci-dessous

    metrics(m).Qspec = Mx.Qspec;
    metrics(m).Qent  = Mx.Qent;
    metrics(m).Qsfm  = Mx.Qsfm;
    metrics(m).Qpksep= Mx.Qpksep;
    metrics(m).Qprom = Mx.Qprom;
    metrics(m).Score = score;
    metrics(m).acf   = r;
    metrics(m).lags  = lags/Fs;           % en secondes
    metrics(m).info  = Mx.debug;
end

% --- tableau + tri ---
Qspec = [metrics.Qspec]';
Qent  = [metrics.Qent]';
Qsfm  = [metrics.Qsfm]';
Qpksep= [metrics.Qpksep]';
Qprom = [metrics.Qprom]';
Score = [metrics.Score]';

idx = (1:M).';
ranked = table(idx, Score, Qspec, Qent, Qsfm, Qpksep, Qprom, ...
    'VariableNames', {'Candidate','Score','Qspec','Qent','Qsfm','Qpksep','Qprom'});
ranked = sortrows(ranked, 'Score', 'descend');

% --- plots (optionnels) ---
if opt.DoPlots
    % 1) Barres des scores tous candidats
    figure('Color','w','Name','ACF quality – Scores');
    bar(ranked.Candidate, ranked.Score); grid on
    xlabel('Candidat'); ylabel('Score ACF (0..1)');
    title('Classement par qualité ACF (plus haut = mieux)');

    % 2) Comparaison des métriques pour TopK
    K = min(opt.TopK, M);
    topIdx = ranked.Candidate(1:K);
    Mmat = [ [metrics(topIdx).Qspec]' [metrics(topIdx).Qent]' [metrics(topIdx).Qsfm]' ...
             [metrics(topIdx).Qpksep]' [metrics(topIdx).Qprom]' ];
    figure('Color','w','Name','ACF quality – Metrics (TopK)');
    bar(Mmat, 'grouped'); grid on
    xticklabels(string(topIdx));
    legend({'Qspec','Qent','Qsfm','Qpksep','Qprom'},'Location','northoutside','Orientation','horizontal');
    xlabel('Candidat (index colonne S)'); ylabel('Score métrique (0..1)');
    title('Comparaison des métriques ACF – TopK');

    % 3) Détail du #1 : ACF + spectre de l’ACF (débaseliné)
    c1 = ranked.Candidate(1);
    r  = metrics(c1).acf;
    tlag = metrics(c1).lags;
    figure('Color','w','Name','ACF quality – Detail Top1');
    tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

    % ACF
    nexttile; plot(tlag, r, 'LineWidth',1.2); grid on
    xlabel('Lag (s)'); ylabel('ACF (norm.)'); title(sprintf('ACF – candidat #%d', c1));

    % Spectre ACF (avec baseline médiane)
    nexttile;
    Rspec = abs(fft(r)).^2; Rspec = Rspec(1:floor(numel(r)/2));
    Lb = max(5, round(0.01*numel(Rspec))); base = movmedian(Rspec,Lb);
    Rs = max(Rspec - base, 0);
    fbin = (0:numel(Rs)-1)';  % en "bins" (pas de Fs ici: ACF de x)
    plot(fbin, 10*log10(Rs+eps), 'LineWidth',1.0); grid on
    xlabel('Bin'); ylabel('ACF-spectrum (dB, baseline-removed)');
    title('Spectre de l’ACF (débaseliné)');
end
end

% ---------- helpers locaux ----------
function [r,lags] = local_acf_fft(x, lagMax, Lmax)
% ACF à partir d'une fenêtre max Lmax, lissée par Hann, normalisée r(0)=1
x = x(:) - mean(x,'omitnan');
N = numel(x);
L = min([N, Lmax]);
x = x(1:L);
w = hann(L,'periodic');
xw = x .* w;
M = 2^nextpow2(2*L-1);
S = fft(xw, M);
r = ifft(abs(S).^2, 'symmetric');
r = r(1:lagMax+1);
if r(1) ~= 0, r = r / r(1); else, r = r / max(abs(r)+eps); end
lags = 0:lagMax;
end

function M = acf_quality_metrics2(r, Opt)
% Qualité d'une ACF (sans Qacf), robustifiée.
% Sortie M: struct avec Qspec, Qent, Qsfm, Qpksep, Qprom (+ debug)

if nargin<2, Opt = struct; end
DEF = struct('IgnoreDC', true, 'BaselineFrac', 0.01, 'MedianLenMin', 5, ...
             'Taper', true, 'Win', 'tukey', 'TukeyAlpha', 0.1, ...
             'MinPeakDistBins', 2);
fn = fieldnames(DEF);
for k=1:numel(fn), if ~isfield(Opt,fn{k}), Opt.(fn{k}) = DEF.(fn{k}); end, end

r = r(:);
r = r - mean(r,'omitnan');
if ~isempty(r) && r(1) ~= 0
    r = r / max(r(1), eps);
else
    r = r / max(abs(r)+eps);
end
N = numel(r);
if N < 16
    M = empty_out_acf_metrics(); return;
end

if Opt.Taper
    switch lower(Opt.Win)
        case 'hann'
            w = hann(N);
        otherwise
            w = tukeywin(N, Opt.TukeyAlpha);
    end
    r = r .* w;
end

Rspec = abs(fft(r)).^2;
Rspec = Rspec(1:floor(N/2));
if Opt.IgnoreDC && ~isempty(Rspec)
    if numel(Rspec) >= 4
        Rspec(1) = median(Rspec(2:min(8,end)));
    else
        Rspec(1) = 0;
    end
end

Lb = max(Opt.MedianLenMin, round(Opt.BaselineFrac * numel(Rspec)));
base = movmedian(Rspec, Lb, 'omitnan');
Rs = max(Rspec - base, 0);

S = sum(Rs);
if S <= 0
    M = empty_out_acf_metrics(); return;
end

Qspec = max(Rs) / S;

ps = Rs / S; ps = ps(ps>0);
H  = -sum(ps .* log(ps));
Hmax = log(max(numel(ps),1));
Qent = 1 - H / max(Hmax, eps);

AF = mean(Rs);
GF = exp(mean(log(Rs + eps)));
Qsfm = 1 - min(GF/AF, 1);

[pkVals, pkLocs] = findpeaks(Rs, 'MinPeakDistance', Opt.MinPeakDistBins);
pkVals = sort(pkVals, 'descend');
if numel(pkVals) >= 2
    p1 = pkVals(1); p2 = pkVals(2);
    Qpksep = (p1 - p2) / max(p1 + p2, eps);
elseif numel(pkVals) == 1
    Qpksep = 1;
else
    Qpksep = 0;
end

if ~isempty(pkLocs)
    [~, imax] = max(Rs);
    [~, ~, ~, prom] = findpeaks(Rs, 'MinPeakDistance', Opt.MinPeakDistBins);
    if ~isempty(prom), pmax = max(prom); else, pmax = Rs(imax); end
    Qprom = pmax / max(pmax + median(Rs), eps);
else
    Qprom = 0;
end

M = struct('Qspec', Qspec, ...
           'Qent',  Qent, ...
           'Qsfm',  Qsfm, ...
           'Qpksep',Qpksep, ...
           'Qprom', Qprom, ...
           'peaks', struct('locs', pkLocs, 'vals', pkVals), ...
           'debug', struct('Rspec', Rspec, 'Rs', Rs, 'baselineLen', Lb));
end

function M = empty_out_acf_metrics()
M = struct('Qspec',0,'Qent',0,'Qsfm',0,'Qpksep',0,'Qprom',0, ...
           'peaks',struct('locs',[],'vals',[]), ...
           'debug',struct('Rspec',[],'Rs',[],'baselineLen',[]));
end


function s = acf_quality_score(M, w)
% Combinaison linéaire de métriques en [0,1]
if nargin<2, w = [0.25 0.25 0.20 0.15 0.15]; end
v = [M.Qspec, M.Qent, M.Qsfm, M.Qpksep, M.Qprom];
s = sum(w(:)'.*v);
end

