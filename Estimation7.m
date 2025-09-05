% close all;
% clear;
FsOr = 12000;
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


array   =   array1A;
freqRot = freqRot1 ;

% test = load('C:/Users/GJ2640/Downloads/KI21/KI21/N15_M07_F10_KI21_1.mat'); % 1500
% FsOr=64000;
% freqRot = 1500/60;
% array=test.N15_M07_F10_KI21_1.Y(7).Data(:);
% sound(array,FsOr);
%%


sizeSlide=8192*2;
n=floor(numel(array)/sizeSlide);

for slide_idx = 1:n
    
    accOrigi=array((sizeSlide*(slide_idx-1))+1:(sizeSlide*slide_idx));
    
    dec = 8;
    accTest = decimate(accOrigi, 8, 1, 'fir');
    Fs = FsOr /dec;
    

    dec = 2;
    acc = decimate(accOrigi, dec, 1, 'fir');
    Fs = FsOr /dec;
    plotFFT(acc,Fs);
    
    signalEnvelope = getEnv(acc, Fs, 1300);  %   extracts the envelope
    [AcfIn] = nsdf_mcleod_fast(signalEnvelope, Fs, 0, 160, numel(signalEnvelope));

%     N=numel(nsdf);
%     [nsdf1] = nsdf_mcleod_fast(signalEnvelope(1:N/2), Fs, 0.5, 120, N);
%     [nsdf2] = nsdf_mcleod_fast(signalEnvelope(N/4+1:N/4+N/2), Fs, 0.5, 120, N);
%     [nsdf3] = nsdf_mcleod_fast(signalEnvelope(N/2+1:end), Fs, 0.5, 120, N);
%     AcfIn = nsdf1.*nsdf2.*nsdf3;
    
    
    %% parametres
    
    
    MaxLagSec =   1.0; % analyser l'ACF jusqu'à 1 s
    TopPeaks =   6; % # pics max pour la médiane des intervalles
    Prominence = 0.05;  % proéminence min (fraction du max hors lag 0)
    SmoothMs =   4.0;   % lissage avant détection de pics (ms)
    PsdFmax  =   300;
    fmin        = 0.5;     % bande de recherche (Hz)
    fmax        = 120;     % bande de recherche (Hz)
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
    
    %% cepstrum
    [ c_q, q_axis] = coarse_cepstrum_f0(acc, Fs, [0.5 120], 'Win','hann','Dither',1e-3);
    c_q=c_q(1:min(numel(acc),Fs/fmin));
    q_axis=q_axis(1:min(numel(acc),Fs/fmin));
    f_axis = 1 ./ q_axis;
    mask = (f_axis > fmax);
    c_q(mask) = 0;
    [~,idx]=max(c_q);
    candf_Cepstrum = 1/q_axis(idx)*4;

    %% Estimation de période par médiane des différences successives
    [pkVal, pkLoc] = findpeaks(AcfIn, 'MinPeakProminence', prom, 'MinPeakDistance', minPkDist);   % pkLoc = lags, pkVal = amplitudes
    medianVal = median(pkVal);
    mask = (pkVal > medianVal);
    pkVal=pkVal(mask);
    pkLoc=pkLoc(mask);
    L_best = median(diff(pkLoc));
    [~, iNear] = min(abs(double(pkLoc) - L_best));
    L = pkLoc(iNear);
    [T_hat,L_ref] =interp_parabolic_acf(AcfIn, L, Fs, tauMaxSamp);

    candf_mediane=Fs/L_ref;

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
    fcandstring = ["Cepst  ","yin    ", "comb   ", "hps    ","mediane"];
    fcands = [candf_Cepstrum,candf_yin, candf_comb, candf_hps,candf_mediane];
    


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
    vals = [freqRot,candf_Cepstrum,candf_yin, candf_comb, candf_hps, f_twm, candf_mediane];
    names= {'real','Cpst','YIN','Comb','HPS','TWM','MEDIANE'};
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

% function v = safe_interp(x, y, q)
% % interpolation linéaire sûre aux indices fractionnaires (en échantillons)
% % x : 0..Nlag-1, y : même taille, q : scalaire (ou vecteur) en "samples"
% q = q(:);
% v = nan(size(q));
% Nlag = numel(x);
% for i=1:numel(q)
%     if q(i) < 1 || q(i) > Nlag-2
%         v(i) = NaN; %#ok<*AGROW>
%     else
%         i0 = floor(q(i));
%         a  = q(i) - i0;
%         v(i) = (1-a)*y(i0+1) + a*y(i0+2);
%     end
% end
% if isscalar(v), v=v(1); end
% end



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



function out = tern(c,a,b), if c, out=a; else, out=b; end, end

function dp = interp_dprime(dprime, q)
Nlag = numel(dprime);
if q<1, q=1; end, if q>Nlag-1, q=Nlag-1; end
i0 = floor(q); a = q - i0;
dp = (1-a)*dprime(i0) + a*dprime(i0+1);
end

function [fr_hat, T_hat, DEC] = anti_half_acf_unified(r_in, Fs, T_hat_in, tauMaxSamp, K, lambdaEven, gammaHalf, opts)
% ANTI_HALF_ACF_UNIFIED
%  Décide entre T, T/2 (et optionnellement 2T) en combinant :
%   - ratio ACF : R(T/2) vs R(T)
%   - peigne odd/even ACF
%   - (option) critère YIN : d'(T/2) << d'(T)
%
% In:
%   r_in        : ACF (lags >=0). Peut être brute; on normalise R(0)=1.
%   Fs          : Hz
%   T_hat_in    : période candidate (secondes)
%   tauMaxSamp  : max lag utile (échantillons)
%   K           : nb harmoniques pour peigne odd/even (2..6)
%   lambdaEven  : pondération "even" (0.8..1.0)
%   gammaHalf   : seuil ratio R(T/2) > gamma * R(T) (1.05..1.2)
%   opts (struct, facultatif):
%       .UseYIN      (default true)  % utilise dprime_s si fourni
%       .dprime_s    ([])            % CMNDF lissée si dispo
%       .Lref        ([])            % lag frac. correspondant à T_hat_in (échantillons)
%       .Consider2T  (default true)  % essaie aussi 2T
%       .ClampBand   ([fmin fmax])   % clamp final en Hz (ex: [0.5 120])
%
% Out:
%   fr_hat  : fréquence corrigée (Hz)
%   T_hat   : période corrigée (s)
%   DEC     : diagnostics (struct)

if nargin<8, opts = struct(); end
if ~isfield(opts,'UseYIN'),     opts.UseYIN = true; end
if ~isfield(opts,'dprime_s'),   opts.dprime_s = []; end
if ~isfield(opts,'Lref'),       opts.Lref = []; end
if ~isfield(opts,'Consider2T'), opts.Consider2T = true; end
if ~isfield(opts,'ClampBand'),  opts.ClampBand = []; end

% --- normalisation sûre ---
r = r_in(:);
r = r / max(r(1), eps);
lags = (0:numel(r)-1).';

interpR = @(q) safe_interp(lags, r, q);   % linéaire
T   = max(1/Fs, T_hat_in);
Ts  = T * Fs;                              % en échantillons
tauMaxSamp = min(tauMaxSamp, numel(r)-1);

% --- fonctions utilitaires ---
    function S = comb_odd_even(Tsamp)
        Sodd=0; Seven=0;
        for kk=1:K
            q = kk*Tsamp;
            if q<=tauMaxSamp
                val = interpR(q); val = max(val,0);
                if mod(kk,2)==1, Sodd=Sodd+val; else, Seven=Seven+val; end
            end
        end
        S = (Sodd - lambdaEven*Seven) / max(Sodd + lambdaEven*Seven, eps);
    end

% --- scores pour T, T/2, (2T option) ---
RT  = interpR(Ts);
RT2 = interpR(Ts/2);
S_T  = comb_odd_even(Ts);
S_H  = comb_odd_even(Ts/2);

scoreT  = (RT  - gammaHalf*max(RT2,0)) + 0.20*S_T;
scoreH  = (RT2 - gammaHalf*max(RT ,0)) + 0.20*S_H;

% Critère YIN optionnel
yinBoost_T = 0; yinBoost_H = 0;
if opts.UseYIN && ~isempty(opts.dprime_s)
    % si Lref inconnu, approx via T*Fs
    Lref = tern(~isempty(opts.Lref), opts.Lref, Ts);
    dT  = interp_dprime(opts.dprime_s, Lref);
    dH  = interp_dprime(opts.dprime_s, Lref/2);
    if isfinite(dT) && isfinite(dH)
        if dH < 0.85 * dT, yinBoost_H = 0.15; end  % bonus si vallée T/2 nettement meilleure
    end
end
scoreT = scoreT + yinBoost_T;
scoreH = scoreH + yinBoost_H;

% Option 2T
score2 = -Inf; T2 = 2*T; S_2 = NaN; R2 = NaN;
if opts.Consider2T
    R2 = interpR(2*Ts);
    S_2 = comb_odd_even(2*Ts);
    score2 = (R2 - gammaHalf*max(RT,0)) + 0.20*S_2;
end

% Décision
[~,ix] = max([scoreT, scoreH, score2]);
if ix==2
    T_hat = T/2;
elseif ix==3
    T_hat = T*2;
else
    T_hat = T;
end

% Clamp bande fréquence (facultatif)
fr_hat = 1 / max(T_hat, eps);
if ~isempty(opts.ClampBand) && numel(opts.ClampBand)==2
    fr_hat = min(max(fr_hat, opts.ClampBand(1)), opts.ClampBand(2));
    T_hat  = 1 / fr_hat;
end

% Diagnostics
combContrast = comb_odd_even(T_hat*Fs);
DEC = struct('RT',RT,'RT2',RT2,'R2',R2, ...
    'S_T',S_T,'S_H',S_H,'S_2',S_2, ...
    'scoreT',scoreT,'scoreH',scoreH,'score2',score2, ...
    'combContrast',combContrast);
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

function [nsdfW, Vcorr, nsdfStd] = nsdf_mcleod_correntropy_mix(x, Fs, fmin, fmax, sigma, alpha)
% alpha in [0..1] : mix additif final  S = alpha*NSDF_std + (1-alpha)*Vcorr

x = double(x(:));
N = numel(x);
tauMax = min(N-2, ceil(Fs/max(fmin,eps)));
tauMin = max(2, floor(Fs/max(fmax,eps)));

% NSDF standard (pour mix)
x2 = x.^2; E = cumsum(x2); E0 = E(end);
nsdfStd = zeros(tauMax+1,1); nsdfStd(1)=1;

% Correntropy simple
Vcorr = zeros(tauMax+1,1); Vcorr(1)=1;

% NSDF pondérée par correntropy
nsdfW = zeros(tauMax+1,1); nsdfW(1)=1;

for T = tauMin:tauMax
    M = N - T;
    if M<=0, break; end
    x1 = x(1+T:N);
    x0 = x(1:N-T);
    d  = x1 - x0;
    
    % ----- correntropy -----
    w = exp( - (d.^2) / (2*sigma^2) );     % poids gaussiens
    Vcorr(1+T) = mean(w);
    
    % ----- NSDF standard -----
    R = sum(x1 .* x0);
    denomStd = E0 + E(N-T);
    nsdfStd(1+T) = (2*R) / max(denomStd, eps);
    
    % ----- NSDF pondérée (kernel-weighted) -----
    Rw = sum(w .* (x1 .* x0));
    denomW = sum(w .* (x1.^2 + x0.^2));
    nsdfW(1+T) = (2*Rw) / max(denomW, eps);
end

% Option : fusion additive (si tu veux un seul score robuste)
if nargin >= 6 && ~isempty(alpha)
    nsdfW = alpha*nsdfStd + (1-alpha)*Vcorr;  % score hybride
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

function w = weights_K(K, mode)
switch lower(string(mode))
    case "equal", w = ones(1,K);
    otherwise,    w = 1./(1:K);
end
w = w / sum(w);
end
function [scoreComb, Sodd, Seven] = comb_score(r, T_samp, K, wk, lambdaEven, useLog)
tauMax = numel(r)-1; interpR = @(q) safe_interp((0:tauMax).', r, q);
Sodd=0; Seven=0;
for k=1:K
    q = k*T_samp; if q>tauMax, break; end
    val = interpR(q); val = max(val,0);
    if useLog, val = log(1e-3 + val); end
    if mod(k,2)==1, Sodd = Sodd + wk(k)*val; else, Seven = Seven + wk(k)*val; end
end
scoreComb = (Sodd - lambdaEven*Seven) / max(abs(Sodd)+lambdaEven*abs(Seven), eps);
end

function s = odd_even_spec(Spec, f0, K, wk, lambdaEven)
Sodd=0; Seven=0;
for k=1:K
    fk = k*f0;
    if mod(k,2)==1, Sodd = Sodd + wk(k)*Spec(fk);
    else            Seven= Seven+ wk(k)*Spec(fk);
    end
end
s = (Sodd - lambdaEven*Seven) / max(Sodd + lambdaEven*Seven, eps);
end

function E = band_energy(F, Pxx, f0, bw)
if f0<=0 || f0>F(end), E=0; return; end
m = (F >= (f0-bw)) & (F <= (f0+bw));
if ~any(m), E=0; return; end
E = trapz(F(m), Pxx(m));
end

function [L0, S] = pickL_from_acf_harmonics(AcfIn, Lref, Kharm, tauMaxSamp, useInterp)
% PICKL_FROM_ACF_HARMONICS
%  Compare la somme des harmoniques de l'ACF aux multiples de Lref:
%  Lref, Lref/2, 2*Lref, puis choisit celui qui maximise la somme.
%
% In:
%   AcfIn        : ACF (lags >= 0), vecteur ligne/colonne
%   Lref         : lag de référence (en échantillons, flottant possible)
%   Kharm        : nombre d'harmoniques à sommer (ex: 6)
%   tauMaxSamp   : lag max exploitable (<= numel(AcfIn)-1)
%   useInterp    : true -> interp linéaire (plus précis), false -> round index
%
% Out:
%   L0 : lag choisi (Lref, Lref/2 ou 2*Lref)
%   S  : struct avec les sommes pour debug (S.base, S.half, S.double)

if nargin < 5, useInterp = true; end

AcfIn = AcfIn(:);                     % colonne
nLag  = numel(AcfIn) - 1;             % lags valides: [0 .. nLag]
tauMaxSamp = min(max(1, floor(tauMaxSamp)), nLag);

H   = (1:Kharm);                      % harmoniques 1..K
Ls  = [Lref, Lref/2, 2*Lref];         % candidats
Sarr = zeros(1,3);

if useInterp
    % --- Interpolation linéaire sur ACF ---
    x = (0:nLag).';                   % grille des lags
    for j = 1:3
        q = H .* Ls(j);               % positions harmoniques (flottantes)
        m = (q <= tauMaxSamp);        % masque bornes
        if any(m)
            % extrapolation=0 pour ignorer hors-borne
            v = interp1(x, AcfIn, q(m), 'linear', 0);
            % clamp négatifs légers si besoin (optionnel) :
            % v = max(v, 0);
            Sarr(j) = sum(v);
        end
    end
else
    % --- Version indices entiers (plus rapide, moins précise) ---
    for j = 1:3
        idx = round(H .* Ls(j));
        m = (idx >= 1) & (idx <= tauMaxSamp);
        if any(m)
            Sarr(j) = sum(AcfIn(idx(m)));
        end
    end
end

% Décision
[~, k] = max(Sarr);
L0 = Ls(k);

% Debug / diagnostics
S = struct('base', Sarr(1), 'half', Sarr(2), 'double', Sarr(3));
end

function [L_best, score_best, keptIdx_best, keptVal_best, H_best, diag] = ...
    pick_f0_by_self_harmonics(peakIdx, peakVal, K, varargin)
% PICK_F0_BY_SELF_HARMONICS
%  Estime le lag fondamental L_best (en échantillons) en testant chaque pic
%  comme candidat fondamental Li = peakIdx(i), et en ne gardant que ses K
%  harmoniques (h=1..K) présents parmi les pics détectés, avec tolérance.
%  Pour chaque h, on garde au plus un pic (le plus fort), on somme (pondération
%  optionnelle), et on retient le Li qui maximise le score.
%
% In:
%  peakIdx : [N×1] lags (échantillons) des pics (ex: de l'ACF)
%  peakVal : [N×1] amplitudes/poids des pics
%  K       : nombre d’harmoniques considérées (h=1..K)
%
% Options (name/value):
%  'RelTol'     (0.06)   : tolérance relative ±RelTol*Li (si AbsTol vide)
%  'AbsTol'     ([])     : tolérance absolue (échantillons) (prioritaire)
%  'Weights'    ([])     : pondération w(h) (def: 1./sqrt(h))
%                          - handle @(h)->vector, ou vecteur 1..K
%  'nTopCand'   (30)     : nb de pics candidats (les plus forts)
%  'MinHarmHit' (2)      : nb min d’harmoniques trouvées pour accepter Li
%
% Out:
%  L_best       : lag fondamental estimé (échantillons)
%  score_best   : score du meilleur candidat
%  keptIdx_best : lags des pics retenus pour L_best (un max par h)
%  keptVal_best : amplitudes correspondantes
%  H_best       : numéros d’harmoniques associés à keptIdx_best
%  diag         : diagnostics (scores par candidat, mapping, etc.)
%
% Notes:
%  - Complexité ~ O(nTopCand * N). Très rapide pour N ~ 100–300.
%  - Si tu veux un poids uniforme, passe 'Weights', @(h) ones(size(h)).

% ---------- parsing options ----------
p = inputParser;
addParameter(p,'RelTol',0.06);
addParameter(p,'AbsTol',[]);
addParameter(p,'Weights',[]);
addParameter(p,'nTopCand',30);
addParameter(p,'MinHarmHit',2);
parse(p,varargin{:});
RelTol     = p.Results.RelTol;
AbsTol     = p.Results.AbsTol;
Weights    = p.Results.Weights;
nTopCand   = p.Results.nTopCand;
MinHarmHit = p.Results.MinHarmHit;

% ---------- préparation ----------
peakIdx = peakIdx(:);
peakVal = peakVal(:);
N = numel(peakIdx);
assert(N==numel(peakVal),'peakIdx/peakVal size mismatch');
assert(K>=1 && floor(K)==K,'K doit être un entier >=1');

% candidats = les nTopCand pics les plus forts
[~, ordByVal] = sort(peakVal, 'descend');
candOrd = ordByVal(1:min(nTopCand, N));
Li_list = double(peakIdx(candOrd));   % candidats fondamental Li

% pondération w(h)
h1K = (1:K).';
if isempty(Weights)
    w = 1 ./ sqrt(h1K);    % défaut doux
elseif isa(Weights,'function_handle')
    w = Weights(h1K);
else
    w = Weights(:);
    if numel(w) < K, error('Weights length < K'); end
    w = w(1:K);
end

% ---------- boucle candidats ----------
scores_all   = -inf(numel(Li_list),1);
hits_all     = zeros(numel(Li_list),1);
keptMapsCell = cell(numel(Li_list),1);

L_best       = NaN;
score_best   = -Inf;
keptIdx_best = [];
keptVal_best = [];
H_best       = [];

for ic = 1:numel(Li_list)
    Li = Li_list(ic);
    
    % tolérance
    if isempty(AbsTol)
        tol = RelTol * max(Li, eps);
    else
        tol = AbsTol;
    end
    
    % attribution harmonique h ≈ round(idx/Li)
    ratio = double(peakIdx) ./ Li;
    hCand = round(ratio);
    delta = abs(double(peakIdx) - hCand * Li);
    
    % garder seulement h dans [1..K] et delta <= tol
    valid = (hCand >= 1) & (hCand <= K) & (delta <= tol);
    if ~any(valid), continue; end
    
    vIdx = peakIdx(valid);
    vVal = peakVal(valid);
    vH   = hCand(valid);
    % (optionnel : vDel = delta(valid);)
    
    % un pic max par harmonique
    [uniqH, ~, grp] = unique(vH, 'stable');
    bestLag = zeros(numel(uniqH),1);
    bestVal = zeros(numel(uniqH),1);
    for u = 1:numel(uniqH)
        sel  = (grp == u);
        tVal = vVal(sel);
        tIdx = vIdx(sel);
        [bestVal(u), k] = max(tVal);
        bestLag(u)      = tIdx(k);
    end
    
    % score
    w_used = w(uniqH);
    score  = sum(w_used .* bestVal);
    hits   = numel(uniqH);
    
    scores_all(ic)   = score;
    hits_all(ic)     = hits;
    keptMapsCell{ic} = struct('H',uniqH,'Idx',bestLag,'Val',bestVal);
    
    if (hits >= MinHarmHit) && (score > score_best)
        score_best   = score;
        L_best       = Li;
        keptIdx_best = bestLag;
        keptVal_best = bestVal;
        H_best       = uniqH;
    end
end

% ---------- fallback si aucun candidat ne satisfait MinHarmHit ----------
if isnan(L_best)
    [score_best, ib] = max(scores_all);
    if isfinite(score_best)
        L_best = Li_list(ib);
        m      = keptMapsCell{ib};
        if ~isempty(m)
            keptIdx_best = m.Idx;
            keptVal_best = m.Val;
            H_best       = m.H;
        else
            keptIdx_best = [];
            keptVal_best = [];
            H_best       = [];
        end
    else
        % dernier repli : choisir le pic d'amplitude max
        [score_best, iMax] = max(peakVal);
        L_best       = double(peakIdx(iMax));
        keptIdx_best = peakIdx(iMax);
        keptVal_best = peakVal(iMax);
        H_best       = 1;
    end
end

% ---------- diagnostics ----------
diag = struct();
diag.candidates  = Li_list;
diag.scores      = scores_all;
diag.hits        = hits_all;
diag.maps        = keptMapsCell;
diag.params      = struct('K',K,'RelTol',RelTol,'AbsTol',AbsTol, ...
    'nTopCand',nTopCand,'MinHarmHit',MinHarmHit);
end

function [L_refined, Lh, wh] = refine_L_from_harmonics(acf, keptIdx, keptVal, H, huberK)
% keptIdx/keptVal/H: from your harmonic selector (one per harmonic)
% huberK: Huber tuning (in samples), e.g. 0.4–0.8; [] -> no robust step

acf = acf(:);
n = numel(acf);
epsDen = 1e-12;

% 1) Parabolic interpolation per harmonic
Lh = nan(numel(H),1);
wh = zeros(numel(H),1);
for t=1:numel(H)
    i = keptIdx(t);
    if i<=1 || i>=n-1, continue; end
    y1 = acf(i-1); y2 = acf(i); y3 = acf(i+1);
    den = (y1 - 2*y2 + y3);
    if abs(den) < 1e-9, continue; end
    delta = 0.5*(y1 - y3)/(den);     % sub-sample shift
    ell   = i + delta;               % refined lag at harmonic
    Lh(t) = ell / H(t);
    
    % 2) weight = amplitude * curvature proxy
    curv  = abs(den);
    edge  = abs(y1 - y3) + abs(y2);  % helps demote flat/ambiguous peaks
    w0    = max(keptVal(t),0);
    wh(t) = w0 * (curv / (edge + epsDen));
end

% remove invalids
v = isfinite(Lh) & (wh>0);
Lh = Lh(v); wh = wh(v);
if isempty(Lh)
    L_refined = NaN; return;
end

% 3) Robust fuse (weighted)
L0 = sum(wh.*Lh)/sum(wh);   % weighted mean init

if ~isempty(huberK) && isfinite(huberK) && huberK>0
    r = Lh - L0;
    a = abs(r);
    % Huber weights
    hub = ones(size(a));
    m = a>huberK;
    hub(m) = huberK ./ a(m);
    wR = wh .* hub;
    L_refined = sum(wR.*Lh)/sum(wR);
else
    % simple weighted mean
    L_refined = L0;
end
end

function [f_hat, diag] = fuse_f0_candidates_fast(fcands, acf_s, Fs, varargin)
% FUSE_F0_CANDIDATES_FAST
%  Fusionne des estimateurs de f0 {fcands} en 2 étapes légères:
%   (i) Regroupement par octave sur log2(f), (ii) médiane pondérée.
%  Option: mini anti-octave via R(T) vs R(T/2) (coût O(1)).
%
% In:
%   fcands : vecteur des f0 candidats (Hz), NaN/Inf/<=0 ignorés
%   acf_s  : ACF(lags>=0) normalisée (R(0)=1), [] si indisponible
%   Fs     : Hz
%
% Options (name/value):
%   'Weights'    : poids par candidat (par ex. courbure parab., NSDF score). def = []
%   'OctaveTol'  : tolérance sur log2 pour regrouper (def 0.20 -> ±20%)
%   'AntiHalf'   : true/false, active test ACF R(T) vs R(T/2) (def true)
%   'GammaHalf'  : seuil pour anti-half (def 1.05)
%
% Out:
%   f_hat : fréquence finale (Hz)
%   diag  : struct avec clusters, scores, etc.

p = inputParser;
addParameter(p,'Weights',[]);
addParameter(p,'OctaveTol',0.20);
addParameter(p,'AntiHalf',true);
addParameter(p,'GammaHalf',1.05);
parse(p,varargin{:});
w         = p.Results.Weights;
octTol    = p.Results.OctaveTol;
doAnti    = p.Results.AntiHalf;
gammaHalf = p.Results.GammaHalf;

% 0) Nettoyage
fc = fcands(:);
fc = fc(isfinite(fc) & fc>0);
if isempty(fc), f_hat = NaN; diag=[]; return; end
if isempty(w),  w = ones(size(fc)); else, w = w(:); w = w(1:numel(fc)); end

% 1) Consensus d’octave (clustering 1D sur log2)
lf = log2(fc);
% centre initial: médiane pondérée
lf0 = weightedMedian(lf, w);
% normalise chaque f par un décalage d’octave entier vers lf0
k = round(lf - lf0);                % entier d’octave
lfn = lf - k;                       % replie près de lf0 (même bande)
% clusterise par binning fin
[lfn_sorted, ord] = sort(lfn);
w_sorted = w(ord);
clusters = {}; scores = [];
start = 1;
for i = 2:numel(lfn_sorted)+1
    if i==numel(lfn_sorted)+1 || abs(lfn_sorted(i)-lfn_sorted(start))>octTol
        clusters{end+1} = ord(start:i-1); %#ok<AGROW>
        scores(end+1)   = sum(w_sorted(start:i-1)); %#ok<AGROW>
        start = i;
    end
end
[~, ibest] = max(scores);
idxBest = clusters{ibest};
% fréquence fusionnée (log-médiane pondérée)
lf_best = weightedMedian(lf(idxBest), w(idxBest));
f0 = 2.^lf_best;

% 2) Anti-octave ultra-léger (optionnel, si ACF dispo)
if doAnti && ~isempty(acf_s) && numel(acf_s)>=3
    T0 = 1/max(f0, eps);
    L0 = T0 * Fs;
    % interpolation linéaire très simple
    interpR = @(q) simple_interp((0:numel(acf_s)-1).', acf_s(:), q);
    RT = interpR(L0);
    RH = interpR(L0/2);
    if (RH > gammaHalf*RT)
        f0 = 2*f0;  % T/2 => double fréquence
    end
end

f_hat = f0;

% diag
diag = struct();
diag.fc_clean   = fc;
diag.weights    = w;
diag.lf         = lf;
diag.lf_ref     = lf0;
diag.lf_folded  = lfn;
diag.clusters   = clusters;
diag.clusterScores = scores;
diag.f_beforeAnti  = 2.^lf_best;

end

% --- utilitaires compacts ---
function m = weightedMedian(x, w)
[xs,ord] = sort(x(:)); ws = w(ord);
c = cumsum(ws)/sum(ws);
i = find(c>=0.5,1,'first');
m = xs(i);
end

function v = simple_interp(x, y, q)
% linéaire, hors-borne -> clamp
q = q(:);
q = max(min(q, x(end)), x(1));
i = floor(q); a = q - i;
i = max(1, min(numel(x)-1, i+1)); % 1-based
v = (1-a).*y(i) + a.*y(i+1);
end


function [L0, S] = pickL_from_acf_harmonics_win(AcfIn, Lref, Kharm, tauMaxSamp, Fs, fmax, gammaHalf, lambdaEven, windowSamp)
% S = struct('base',   Lref);
% L0 = Lref;
% return 
% PICKL_FROM_ACF_HARMONICS_WIN
%  Compare la somme des pics ACF aux harmoniques de Lref avec une recherche locale
%  dans une fenêtre autour de chaque harmonique, puis choisit Lref, Lref/2 ou 2*Lref.
%
% In:
%   AcfIn        : ACF (lags >= 0), vecteur ligne/colonne
%   Lref         : lag de référence (échantillons, flottant OK)
%   Kharm        : nb d'harmoniques à sommer (ex: 6)
%   tauMaxSamp   : lag max exploitable (<= numel(AcfIn)-1)
%   Fs           : Hz (utilisé pour estimer la fenêtre si windowSamp vide)
%   fmax         : Hz (fréq maxi attendue; sert à fixer la fenêtre si windowSamp vide)
%   gammaHalf    : pénalité T/2 (ex: 1.05)  => sumHalf = sumHalf / gammaHalf
%   lambdaEven   : pénalité 2T  (ex: 1.10)  => sumDouble = sumDouble / lambdaEven
%   windowSamp   : (optionnel) demi-largeur fenêtre en échantillons (par défaut ~Fs/(2*fmax))
%
% Out:
%   L0 : lag choisi (Lref, Lref/2 ou 2*Lref)
%   S  : struct diagnostics:
%        .base/.half/.double  (sommes après pénalités)
%        .base_raw/.half_raw/.double_raw (avant pénalités)
%        .pos_base/.pos_half/.pos_double (indices retenus par harmonique)
%        .win = windowSamp

    AcfIn = AcfIn(:);
    nLag  = numel(AcfIn) - 1;
    tauMaxSamp = min(max(1, floor(tauMaxSamp)), nLag);

    % Définir la fenêtre si non fournie
%     if nargin < 9 || isempty(windowSamp)
%         if nargin >= 6 && ~isempty(Fs) && ~isempty(fmax) && isfinite(Fs) && isfinite(fmax) && fmax>0
%             windowSamp = max(1, round(Fs/(2*fmax)));  % ~ ce que tu faisais
%         else
            windowSamp = max(1, round(0.05 * Lref));  % fallback ~5% de Lref
%         end
%     end
    windowSamp = max(1, round(windowSamp));

    % Helper: somme des maxima locaux autour des harmoniques d'un L
    function [s, pos_list] = sum_harmonics_with_window(Lcand)
        s = 0;
        pos_list = nan(Kharm,1);
        for h = 1:Kharm
            center = round(h * Lcand);
            if center < 1, continue; end
            if center > tauMaxSamp, break; end
            a = max(1, center - windowSamp);
            b = min(tauMaxSamp, center + windowSamp);
            % pic local dans la fenêtre
            [~,rel] = max(AcfIn(a:b));
            pos = a + rel - 1;
            pos_list(h) = pos;
            % on peut clamp à >=0 si souhaité (décommenter)
            % s = s + max(AcfIn(pos), 0);
            s = s + AcfIn(pos);
        end
    end

    % Candidats : L, L/2, 2L
    [sumB_raw, posB] = sum_harmonics_with_window(Lref);
    [sumH_raw, posH] = sum_harmonics_with_window(Lref/2);
    [sumD_raw, posD] = sum_harmonics_with_window(Lref*2);

    % Pénalités (comme ton code)
    sumB = sumB_raw;
    sumH = sumH_raw / max(gammaHalf, eps);
    sumD = sumD_raw / max(lambdaEven, eps);

    % Choix
    [~, k] = max([sumB, sumH, sumD]);
    if k == 2
        L0 = Lref/2;
    elseif k == 3
        L0 = Lref*2;
    else
        L0 = Lref;
    end

    % Diagnostics
    S = struct( ...
        'base',   sumB, 'half',   sumH, 'double',   sumD, ...
        'base_raw',sumB_raw,'half_raw',sumH_raw,'double_raw',sumD_raw, ...
        'pos_base',posB,'pos_half',posH,'pos_double',posD, ...
        'win',windowSamp);

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

function [M, Fs2, alpha] = choose_decimation_safe(Fs, fmax, varargin)
% CHOOSE_DECIMATION_SAFE — borne M, choisit un facteur dans un set autorisé
%  Fs2 = Fs/M >= alpha*fmax, avec M <= Mmax et M ∈ allowedMs.

p = inputParser;
addParameter(p,'Alpha',10);             % échant./période visés
addParameter(p,'Mmax',8);               % borne supérieure raisonnable
addParameter(p,'Allowed',[2 3 4 5 6 8 10 12 14 16]);% facteurs autorisés
parse(p,varargin{:});
alpha     = p.Results.Alpha;
Mmax      = p.Results.Mmax;
allowedMs = sort(unique(p.Results.Allowed));

if ~isfinite(fmax) || fmax<=0
    M = 1; Fs2 = Fs; return;
end

% M cible naïf
M_target = floor(Fs / (alpha*fmax));
M_target = max(1, M_target);

% borne supérieure
M_cap = min(M_target, Mmax);

% choisir le plus grand M autorisé <= M_cap
M_opts = allowedMs(allowedMs <= M_cap);
if isempty(M_opts)
    M = 1;
else
    M = M_opts(end);
end

Fs2 = Fs / M;

% garde-fou: si Fs2 < alpha*fmax (rare avec floor), on réduit encore M
while Fs2 < alpha*fmax && M > 1
    % prendre l'option autorisée suivante plus petite
    prev = allowedMs(allowedMs < M);
    if isempty(prev), M = 1; break; end
    M = prev(end);
    Fs2 = Fs / M;
end
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


function [f_pick, S, details] = decide_T_vs_half(x, Fs, Lcand, r_acf, K)
% DECIDE_T_VS_HALF  -- disambiguation f0 vs 2f0 (T vs T/2)
% In:
%   x      : signal (col ou ligne)
%   Fs     : Hz
%   Lcand  : lag candidat T (en échantillons, flottant possible)
%   r_acf  : ACF(lags>=0), normalisée: r(0)=1 (vector)
%   K      : nb d’harmoniques pour peigne ACF (ex 6)
% Out:
%   f_pick : fréquence choisie (Hz): soit Fs/Lcand, soit 2*Fs/Lcand
%   S      : struct avec scores/votes
%   details: struct diag

x = x(:); r_acf = r_acf(:);
N  = numel(x);
tauMax = numel(r_acf)-1;

% --- utilitaires ---
lags = (0:tauMax).';
interpR = @(q) simple_interp(lags, r_acf, q);
function v = simple_interp(xi, yi, q)
    q = max(min(q, xi(end)), xi(1));
    i = floor(q); a = q - i;
    i = max(1, min(numel(xi)-1, i+1));
    v = (1-a).*yi(i) + a.*yi(i+1);
end

% 0) Pré-calculs
LT   = max(Lcand, 1);     % T (samples)
LH   = LT/2;              % T/2 (samples)
f_T  = Fs / LT;
f_H  = 2 * f_T;

% 1) ACF signée: R(T) vs R(T/2) (avec signe)
RT = interpR(LT);
RH = interpR(LH);
score_acf = RT - abs(RH);      % si RH (en valeur absolue) rivalise, pénalise T
vote_T_acf = (RT > abs(RH));
vote_H_acf = (abs(RH) > RT && RH > 0);  % 2f0 plausible seulement si RH positif

% 2) Peigne odd/even sur ACF
if nargin<5 || isempty(K), K = 6; end
tauMax  = numel(r_acf)-1;
S_T  = comb_odd_even_acf(LT, K, 1.0, interpR, tauMax);     % grand => impairs > pairs -> en faveur de T
S_H  = comb_odd_even_acf(LH, K, 1.0, interpR, tauMax);    % grand => impairs > pairs à T/2 (peu probable si 2f0 réel)
score_comb = S_T - S_H;
vote_T_comb = (S_T > S_H);
vote_H_comb = (S_H > S_T && RH>0);

% 3) Half-wave flip dans le temps: corr(x, -x décalé de T/2)
rho_flip = halfwave_flip_corr(x, LH);
% si rho_flip est grand et positif -> T (car x(t+T/2) ≈ -x(t))
vote_T_flip = (rho_flip > 0.25);  % seuil doux
vote_H_flip = (rho_flip < 0.00 && RH>0); % si corr≈0 et RH positif, ça n’aide pas T

% 4) Ancrage cepstre (optionnel): si tu as Tcep, ajoute un petit vote
vote_T_cep = false; vote_H_cep = false; score_cep = 0;
% (tu peux passer Tcep via details et activer un vote faible)

% ---- Pondération & décision ----
w_acf  = 0.45;    % ACF signée
w_comb = 0.20;    % peigne odd/even
w_flip = 0.35;    % temps (flip)
score_T = w_acf*(vote_T_acf) + w_comb*(vote_T_comb) + w_flip*(vote_T_flip);
score_H = w_acf*(vote_H_acf) + w_comb*(vote_H_comb) + w_flip*(vote_H_flip);

if score_H > score_T
    f_pick = f_H;     % T/2 -> 2*f0
else
    f_pick = f_T;     % garder T
end

% sorties diag
S = struct('score_T',score_T,'score_H',score_H, ...
           'score_acf',score_acf,'score_comb',score_comb,'rho_flip',rho_flip, ...
           'RT',RT,'RH',RH);
details = struct('LT',LT,'LH',LH,'f_T',f_T,'f_H',f_H,'S_T',S_T,'S_H',S_H);
end

% function S = comb_odd_even_acf(Ts, K, lambdaEven, interpR, tauMax)
%     Sodd=0; Seven=0;
%     for h=1:K
%         q = h*Ts;
%         if q>=1 && q<=tauMax
%             r = interpR(q);       % <--- ici on utilise ton interp linear
%             r = max(r,0);         % on évite de récompenser les valeurs négatives
%             if mod(h,2)==1
%                 Sodd = Sodd + r;
%             else
%                 Seven = Seven + r;
%             end
%         end
%     end
%     S = (Sodd - lambdaEven*Seven) / max(Sodd + lambdaEven*Seven, eps);
% end

function rho = halfwave_flip_corr(x, LH)
    % corrélation normalisée entre x[n] et -x[n+LH]
    % LH peut être non entier -> on interpole
    N = numel(x);
    nmax = floor(N - LH - 1);
    if nmax < 10, rho = 0; return; end
    n = (1:nmax).';
    x1 = x(n);
    % interp linéaire de x à demi-période
    idx = n + LH;
    i0 = floor(idx); a = idx - i0;
    i0 = max(1, min(N-1, i0));
    x2 = (1-a).*x(i0) + a.*x(i0+1);
    num = sum( x1 .* (-x2) );
    den = sqrt(sum(x1.^2) * sum(x2.^2)) + eps;
    rho = num / den;   % ~ +1 si T (demi-ondes opposées), ~0 ou négatif si T/2
end

function [f_pick, S, details] = decide_T_T2_Tover2(x, Fs, Lcand, r_acf, K, varargin)
% Choix parmi {T/2, T, 2T} avec flip robuste + biais anti-2T
p = inputParser;
addParameter(p,'GammaLo',1.10);    % pénalise R(L/2) si positif
addParameter(p,'GammaHi',0.60);    % pénalise R(2L)
addParameter(p,'LambdaEven',1.00); % peigne
addParameter(p,'Wacf',0.45);
addParameter(p,'Wcomb',0.25);
addParameter(p,'Wflip',0.30);
addParameter(p,'TieMargin',0.05);  % marge "quasi-égalité"
addParameter(p,'MinFlipSamples',64); % min points pour calculer flip
parse(p,varargin{:});
gLo=p.Results.GammaLo; gHi=p.Results.GammaHi;
lambdaEven=p.Results.LambdaEven;
w_acf=p.Results.Wacf; w_comb=p.Results.Wcomb; w_flip=p.Results.Wflip;
tieMargin=p.Results.TieMargin; minFlipN=p.Results.MinFlipSamples;

x = x(:);
x = x - mean(x); s=std(x); if s>0, x=x/s; end           % <-- important
r_acf = r_acf(:);
tauMax = numel(r_acf)-1; lags=(0:tauMax).';
interpR = @(q) safe_interp(lags, r_acf, q);

if nargin<5 || isempty(K), K=6; end
L_H = max(Lcand/2,1);  L_T = max(Lcand,1);  L_2 = max(2*Lcand,1);
L_all = [L_H, L_T, L_2]; names = {'T/2','T','2T'};

S_acf=zeros(1,3); S_comb=zeros(1,3); S_flip=zeros(1,3);

for ic=1:3
    Lc=L_all(ic); if Lc>tauMax, continue; end
    % ---- ACF: favorise un "vrai" pic isolé et pénalise voisins d'octave
    R_c  = interpR(Lc);
    R_lo = interpR(max(Lc/2,1));
    R_hi = interpR(min(2*Lc,tauMax));
    % combinaison différence + ratio (plus stable)
    diffTerm = R_c - gLo*max(R_lo,0) - gHi*abs(R_hi);
    ratioTerm = log( (R_c+eps) / (max(abs([R_lo R_hi]))+eps) );
    S_acf(ic) = 0.6*diffTerm + 0.4*ratioTerm;

    % ---- Peigne odd/even ACF (valeurs négatives clampées à 0)
    S_comb(ic) = comb_odd_even_acf(Lc, K, lambdaEven, interpR, tauMax);

    % ---- Flip: corr(x, x décalé de Lc/2) et corr(x, -x décalé de Lc/2)
    Lh = max(Lc/2,1);
    [rho_pos, rho_neg, okFlip] = flip_corr_pos_neg(x, Lh, minFlipN);
    if okFlip
        % Score >0 attendu pour T (demi-ondes opposées -> rho_neg négatif fort)
        % et <0 pour 2T (rho_pos positif fort).
        S_flip(ic) = (-rho_neg) - (rho_pos);
    else
        S_flip(ic) = 0;   % on neutralise le flip si pas assez d'échantillons
    end
end

% ---- Score global
S_tot = w_acf*S_acf + w_comb*S_comb + w_flip*S_flip;

% ---- Choix + tie-break "plus petit L si quasi-égalité"
[best, ix] = max(S_tot);
[~, secondIdx] = max(S_tot .* (1 - (1:3==ix)));
if (best - S_tot(secondIdx)) < tieMargin
    % petite préférence structurée: L plus petit gagne en cas d'ex-aequo
    [~, ix] = min(L_all);
end

% fréquence finale
f_pick = Fs / L_all(ix);

% sorties
S = struct('names',{names},'L',L_all,'S_acf',S_acf,'S_comb',S_comb,'S_flip',S_flip,'S_tot',S_tot);
details = struct('RT',interpR(L_T),'RH',interpR(max(L_T/2,1)),'R2',interpR(min(2*L_T,tauMax)));
end

% ---------- helpers ----------
function v = safe_interp(x, y, q)
    q = min(max(q, x(1)), x(end));
    i = floor(q); a = q - i;
    i = max(1, min(numel(x)-1, i+1));
    v = (1-a).*y(i) + a.*y(i+1);
end

function S = comb_odd_even_acf(Ts, K, lambdaEven, interpR, tauMax)
    Sodd=0; Seven=0;
    for h=1:K
        q=h*Ts;
        if q>=1 && q<=tauMax
            r=interpR(q); r=max(r,0);
            if mod(h,2)==1, Sodd=Sodd+r; else, Seven=Seven+r; end
        end
    end
    denom = max(Sodd + lambdaEven*Seven, eps);
    S = (Sodd - lambdaEven*Seven) / denom;
end

function [rho_pos, rho_neg, ok] = flip_corr_pos_neg(x, Lh, minN)
    N = numel(x);
    nmax = floor(N - Lh - 1);
    if nmax < minN, rho_pos=0; rho_neg=0; ok=false; return; end
    n=(1:nmax).';
    x1=x(n);
    idx=n+Lh; i0=floor(idx); a=idx-i0; i0=max(1,min(N-1,i0));
    x2=(1-a).*x(i0)+a.*x(i0+1);
    den = sqrt(sum(x1.^2)*sum(x2.^2)) + eps;
    rho_pos = sum(x1.*x2)/den;       % corr avec +x
    rho_neg = sum(x1.*(-x2))/den;    % corr avec -x
    ok=true;
end

function score = srh_score(x, Fs, f0, K, varargin)
% SRH_SCORE  — somme énergie aux harmoniques − alpha*énergie inter-harmoniques
% Options: 'Alpha',1.0, 'BWrel',0.03, 'Nfft',[], 'DoLPC',true, 'UseLogMag',false
p = inputParser;
addParameter(p,'Alpha',1.0);
addParameter(p,'BWrel',0.03);   % demi-bande relative: bwHz = max(BWrel*fh, 2*Fbin)
addParameter(p,'Nfft',[]);
addParameter(p,'DoLPC',true);
addParameter(p,'UseLogMag',false);
parse(p,varargin{:});
alpha = p.Results.Alpha; BWrel = p.Results.BWrel; Nfft = p.Results.Nfft;
DoLPC = p.Results.DoLPC; UseLogMag = p.Results.UseLogMag;

if ~isfinite(f0) || f0<=0, score = -Inf; return; end

x = x(:); x = x - mean(x); s = std(x); if s>0, x = x/s; end

% 1) pré-blanchiment (optionnel)
if DoLPC
    pL = 10;
    a  = lpc(x, pL);
    r  = filter(a, 1, x);           % inverse filtering: approx. excitation
else
    r  = x;
end

% 2) spectre 1-côté
if isempty(Nfft), Nfft = 2^nextpow2(numel(r)*2); end
w   = hann(numel(r),'periodic');
R   = abs(fft(w.*r, Nfft));
R   = R(1:floor(Nfft/2)+1);         % demi-spectre
F   = (0:floor(Nfft/2))' * (Fs/Nfft);
Fbin= Fs/Nfft;

if UseLogMag
    R = log(R + eps);               % option: compression dynamique
end

% 3) SRH
Eh = 0; Ei = 0;
for h = 1:K
    fh = h*f0;          if fh >= Fs/2, break; end
    fi = (h+0.5)*f0;    if fi >= Fs/2, fi = NaN; end

    % demi-bande en Hz: proportionnelle à fh, jamais < 2 bins
    bwHz = max(BWrel*fh, 2*Fbin);

    % fenêtres d’intégration
    idxh = abs(F - fh) <= bwHz;
    Eh   = Eh + sum(R(idxh));

    if ~isnan(fi)
        idxi = abs(F - fi) <= bwHz;
        Ei   = Ei + sum(R(idxi));
    end
end

score = Eh - alpha*Ei;              % --> A MAXIMISER
end


% function score = srh_score(x, Fs, f0, K)
% % 1) pré-blanchiment (résiduel LPC court) pour annuler l’enveloppe
% p = 10;
% a = lpc(x, p);                  % coefficients LPC
% r = filter(a, 1, x);            % résiduel (approx “flat”)
% % 2) spectre
% N = 2^nextpow2(numel(r));
% R = abs(fft(hann(numel(r),'periodic').*r, N));
% F = (0:N-1)*(Fs/N);
% % 3) somme harmonique - inter-harmonique
% bw = max(1, round(0.015*Fs));   % demi-fenêtre en bins (~1.5% Fs) à ajuster
% score = 0;
% for h = 1:K
%     fh = h*f0;  
%     if fh>Fs/2, 
%         break; 
%     end
%     fi = (h+0.5)*f0;            % inter-harmonique
%     % fenêtres d’intégration
%     idxh = abs(F - fh) <= (bw*Fs/N);
%     idxi = abs(F - fi) <= (bw*Fs/N);
%     score = score + sum(R(idxh)) - sum(R(idxi));
% end
% end


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

function [newVal, isCorrection, scores, multGrid] = searchSRH(acc, Fs, candf, Kharm_KSRH, varargin)
% SEARCHSRH  — choisit entre {f0/8, f0/4, f0/2, f0, 2f0, 4f0, 8f0} via SRH
% Options: 'Alpha',1.0, 'BWrel',0.03, 'Nfft',[], 'DoLPC',true
p = inputParser;
addParameter(p,'Alpha',1.0);
addParameter(p,'BWrel',0.03);      % demi-bande relative autour de chaque (inter)harmonique (ex: 3%)
addParameter(p,'Nfft',[]);
addParameter(p,'DoLPC',true);
parse(p,varargin{:});
alpha = p.Results.Alpha; BWrel = p.Results.BWrel; Nfft = p.Results.Nfft; DoLPC = p.Results.DoLPC;

multGrid = [1/8, 1/4, 1/2, 1, 2, 4, 8];
scores   = zeros(numel(multGrid),1);

for i = 1:numel(multGrid)
    ftest = candf * multGrid(i);
    scores(i) = srh_score(acc, Fs, ftest, Kharm_KSRH, ...
                          'Alpha',alpha, 'BWrel',BWrel, 'Nfft',Nfft, 'DoLPC',DoLPC);
end

% SRH = score à MAXIMISER
[~, idxMax] = max(scores);
newVal       = candf * multGrid(idxMax);
isCorrection = (multGrid(idxMax) ~= 1);
end

% function [newVal , iscorrection] = searchSRH(acc, Fs, candf, Kharm_KSRH)
% scores=zeros(8,1);
% divArray = [1/8,1/4,1/2,1,2,4,8];
% for i = 1:numel(divArray)
%     scores(i) = srh_score(acc, Fs, candf*divArray(i), Kharm_KSRH);
% end
% [~,idxMin] = min(scores);
% iscorrection = false;
% if divArray(idxMin) ~= 1
%     iscorrection = true;
% end
% newVal = candf*divArray(idxMin);
% end

function [fhat, out] = cand_hps_from_psd_robust(r, Fs, fmin, fmax, dfHz, K, varargin)
% PSD via ACF (Wiener–Khinchin) + Harmonic Sum robuste (anti demi-tour).
% Choisit f0 en maximisant un score pénalisé par l'énergie à f0/2 et 2f0,
% et en privilégiant les peignes où les harmoniques impairs dominent.
%
% Name-Value (optionnel)
%   'WeightMode'   : '1/k' (def) ou 'equal'
%   'LambdaEven'   : 0.8   (pondération des harmonique pairs dans le contraste)
%   'AlphaHalf'    : 0.8   (poids de la pénalité liée à f0/2)
%   'AlphaDouble'  : 0.6   (poids de la pénalité liée à 2*f0)
%   'MinHarmonics' : 3     (au moins 3 harmoniques dans la bande)
%   'UseLog'       : true  (somme des log(ε+P) pour plus de robustesse)
%   'Eps'          : 1e-3  (ε pour les logs)
%   'Parabolic'    : true  (raffinement sub-bin sur la grille f0)
%
% Out:
%   fhat : Hz
%   out  : struct diag (F, Pxx, fgrid, S, Spen, combContrast, idx, ...)

% ---- options
p = inputParser;
p.addParameter('WeightMode','1/k');
p.addParameter('LambdaEven',0.8);
p.addParameter('AlphaHalf',0.8);
p.addParameter('AlphaDouble',0.6);
p.addParameter('MinHarmonics',3);
p.addParameter('UseLog',true);
p.addParameter('Eps',1e-3);
p.addParameter('Parabolic',true);
p.parse(varargin{:});
P = p.Results;

% ---- PSD via ACF
[F, Pxx] = psd_from_acf(r, Fs);     % suppose r déjà normalisée ou pas, OK pour scoring
% coupe au-dessus de fmax
Pxx(F > fmax+1) = 0;

% ---- grille f0
fgrid = (fmin:dfHz:fmax).';
S     = -inf(numel(fgrid),1);
Sodd  = zeros(size(S));  Seven = zeros(size(S));
cntH  = zeros(size(S));
Spec  = @(f) interp1(F, Pxx, f, 'linear', 0);

% poids harmoniques
Kmax = K;
switch lower(string(P.WeightMode))
    case "equal", wk = ones(1,Kmax);
    otherwise,    wk = 1./(1:Kmax);   % défaut : 1/k
end
wk = wk / sum(wk);

% ---- score de peigne + séparation impairs/pairs
for i=1:numel(fgrid)
    f0 = fgrid(i);
    s = 0; sodd = 0; seven = 0; used = 0;
    for k=1:Kmax
        fk = k*f0;
        if fk > F(end), break; end
        val = Spec(fk);
        if P.UseLog, val = log(P.Eps + val); end
        s = s + wk(k)*val;
        if mod(k,2)==1, sodd = sodd + wk(k)*val; else, seven = seven + wk(k)*val; end
        used = used + 1;
    end
    if used >= P.MinHarmonics
        S(i)    = s;
        Sodd(i) = sodd;
        Seven(i)= seven;
        cntH(i) = used;
    end
end

% ---- pénalités f/2 et 2f (anti sous-/sur-harmonique)
% on recalcule un score HPS aux fréquences voisines f/2 et 2f avec mêmes poids
S_half   = local_hps_at(fgrid/2,   F, Pxx, wk, Kmax, P.UseLog, P.Eps);
S_double = local_hps_at(2*fgrid,   F, Pxx, wk, Kmax, P.UseLog, P.Eps);

% contraste odd/even (sur le peigne du candidat)
combContrast = (Sodd - P.LambdaEven*Seven) ./ max(Sodd + P.LambdaEven*Seven, eps);

% score pénalisé
Spen = S - P.AlphaHalf  * max(S_half,  0) ...
         - P.AlphaDouble* max(S_double,0) ...
         + 0.10 * combContrast;    % petit bonus si impairs > pairs

% garde-fous: invalide si cntH trop petit
Spen(cntH < P.MinHarmonics) = -inf;

% ---- choix du candidat
[~, ix] = max(Spen);
f0 = fgrid(ix);

% ---- raffinement parabolique optionnel sur Spen
if P.Parabolic && ix>1 && ix<numel(Spen)
    y1=Spen(ix-1); y2=Spen(ix); y3=Spen(ix+1);
    denom = (y1 - 2*y2 + y3);
    if isfinite(denom) && denom < 0
        delta = 0.5*(y1 - y3) / denom;
        delta = max(min(delta, 0.5), -0.5);
    else
        delta = 0;
    end
    f0 = fgrid(ix) + delta*dfHz;
end

fhat = f0;

% ---- diagnostics
out = struct();
out.F = F; out.Pxx = Pxx;
out.fgrid = fgrid; out.S = S; out.Spen = Spen;
out.S_half = S_half; out.S_double = S_double;
out.combContrast = combContrast;
out.idx = ix; out.f_cand = fgrid(ix);

end

% ===== helpers =====
function Sg = local_hps_at(fvec, F, Pxx, wk, K, useLog, epslog)
Spec = @(f) interp1(F, Pxx, f, 'linear', 0);
Sg = -inf(numel(fvec),1);
for i=1:numel(fvec)
    f0 = fvec(i);
    if f0 <= 0 || f0 > F(end), Sg(i) = -inf; continue; end
    s=0; used=0;
    for k=1:K
        fk = k*f0;
        if fk > F(end), break; end
        v = Spec(fk);
        if useLog, v = log(epslog + v); end
        s = s + wk(k)*v;
        used = used + 1;
    end
    if used>=2, Sg(i)=s; end
end
end

function BEST = auto_select_band_for_envelope(x, Fs, opts)
% Auto-sélection de bande pour extraction d’enveloppe (AM autour d’une porteuse).
% 1) trouve la porteuse f_c et l’offset Δ (sidebands) dans le spectre
% 2) scanne des bandes autour de f_c pour maximiser un score de périodicité
%
% Retourne BEST avec :
%   .fc, .Delta, .band = [f1 f2], .score
%   .fr_est (Hz), .envFs, .diag (structs)

if nargin<3, opts = struct; end
% ---- défauts
DEF = struct('searchRange',[200 5000], 'nfft', 2^15, ...
             'DeltaScan',[5 150], 'DeltaStep', 0.5, ...
             'Kmax', 3, 'Bmult', [1.2 2 3 4], ...
             'envFs', 500, 'lpMult', 2.8, ...
             'yinThresh', 0.15, 'Kharm', 5, 'lambdaEven', 0.8, 'gammaHalf', 1.12);
fn=fieldnames(DEF); for k=1:numel(fn), if ~isfield(opts,fn{k}), opts.(fn{k})=DEF.(fn{k}); end, end

x = x(:) - mean(x(:));

% ===== 1) Spectre & détection porteuse =====
nfft = opts.nfft;
X = fft(x, nfft);
F = (0:nfft-1).' * (Fs/nfft);
P = abs(X).^2; P = P(1:floor(nfft/2)+1); F = F(1:numel(P));

% limite de recherche
m = (F>=opts.searchRange(1) & F<=opts.searchRange(2));
Fsr = F(m); Psr = P(m);

% pic principal
[~,iMax] = max(Psr);
fc = Fsr(iMax);

% ===== 2) Scan des Δ par symétrie de sidebands =====
DeltaVals = opts.DeltaScan(1):opts.DeltaStep:opts.DeltaScan(2);
Spec = @(f) interp1(F, P, f, 'linear', 0);
scoreDelta = zeros(numel(DeltaVals),1);
for i=1:numel(DeltaVals)
    d = DeltaVals(i);
    s = 0;
    for k=1:opts.Kmax
        s = s + min(Spec(fc-k*d), Spec(fc+k*d));  % symétrie
    end
    scoreDelta(i) = s;
end
[~,iBest] = max(scoreDelta);
Delta = DeltaVals(iBest);
end

function BEST = pick_best_carrier(x, Fs, Opt)
% PICK_BEST_CARRIER  Trouve la "meilleure" porteuse (résonance) dans le spectre.
% Combine proéminence, facteur Q et symétrie de sidebands (AM).
%
% In:
%   x   : signal temporel (vecteur)
%   Fs  : Hz
%   Opt : struct optionnel
%     .rangeHz      = [200 5000]   % bande de recherche des porteuses
%     .nfft         = 2^15
%     .winSec       = 1.0          % fenêtre Welch
%     .ovlp         = 0.5          % recouvrement Welch
%     .DeltaScan    = [5 150]      % plage Δ (Hz) à tester
%     .DeltaStep    = 0.5
%     .Kharm        = 3            % nb d’ordres de sidebands
%     .minProm      = 6            % dB, proéminence min des pics candidats
%     .maxPeaks     = 20           % nb max de pics à évaluer
%     .wProm        = 0.35         % poids proéminence
%     .wQ           = 0.25         % poids Q
%     .wSB          = 0.40         % poids sidebands
%     .useLogSpec   = true         % travailler sur log10(P) pour robustesse
%
% Out:
%   BEST struct:
%     .fc           : Hz (porteuse)
%     .Delta        : Hz (espacement sidebands estimé)
%     .band         : [f1 f2] Hz bande utile autour de fc (≈ ±3Δ)
%     .score        : score total normalisé
%     .peakProm_dB  : proéminence (dB)
%     .Q            : facteur Q estimé
%     .SBscore      : score sidebands
%     .F, .P        : spectre (pour plots)
%     .diag         : diagnostics (candidats, scores par Δ, etc.)

if nargin<3, Opt = struct; end
DEF = struct('rangeHz',[200 5000],'nfft',2^15,'winSec',1.0,'ovlp',0.5, ...
             'DeltaScan',[1 120],'DeltaStep',0.5,'Kharm',3,'minProm',6, ...
             'maxPeaks',20,'wProm',0.35,'wQ',0.25,'wSB',0.40,'useLogSpec',true);
fn=fieldnames(DEF); for k=1:numel(fn), if ~isfield(Opt,fn{k}), Opt.(fn{k})=DEF.(fn{k}); end, end

x = x(:) - mean(x(:));

% ===== 1) PSD (Welch) =====
Lwin = max(128, round(min(Opt.winSec*Fs,numel(x))));
Lhop = max(1, round(Lwin*(1-Opt.ovlp)));
win  = hann(Lwin,'periodic');
[Pxx, F] = pwelch(x, win, Lwin-Lhop, Opt.nfft, Fs, 'onesided');

% Nettoyage bande
m = (F>=Opt.rangeHz(1) & F<=Opt.rangeHz(2));
Fsr = F(m); Psr = Pxx(m);
if Opt.useLogSpec
    S = 10*log10(Psr + eps);   % dB
else
    S = Psr;
end

% ===== 2) Pics candidats (proéminence) =====
[pks, locs, w, prom] = findpeaks(S, Fsr, 'MinPeakProminence', Opt.minProm);
% Garder maxPeaks plus proéminents
if numel(pks) > Opt.maxPeaks
    [~,ord] = maxk(prom, Opt.maxPeaks);
    pks  = pks(ord); locs = locs(ord); w = w(ord); prom = prom(ord);
end

% Garde-fou si aucun pic:
if isempty(locs)
    BEST = struct('fc',NaN,'Delta',NaN,'band',[NaN NaN],'score',0, ...
                  'peakProm_dB',NaN,'Q',NaN,'SBscore',0,'F',F,'P',Pxx,'diag',struct());
    return;
end

% ===== 3) Score sidebands pour chaque candidat =====
DeltaVals = Opt.DeltaScan(1):Opt.DeltaStep:Opt.DeltaScan(2);
Spec = @(f) interp1(F, (Opt.useLogSpec*(10*log10(Pxx+eps)) + ~Opt.useLogSpec*Pxx), f, 'linear', -Inf);

nC = numel(locs);
SB   = zeros(nC,1);
Qfac = zeros(nC,1);
PromN = zeros(nC,1);
bestDelta = zeros(nC,1);

for i=1:nC
    fc = locs(i);

    % --- facteur Q : largeur -3dB locale (sur S en dB)
    Qfac(i) = local_Q_estimate(Fsr, S, fc);

    % --- proéminence normalisée (vers 0..1)
    PromN(i) = 1 - 10.^(-prom(i)/20);  % ~ 0 si petite prom, -> 1 si grande

    % --- scan Δ : symétrie min() des paires ±kΔ
    sDelta = zeros(numel(DeltaVals),1);
    for j=1:numel(DeltaVals)
        d = DeltaVals(j); s = 0; used=0;
        for k=1:Opt.Kharm
            fL = fc - k*d; fH = fc + k*d;
            if fL < F(2) || fH > F(end), break; end
            vL = Spec(fL); vR = Spec(fH);
            s  = s + min(vL, vR); used = used + 1;
        end
        % on exige au moins 2 paires
        if used>=2, sDelta(j)=s; else, sDelta(j)=-Inf; end
    end
    [sBest, jBest] = max(sDelta);
    SB(i) = sBest;
    bestDelta(i) = DeltaVals(jBest);
end

% Normaliser SB et Q sur [0,1] (robuste)
SBn = robust_unit(SB);
Qn  = robust_unit(Qfac);

% Score total et meilleur
Score = Opt.wProm*PromN(:) + Opt.wQ*Qn(:) + Opt.wSB*SBn(:);
[~, ibest] = max(Score);

% ===== 4) Sortie =====
fc  = locs(ibest);
dlt = bestDelta(ibest);
band = [max(10, fc - 3*dlt), min(Fs/2-50, fc + 3*dlt)];  % ±3Δ par défaut

BEST = struct();
BEST.fc           = fc;
BEST.Delta        = dlt;
BEST.band         = band;
BEST.score        = Score(ibest);
BEST.peakProm_dB  = prom(ibest);
BEST.Q            = Qfac(ibest);
BEST.SBscore      = SB(ibest);
BEST.F            = F; 
BEST.P            = Pxx;
BEST.diag = struct('candidates',locs,'prom_dB',prom,'Q',Qfac, ...
                   'SB',SB,'Score',Score,'DeltaVals',DeltaVals);

end

% ===== helpers =====

function Q = local_Q_estimate(F, Sdb, fc)
% est. Q autour de fc en cherchant -3 dB
[~,i0] = min(abs(F - fc));
peak = Sdb(i0);
th = peak - 3;
% gauche
iL = i0;
while iL>1 && Sdb(iL) > th, iL=iL-1; end
fL = F(iL);
% droite
iR = i0;
while iR<numel(F) && Sdb(iR) > th, iR=iR+1; end
fR = F(iR);
BW = max(1e-6, fR - fL);
Q  = fc / BW;
% clamp raisonnable
Q = max(0, min(Q, 1e4));
end

function y = robust_unit(x)
% normalisation robuste [0..1] en coupant les extrêmes (5e–95e pct)
x = x(:);
lo = prctile(x(isfinite(x)),5);
hi = prctile(x(isfinite(x)),95);
if hi<=lo, y = 0*x; return; end
y = (x - lo) / (hi - lo);
y = max(0, min(1, y));
end

