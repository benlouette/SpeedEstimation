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
array   =   array2A;
freqRot = freqRot3 ;
ISPLOTING = false;

sizeSlide=8192*2;
n=floor(numel(array)/sizeSlide);

for slide_idx = 1:n
    
    accOrigi=array((sizeSlide*(slide_idx-1))+1:(sizeSlide*slide_idx));
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
    if(ISPLOTING)
        plot_pick_best_comb(res, Opt, 10);
    end
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
    if ISPLOTING
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
    
    [~,idxBestScore] = max(Score);
    
    accFir = S_fft(:,idxBestScore);
    
    signalEnvelope = getEnv(accFir, Fs, info_fft.fc_Hz(idxBestScore));% res.best.fc BEST.band(1)-50);  %   extracts the envelope
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
    
    Opt = struct('SmoothMs',4.0,'PromRel',0.2,'Parabolic',true,'AntiHalf',false);
    [OUT] = period_from_acf_median_robust(AcfIn, Fs, tauMinSamp, tauMaxSamp, Opt);
    candf_mediane = OUT.f_hat;
    
    %% Estimateur YIN/CMNDF
    [L, idxRange, dprime] = yin_pick_L_from_acf_cmndf(AcfIn, Fs, tauMinSamp, tauMaxSamp, Threshold, true, 0.10);
    
    %     [L, idxRange, dprime] = yin_pick_L_from_acf_cmndf_robust(AcfIn, Fs, tauMinSamp, tauMaxSamp);
    %     % On choisit le pic le plus proche de T_samp et on récupère sa position entière L.
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
    [candf_hps, Lref] = cand_hps_from_psd_robust(AcfIn, Fs, fmin, fmax, HPS_dfHz, HPS_K);
    candf_hps = Fs/Lref;
    %% Decide

    fcands = [candf_yin, candf_comb, candf_hps,candf_mediane];
%     [f_twm, info] = decide_f0_twm_cepst(AcfIn, Fs, fcands, candf_Cepstrum, ...
%     'RelPerturb',-0.03:0.005:0.03, ...
%     'MaxPeaks',50, 'Kmax',18, 'wMode','1/k', ...
%     'PriorLambda',0.6, 'PriorCents',35, ...
%     'RefinePct',0.04, 'RefineStep',0.001, ...
%     'SnapTolRel',0.04);
    f_med = median_harmonic_aligned(fcands,'MaxHarm',4,'TolRel',0.06);

% ---------- 5) Cepstrum reconciliation (only if already close)
f_snap = f_med; snapped=false; snap_m=1;
SnapTolRel = 0.06;
if isfinite(candf_Cepstrum) && candf_Cepstrum>0
    M = [0.5 1 2 3 4].';
    f_targets = candf_Cepstrum .* M;
    relerr = abs(f_med - f_targets) ./ max(f_med, eps);
    [emin, idxm] = min(relerr);
    if emin <= SnapTolRel
        candf_Cepstrum = f_targets(idxm);
        snapped = true;
        snap_m = M(idxm);
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
    title(sprintf('ACF %d + consensus final ',idxBestScore));
    
    % --- subplot 2 : barres candidats ---
    subplot(3,1,3);
    vals = [freqRot, candf_Cepstrum,candf_yin, candf_comb, candf_hps,f_med, candf_mediane];
    names= {'real','Cpst','YIN','Comb','HPS','MED','MEDIANE'};
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
    
    % --- Q "par rapport au bruit" (BW mesurée au niveau bruit + τ dB) ---
    tau_dB = 3;            % seuil au-dessus du bruit (à ajuster: 2..6 dB)
    winHz  = 50;           % demi-fenêtre pour estimer le bruit (à ajuster)
    excHz  = 5;            % exclure ±excHz autour de la porteuse (évite de "voir" le pic)
    
    [Qnoise, fL, fR, noise_dB, prom_dB] = q_from_noise(Fsr, Pdbsr, fc, winHz, excHz, tau_dB, df);
    
    % ---- métrique combinée "hauteur × finesse"
    tau_prom = 3;                              % dB (à ajuster 2..6)
    SharpPW  = Qnoise * max(0, prom_dB - tau_prom);
    % (option) compression douce avant normalisation globale
    SharpPW_log = log1p(SharpPW);              % évite outliers
    Qval(i) = SharpPW_log;
    
    % --- recherche Δ par corrélation peigne ---
    %     [Delta(i), SBcomb(i), Cov(i), Asym(i), NoiseL(i)] = ...
    %         best_delta_comb(Fsr, Pdbsr, fc, Opt.K, Opt.DeltaHz, df, Opt.etaBox, Opt.thrSNR_dB, Opt.Delta0_small, F2i);
    [Delta(i), SBcomb(i), Cov(i), Asym(i), NoiseL(i), scanSym] = ...
        best_delta_sym(Fsr, Pdbsr, fc, Opt.K, Opt.DeltaHz, df, Opt.etaBox, Opt.thrSNR_dB, F2i);
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



function [Qnoise, fL, fR, noise_dB,prom_dB] = q_from_noise(F, Pdb, fc, winHz, excHz, tau_dB, df)
% Qnoise = fc / BW_noise, où BW_noise est la largeur du pic mesurée
% au niveau (noise_dB + tau_dB). Interpolation linéaire en dB aux croisements.

if nargin < 7 || isempty(df), df = mean(diff(F)); end
[~, i0] = min(abs(F - fc));
% 1) recentrage sur le max local (±3 bins)
w = 3;
iA = max(2, i0-w);
iB = min(numel(Pdb)-1, i0+w);
[~, krel] = max(Pdb(iA:iB));
i0 = iA + krel - 1;
pk_db = Pdb(i0);
% ---- 1) Bruit local (médiane) dans [fc-winHz, fc+winHz], en excluant ±excHz autour du pic
% iA = max(1, i0 - round(winHz/df));
% iB = min(numel(F), i0 + round(winHz/df));

segP = Pdb(iA:iB);
segF = F(iA:iB);
mask = abs(segF - fc) >= excHz;           % enlève la porteuse ±excHz
if any(mask)
    noise_dB = median(segP(mask), 'omitnan');
else
    noise_dB = median(Pdb, 'omitnan');    % repli global si fenêtre vide
end

% ---- 2) Seuil relatif au bruit
level = noise_dB + tau_dB;

% ---- 3) Recentrage doux: prends le max local proche de fc pour pointer le pic
w = 3; % ±3 bins
iLmtA = max(2, i0-w);
iLmtB = min(numel(Pdb)-1, i0+w);
[~, krel] = max(Pdb(iLmtA:iLmtB));
i0 = iLmtA + krel - 1;

% ---- 4) Croisement gauche (descente en dessous du level)
iL = i0;
while iL > 1 && Pdb(iL) > level, iL = iL - 1; end
if iL == i0 || iL <= 1
    fL = F(max(1,iL));
else
    % interpolation linéaire en dB entre (iL, iL+1)
    p1 = Pdb(iL);   f1 = F(iL);
    p2 = Pdb(iL+1); f2 = F(iL+1);
    alpha = (level - p1) / max(p2 - p1, eps);
    alpha = max(0, min(1, alpha));
    fL = f1 + alpha*(f2 - f1);
end

% ---- 5) Croisement droit
iR = i0;
while iR < numel(Pdb) && Pdb(iR) > level, iR = iR + 1; end
if iR == i0 || iR >= numel(Pdb)
    fR = F(min(numel(F), iR));
else
    p1 = Pdb(iR-1); f1 = F(iR-1);
    p2 = Pdb(iR);   f2 = F(iR);
    alpha = (level - p1) / max(p2 - p1, eps);
    alpha = max(0, min(1, alpha));
    fR = f1 + alpha*(f2 - f1);
end

% ---- 6) Largeur et Q (garde-fous)
BW_noise = max(2*df, fR - fL);        % au moins ~2 bins pour éviter Q infinis
Qnoise   = fc / max(BW_noise, 1e-6);
prom_dB = max(0, pk_db - noise_dB);       % pk_db = Pdbsr(i0_recentré)
end

function [DeltaBest, ScoreBest, CovBest, AsymBest, noiseMed, scan] = ...
    best_delta_sym(F, Pdb, fc, K, DeltaHz, df, eta, thrSNR, F2i)

% Grille Δ
dGrid = linspace(DeltaHz(1), DeltaHz(2), max(41, ceil((DeltaHz(2)-DeltaHz(1))/1))); % pas ≈ 1 Hz

nD   = numel(dGrid);
Score= -inf(nD,1);
Cov  = zeros(nD,1);
Asym = zeros(nD,1);

% Pré-w: poids croissant avec l'ordre (léger)
w_k = @(k) (1 + 0.15*(k-1));

noiseMed = median(Pdb,'omitnan');  % default (au cas où)

for j = 1:nD
    d = dGrid(j);
    boxHz = max(2*df, eta*d);
    
    % --- bruit local hors boîtes (autour de [fc-3d, fc+3d])
    iN1 = F2i(fc - 3*d - 2*boxHz);  iN1 = max(1, min(iN1, numel(F)));
    iN2 = F2i(fc + 3*d + 2*boxHz);  iN2 = max(1, min(iN2, numel(F)));
    if iN2 < iN1, [iN1,iN2] = deal(iN2,iN1); end
    
    mask = true(iN2-iN1+1,1);
    
    % masque porteuse
    iC1 = F2i(fc - boxHz); iC2 = F2i(fc + boxHz);
    iC1 = max(iN1, min(iC1, numel(F)));
    iC2 = max(1,   min(iC2, numel(F)));
    iC1 = max(iN1, iC1); iC2 = min(iN2, iC2);
    if iC2 >= iC1, mask((iC1-iN1+1):(iC2-iN1+1)) = false; end
    
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
        if iL2 >= iL1, mask((iL1-iN1+1):(iL2-iN1+1)) = false; end
        if iR2 >= iR1, mask((iR1-iN1+1):(iR2-iN1+1)) = false; end
    end
    
    seg = Pdb(iN1:iN2);
    if any(mask)
        noiseMed = median(seg(mask), 'omitnan');
    else
        noiseMed = median(Pdb, 'omitnan');
    end
    
    % --- accumulation symétrique
    Wsum=0; EsymSum=0; Asum=0; used=0; covered=0;
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
        SNRL = max(0, EL - noiseMed);
        SNRR = max(0, ER - noiseMed);
        
        % --- symétrie stricte
        % option A (simple, robuste) :
        symk = min(SNRL, SNRR);
        % option B (plus sensible à un déséquilibre) :
        % symk = 2 / ( 1/max(SNRL,eps) + 1/max(SNRR,eps) );  % moyenne harmonique
        
        wk = w_k(k);
        EsymSum = EsymSum + wk * symk;
        Wsum    = Wsum + wk;
        
        used = used + 1;
        if SNRL>thrSNR && SNRR>thrSNR, covered = covered + 1; end
        Asum = Asum + abs(SNRL - SNRR);
    end
    
    if used >= 1 && Wsum>0
        Score(j) = EsymSum / Wsum;  % pur score symétrique
        Cov(j)   = covered / used;
        Asym(j)  = Asum / used;
    end
end

% --- meilleur Δ sur la grille
[ScoreBest, jbest] = max(Score);
DeltaBest = dGrid(jbest);
CovBest   = Cov(jbest);
AsymBest  = Asym(jbest);

% --- raffinement parabolique sur (Δ, Score)
if jbest>1 && jbest<nD && all(isfinite(Score([jbest-1 jbest jbest+1])))
    d1 = dGrid(jbest-1); s1 = Score(jbest-1);
    d2 = dGrid(jbest);   s2 = Score(jbest);
    d3 = dGrid(jbest+1); s3 = Score(jbest+1);
    % parabole passant par (d1,s1),(d2,s2),(d3,s3)
    denom = (d1-d2)*(d1-d3)*(d2-d3);
    if abs(denom) > eps
        A = (d3*(s2-s1)+d2*(s1-s3)+d1*(s3-s2)) / denom;
        B = (d3^2*(s1-s2)+d2^2*(s3-s1)+d1^2*(s2-s3)) / denom;
        d_refine = -B/(2*A);
        if isfinite(d_refine) && d_refine>=min(dGrid) && d_refine<=max(dGrid)
            DeltaBest = d_refine;   % Δ affiné
        end
    end
end

% (option) retourner la courbe pour debug
scan = struct('Delta', dGrid(:), 'ScoreSym', Score(:), 'Coverage', Cov(:), 'Asym', Asym(:));
end

function [Lbest, Sbest] = scorePeigneGridRefine(acf_s, tauMinSamp, tauMaxSamp, opts)
% SCOREPEIGNEGRIDREFINE (version robuste fondamental)
%  - Score peigne log sur ACF lissée/normalisée positive
%  - Grille grossière (>=1)
%  - Pénalisation implicite des harmoniques manquants (log(eps))
%  - Sélection du fondamental par peigne de multiples (1/m)
%  - Raffinement local (pas 1) + interpolation parabolique

if nargin<4, opts = struct(); end
if ~isfield(opts,'Kharm'),           opts.Kharm = 6; end
if ~isfield(opts,'WeightMode'),      opts.WeightMode = '1/k'; end   % ou 'equal'
if ~isfield(opts,'Eps'),             opts.Eps = 1e-3; end
if ~isfield(opts,'gridStepSamp'),    opts.gridStepSamp = 2; end
if ~isfield(opts,'needAtLeastHarm'), opts.needAtLeastHarm = 2; end
if ~isfield(opts,'FundMaxMult'),     opts.FundMaxMult = 6; end      % nb de multiples pour le vote
if ~isfield(opts,'RefineHalfWin'),   opts.RefineHalfWin = 6; end    % demi-fenêtre de raffinement en échantillons (>=2)

acf_s = acf_s(:);
Nlag  = numel(acf_s);
tauMaxSamp = min(tauMaxSamp, Nlag-1);

% ---- Poids harmonique ----
K  = opts.Kharm;
if strcmpi(opts.WeightMode,'equal')
    wk = ones(1,K);
else
    wk = 1./(1:K);    % par défaut: 1/k
end
wk = wk / sum(wk);    % normalisation globale (constante), PAS par kmax

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
    kmax = min(K, floor(tauMaxSamp / max(Lg,1)));
    if kmax < opts.needAtLeastHarm, continue; end
    
    % on construit un vecteur 'vals_full' de longueur K (les manquants valent 0)
    vals_full = zeros(K,1);
    idx = (1:kmax).*Lg + 1;     % +1 car r(1) = lag0
    vals_full(1:kmax) = rpos(idx);
    
    % score = somme pondérée des logs, sans renormaliser par kmax
    ScoreG(s) = sum( wk(:) .* logeps(vals_full(:)) );  % /sum(wk) optionnel mais constant
end

% --- 1) Choisir le fondamental par vote sur les multiples
[L0, ~] = pickFundamentalByMultiples(tauGrid, ScoreG, opts.FundMaxMult);

% --- 2) Raffinement local (pas 1) autour de L0
halfW = max(2, round(opts.RefineHalfWin));
Lmin  = max(tauMinSamp, L0 - halfW);
Lmax  = min(tauMaxSamp, L0 + halfW);
Lloc  = (Lmin:1:Lmax).';
ScoreLoc = arrayfun(@(L) scoreOneLag(rpos, L, K, wk, logeps, tauMaxSamp, opts.needAtLeastHarm), Lloc);

% garder le meilleur entier local
[SlocBest, iLoc] = max(ScoreLoc);
Lbest = Lloc(iLoc);
Sbest = SlocBest;

% --- 3) Interpolation parabolique sur (Lbest-1, Lbest, Lbest+1) pour précision sub-échantillon
if iLoc>1 && iLoc<numel(Lloc) && isfinite(ScoreLoc(iLoc-1)) && isfinite(ScoreLoc(iLoc+1))
    f1 = Lloc(iLoc-1); s1 = ScoreLoc(iLoc-1);
    f2 = Lloc(iLoc);   s2 = ScoreLoc(iLoc);
    f3 = Lloc(iLoc+1); s3 = ScoreLoc(iLoc+1);
    denom = (f1-f2)*(f1-f3)*(f2-f3);
    if abs(denom) > eps
        A = (f3*(s2-s1)+f2*(s1-s3)+f1*(s3-s2)) / denom;
        B = (f3^2*(s1-s2)+f2^2*(s3-s1)+f1^2*(s2-s3)) / denom;
        Lref = -B/(2*A);
        % on ne renvoie que Lbest entier + Sbest demandé; mais tu peux récupérer Lref si besoin
        % (si tu veux sortir Lref aussi: change la signature)
    end
end
end

% ---------- helpers ----------

function S = scoreOneLag(rpos, L, K, wk, logeps, tauMaxSamp, needAtLeastHarm)
kmax = min(K, floor(tauMaxSamp / max(L,1)));
if kmax < needAtLeastHarm
    S = -inf; return;
end
vals_full = zeros(K,1);
idx = (1:kmax).*L + 1;
vals_full(1:kmax) = rpos(idx);
S = sum( wk(:) .* logeps(vals_full(:)) );
end

function [Lbest, Scomb] = pickFundamentalByMultiples(tauGrid, ScoreG, M)
% Vote fondamental : moyenne pondérée des scores aux multiples (poids 1/m)
if nargin<3, M=6; end
tauGrid = tauGrid(:); ScoreG = ScoreG(:);
S = ScoreG; S(~isfinite(S)) = -inf;

w = 1./(1:M); w = w(:);
Scomb = -inf(size(tauGrid));
for i = 1:numel(tauGrid)
    L = double(tauGrid(i));
    if ~isfinite(S(i)), continue; end
    vals = []; ww = [];
    for m = 1:M
        pos = L*m;
        if pos > double(tauGrid(end)), break; end
        [~, j] = min(abs(double(tauGrid) - pos));
        if isfinite(S(j))
            vals(end+1,1) = S(j); %#ok<AGROW>
            ww(end+1,1)   = w(m);
        end
    end
    if ~isempty(vals)
        Scomb(i) = sum(ww.*vals)/sum(ww);   % moyenne pondérée (évite de favoriser trop les tout petits L par la somme)
    end
end
[~, ib] = max(Scomb);
Lbest = tauGrid(ib);
end


function [fhat, L0] = cand_hps_from_psd_robust(r, Fs, fmin, fmax, dfHz, K, Opt)
% Robust HPS on PSD (from ACF) with local-SNR, anchor on k=1..2, and anti-subharmonic bias.
% Defaults:
%   Opt.bwHz = 0.8;        % averaging half-band around fk (≈ peak width)
%   Opt.noiseRingHz = 5;   % ring to estimate noise (median) around fk
%   Opt.tau_dB = 0;        % SNR threshold in dB (0..3 typical if you work in dB)
%   Opt.p = 1.0;           % weights 1/k^p
%   Opt.alpha_anchor = 0.5;% mixes overall sum with anchor fraction
%   Opt.require12 = [true true]; % require k=1 and k=2 above threshold
%   Opt.antiSub = 0.5;     % penalize when SNR1 << SNR2 (0..1)
%   Opt.fminHard = max(3,fmin); % hard lower bound to avoid LF artefacts

if nargin<7, Opt = struct; end
if ~isfield(Opt,'bwHz'),          Opt.bwHz = 0.8; end
if ~isfield(Opt,'noiseRingHz'),   Opt.noiseRingHz = 5; end
if ~isfield(Opt,'tau_dB'),        Opt.tau_dB = 0; end
if ~isfield(Opt,'p'),             Opt.p = 1.0; end
if ~isfield(Opt,'alpha_anchor'),  Opt.alpha_anchor = 0.5; end
if ~isfield(Opt,'require12'),     Opt.require12 = [true true]; end
if ~isfield(Opt,'antiSub'),       Opt.antiSub = 0.5; end
if ~isfield(Opt,'fminHard'),      Opt.fminHard = max(3, fmin); end

% --- PSD via Wiener–Khinchin
[F, Pxx] = psd_from_acf(r, Fs);       % your helper (linear power)
Pxx = max(Pxx, 0);
% zero anything above fmax-1 as you did
Pxx(F >= fmax-1) = 0;
% (optional) also zero ultra-low freq region below ~fminHard-1
Pxx(F <= max(0, Opt.fminHard-1)) = 0;

% local SNR helper at frequency f (linear scale)
df = mean(diff(F));
specBand = @(f,halfHz) mean_at_band(F,Pxx,f,halfHz);
specNoise= @(f,ringHz,halfHz) median_ring(F,Pxx,f,ringHz,halfHz);

fgrid = (max(fmin,Opt.fminHard):dfHz:fmax).';
S = -inf(numel(fgrid),1);

for i = 1:numel(fgrid)
    f0 = fgrid(i);
    snr = zeros(K,1); useK = 0;
    for k=1:K
        fk = k*f0;
        if fk > F(end), break; end
        band = specBand(fk, Opt.bwHz/2);
        noise= specNoise(fk, Opt.noiseRingHz, Opt.bwHz/2);
        snr_lin = max(0, band - noise);
        % work in dB if you prefer thresholds in dB:
        snr_db  = 10*log10(max(snr_lin, eps));
        snr(k)  = max(0, snr_db - Opt.tau_dB);   % thresholded SNR(dB)
        useK    = useK + 1;
    end
    if useK < 2, continue; end
    
    snr = snr(1:useK);
    w   = (1./((1:useK)'.^Opt.p));  w = w/sum(w);
    
    % anchors on k=1 and k=2
    if Opt.require12(1) && snr(1)<=0,  continue; end
    if useK>=2 && Opt.require12(2) && snr(2)<=0, continue; end
    
    % base score = weighted average of SNRs
    base = sum(w .* snr);
    
    % anchor fraction: ensure k=1 contributes
    frac1 = snr(1) / max(sum(snr), eps);
    score = base * (Opt.alpha_anchor + (1-Opt.alpha_anchor)*frac1);
    
    % anti-subharmonic: penalize if snr1 << snr2
    if useK>=2
        ratio12 = snr(1) / max(snr(2), eps);
        pen = min(1, ratio12);           % in [0..1]
        score = score * (Opt.antiSub*pen + (1-Opt.antiSub));
    end
    
    S(i) = score;
end

% argmax + parabolic refine (as you did)
[~,ix] = max(S);
if ix>1 && ix<numel(S)
    y1=S(ix-1); y2=S(ix); y3=S(ix+1);
    delta = 0.5*(y1-y3)/max(y1-2*y2+y3, eps);
    delta = max(min(delta,0.5),-0.5);
else
    delta = 0;
end
fhat = fgrid(ix) + delta*dfHz;
L0 = Fs / fhat;
end

% --------- helpers ----------
function m = mean_at_band(F,P,f,halfHz)
if f< F(1) || f>F(end), m=0; return; end
df = mean(diff(F));
h  = max(1, round(halfHz/df));
[~,i0] = min(abs(F-f));
i1 = max(1, i0-h); i2 = min(numel(F), i0+h);
m  = mean(P(i1:i2), 'omitnan');
end
function med = median_ring(F,P,f,ringHz,halfHz)
df = mean(diff(F));
[~,i0]= min(abs(F-f));
h    = max(1, round(halfHz/df));
r    = max(1, round(ringHz/df));
iL1  = max(1, i0-r); iL2 = max(1, i0-h-1);
iR1  = min(numel(F), i0+h+1); iR2 = min(numel(F), i0+r);
seg = [];
if iL2>=iL1, seg = [seg; P(iL1:iL2)]; end
if iR2>=iR1, seg = [seg; P(iR1:iR2)]; end
if isempty(seg), med = median(P,'omitnan'); else, med = median(seg,'omitnan'); end
end

function [L, idxRange, dprime] = yin_pick_L_from_acf_cmndf_robust(r_ns, Fs, tauMinSamp, tauMaxSamp, opts)
% YIN/CMNDF robuste sans seuil absolu : sélection par score de vallées.
% opts (facultatif):
%   .Eps         (1e-12)   % pour CMNDF
%   .ValleyTol   (0.05)    % tolérance pour mesurer largeur de vallée
%   .PrefShort   (true)    % préférer le plus petit lag si ~équivalent
%   .KeepFrac    (0.10)    % tolérance (10%) pour garder candidats ~meilleurs
%   .PenTau      (0.25)    % pénalité ~ 1/(1+PenTau*(tau/tauMax))
%   .HarmK       (3)       % nb d’harmoniques pour cohérence (2..3)
%   .HarmTol     (0.05)    % d'(m*tau) < d'(tau) - HarmTol  => suspicion sous-harmonique
%
% Sorties:
%   L        : lag entier choisi
%   idxRange : fenêtre locale (utile pour ta parabole)
%   dprime   : CMNDF (pour debug)

if nargin<5, opts = struct; end
if ~isfield(opts,'Eps'),       opts.Eps = 1e-12; end
if ~isfield(opts,'ValleyTol'), opts.ValleyTol = 0.05; end
if ~isfield(opts,'PrefShort'), opts.PrefShort = true; end
if ~isfield(opts,'KeepFrac'),  opts.KeepFrac = 0.10; end
if ~isfield(opts,'PenTau'),    opts.PenTau = 0.25; end
if ~isfield(opts,'HarmK'),     opts.HarmK = 5; end
if ~isfield(opts,'HarmTol'),   opts.HarmTol = 0.05; end

% --- garde-fous ---
r_ns = r_ns(:);
Nlag  = numel(r_ns);
tauMinSamp = max(2, floor(tauMinSamp));
tauMaxSamp = min(max(tauMinSamp+2, floor(tauMaxSamp)), Nlag-1);

% --- CMNDF ---
tauIdx = (1:tauMaxSamp).';
d      = 2.0 * (1 - r_ns(tauIdx+1));
cum    = cumsum(d);
dprime = d .* (tauIdx ./ max(cum, opts.Eps));
dprime(1:tauMinSamp-1) = +Inf;

% --- détecter les minima locaux de dprime ---
isMin = false(size(dprime));
for i = (tauMinSamp):(tauMaxSamp-1)
    if isfinite(dprime(i)) && dprime(i) <= dprime(i-1) && dprime(i) < dprime(i+1)
        isMin(i) = true;
    end
end
mins = find(isMin);
if isempty(mins)
    % fallback: min global dans la bande
    [~,L] = min(dprime);
    if L<tauMinSamp || L>tauMaxSamp
        [~,L] = min(dprime(tauMinSamp:tauMaxSamp)); L = L + tauMinSamp - 1;
    end
    idxRange = [max(tauMinSamp, L-3), min(tauMaxSamp, L+3)];
    return
end

% --- scorer chaque vallée ---
scores = zeros(numel(mins),1);
widths = zeros(numel(mins),1);
for k = 1:numel(mins)
    i = mins(k);             % indice tau
    tau = i;
    dval = dprime(i);
    
    % largeur de vallée: d' <= dval + tol autour de i
    tol  = opts.ValleyTol;
    L1 = i; while L1>tauMinSamp && dprime(L1) <= dval + tol, L1 = L1-1; end
    R1 = i; while R1<tauMaxSamp && dprime(R1) <= dval + tol, R1 = R1+1; end
    widths(k) = (R1 - L1 - 1);
    
    % force ACF au lag (pic corrélé)
    rpk = max(0, r_ns(tau+1));  % r_ns(1)=1 donc tau+1 ok
    
    % cohérence harmonique: punir si d'(m*tau) << d'(tau)
    harmOK = 1.0;
    for m = 2:opts.HarmK
        j = m*tau;
        if j <= tauMaxSamp && isfinite(dprime(j))
            if dprime(j) < dval - opts.HarmTol
                harmOK = harmOK * 0.7;  % pénalise sous-harmonique
            end
        end
    end
    
    % pénalité τ long (préférence lags courts si proche)
    penTau = 1.0 / (1.0 + opts.PenTau * (tau / tauMaxSamp));
    
    % score composite (exposants doux)
    scores(k) = (1 - dval)^(1.0) * (0.5 + 0.5*rpk) * (1 + 0.1*widths(k)) * harmOK * penTau;
end

% --- choisir meilleur, puis préférer le plus petit dans la marge ---
[smax, imax] = max(scores);
bestIdx = mins(imax);
if opts.PrefShort
    keep = mins(scores >= (1 - opts.KeepFrac) * smax);
    L = min(keep);
else
    L = bestIdx;
end

% --- fenêtre locale pour ta parabole d’affinage en aval ---
win = max(3, round(0.1*L));
i1  = max(tauMinSamp, L - win);
i2  = min(tauMaxSamp, L + win);
idxRange = [i1 i2];
% --- Post-fix: promote submultiples if comparable (avoid picking m*T)
betaRel   = 0.20;       % tolerance: accept submultiple if score >= (1-beta)*score_best
harmBoost = 0.10;       % small bonus favoring shorter period (more cycles in window)
candDiv   = [2 3 4];    % test L/2, L/3, L/4

L_ref = L;
S_ref = smax;           % best score at L (from your loop)

bestL = L_ref; bestS = S_ref;
for d = candDiv
    Lc = round(L_ref / d);
    if Lc < tauMinSamp || Lc > tauMaxSamp, continue; end
    
    % recompute valley score at Lc (same recipe as in the main loop)
    i = Lc; dval = dprime(i);
    tol = opts.ValleyTol;
    L1=i; while L1>tauMinSamp && dprime(L1) <= dval+tol, L1=L1-1; end
    R1=i; while R1<tauMaxSamp && dprime(R1) <= dval+tol, R1=R1+1; end
    widthLc = (R1 - L1 - 1);
    rpkLc   = max(0, r_ns(i+1));
    
    harmOK  = 1.0;
    for m = 2:opts.HarmK
        j = m*i;
        if j <= tauMaxSamp && isfinite(dprime(j))
            if dprime(j) < dval - opts.HarmTol
                harmOK = harmOK * 0.7;
            end
        end
    end
    
    penTau  = 1.0 / (1.0 + opts.PenTau * (i / tauMaxSamp));
    S_cand  = (1 - dval) * (0.5 + 0.5*rpkLc) * (1 + 0.1*widthLc) * harmOK * penTau;
    
    % small bias for shorter, plausible fundamentals
    S_adj = S_cand * (1 + harmBoost);
    
    if S_adj >= (1 - betaRel) * S_ref
        bestL = i; bestS = S_cand;
        % keep looping to allow L/3 to beat L/2 if both qualify; we take smallest later
    end
end

% prefer the smallest qualifying lag
if bestL ~= L
    L = bestL;
    win = max(3, round(0.1*L));
    i1  = max(tauMinSamp, L - win);
    i2  = min(tauMaxSamp, L + win);
    idxRange = [i1 i2];
end
end
function [f_final, out] = decide_f0_twm_cepst( ...
    x, Fs, fcands, f_cepst, varargin)
% Robust TWM decision with prior-to-seeds, local refine, guards, and cepstrum snap.

% ---------- parse
p = inputParser;
% Grid expansion
addParameter(p,'Octaves',[0.5 1 2]);
addParameter(p,'RelPerturb',-0.03:0.005:0.03);
addParameter(p,'fminmax',[0.5 200]);
addParameter(p,'MaxGrid',250);
addParameter(p,'Nfft',[]);
% TWM params
addParameter(p,'MaxPeaks',50);
addParameter(p,'Kmax',18);
addParameter(p,'rho',0.33);
addParameter(p,'q',1.4);
addParameter(p,'wMode','1/k');
% Guards
addParameter(p,'GapRel',0.05);
addParameter(p,'NearRel',0.03);
addParameter(p,'PromoteRel',0.08);
% Anchor check
addParameter(p,'AnchorLambda',0.15);
addParameter(p,'AnchorBwHz',1.0);
addParameter(p,'AnchorRingHz',5.0);
% Prior to seeds (in cents on log-f)
addParameter(p,'PriorLambda',0.6);    % strength of prior (0 = off)
addParameter(p,'PriorCents',35);      % width (1 semitone = 100 cents)
% Local refinement around TWM winner
addParameter(p,'RefinePct',0.04);     % ±4% window
addParameter(p,'RefineStep',0.001);   % 0.1% step
% Cepstral snap
addParameter(p,'SnapTolRel',0.04);


parse(p,varargin{:});
opt = p.Results;

% ---------- 1) TWM with octave grid
[f_twm_raw, costs_raw, candGrid, infoGrid] = twm_select_with_octave_grid( ...
    x, Fs, fcands, ...
    'Octaves',opt.Octaves, 'RelPerturb',opt.RelPerturb, ...
    'fminmax',opt.fminmax, 'MaxGrid',opt.MaxGrid, ...
    'Nfft',opt.Nfft, 'MaxPeaks',opt.MaxPeaks, 'Kmax',opt.Kmax, ...
    'rho',opt.rho, 'q',opt.q, 'wMode',opt.wMode);

% Spectrum once for anchor scoring
[Fspec, Pspec] = local_spectrum_periodogram(x, Fs, opt.Nfft);

% ---------- 2) Anchor-adjusted cost
costs_adj = costs_raw;
for i = 1:numel(candGrid)
    aScore = anchor_score(Fspec,Pspec,candGrid(i),opt.AnchorBwHz,opt.AnchorRingHz); % 0..1
    costs_adj(i) = costs_adj(i) * (1 + opt.AnchorLambda*(1 - aScore));
end

% ---------- 3) Prior to seeds (cents penalty)
if ~isempty(fcands) && opt.PriorLambda > 0
    fc = fcands(:); fc = fc(isfinite(fc) & fc>0);
    if ~isempty(fc)
        prior = zeros(size(candGrid));
        for i = 1:numel(candGrid)
            f0 = candGrid(i);
            cents = 1200*min(abs(log2(f0./fc)));  % distance to nearest seed (in cents)
            prior(i) = opt.PriorLambda * (cents/opt.PriorCents)^2; % quadratic penalty
        end
        costs_adj = costs_adj .* (1 + prior);
    end
end

% pick best with guards (gap tie-break)
[cs, ord] = sort(costs_adj,'ascend');
bestIdx = ord(1); f_twm = candGrid(bestIdx);
if numel(cs) >= 2
    gap = (cs(2)-cs(1)) / max(cs(2), eps);
else
    gap = inf;
end
if gap < opt.GapRel
    near = candGrid(costs_adj <= cs(1)*(1+opt.NearRel));
    if ~isempty(near), f_twm = min(near); end    % prefer smallest in near-tie
end

% submultiple promotion: f/2, f/3
promote = [0.5 1/3];
candProm = f_twm .* promote;
candProm = candProm(candProm>=opt.fminmax(1) & candProm<=opt.fminmax(2));
if ~isempty(candProm)
    [~, cProm, ~] = twm_select_f0(x, Fs, [f_twm candProm], ...
        'Nfft',opt.Nfft, 'MaxPeaks',opt.MaxPeaks, 'Kmax',opt.Kmax, ...
        'rho',opt.rho, 'q',opt.q, 'wMode',opt.wMode);
    % anchor+prior adjust promoted costs
    for j = 1:numel(candProm)
        aScore = anchor_score(Fspec,Pspec,candProm(j),opt.AnchorBwHz,opt.AnchorRingHz);
        adj = (1 + opt.AnchorLambda*(1 - aScore));
        % prior
        if ~isempty(fcands) && opt.PriorLambda > 0
            cents = 1200*min(abs(log2(candProm(j)./fcands(:))));
            prior = opt.PriorLambda * (min(cents)/opt.PriorCents)^2;
            adj = adj * (1 + prior);
        end
        cProm(1+j) = cProm(1+j) * adj;
    end
    if any(cProm(2:end) <= cs(1)*(1+opt.PromoteRel))
        f_twm = min([f_twm, candProm(cProm(2:end) <= cs(1)*(1+opt.PromoteRel))]);
    end
end

% ---------- 4) Local refinement (±RefinePct around f_twm)
if isfinite(f_twm) && f_twm>0
    sLo = 1 - opt.RefinePct; sHi = 1 + opt.RefinePct;
    sGrid = sLo:opt.RefineStep:sHi;
    fLoc  = f_twm .* sGrid;
    fLoc  = fLoc(fLoc>=opt.fminmax(1) & fLoc<=opt.fminmax(2));
    [~, cLoc, ~] = twm_select_f0(x, Fs, fLoc, ...
        'Nfft',opt.Nfft, 'MaxPeaks',opt.MaxPeaks, 'Kmax',opt.Kmax, ...
        'rho',opt.rho, 'q',opt.q, 'wMode',opt.wMode);
    % anchor+prior adjust local costs
    for i = 1:numel(fLoc)
        aScore = anchor_score(Fspec,Pspec,fLoc(i),opt.AnchorBwHz,opt.AnchorRingHz);
        adj = (1 + opt.AnchorLambda*(1 - aScore));
        if ~isempty(fcands) && opt.PriorLambda > 0
            cents = 1200*min(abs(log2(fLoc(i)./fcands(:))));
            prior = opt.PriorLambda * (min(cents)/opt.PriorCents)^2;
            adj = adj * (1 + prior);
        end
        cLoc(i) = cLoc(i) * adj;
    end
    [~, iBestLoc] = min(cLoc);
    f_twm = fLoc(iBestLoc);
end



f_final = f_twm;

% ---------- outputs
out = struct();
out.f_twm_raw     = f_twm_raw;
out.f_twm_guarded = f_twm;
out.costs_raw     = costs_raw;
out.costs_adj     = costs_adj;
out.candGrid      = candGrid;
out.infoGrid      = infoGrid;
out.params        = opt;
end

% ===== helpers =====
function [F, P] = local_spectrum_periodogram(x, Fs, Nfft)
x = x(:); x = x - mean(x,'omitnan');
if isempty(Nfft), Nfft = 2^nextpow2(max(4096, numel(x))); end
Nw = min(numel(x), Nfft);
w  = hann(Nw,'periodic');
X  = fft([x(1:Nw).*w; zeros(Nfft-Nw,1)], Nfft);
P  = (abs(X).^2) / (sum(w.^2)+eps);
F  = (0:Nfft-1).' * (Fs/Nfft);
keep = (F <= Fs/2); F=F(keep); P=P(keep);
end

function a = anchor_score(F,P,f0,bwHz,ringHz)
snr1 = local_snr(F,P,f0,  bwHz,ringHz);
snr2 = local_snr(F,P,2*f0,bwHz,ringHz);
map = @(s) min(1, max(0, (10*log10(max(s,eps)) + 3) / 15)); % -3..+12 dB -> 0..1
a = 0.6*map(snr1) + 0.4*map(snr2);
end

function s = local_snr(F,P,f, bwHz, ringHz)
if f <= F(1) || f >= F(end), s = 0; return; end
df = mean(diff(F));
h  = max(1, round((bwHz/2)/df));
r  = max(1, round(ringHz/df));
[~,i0] = min(abs(F - f));
i1 = max(1, i0 - h); i2 = min(numel(F), i0 + h);
band = mean(P(i1:i2), 'omitnan');
jL1 = max(1, i0 - r); jL2 = max(1, i0 - h - 1);
jR1 = min(numel(F), i0 + h + 1); jR2 = min(numel(F), i0 + r);
seg = [];
if jL2 >= jL1, seg = [seg; P(jL1:jL2)]; end
if jR2 >= jR1, seg = [seg; P(jR1:jR2)]; end
noise = median(seg, 'omitnan'); if isempty(seg), noise = median(P,'omitnan'); end
s = max(0, band - noise);
end

function f_med = median_harmonic_aligned(f_cands, varargin)
% MEDIAN_HARMONIC_ALIGNED
%   Aligne un ensemble de fréquences candidates (Hz) sur une fondamentale
%   commune en corrigeant les octaves/harmoniques, puis retourne la médiane.
%
% Usage:
%   f_med = median_harmonic_aligned(f_cands)
%   f_med = median_harmonic_aligned(f_cands, 'MaxHarm',4,'TolRel',0.05)
%
% In:
%   f_cands : vecteur de fréquences (Hz)
% Options:
%   'MaxHarm' (4)  : nombre max d’harmoniques testés (÷2, ÷3, ÷4)
%   'TolRel' (0.05): tolérance relative (5%) pour décider si f/k ≈ ref
%
% Out:
%   f_med   : médiane des valeurs alignées (Hz)

p = inputParser;
addParameter(p,'MaxHarm',4);
addParameter(p,'TolRel',0.05);
parse(p,varargin{:});
MaxH = p.Results.MaxHarm;
TolR = p.Results.TolRel;

f_cands = f_cands(:);
f_cands = f_cands(isfinite(f_cands) & f_cands>0);
if isempty(f_cands)
    f_med = NaN; return;
end

% 1) fréquence de référence brute
f_ref = median(f_cands);

% 2) ramener chaque candidat au plus proche sous-multiple (1..MaxHarm)
f_aligned = nan(size(f_cands));
for i=1:numel(f_cands)
    f = f_cands(i);
    best = f;
    for k=1:MaxH
        f_try = f / k;
        if abs(f_try - f_ref)/f_ref <= TolR
            best = f_try; break;  % on garde le premier qui colle
        end
    end
    f_aligned(i) = best;
end

% 3) médiane des valeurs alignées
f_med = median(f_aligned);
end


