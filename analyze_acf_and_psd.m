function R = analyze_acf_and_psd(acf, Fs, varargin)
% ANALYZE_ACF_AND_PSD  Analyse un ACF pour extraire la période/rotation + PSD.
%
% R = analyze_acf_and_psd(acf, Fs, 'Param', value, ...)
%
% Obligatoire
%   acf : vecteur ACF (lags >= 0), taille N
%   Fs  : fréquence d'échantillonnage (Hz)
%
% Paramètres (optionnels)
%   'MaxLagSec'  : (default 1.0)  durée max analysée sur l'ACF (s)
%   'TopPeaks'   : (default 15)   nb max de pics utilisés pour la médiane des périodes
%   'Prominence' : (default 0.05) proéminence min (fraction du max ACF norm.)
%   'SmoothMs'   : (default 2)    lissage (ms) avant détection de pics
%   'PsdFmax'    : (default 300)  borne d'affichage de la PSD (Hz)
%   'ExportPrefix': (default '')  préfixe fichiers d'export (PNG/CSV); vide => pas d'export
%
% Sortie (struct R)
%   R.period_s, R.f_hz, R.rpm, R.omega
%   R.peaks_samples, R.peaks_time_s, R.peaks_values
%   R.psd_F, R.psd_P, R.psd_top_freqs
%   R.params : paramètres effectifs

% ---------- Params ----------
p = inputParser;
p.addParameter('MaxLagSec',   1.0,   @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('TopPeaks',    15,    @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('Prominence',  0.05,  @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('SmoothMs',    2,     @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('PsdFmax',     300,   @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('ExportPrefix','',    @(x)ischar(x)||isstring(x));
p.parse(varargin{:});
P = p.Results;

acf = acf(:);
N   = numel(acf);
lags = (0:N-1).';% [0 ,1 ,2 ... ,5999]
tau  = lags / Fs;% [0 ,0.0002 ,0.0003,... ,0.999] 


acf_n = acf;


% Limite en lag (échantillons)
maxLagSamp = min(N-1, round(P.MaxLagSec * Fs));%5999 le dernier

% ---------- Détection de pics ACF ----------
% wSamp représente la taille de la fenêtre de lissage (en nombre d’échantillons) appliquée sur l’ACF.
wSamp = max(1, round(P.SmoothMs * 1e-3 * Fs)); % P.SmoothMs * 1e-3 * Fs : convertit la durée de lissage en nombre d’échantillons.
acf_s = movavg(acf_n, wSamp);

% seuil de proéminence (relatif au max hors lag 0)
base = max(acf_s(2:min(end,maxLagSamp+1)));
prom = P.Prominence * base;

[pkVal, pkLoc] = findpeaks_safe(acf_s(1:maxLagSamp+1), ...
    'MinPeakProminence', prom);

% Retire lag 0 si présent
pkMask = pkLoc > 1; 
pkLoc  = pkLoc(pkMask);
pkVal  = pkVal(pkMask);

% Garde au plus TopPeaks (les plus forts au début)
if numel(pkLoc) > P.TopPeaks
    [~, idx] = maxk(pkVal, P.TopPeaks);
    idx = sort(idx);
    pkLoc = pkLoc(idx);
    pkVal = pkVal(idx);
end

% Estimation de période par médiane des différences successives
period_samples = NaN;  fr_est = NaN;
if numel(pkLoc) >= 2
    d = diff(pkLoc);
    period_samples = median(double(d));
    fr_est = Fs / period_samples;
elseif numel(pkLoc) == 1
    period_samples = double(pkLoc);
    fr_est = Fs / period_samples;
end

% ---------- PSD via Wiener–Khinchin ----------
[Freq, Pxx] = psd_from_acf(acf_n, Fs);

% Pics PSD sous PsdFmax
Fmask = (Freq >= 0.1) & (Freq <= P.PsdFmax);
[Ftops, ~] = local_peaks(Freq(Fmask), Pxx(Fmask), 8); % top 8

% ---------- Figures ----------
figure('Name','ACF analysée'); 
plot(tau, acf_n, 'b'); grid on; hold on;
stem(tau(pkLoc), acf_n(pkLoc), 'filled', 'Color',[0.85 0.33 0.10]);
xlabel('\tau (s)'); ylabel('ACF (norm.)'); xlim([0, P.MaxLagSec]);
title('ACF normalisée & pics');
if isfinite(fr_est)
    T_est = 1/fr_est;
    xline(T_est, '--', sprintf('T≈%.4f s  (f≈%.2f Hz)', T_est, fr_est), ...
        'Color',[0.2 0.2 0.2]);
end
legend('ACF','Pics','Location','best');

figure('Name','PSD (Wiener–Khinchin)');
plot(Freq, Pxx, 'b'); grid on; hold on;
for k=1:numel(Ftops)
    xline(Ftops(k), 'r--', sprintf('%.2f Hz', Ftops(k)));
end
xlabel('Fréquence (Hz)'); ylabel('Puissance (u.a.)');
xlim([0, P.PsdFmax]);
title('PSD estimée à partir de l''ACF');

% ---------- Export (optionnel) ----------
if ~isempty(P.ExportPrefix)
    prefix = char(P.ExportPrefix);
    % PNG
    figs = findobj('Type','figure');
    for i=1:numel(figs)
        figname = sprintf('%s_fig%d.png', prefix, i);
        saveas(figs(i), figname);
    end
    % CSV pics ACF
    T = table((pkLoc-1)/Fs, acf_n(pkLoc), 'VariableNames',{'tau_s','acf_norm'});
    writetable(T, sprintf('%s_peaks_acf.csv', prefix));
    % CSV PSD (tronquée à F<=PsdFmax)
    Tpsd = table(Freq(Fmask), Pxx(Fmask), 'VariableNames',{'F_Hz','P'});
    writetable(Tpsd, sprintf('%s_psd.csv', prefix));
end

% ---------- Sorties ----------
R = struct();
R.period_s   = 1/fr_est;
R.f_hz       = fr_est;
R.rpm        = fr_est * 60;
R.omega      = 2*pi*fr_est;
R.peaks_samples = pkLoc-1;
R.peaks_time_s  = (pkLoc-1)/Fs;
R.peaks_values  = acf_n(pkLoc);
R.psd_F      = Freq;
R.psd_P      = Pxx;
R.psd_top_freqs = Ftops;
R.params     = P;

% --------- sous-fonctions locales ---------
function y = movavg(x,w)
    if w<=1, y=x; return; end
    k = ones(w,1)/w;
    y = filtfilt(k,1,double(x));
end

function [f_peaks, v_peaks] = local_peaks(F, Pxx, K)
    % petite détection de pics sans toolbox
    idx = [];
    for i=2:numel(F)-1
        if Pxx(i)>Pxx(i-1) && Pxx(i)>Pxx(i+1)
            idx(end+1)=i; %#ok<AGROW>
        end
    end
    if isempty(idx)
        f_peaks = []; v_peaks = [];
        return;
    end
    [~, order] = maxk(Pxx(idx), min(K, numel(idx)));
    idx = sort(idx(order));
    f_peaks = F(idx); v_peaks = Pxx(idx);
end

end % ===== fin fonction principale =====


% ===== findpeaks_safe : utilise findpeaks si dispo, sinon fallback =====
function [pks, locs] = findpeaks_safe(x, varargin)
try
    [pks, locs] = findpeaks(x, varargin{:});
catch
    % Fallback simple: pics stricts + seuil de proéminence (approx.)
    p = inputParser; p.addParameter('MinPeakProminence', 0); p.parse(varargin{:});
    prom = p.Results.MinPeakProminence;
    locs = [];
    for i=2:numel(x)-1
        if x(i)>x(i-1) && x(i)>x(i+1) && x(i)>prom
            locs(end+1)=i; %#ok<AGROW>
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
N = numel(acf_sym);
n = (0:N-1).';
w = 0.5 - 0.5*cos(2*pi*n/(N-1)); % fenetre de Hann
R = acf_sym .* w; % appliques cette fenêtre à l’ACF avant de faire la FFT :
% PSD = FFT(R) réelle (procès réel)
S = fft(R);
Pxx = real(S(1:floor(N/2)+1));
F   = (0:floor(N/2)).' * (Fs/N);
end
