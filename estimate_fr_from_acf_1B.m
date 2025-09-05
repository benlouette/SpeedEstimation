function [fr_hat, OUT] = estimate_fr_from_acf_1B(acf, Fs, P)
% Estimation de la fréquence de rotation depuis l'ACF avec YIN/CMNDF (solution 1B).
% - d(tau) = 2*(R(0) - R(tau)) dérivé de l'ACF
% - CMNDF (YIN) : d'(tau) = d(tau) / ( (1/tau) * sum_{j=1..tau} d(j) )
% - Règle du seuil absolu (0.1–0.2) => première période fiable
% - Raffinement par interpolation parabolique autour du minimum local
% - Anti demi-tour (ACF): compare d'(T) vs d'(T/2), R(T) vs R(T/2),
%   et peigne impairs/pairs sur l'ACF.
%
% In:
%   acf : vecteur ACF (lags >= 0, acf(1) = R(0))
%   Fs  : fréquence d'échantillonnage (Hz)
%   P   : struct de paramètres (voir demo)
%
% Out:
%   fr_hat : fréquence estimée (Hz)
%   OUT    : diagnostics (struct)

acf = acf(:);
N   = numel(acf);
lags = (0:N-1).';
tau  = lags / Fs;

% ---------- paramètres par défaut ----------
def = struct('MaxLagSec',1.0,'fmin',0.5,'fmax',120,'Threshold',0.15, ...
             'SmoothMs',1.0,'Kharm',4,'lambdaEven',0.8,'gammaHalf',1.10,'Plot',true);
fn = fieldnames(def);
for k=1:numel(fn)
    if ~isfield(P, fn{k}), P.(fn{k}) = def.(fn{k}); end
end

% ---------- normalisation ACF ----------
% R0 = acf(1);
% r  = (R0~=0) * (acf / R0) + (R0==0) * acf;
r = acf;

% ---------- bornes de recherche ----------
tauMinSamp = max(2, floor(Fs / P.fmax)); % >= 2 pour parabolique
tauMaxSamp = min(N-2, ceil(Fs / P.fmin));
maxLagSamp = min(N-1, round(P.MaxLagSec*Fs));
tauMaxSamp = min(tauMaxSamp, maxLagSamp);
if tauMinSamp+2 >= tauMaxSamp
    error('Plage de lag invalide: ajuster fmin/fmax/MaxLagSec.');
end

%Construction de la NDF (difference function) à partir de l’ACF
% ---------- YIN / CMNDF à partir de l'ACF ----------
% d(tau) = 2*(R0 - R(tau)) ; avec r=R/R0 => d ~ 2*(1 - r)
tauIdx = (1:tauMaxSamp).'; % tauIdx est le vecteur des lags candidats 𝜏=1,2,…,𝜏max(en échantillons).
d = 2.0 * (1 - r(tauIdx+1)); % ignore le facteur R0 (s'annule dans CMNDF)
% CMNDF : d'(tau) = d(tau) * tau / sum_{j=1..tau} d(j)
cum = cumsum(d); %CMNDF (Cumulative Mean Normalized Difference Function)
dprime = d .* (tauIdx ./ max(cum, eps));

% Lissage léger (en échantillons)
w = max(1, round(P.SmoothMs*1e-3*Fs));
if w>1, dprime_s = movavg(dprime, w); else, dprime_s = dprime; end

% Masque hors bande (protéger contre f trop haut/bas)
mask = true(size(dprime_s));
mask(1:tauMinSamp-1) = false;  % on ne regarde pas < tauMin
% (tauMaxSamp est déjà limite supérieure)

% ---------- règle du seuil absolu ----------
th = P.Threshold;
idxCross = find(dprime_s(mask) < th, 1, 'first');
if ~isempty(idxCross)
    % réindexer vers indices absolus
    firstIdx = find(mask,1,'first') + idxCross - 1;
    % chercher un minimum local autour (fenêtre proportionnelle à tau)
    win = max(3, round(0.1 * firstIdx)); 
    i1 = max(tauMinSamp, firstIdx - win);
    i2 = min(tauMaxSamp, firstIdx + win);
    [~,rel] = min(dprime_s(i1:i2));
    L = i1 + rel - 1;
else
    % pas de crossing: prendre le minimum global dans la bande
    [~, L] = min( replace_inf(~mask, dprime_s, +Inf) );
    if L<tauMinSamp || L>tauMaxSamp
        % secours : contraindre à la bande
        [~, L] = min(dprime_s(tauMinSamp:tauMaxSamp));
        L = L + tauMinSamp - 1;
    end
end

% ---------- interpolation parabolique autour du minimum ----------
if L<2 || L>tauMaxSamp-1
    Lref = L;
else
    y1 = dprime_s(L-1); y2 = dprime_s(L); y3 = dprime_s(L+1);
    delta = 0.5*(y1 - y3) / max(y1 - 2*y2 + y3, eps);
    delta = max(min(delta, 0.5), -0.5);
    Lref = L + delta;
end
T_hat = Lref / Fs;   % secondes

% ---------- anti-demi-tour (ACF) ----------
interp = @(q) safe_interp(lags, r, q);  % interp sur ACF normalisée
T  = T_hat * Fs;   % en échantillons (fractionnaire)
R_T   = interp(T);
R_T2  = interp(T/2);

% Peigne odd/even sur ACF
K = P.Kharm;
S_odd = 0; S_even = 0;
for k=1:K
    q = k*T;
    if q <= tauMaxSamp
        if mod(k,2)==1, S_odd  = S_odd  + interp(q);
        else            S_even = S_even + interp(q);
        end
    end
end
combContrast = (S_odd - P.lambdaEven*S_even) / max(S_odd + P.lambdaEven*S_even, eps);

% Critère demi-tour : d'(T/2) << d'(T) OU (R(T/2) > gamma*R(T) ET pairs>impairs)
dprime_T  = interp_dprime(dprime_s, Lref);
dprime_T2 = interp_dprime(dprime_s, Lref/2);
isHalf = (dprime_T2 < 0.85*dprime_T) || ((R_T2 > P.gammaHalf*R_T) && (S_even > S_odd));

if isHalf
    T_hat = T_hat/2;
    % recompute diagnostics at corrected T
    T  = T_hat * Fs;
    R_T  = interp(T);
    dprime_T = interp_dprime(dprime_s, T);
end

fr_hat = 1 / max(T_hat, eps);

% ---------- sorties & tracés ----------
OUT = struct();
OUT.dprime = dprime_s; OUT.tauIdx = tauIdx; OUT.mask = mask;
OUT.L = L; OUT.Lref = Lref; OUT.T_hat = T_hat;
OUT.dprime_T = dprime_T;
OUT.periodicity = R_T;
OUT.R_T = R_T; OUT.R_T2 = R_T2;
OUT.combContrast = combContrast; OUT.isHalf = isHalf;
OUT.params = P;

if P.Plot
    % d' (YIN)
    figure('Name','YIN / CMNDF (d'')');
    tt = tauIdx/Fs;
    plot(tt, dprime_s, 'b'); grid on; hold on;
    yline(th,'--','Threshold'); 
    xline(T_hat, '--', sprintf('T\\_hat=%.4f s  (f=%.2f Hz)', T_hat, fr_hat), 'Color',[0.2 0.2 0.2]);
    xlabel('\tau (s)'); ylabel('d''(\tau)'); xlim([0, P.MaxLagSec]);
    title('YIN / CMNDF à partir de l''ACF');

    % ACF avec T, T/2
    figure('Name','ACF & périodes');
    plot(tau, r, 'b'); grid on; hold on;
    xline(T_hat, '--', 'T'); xline(T_hat/2, ':', 'T/2');
    xlabel('\tau (s)'); ylabel('ACF (norm.)'); xlim([0, P.MaxLagSec]);
    title('ACF normalisée & marquage périodicités');

    % Peigne odd/even
    figure('Name','Peigne odd/even sur ACF');
    tmarks = (1:K)*T_hat; yodd = nan(size(tmarks)); yeven = nan(size(tmarks));
    for k=1:K
        if mod(k,2)==1, yodd(k) = interp(k*T_hat*Fs);
        else            yeven(k)= interp(k*T_hat*Fs);
        end
    end
    stem(tmarks, yodd, 'g','filled'); hold on; stem(tmarks, yeven,'r','filled'); grid on;
    xlabel('k*T (s)'); ylabel('R(kT)/R(0)'); legend('impairs','pairs','Location','best');
    title(sprintf('Contraste (odd-even)=%.3f  %s', combContrast, tern(isHalf,'→ demi-tour corrigé','')));
end
end

% ===== helpers =====

function y = movavg(x,w)
if w<=1, y=x; return; end
k = ones(w,1)/w;
y = filtfilt(k,1,double(x));
end

function v = safe_interp(x, y, q)
% interpolation linéaire sûre aux indices fractionnaires
q = q(:); v = nan(size(q)); N = numel(x);
for i=1:numel(q)
    if q(i) < 1 || q(i) > N-2
        v(i) = NaN;
    else
        i0 = floor(q(i));
        a  = q(i) - i0;
        v(i) = (1-a)*y(i0+1) + a*y(i0+2);
    end
end
if isscalar(v), v=v(1); end
end

function z = tern(c,a,b)
if c, z=a; else, z=b; end
end

function dp = interp_dprime(dprime, q)
% Interp parabolique de d' au voisinage d'un index fractionnaire q
N = numel(dprime);
if q<2 || q>N-1
    % linéaire si trop près du bord
    i0 = max(1, min(N-1, floor(q)));
    a  = q - i0;
    dp = (1-a)*dprime(i0) + a*dprime(i0+1);
    return;
end
L = round(q);
y1 = dprime(L-1); y2 = dprime(L); y3 = dprime(L+1);
delta = 0.5*(y1 - y3) / max(y1 - 2*y2 + y3, eps);
delta = max(min(delta, 0.5), -0.5);
Lref  = L + delta;
i0 = floor(Lref); a = Lref - i0;
dp = (1-a)*dprime(i0) + a*dprime(i0+1);
end

function y = replace_inf(mask, x, val)
y = x;
y(mask) = val;
end
