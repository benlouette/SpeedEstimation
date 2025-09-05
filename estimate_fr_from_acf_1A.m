function [fr_hat, OUT] = estimate_fr_from_acf_1A(acf, Fs, P)
% Estimation tachless depuis l'ACF (solution 1A):
%  - détection de pics ACF
%  - médiane des intervalles entre pics -> période T
%  - interpolation parabolique locale
%  - anti-demi-tour: tests R(T) vs R(T/2) + peigne odd/even
%
% In:
%   acf  : vecteur ACF pour lags >= 0 (acf(1) = R(0))
%   Fs   : fréquence d'échantillonnage (Hz)
%   P    : struct de paramètres (voir demo)
% Out:
%   fr_hat : estimation de la fréquence (Hz)
%   OUT    : diagnostics (struct)

acf = acf(:);
N   = numel(acf);
lags = (0:N-1).';
tau  = lags / Fs;

% ---------- paramètres par défaut ----------
def = struct('MaxLagSec',1.0,'SmoothMs',2,'Prominence',0.05,'TopPeaks',15, ...
             'fmin',0.5,'fmax',120,'Kharm',4,'lambdaEven',0.8,'gammaHalf',1.15,'Plot',true);
fn = fieldnames(def);
for k=1:numel(fn)
    if ~isfield(P, fn{k}), P.(fn{k}) = def.(fn{k}); end
end

% ---------- normalisation ----------
% R0 = acf(1);
% if R0 ~= 0, r = acf / R0; else, r = acf; end
r = acf;
% ---------- bornes en lag ----------
tauMinSamp = max(2, floor(Fs / P.fmax)); % >= 2 pour parabolique on ne cherche pas de période plus courte que 1/fmax, et jamais en dessous de 2 échantillons.
tauMaxSamp = min(N-2, ceil(Fs / P.fmin));% on ne cherche pas de période plus grande que 1/fmin, ni au-delà de ce que l’ACF permet.
maxLagSamp = min(N-1, round(P.MaxLagSec*Fs));% borne maximale imposée par un paramètre utilisateur (MaxLagSec).
tauMaxSamp = min(tauMaxSamp, maxLagSamp);% la borne maximale effective est la plus restrictive des deux.
if tauMinSamp+2 >= tauMaxSamp
    error('Plage de lag invalide: ajuster fmin/fmax/MaxLagSec.');
end

% ---------- lissage et détection de pics ----------
w = max(1, round(P.SmoothMs*1e-3*Fs));% fenêtre max pour lissage
rs = movavg(r, w); % lissage
% proéminence relative
base = max(rs(2:tauMaxSamp+1));
prom = P.Prominence * max(base, eps);

% impose une distance min entre pics pour éviter les doublons (≈ T/3 à fmax)
minPkDist = max(2, round(0.33 * Fs / P.fmax));

[pkVal, pkLoc] = findpeaks_safe(rs(1:tauMaxSamp+1), ...
    'MinPeakProminence', prom, 'MinPeakDistance', minPkDist);

% retire lag 0 si présent
mask = pkLoc > 1;
pkLoc = pkLoc(mask);
pkVal = pkVal(mask);

% garde au plus TopPeaks (les plus forts)
if numel(pkLoc) > P.TopPeaks
    [~, idx] = maxk(pkVal, P.TopPeaks);
    idx = sort(idx);
    pkLoc = pkLoc(idx);
    pkVal = pkVal(idx);
end

if isempty(pkLoc)
    fr_hat = 0; OUT = struct('reason','no_peaks'); return;
end

% ---------- période robuste: médiane des intervalles ----------
d = diff(pkLoc); % en échantillons
if isempty(d)
    T_samp = double(pkLoc(1));  % si 1 seul pic exploitable
else
    T_samp = median(double(d));
end

% ---------- interpolation parabolique autour du 1er pic près de T ----------
% cherche le pic le plus proche de T_samp
%On choisit le pic le plus proche de T_samp et on récupère sa position entière L.
[~, iNear] = min(abs(double(pkLoc) - T_samp));
L = pkLoc(iNear);
if L<2 || L>tauMaxSamp-1
    % Il faut 3 points consécutifs pour une parabole : L-1, L, L+1.
    % Si L est trop près d’un bord (<2 ou > tauMaxSamp-1), on n’interpole pas → on renvoie une estimation de repli (T_samp/Fs, en secondes).
    % bord → pas d'interpolation possible
    T_hat = T_samp / Fs;
else
    % parabolique pour sous-échantillon
    y1 = rs(L-1); y2 = rs(L); y3 = rs(L+1);
    if y2 < max(y1,y3) || (y1 - 2*y2 + y3) >= 0
        T_hat = T_samp / Fs; % pas d’interpolation si pas concave
    else
        delta = 0.5*(y1 - y3) / max(y1 - 2*y2 + y3, eps);
        delta = max(min(delta, 0.5), -0.5);
        Lref  = L + delta;
        T_hat = Lref / Fs;
    end
end

% ---------- anti-demi-tour : tests R(T) vs R(T/2) + peigne odd/even ----------
% interpolation linéaire de r aux positions fractionnaires
% Si le signal a une vraie période T, on s’attend à ce que les impairs (T, 3T, 5T, …) soient globalement plus forts.
% Si c’est plutôt 𝑇/2, les pairs (2×(T/2)=T, 4×(T/2)=2T, …) dominent.

% P.Kharm : nombre d’harmoniques considérées dans le peigne.
% Trop petit → décision bruitée ; trop grand → sensible aux pics lointains/noise.
% Souvent 4–8 est un bon compromis.
% P.lambdaEven : pondération des pairs.
% Si tes signaux ont naturellement des pairs plus faibles/forts, ajuste pour équilibrer.
% P.gammaHalf : seuil “à partir de quand” T/2 est crédible.Plus bas → bascule plus facilement à T/2.Plus haut → plus conservateur (reste à T).
% tauMaxSamp : borne max en lag (évite de sommer hors zone fiable de l’ACF, cf. troncage/fenêtrage).

interp = @(q) safe_interp(lags, rs, q);
T  = T_hat * Fs;   % en échantillons (fractionnaire)
R_T   = interp(T); % force du pic à la période estimée.
R_T2  = interp(T/2); % “tentation demi-tour”.
R_2T  = interp(2*T); %utile pour vérifier la cohérence (harmoniques).

% Peigne odd/even
K = P.Kharm;
S_odd = 0; S_even = 0;
for k=1:K
    q = k*T;
    if q<=tauMaxSamp
        if mod(k,2)==1, S_odd  = S_odd  + interp(q);
        else             S_even = S_even + interp(q);
        end
    end
end
combContrast = (S_odd - P.lambdaEven*S_even) / max(S_odd + P.lambdaEven*S_even, eps);

% Décision demi-tour
isHalf = (R_T2 > P.gammaHalf*R_T) && (S_even > S_odd);
if isHalf
    T_hat = T_hat/2;
end

fr_hat = 1 / max(T_hat, eps);

% ---------- sorties & tracés ----------
OUT = struct();
OUT.peaks_samples = pkLoc-1;
OUT.peaks_time_s  = (pkLoc-1)/Fs;
OUT.peaks_values  = pkVal;
OUT.periodicity   = interp(T_hat*Fs); % ≈ R(T)/R(0)
OUT.combContrast  = combContrast;
OUT.R_T  = R_T; OUT.R_T2 = R_T2; OUT.R_2T = R_2T;
OUT.isHalf = isHalf;
OUT.params = P;

if P.Plot
    figure('Name','ACF 1A');
    plot(tau, r, 'b'); grid on; hold on;
    stem(tau(pkLoc), r(pkLoc), 'filled', 'Color',[0.85 0.33 0.10]);
    xline(T_hat, '--', sprintf('T\\_hat=%.4f s  (f=%.2f Hz)', T_hat, fr_hat), 'Color',[0.2 0.2 0.2]);
    if ~isHalf, xline(T_hat/2, ':', 'T/2'); else, xline(T_hat*2, ':', '2T'); end
    xlabel('\tau (s)'); ylabel('ACF (norm.)'); xlim([0, P.MaxLagSec]);
    title('ACF normalisée, pics et T estimé');

    % Visualisation peigne odd/even
    figure('Name','Peigne odd/even sur ACF');
    tmarks = (1:K)*T_hat;
    yodd = nan(size(tmarks)); yeven = nan(size(tmarks));
    for k=1:K
        if mod(k,2)==1, yodd(k) = interp(k*T_hat*Fs);
        else            yeven(k)= interp(k*T_hat*Fs);
        end
    end
    stem(tmarks, yodd, 'g','filled'); hold on; stem(tmarks, yeven,'r','filled'); grid on;
    xlabel('k*T (s)'); ylabel('R(kT)/R(0)'); legend('impairs','pairs');
    title(sprintf('Contraste peigne (odd-even) = %.2f  %s', ...
        combContrast, tern(isHalf,'→ demi-tour détecté','')));
end
end

% ===== helpers =====

function y = movavg(x,w)
if w<=1, y=x; return; end
k = ones(w,1)/w;
y = filtfilt(k,1,double(x));
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

function v = safe_interp(x, y, q)
% interpolation linéaire sûre aux indices fractionnaires (en échantillons)
% x : 0..N-1, y : même taille, q : scalaire (ou vecteur) en "samples"
q = q(:);
v = nan(size(q));
N = numel(x);
for i=1:numel(q)
    if q(i) < 1 || q(i) > N-2
        v(i) = NaN; %#ok<*AGROW>
    else
        i0 = floor(q(i));
        a  = q(i) - i0;
        v(i) = (1-a)*y(i0+1) + a*y(i0+2);
    end
end
if isscalar(v), v=v(1); end
end

function s = tern(cond, a, b)
if cond, s=a; else, s=b; end
end
