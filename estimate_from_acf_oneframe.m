function [f_hat, f_yin, f_peigne, OUT] = estimate_from_acf_oneframe(acf, Fs, P)
% ESTIMATE_FROM_ACF_ONEFRAME
% Estimation robuste de la fondamentale à partir d’UNE ACF (lags >= 0).
%   • Mesure 1 : YIN/CMNDF appliqué à l’ACF (on cherche un MIN de d'(τ))
%   • Mesure 2 : Peigne log (score multi-harmonique ; on cherche un MAX)
%   • Anti-demi-tour YIN : test local T vs 2T sur CMNDF + ACF
%   • Recoupement YIN ↔ Peigne : correction octave/½-octave selon cohérence
%
% SORTIE :
%   f_hat    : fréquence finale (Hz) — par défaut, la version PEIGNE corrigée
%   f_yin    : fréquence estimée par YIN/CMNDF (Hz) après corrections
%   f_peigne : fréquence estimée par peigne log (Hz) après corrections
%   OUT      : diagnostics complets (structures/curves/indices)
%
% ENTREE :
%   acf : vecteur ACF d’une frame, lags >= 0, avec acf(1) = R(0) (lag 0)
%   Fs  : fréquence d’échantillonnage du signal d’origine [Hz]
%   P   : (optionnel) params. Tous ont des défauts raisonnables :
%         Bande/échelle :
%           .fmin (0.5) .fmax (120) .maxLagSec (1.0)
%         Lissages / robustesse :
%           .smoothMs (1.0)            % lissage léger (ACF)
%           .yinThresh (0.15)          % seuil absolu YIN
%         Peigne :
%           .Kharm (5)                 % # harmoniques max
%           .WeightMode ('1/k'|'equal')
%           .Eps (1e-3)                % log(eps + r+)
%           .gridStepSamp (1)          % pas de grille (échantillons)
%         Corrections :
%           .useAntiHalfYin (true)     % T vs 2T autour du choix YIN
%           .alphaYin (0.85)           % d'(2T) < alpha*d'(T) -> 2T meilleur
%           .betaYin  (1.10)           % R(2T)  > beta*R(T)
%           .useCrossWithComb (true)   % recouper YIN avec le peigne
%           .crossTolMain (0.20)       % tolérance relative L_yin ~ L_comb
%           .crossTolHalf (0.15)       % tol. pour décider ×2 / ÷2
%
% RECOMMANDE :
%   — Laisse les défauts, commence par comparer f_yin vs f_peigne vs f_hat.
%   — Si YIN tombe souvent à F/2 : baisse .yinThresh (p.ex. 0.10–0.12),
%     garde .useAntiHalfYin = true, et .useCrossWithComb = true.

% ------------------ Défauts ------------------
DEF = struct( ...
    'fmin',0.5,'fmax',120,'maxLagSec',1.0, ...
    'smoothMs',1.0,'yinThresh',0.15, ...
    'Kharm',5,'WeightMode','1/k','Eps',1e-3,'gridStepSamp',1, ...
    'useAntiHalfYin',true,'alphaYin',0.85,'betaYin',1.10, ...
    'useCrossWithComb',true,'crossTolMain',0.20,'crossTolHalf',0.15);
if nargin<3, P = struct(); end
fn = fieldnames(DEF);
for k=1:numel(fn), if ~isfield(P,fn{k}), P.(fn{k}) = DEF.(fn{k}); end, end

% ------------------ Mise en forme ------------------
r = acf(:);                  % colonne
Nlag = numel(r);
if Nlag < 5, error('ACF trop courte.'); end
if r(1) ~= 0, r = r / r(1); end  % normalise : r(0)=1 (r(1) MATLAB)

% Lissage doux (stabilise évaluations à kT)
w = max(1, round(P.smoothMs * 1e-3 * Fs));
r_eval = (w>1) * movavg(r,w) + (w<=1) * r;

% Bornes lags utiles (en échantillons)
tauMin = max(2, floor(Fs / P.fmax));           % >=2 pour parabolique
tauMax = min(Nlag-2, ceil(Fs / P.fmin));
tauMax = min(tauMax, min(Nlag-1, round(P.maxLagSec*Fs)));
if tauMin+2 >= tauMax
    [f_hat,f_yin,f_peigne,OUT] = deal(0,0,0, struct('reason','invalid_range'));
    return;
end

% =====================================================================
% 1) YIN/CMNDF sur ACF — on cherche un MINIMUM de d'(τ)
% =====================================================================
% d(τ) ≈ 2*(1 - r(τ)), pour τ = 1..tauMax (note : r_eval(1)=lag 0)
d = 2 * (1 - r_eval(2:tauMax+1));
cum = cumsum(d);
tIdx = (1:tauMax).';
dprime = d .* (tIdx ./ max(cum, eps));  % CMNDF
% masque < tauMin
dprime(1:tauMin-1) = +Inf;

% seuil absolu YIN : premier crossing
idx = find(dprime < P.yinThresh, 1, 'first');
if ~isempty(idx)
    Lc  = idx;
    win = max(3, round(0.1 * Lc));      % ~10% de τ
    i1  = max(tauMin, Lc - win);
    i2  = min(tauMax, Lc + win);
    [~,rel] = min(dprime(i1:i2)); L = i1 + rel - 1;
else
    [~,L] = min(dprime);
    if L<tauMin || L>tauMax, [~,L]=min(dprime(tauMin:tauMax)); L=L+tauMin-1; end
end

% interpolation parabolique autour du minimum (sub-échantillon)
if L>=2 && L<=tauMax-1
    y1=dprime(L-1); y2=dprime(L); y3=dprime(L+1);
    delta = 0.5*(y1 - y3) / max(y1 - 2*y2 + y3, eps);
    delta = max(min(delta, 0.5), -0.5);
    Lref  = L + delta;
else
    Lref = L;
end

% -------- Anti-demi-tour LOCAL pour YIN (option) --------
% idée : si 2T a un d' bien meilleur ET une corrélation R plus forte, doubler.
if P.useAntiHalfYin
    L1 = Lref;
    L2 = 2*L1;
    if L2 <= tauMax
        d1 = interp_lin(dprime, L1);
        d2 = interp_lin(dprime, L2);
        R  = @(q) safe_interp((0:Nlag-1).', r_eval, q);
        R1 = R(L1);
        R2 = R(L2);
        if (d2 < P.alphaYin * d1) && (R2 > P.betaYin * R1)
            Lref = L2;   % corrige YIN: T/2 -> T
        end
    end
end
T_yin = Lref / Fs; f_yin = 1 / max(T_yin, eps);
nsdf_like = max(0, min(1, 1 - interp_lin(dprime, Lref)/2));  % indicateur qualité

% =====================================================================
% 2) Peigne log — on cherche un MAXIMUM de Score(L)
% =====================================================================
gridStep = max(1, round(P.gridStepSamp));
tauGrid  = (tauMin:gridStep:tauMax).';   Sg = numel(tauGrid);

% Poids harmoniques (normalisés)
K = P.Kharm;
wk = strcmpi(P.WeightMode,'equal') * ones(1,K) + ...
     ~strcmpi(P.WeightMode,'equal') * (1./(1:K));
wk = wk / sum(wk);

rpos   = max(r_eval, 0);               % peigne “positif” (enveloppe/AM)
logeps = @(x) log(P.Eps + x);
Score  = -inf(Sg,1);

for s=1:Sg
    Lg   = tauGrid(s);
    kmax = min(K, floor(tauMax / Lg));     % nb d'harmoniques disponibles
    if kmax < 2, continue; end             % demande >=2 harmoniques
    used = 1:kmax;
    wloc = wk(used); wloc = wloc / sum(wloc);    % renormalise localement
    idx  = used .* Lg + 1;                       % +1 -> lag0
    vals = rpos(idx);
    Score(s) = dot(wloc(:), logeps(vals(:)));
end

% max du score + raffinement parabolique
[~, imax] = max(Score);
L0 = tauGrid(imax);
if imax>1 && imax<Sg
    y1=Score(imax-1); y2=Score(imax); y3=Score(imax+1);
    delta = 0.5*(y1 - y3) / max(y1 - 2*y2 + y3, eps);
    delta = max(min(delta, 0.5), -0.5);
    L_comb = L0 + delta*gridStep;
else
    L_comb = L0;
end
T_comb = L_comb / Fs; f_peigne = 1 / max(T_comb, eps);

% -------- Recoupement YIN ↔ Peigne (option) --------
% règle : si L_yin est “à la moitié” de L_comb, doubler L_yin ; si à “deux fois”,
% diviser ; sinon on laisse tel quel. On corrige ensuite le PEIGNE si besoin.
if P.useCrossWithComb
    rel = @(a,b) abs(a-b)/b;
    % tente de corriger YIN d'abord (si écart significatif)
    if rel(Lref, L_comb) > P.crossTolMain
        if rel(2*Lref, L_comb) < P.crossTolHalf
            Lref = 2*Lref; T_yin = Lref/Fs; f_yin = 1/max(T_yin,eps);
        elseif rel(Lref, 2*L_comb) < P.crossTolHalf
            Lref = Lref/2;  T_yin = Lref/Fs; f_yin = 1/max(T_yin,eps);
        end
    end
    % possibilité inverse : si le peigne est manifestement “à ×2 ou ÷2” de YIN
    if rel(L_comb, Lref) > P.crossTolMain
        if rel(2*L_comb, Lref) < P.crossTolHalf
            L_comb = 2*L_comb; T_comb = L_comb/Fs; f_peigne = 1/max(T_comb,eps);
        elseif rel(L_comb, 2*Lref) < P.crossTolHalf
            L_comb = L_comb/2; T_comb = L_comb/Fs; f_peigne = 1/max(T_comb,eps);
        end
    end
end

% -------- Choix final --------
% Stratégie simple : peigne (souvent le plus robuste aux ambiguïtés d’octave).
f_hat = f_peigne;

% ------------------ Diagnostics ------------------
OUT = struct();
OUT.params      = P;
OUT.r_eval      = r_eval;
OUT.dprime      = dprime;
OUT.nsdf_like   = nsdf_like;

OUT.tauMin      = tauMin;
OUT.tauMax      = tauMax;
OUT.grid        = tauGrid;
OUT.Score       = Score;

OUT.L_yin       = Lref;
OUT.T_yin       = T_yin;

OUT.L_peigne    = L_comb;
OUT.T_peigne    = T_comb;
end

% ===================== Helpers =====================
function y = movavg(x,w)
% Moyenne glissante (zéro-phase si filtfilt est dispo)
if w<=1, y=x; return; end
k = ones(w,1)/w;
try,    y = filtfilt(k,1,double(x));
catch,  y = conv(double(x),k,'same');
end
y = cast(y, 'like', x);
end

function v = interp_lin(x, q)
% Interpolation linéaire sur indices 1..N (q réel)
N = numel(x);
if q<=1, v=x(1);  return; end
if q>=N, v=x(end); return; end
i0 = floor(q); a = q - i0;
v = (1-a)*x(i0) + a*x(i0+1);
end

function v = safe_interp(lags, y, q)
% Interp. linéaire sur axe "lags = 0..N-1" (q en échantillons)
N = numel(y);
if q<=0, v=y(1); return; end
if q>=N-1, v=y(end); return; end
i0 = floor(q); a = q - i0;
v  = (1-a)*y(i0+1) + a*y(i0+2);   % +1 car y(1) = lag0
end
