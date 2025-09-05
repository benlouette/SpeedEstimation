function [fr_hat, OUT] = reco_mcu_acf(acf, Fs, P)
% RECO_MCU_ACF — Estimation de la fréquence de rotation depuis UNE ACF (MCU-friendly)
%   • Zéro FFT/log — uniquement ACF + arithmétique légère
%   • Pics ACF -> médiane des premiers intervalles -> parabolique local
%   • Anti-½ complet : compare T, T/2 et 2T (linéaire pour la décision)
%   • Ceintures de sécurité : bornes période/fréquence, clamp final
%
% In :
%   acf : vecteur ACF (lags >= 0), acf(1)=R(0) (pas besoin de normaliser)
%   Fs  : fréquence d’échantillonnage (Hz)
%   P   : struct paramètres (tous optionnels)
%
% Params (défauts pensés “terrain”)
DEF = struct( ...
    'fmin',0.5, ...           % Hz, borne basse plausible
    'fmax',120, ...           % Hz, borne haute plausible
    'maxLagSec',1.0, ...      % s, portée utile max de l’ACF
    'smoothN',5, ...          % taille moyenne mobile (1=off) — léger
    'minProm',0.06, ...       % proéminence relative de base (0.03..0.1)
    'minPkDist',[], ...       % si vide -> 0.8*Fs/fmax (cohérent avec Tmin)
    'useAdaptiveProm',true, ... % proéminence adaptative sur l’ACF hors lag0
    'M_firstIntervals',4, ... % médiane des M premiers intervalles (2..4)
    'parabolic',true, ...     % raffinement sub-échantillon local
    'Kharm',3, ...            % nb. harmoniques pour peigne ACF (2..4)
    'gammaHalf',1.15, ...     % seuil R(T/2) vs R(T)
    'lambdaEven',0.85, ...    % pénalisation des pairs
    'roundMult',true, ...     % nearest-neighbour pour évaluations NON critiques
    'linInterpInDecision',true, ... % linéaire pour anti-½/2T (recommandé)
    'considerDouble',true, ... % teste 2T (pour corriger un harmonique)
    'Plot',false ...          % graphiques de debug
);
if nargin<3, P = struct(); end
fn = fieldnames(DEF);
for k=1:numel(fn), if ~isfield(P,fn{k}), P.(fn{k})=DEF.(fn{k}); end, end

% ---- Préparation & bornes utiles ----
r = acf(:);
R0 = max(r(1), eps);
r  = r / R0;                           % normalisation R(0)=1 (sûre)
N  = numel(r);
tauMin = max(2, floor(Fs/P.fmax));     % >=2 pour parabolique
tauMax = min(N-2, ceil(Fs/P.fmin));
tauMax = min(tauMax, round(P.maxLagSec*Fs));
if tauMin+2 >= tauMax
    fr_hat = 0;
    OUT = struct('reason','invalid_range','tauMin',tauMin,'tauMax',tauMax);
    return;
end

% ---- Lissage très léger (stabilise les pics) ----
if P.smoothN > 1
    k = ones(P.smoothN,1)/P.smoothN;
    r = conv(double(r), k, 'same');
end

% ---- Détection de pics ACF (simple mais robuste) ----
% Proéminence relative de base, avec option d’adaptation (évite micro-pics)
acf_hb = r(2:tauMax+1);
base   = max(acf_hb);
promRel = P.minProm * max(base, eps);
if P.useAdaptiveProm
    medp = median(max(acf_hb,0));
    promRel = max(promRel, 0.5*medp);  % adapt simple
end
% Distance mini entre pics (évite doublons/sous-pics)
if isempty(P.minPkDist)
    minDist = max(2, round(0.8 * Fs / P.fmax));  % ~0.8*Tmin
else
    minDist = P.minPkDist;
end

[pkVal, pkLoc] = findpeaks_simple(r(1:tauMax+1), promRel, minDist);
% retirer lag0 s’il est dedans
mask = pkLoc > 1; pkLoc = pkLoc(mask); pkVal = pkVal(mask);
if isempty(pkLoc)
    fr_hat = 0; OUT = struct('reason','no_peaks'); return;
end

% ---- Période robuste = médiane des premiers intervalles ----
d = diff(pkLoc);
if isempty(d)
    T0 = double(pkLoc(1));       % fallback : premier pic
else
    M  = min(P.M_firstIntervals, numel(d));
    T0 = median(double(d(1:M))); % médiane localisée = robuste aux harmoniques
end

% ---- Clamp période dans [Tmin..Tmax] ----
Tmin = Fs / P.fmax;
Tmax = Fs / P.fmin;
T0   = max(Tmin, min(T0, Tmax));

% ---- Raffinement parabolique local (option) ----
if P.parabolic
    % choisir le pic le plus proche de T0
    [~, iNear] = min(abs(double(pkLoc) - T0));
    L = pkLoc(iNear);
    if L>=2 && L<=tauMax-1
        y1 = r(L-1); y2 = r(L); y3 = r(L+1);
        denom = (y1 - 2*y2 + y3);
        if denom > 0             % convexité OK
            delta = 0.5*(y1 - y3)/denom;
            delta = max(min(delta, 0.5), -0.5);
            T0 = L + delta;
        end
    end
end

% ---- Anti-demi-tour complet : comparer T, T/2, 2T ----
% Pour la décision, on force l’interpolation linéaire (plus fiable)
linInterp = P.linInterpInDecision;
interpR = @(q) interp_r(r, q, ~linInterp);  % ~linInterp => roundMult

Tcand = T0;                % en échantillons (fractionnaire OK)
% clamp T/2 et 2T aussi pour éviter hors-bande
Tcand = max(Tmin, min(Tcand, Tmax));
Tc_H  = max(Tmin, min(Tcand/2, Tmax));   % half
Tc_2  = max(Tmin, min(2*Tcand, Tmax));   % double

% Score pour T
R_T  = interpR(Tcand);
R_T2 = interpR(Tcand/2);
S_T  = peigne_score_acf(r, Tcand,  P.Kharm, P.lambdaEven, linInterp);
% Score pour T/2
R_H  = interpR(Tc_H);
S_H  = peigne_score_acf(r, Tc_H,   P.Kharm, P.lambdaEven, linInterp);
% Score pour 2T (option)
if P.considerDouble
    R_2  = interpR(Tc_2);
    S_2  = peigne_score_acf(r, Tc_2,  P.Kharm, P.lambdaEven, linInterp);
else
    R_2 = -Inf; S_2 = -Inf;
end

% Combinaison de votes (léger, sans log/FFT)
%   ratio : (R(T) vs R(T/2)), peigne : odd/even ACF
score_T  = (R_T  - P.gammaHalf*max(R_T2,0)) + 0.20*(S_T);
score_H  = (R_H  - P.gammaHalf*max(R_T ,0)) + 0.20*(S_H);
score_2T = (R_2  - P.gammaHalf*max(R_T ,0)) + 0.20*(S_2);

[~, which] = max([score_T, score_H, score_2T]);
if which==2
    T_hat_samp = Tc_H;
elseif which==3
    T_hat_samp = Tc_2;
else
    T_hat_samp = Tcand;
end
T_hat_samp = max(Tmin, min(T_hat_samp, Tmax));
T_hat = T_hat_samp / Fs;

% ---- Fréquence brute + clamp final [fmin..fmax] ----
fr_hat = 1 / max(T_hat, eps);
fr_hat = max(P.fmin, min(fr_hat, P.fmax));

% ---- Qualité simple (0..1) pour exploitation aval ----
% q_peak   : R(T)/R(0)   (0..1)
% q_spread : régularité des intervalles (MAD/med) -> 1 - spread
% q_comb   : contraste odd/even map (-1..1)->(0..1)
q_peak = max(0, min(1, interpR(T_hat_samp)));
if numel(d) >= 2
    medd = median(double(d(1:min(numel(d), P.M_firstIntervals))));
    madD = median(abs(double(d(1:min(numel(d), P.M_firstIntervals))) - medd));
    spread = madD / max(medd, eps);
else
    spread = 0.5;  % conservateur
end
q_spread = max(0, 1 - spread/0.35);
combContrast = peigne_contrast_only(r, T_hat_samp, P.Kharm, P.lambdaEven, linInterp);
q_comb = 0.5*(combContrast+1);

q = min(1, max(0, 0.55*q_peak + 0.25*q_spread + 0.20*q_comb));

% ---- Sorties ----
OUT = struct();
OUT.T_hat = T_hat;
OUT.T_samp = T_hat_samp;
OUT.peaks_samples = pkLoc-1;
OUT.peaks_values  = pkVal;
OUT.q_peak = q_peak;
OUT.q_spread = q_spread;
OUT.combContrast = combContrast;
OUT.q = q;
OUT.params = P;

% ---- Plots de debug (optionnels) ----
if P.Plot
    tau = (0:numel(r)-1).'/Fs;
    figure('Name','ACF MCU (durcie)','Color','w');
    plot(tau, r, 'b'); grid on; hold on;
    stem(tau(pkLoc), r(pkLoc), 'filled', 'Color',[0.85 0.33 0.10], 'Marker','none');
    xline(T_hat, '--', sprintf('T=%.4f s (%.2f Hz)', T_hat, fr_hat), 'Color',[0.2 0.2 0.2]);
    xline(T_hat/2, ':', 'T/2');
    xlabel('\tau (s)'); ylabel('ACF norm.');
    title(sprintf('MCU reco — q=%.2f  odd-even=%.2f', q, combContrast));

    % Odd vs Even bar
    figure('Name','Odd vs Even (ACF only)','Color','w');
    [Sodd, Seven] = peigne_odd_even(r, T_hat_samp, P.Kharm, linInterp);
    bar([Sodd, Seven]); grid on; set(gca,'XTickLabel',{'odd','even'});
    ylabel('Somme R(kT)^+'); title(sprintf('contrast=%.3f', combContrast));
end
end

% ================== Helpers ==================

function [pks, locs] = findpeaks_simple(x, prom, minDist)
% Détection de pics locale, proéminence globale + distance mini
N = numel(x); locs = []; pks = [];
lastLoc = -minDist;
for i=2:N-1
    if x(i)>x(i-1) && x(i)>=x(i+1) && x(i)>prom
        if (i - lastLoc) >= minDist
            locs(end+1) = i; %#ok<AGROW>
            pks(end+1)  = x(i);
            lastLoc = i;
        elseif ~isempty(pks) && x(i) > pks(end)
            % remplace le pic précédent si plus fort mais trop proche
            locs(end) = i; pks(end) = x(i); lastLoc = i;
        end
    end
end
locs = locs(:); pks = pks(:);
end

function v = interp_r(r, q, roundMult)
% Évalue r au lag fractionnaire q (en échantillons).
% MCU: roundMult=true => nearest-neighbour (aucune division), sinon linéaire
N = numel(r);
if q < 1 || q > N-2
    v = NaN; return;
end
if roundMult
    idx = round(q)+1; idx = max(1, min(N, idx));
    v = r(idx);
else
    i0 = floor(q); a = q - i0;
    v = (1-a)*r(i0+1) + a*r(i0+2);  % +1 car r(1)=lag0
end
end

function [Sodd, Seven] = peigne_odd_even(r, T_samp, K, linInterp)
% Somme des R(kT)^+ pour k impairs/pairs (ACF only, sans log)
Sodd=0; Seven=0;
for k=1:K
    q   = k*T_samp;
    val = interp_r(r, q, ~linInterp);   % forcer linéaire si linInterp=true
    if ~isfinite(val), continue; end
    val = max(val, 0);
    if mod(k,2)==1, Sodd=Sodd+val; else, Seven=Seven+val; end
end
end

function contrast = peigne_contrast_only(r, T_samp, K, lambdaEven, linInterp)
[Sodd, Seven] = peigne_odd_even(r, T_samp, K, linInterp);
contrast = (Sodd - lambdaEven*Seven) / max(Sodd + lambdaEven*Seven, eps);
end

function score = peigne_score_acf(r, T_samp, K, lambdaEven, linInterp)
% Score léger : même formule que le contraste
score = peigne_contrast_only(r, T_samp, K, lambdaEven, linInterp);
end
