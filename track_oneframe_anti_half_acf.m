function [f_final, OUT] = track_oneframe_anti_half_acf(acf, Fs, P)
% TRACK_ONEFRAME_ANTI_HALF_ACF
% Estimation robuste de la fréquence fondamentale à partir d'UNE ACF.
%   1) Génère plusieurs candidats (YIN/CMNDF depuis ACF, NSDF≈ACF, peigne log ACF,
%      HPS sur PSD issue de l'ACF) — indépendants et complémentaires.
%   2) Fait une FUSION robuste par médiane sur log-fréquence (consensus).
%   3) Applique un ANTI-DEMI-TOUR (T vs T/2) basé sur 3 votes : peigne ACF
%      odd/even, peigne spectral odd/even, et ratio R(T) vs R(T/2).
%   4) Calcule une QUALITÉ de frame (0..1) à partir de plusieurs indicateurs,
%      utile pour filtrer/rejeter/pondérer en aval (tracking, Kalman, etc.).
%
% Entrées
%   acf : vecteur ACF (lags >= 0). Pas besoin d’être normalisée, on fera R(0)=1.
%   Fs  : fréquence d’échantillonnage (Hz).
%   P   : struct params (tous optionnels) — voir “Défauts” ci-dessous.
%
% Sorties
%   f_final : fréquence estimée (Hz), après fusion et anti-demi-tour.
%   OUT     : struct diagnostics (tous les candidats, scores, qualités, etc.).
%
% Recommandations pratiques
%   • Laisse les défauts, puis ajuste :
%       - P.fmin / P.fmax  : bande plausible (ex. [0.5..120] Hz).
%       - P.SmoothMs       : 0.5–2 ms pour stabiliser l’ACF.
%       - P.NSDF_min       : 0.5–0.8 selon SNR (seuil pic NSDF).
%       - P.Comb_K         : 4–8 harmoniques; WeightMode='1/k' conseillé.
%       - P.SpecBW_Hz      : largeur d’énergie autour des harmoniques (0.5–2 Hz).
%       - Pondérations anti-½ : P.wComb/P.wSpec/P.wRratio (ex. 0.45/0.40/0.15).
%   • Si tu travailles sur une enveloppe (rectif/Hilbert+LPF+decim), garde
%     SmoothMs modéré et fmin/fmax adaptés à la bande.

% -------------------- Paramètres par défaut --------------------
DEF = struct( ...
    'fmin',0.5,'fmax',120,'maxLagSec',1.0, ...
    'SmoothMs',1.0, ...                % lissage léger de l'ACF (ms)
    'YIN_thresh',0.15, ...             % seuil YIN si tu veux l'utiliser ailleurs
    'NSDF_min',0.55, ...               % seuil pic NSDF≈ACF
    'Comb_K',6,'Comb_eps',1e-3,'Comb_weight','1/k', ...
    'HPS_use',true,'HPS_K',4,'HPS_dfHz',0.1, ...
    'lambdaEven',0.8,'gammaHalf',1.12,'UseLogComb',true,'SpecBW_Hz',1.0, ...
    'wComb',0.45,'wSpec',0.40,'wRratio',0.15, ...  % pondérations votes anti-½
    'considerDouble',false, ...         % optionnel : tester aussi 2T
    'Plot',false ...                    % true : montre quelques graphiques
);
if nargin<3, P = struct(); end
fn = fieldnames(DEF);
for k=1:numel(fn), if ~isfield(P,fn{k}), P.(fn{k}) = DEF.(fn{k}); end, end

% -------------------- Préparation ACF --------------------
r = acf(:);
if r(1) ~= 0, r = r / r(1); end               % normalise : R(0)=1
Nlag = numel(r);
% lissage douceur (stabilise les évaluations aux multiples kT)
w = max(1, round(P.SmoothMs*1e-3*Fs));
if w>1, r2 = movavg(r, w); else, r2 = r; end

% bornes lags utilisables (échantillons)
tauMin = max(2, floor(Fs/P.fmax));            % >=2 pour interpolation parabolique
tauMax = min(Nlag-2, ceil(Fs/P.fmin));
tauMax = min(tauMax, min(Nlag-1, round(P.maxLagSec*Fs)));
if tauMin+2 >= tauMax
    f_final = 0;
    OUT = struct('reason','invalid_range','tauMin',tauMin,'tauMax',tauMax);
    return;
end

% -------------------- Poids harmoniques --------------------
Kc    = P.Comb_K;
wkComb= weights_K(Kc, P.Comb_weight);
Kh    = P.HPS_K;

% -------------------- (1) Candidats multiples --------------------
% 1A-bis : YIN/CMNDF-from-ACF (on renvoie la valeur min de d′, utile comme “profondeur”)
[fyin, yin_dp] = cand_yin_acf(r2, Fs, tauMin, tauMax);
% 1C-bis : NSDF approx (ici nsdf≈r2), 1er pic > P.NSDF_min
[fnsdf, nsdf_pk] = cand_nsdf_acf(r2, Fs, tauMin, tauMax, P.NSDF_min);
% 1E : peigne log ACF
fcomb = cand_comb_log_acf(r2, Fs, tauMin, tauMax, Kc, wkComb, P.Comb_eps);
% HPS (PSD via ACF -> Harmonic Sum)
if P.HPS_use
    fhps = cand_hps_from_psd(r2, Fs, P.fmin, P.fmax, P.HPS_dfHz, Kh);
else
    fhps = NaN;
end

% Collecte candidates valides
fcands = [fyin, fnsdf, fcomb, fhps];
fcands = fcands(isfinite(fcands) & fcands>0);

if isempty(fcands)
    f_final = 0;
    OUT = struct('reason','no_candidates','fyin',fyin,'fnsdf',fnsdf,'fcomb',fcomb,'fhps',fhps);
    return;
end

% -------------------- (2) Fusion robuste (consensus) --------------------
% médiane sur log-f : robuste aux outliers multiplicatifs / erreurs d’octave
lf = log(fcands);
f_cons = exp(median(lf));
% cohérence des méthodes (plus c’est petit, mieux elles s’accordent)
disp_log = mad(lf,1);             % MAD sur log-f (échelle robuste)
q_cons   = max(0, 1 - disp_log/0.25);   % heuristique ~1 si accord ±25%

% -------------------- (3) Anti-demi-tour multi-votes --------------------
Tcand = 1/f_cons;
f_corr = anti_half_decide(r2, Fs, Tcand, P);

% -------------------- (4) Qualité de frame (0..1) --------------------
%   - comb_contrast : peigne ACF odd/even (contraste → périodicité harmonique)
%   - q_yin         : “profondeur” YIN (d′ petit → bonne périodicité)
%   - q_nsdf        : hauteur du pic NSDF (≈ r) → [0..1]
comb_contrast = odd_even_contrast_acf(r2, Fs, f_corr, Kc, wkComb, P.lambdaEven, P.UseLogComb);
q_comb = (comb_contrast+1)/2;               % map [-1..1] -> [0..1]
q_yin  = max(0, 1 - yin_dp/0.35);           % d′ ~ [0..1+], 0.35 ≈ heuristique
q_nsdf = max(0, min(1, nsdf_pk));           % déjà ~[0..1]
% qualité finale (poids heuristiques; ajuste selon ton domaine)
q = min(1, max(0, 0.25*q_cons + 0.35*q_yin + 0.25*q_nsdf + 0.15*q_comb));

% -------------------- Sortie --------------------
f_final = f_corr;
OUT = struct();
OUT.fyin = fyin; OUT.fnsdf = fnsdf; OUT.fcomb = fcomb; OUT.fhps = fhps;
OUT.fcands = fcands; OUT.f_cons = f_cons; OUT.disp_log = disp_log; OUT.q_cons = q_cons;
OUT.f_final = f_final; OUT.q = q; OUT.yin_dp = yin_dp; OUT.nsdf_pk = nsdf_pk; OUT.comb_contrast = comb_contrast;
OUT.params = P;

% -------------------- Plots optionnels --------------------
if P.Plot
    % 1) ACF + positions T, T/2
    tau = (0:Nlag-1)'/Fs;
    figure('Name','ACF (one frame)','Color','w');
    plot(tau, r2, 'b-','LineWidth',1.2); grid on; hold on;
    xline(1/f_cons, 'k--', sprintf('Cons ~ %.2f Hz', f_cons));
    xline(1/f_final,'m-', sprintf('Final %.2f Hz', f_final), 'LineWidth',2);
    xline(0.5/f_final,'m:', 'T/2');
    xlabel('\tau (s)'); ylabel('r(\tau) (lissée)'); title('ACF + consensus + final');

    % 2) Barres candidats & quality
    figure('Name','Candidats & Qualité','Color','w');
    vals = [fyin, fnsdf, fcomb, fhps, f_cons, f_final];
    names= {'YIN','NSDF','Comb','HPS','Cons','Final'};
    bar(vals); grid on; set(gca,'XTick',1:numel(vals),'XTickLabel',names);
    ylabel('Hz'); title(sprintf('q=%.2f  (YIN dp=%.2f, NSDF pk=%.2f, comb=%.2f)', ...
          q, yin_dp, nsdf_pk, comb_contrast));
end
end

% ================== CANDIDATS ==================
function [fhat, dp_min] = cand_yin_acf(r, Fs, tauMin, tauMax)
% YIN/CMNDF appliqué à l'ACF : d(τ)=2(1-r), d'(τ)=d(τ)*τ / cumsum(d)
d = 2*(1 - r(2:tauMax+1));
cum = cumsum(d); idx = (1:tauMax).';
dprime = d .* (idx ./ max(cum, eps));
dprime(1:tauMin-1) = +Inf;
[dp_min, L] = min(dprime);
% interpolation parabolique autour du minimum
if L>=2 && L<=tauMax-1
    y1=dprime(L-1); y2=dprime(L); y3=dprime(L+1);
    delta = 0.5*(y1-y3)/max(y1-2*y2+y3, eps); delta = max(min(delta,0.5),-0.5);
else
    delta = 0;
end
Lref = L + delta;
fhat = Fs / Lref;
end

function [fhat, pk] = cand_nsdf_acf(r, Fs, tauMin, tauMax, minPk)
% NSDF approximée par l'ACF normalisée (ici r). On cherche le 1er pic > minPk.
nsdf = r;
search = nsdf(1:tauMax+1); search(1:tauMin-1) = -Inf;
[pks, locs] = findpeaks_safe(search, 'MinPeakHeight', minPk);
if isempty(locs)
    [pks, locs] = max(search);
else
    pks = pks(1); locs = locs(1);
end
% raffinement local + parabole
guard = max(2, round(0.1*locs));
i1 = max(2, locs-guard); i2 = min(tauMax-1, locs+guard);
[~,rel] = max(nsdf(i1:i2)); L = i1+rel-1;
if L>=2 && L<=tauMax-1
    y1=nsdf(L-1); y2=nsdf(L); y3=nsdf(L+1);
    delta = 0.5*(y1-y3)/max(y1-2*y2+y3, eps); delta = max(min(delta,0.5),-0.5);
else
    delta = 0;
end
Lref = L + delta; fhat = Fs / Lref; pk = pks;
end

function fhat = cand_comb_log_acf(r, Fs, tauMin, tauMax, K, wk, epslog)
% Peigne log uniquement sur la partie positive de r : score(L)=Σ wk log(eps+r(kL))
rpos = max(r, 0);
S = -Inf(tauMax,1);
for L=tauMin:tauMax
    s=0; cnt=0;
    for k=1:K
        q = k*L; if q>tauMax, break; end
        s = s + wk(k)*log(epslog + rpos(q+1));  % +1 car r(1)=lag0
        cnt = cnt + 1;
    end
    if cnt>=2, S(L)=s; end
end
[~,L0] = max(S);
% raffinement parabolique (si concave)
if L0>tauMin && L0<tauMax
    y1=S(L0-1); y2=S(L0); y3=S(L0+1);
    if isfinite(y1)&&isfinite(y2)&&isfinite(y3) && (y1-2*y2+y3)<0
        delta = 0.5*(y1-y3)/max(y1-2*y2+y3, eps); delta = max(min(delta,0.5),-0.5);
    else
        delta = 0;
    end
else
    delta = 0;
end
Lref = L0 + delta; fhat = Fs / Lref;
end

function fhat = cand_hps_from_psd(r, Fs, fmin, fmax, dfHz, K)
% PSD via Wiener–Khinchin (ACF -> FFT), puis Harmonic Sum sur une grille f0.
[F,Pxx] = psd_from_acf(r, Fs);
fgrid = (fmin:dfHz:fmax).';
Spec = @(f) interp1(F, Pxx, f, 'linear', 0);
S = zeros(size(fgrid));
for i=1:numel(fgrid)
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
end

% ================== Décision anti-½ ==================
function f_corr = anti_half_decide(r, Fs, T_cand, P)
% Fusion 3 votes : peigne ACF (odd/even), peigne spectral (odd/even), ratio R(T) vs R(T/2)
K = P.Comb_K; wk = weights_K(K, P.Comb_weight);
tauMax = numel(r)-1;  interpR = @(q) safe_interp((0:tauMax).', r, q);

% (1) Peigne ACF
[sComb_T, ~, ~] = comb_score(r, T_cand*Fs, K, wk, P.lambdaEven, P.UseLogComb);
[sComb_H, ~, ~] = comb_score(r, (T_cand*Fs)/2, K, wk, P.lambdaEven, P.UseLogComb);

% (2) Peigne spectral
[Fspec, Pxx] = psd_from_acf(r, Fs);
BW = max(P.SpecBW_Hz, Fspec(2)-Fspec(1)); 
Spec = @(f0) band_energy(Fspec, Pxx, f0, BW);
f0 = 1/T_cand; fH = 2/T_cand;
sSpec_T = odd_even_spec(Spec, f0, K, wk, P.lambdaEven);
sSpec_H = odd_even_spec(Spec, fH, K, wk, P.lambdaEven);

% (3) Ratio ACF
RT = interpR(T_cand*Fs); RH = interpR((T_cand*Fs)/2);
sRatio_T = (RT - P.gammaHalf*RH);
sRatio_H = (RH - P.gammaHalf*RT);
nr = max(abs([sRatio_T sRatio_H])) + eps; sRatio_T = sRatio_T/nr; sRatio_H = sRatio_H/nr;

% Score global & décision
S_T = P.wComb*sComb_T + P.wSpec*sSpec_T + P.wRratio*sRatio_T;
S_H = P.wComb*sComb_H + P.wSpec*sSpec_H + P.wRratio*sRatio_H;
if S_H > S_T
    f_corr = fH;   % utiliser T/2 -> double fréquence
else
    f_corr = f0;   % garder T
end
end

% ================== Scores / utilitaires ==================
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

function c = odd_even_contrast_acf(r, Fs, f0, K, wk, lambdaEven, useLog)
[score,~,~] = comb_score(r, (Fs/f0), K, wk, lambdaEven, useLog);
c = score;
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

function w = weights_K(K, mode)
switch lower(string(mode))
    case "equal", w = ones(1,K);
    otherwise,    w = 1./(1:K);
end
w = w / sum(w);
end

function [F, Pxx] = psd_from_acf(r, Fs)
% PSD via Wiener–Khinchin (ACF -> FFT), ACF supposée normalisée
rsym = [flipud(r(2:end)); r(:)];
N = numel(rsym);
win = 0.5 - 0.5*cos(2*pi*(0:N-1)'/(N-1)); % Hann
S = fft(rsym .* win);
P = real(S(1:floor(N/2)+1));
F = (0:floor(N/2))' * (Fs/N);
Pxx = max(P, 0);
end

function E = band_energy(F, Pxx, f0, bw)
if f0<=0 || f0>F(end), E=0; return; end
m = (F >= (f0-bw)) & (F <= (f0+bw));
if ~any(m), E=0; return; end
E = trapz(F(m), Pxx(m));
end

function y = movavg(x,w)
if w<=1, y=x; return; end
k = ones(w,1)/w;
try, y = filtfilt(k,1,double(x));
catch, y = conv(double(x), k, 'same'); end
y = cast(y,'like',x);
end

function [pks, locs] = findpeaks_safe(x, varargin)
try
    [pks, locs] = findpeaks(x, varargin{:});
catch
    % fallback minimaliste
    p = inputParser; p.addParameter('MinPeakHeight', -Inf); p.parse(varargin{:});
    h = p.Results.MinPeakHeight;
    locs = [];
    for i=2:numel(x)-1
        if x(i)>x(i-1) && x(i)>x(i+1) && x(i)>h
            locs(end+1)=i; %#ok<AGROW>
        end
    end
    pks = x(locs);
end
end

function v = safe_interp(x, y, q)
% Interp linéaire sur axe lags (x = 0..N-1), q en échantillons (réel)
q = q(:); v = nan(size(q)); N = numel(x);
for i=1:numel(q)
    if q(i) < 1 || q(i) > N-2
        v(i) = NaN;
    else
        i0 = floor(q(i)); a = q(i) - i0;
        v(i) = (1-a)*y(i0+1) + a*y(i0+2);  % +1 car y(1)=lag0
    end
end
if isscalar(v), v=v(1); end
end
