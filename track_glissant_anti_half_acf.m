function [tFrames, f_meas, f_kal, FR] = track_glissant_anti_half_acf(acf_frames, Fs, tFrames, P)
% TRACK_GLISSANT_ANTI_HALF_ACF
%  Entrée : acf_frames [Nlag x Nf] (ACF normalisée NON requise), Fs, tFrames (1xNf)
%  Sorties : temps, mesure (fusion + anti-½), Kalman, diagnostics par frame

acf_frames = double(acf_frames);
if isvector(acf_frames), acf_frames = acf_frames(:); end
[Nlag, Nf] = size(acf_frames);

if nargin<3 || isempty(tFrames)
    tFrames = 0:(Nf-1);
else
    tFrames = tFrames(:).';
end

% ----- Defaults -----
D = struct('fmin',0.5,'fmax',120,'maxLagSec',1.0, ...
           'SmoothMs',1.0,'YIN_thresh',0.15,'NSDF_min',0.55, ...
           'Comb_K',6,'Comb_eps',1e-3,'Comb_weight','1/k', ...
           'HPS_use',true,'HPS_K',4,'HPS_dfHz',0.1, ...
           'lambdaEven',0.8,'gammaHalf',1.12,'UseLogComb',true,'SpecBW_Hz',1.0, ...
           'wComb',0.45,'wSpec',0.40,'wRratio',0.15,'considerDouble',false, ...
           'kal_dt',[],'kal_sigma_fdot',2.0,'kal_sigma_meas',0.8,'kal_adapt_from_quality',true, ...
           'Plot',true);
fn = fieldnames(D); for k=1:numel(fn), if ~isfield(P,fn{k}), P.(fn{k}) = D.(fn{k}); end, end

% ----- Pré-calculs -----
tauMin = max(2, floor(Fs/P.fmax));
tauMax = min(Nlag-2, ceil(Fs/P.fmin));
maxLagSamp = min(Nlag-1, round(P.maxLagSec*Fs));
tauMax = min(tauMax, maxLagSamp);
if tauMin+2 >= tauMax
    error('Plage lag invalide: ajuste fmin/fmax/maxLagSec vs Nlag.');
end

% Poids harmoniques
Kc = P.Comb_K;
wkComb = weights_K(Kc, P.Comb_weight);
Kh = P.HPS_K;
wkHps = weights_K(Kh, '1/k');  %#ok<NASGU>

% Sorties
f_meas = nan(1,Nf);
FR = repmat(struct('fyin',nan,'fnsdf',nan,'fcomb',nan,'fhps',nan, ...
                   'f_cons',nan,'f_final',nan,'q',nan, ...
                   'yin_dp',nan,'nsdf_pk',nan,'comb_contrast',nan), 1, Nf);

% ----- Boucle frames -----
for j=1:Nf
    r = acf_frames(:,j);
    r = r / max(r(1), eps);                   % normalisation R(0)=1
    w = max(1, round(P.SmoothMs*1e-3*Fs));    % lissage léger
    if w>1, r2 = movavg(r, w); else, r2 = r; end

    % --- (1) Candidats multiples ---
    % 1A-bis : YIN/CMNDF-from-ACF
    [fyin, yin_dp] = cand_yin_acf(r2, Fs, tauMin, tauMax, P.YIN_thresh);
    % 1C-bis : NSDF approx (ici nsdf≈r2)
    [fnsdf, nsdf_pk] = cand_nsdf_acf(r2, Fs, tauMin, tauMax, P.NSDF_min);
    % 1E : peigne-produit log
    fcomb = cand_comb_log_acf(r2, Fs, tauMin, tauMax, Kc, wkComb, P.Comb_eps);
    % HPS (PSD via ACF -> Harmonic Sum)
    if P.HPS_use
        fhps = cand_hps_from_psd(r2, Fs, P.fmin, P.fmax, P.HPS_dfHz, Kh);
    else
        fhps = NaN;
    end

    fcands = [fyin, fnsdf, fcomb, fhps];
    fcands = fcands(isfinite(fcands) & fcands>0);
    if isempty(fcands)
        f_meas(j) = NaN; FR(j).q = 0; continue;
    end

    % --- (2) Fusion robuste (consensus) ---
    % Utilise la médiane sur l’échelle log-f -> robuste aux outliers
    lf = log(fcands);
    f_cons = exp(median(lf));
    % dispersion (qualité “cohérence” des méthodes)
    disp_log = mad(lf,1);                   % MAD sur log-f
    q_cons   = max(0, 1 - disp_log/0.25);   % ~1 si méthodes d'accord (≈±25%)

    % --- (3) Anti-demi-tour (vote multi-tests) ---
    Tcand = 1/f_cons;
    f_final = anti_half_decide(r2, Fs, Tcand, P);
    
    % --- (4) Qualité frame ---
    % carte simple : mélange d’indicateurs internes
    comb_contrast = odd_even_contrast_acf(r2, Fs, 1/f_final, Kc, wkComb, P.lambdaEven, P.UseLogComb);
    q_comb = (comb_contrast+1)/2;                   % map -> [0..1]
    q_yin  = max(0, 1 - yin_dp/0.35);               % d' petit => bon (≈ heur.)
    q_nsdf = nsdf_pk;                                % déjà ~[0..1]
    % qualité finale bornée [0..1]
    q = min(1, max(0, 0.25*q_cons + 0.35*q_yin + 0.25*q_nsdf + 0.15*q_comb));

    % --- Sauvegarde ---
    f_meas(j) = f_final;
    FR(j).fyin = fyin; FR(j).fnsdf=fnsdf; FR(j).fcomb=fcomb; FR(j).fhps=fhps;
    FR(j).f_cons = f_cons; FR(j).f_final = f_final; FR(j).q = q;
    FR(j).yin_dp = yin_dp; FR(j).nsdf_pk = nsdf_pk; FR(j).comb_contrast = comb_contrast;
end

% ----- Kalman [f; fdot] -----
if isempty(P.kal_dt)
    if numel(tFrames)>=2
        dt = median(diff(tFrames)); if dt<=0, dt = 1; end
    else
        dt = 1;
    end
else
    dt = P.kal_dt;
end
f_kal = kalman_f(f_meas, dt, P.kal_sigma_fdot, P.kal_sigma_meas, P.kal_adapt_from_quality, [FR.q]);

% ----- Plots optionnels -----
if P.Plot && Nf>1
    figure('Name','Candidats & fusion (ex. derniers frames)');
    jj = max(1,Nf-5):Nf;
    for k=1:numel(jj)
        j = jj(k);
        subplot(numel(jj),1,k);
        bar([FR(j).fyin, FR(j).fnsdf, FR(j).fcomb, FR(j).fhps, FR(j).f_cons, FR(j).f_final]);
        set(gca,'XTickLabel',{'YIN','NSDF','Comb','HPS','Cons','Final'}); grid on;
        ylabel('Hz'); title(sprintf('Frame %d  q=%.2f', j, FR(j).q));
    end
end
end

% ================== CANDIDATS ==================
function [fhat, dp_min] = cand_yin_acf(r, Fs, tauMin, tauMax, thr)
% CMNDF from ACF : d = 2*(1-r), d' = d * tau / cumsum(d)
d = 2*(1 - r(2:tauMax+1));
cum = cumsum(d); idx = (1:tauMax).';
dprime = d .* (idx ./ max(cum,eps));
dprime(1:tauMin-1) = +Inf;
[dp_min, L] = min(dprime);
% interpolation parabolique
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
% NSDF approx : nsdf ~ r (ACF normalisée), on cherche 1ère crête > minPk
nsdf = r; 
search = nsdf(1:tauMax+1); search(1:tauMin-1) = -Inf;
[pkv, pkl] = findpeaks_safe(search, 'MinPeakHeight', minPk);
if isempty(pkl)
    [pkv, pkl] = max(search);
else
    pkv = pkv(1); pkl = pkl(1);
end
% raffinement
guard = max(2, round(0.1*pkl));
i1 = max(2, pkl-guard); i2 = min(tauMax-1, pkl+guard);
[~,rel] = max(nsdf(i1:i2)); L = i1+rel-1;
if L>=2 && L<=tauMax-1
    y1=nsdf(L-1); y2=nsdf(L); y3=nsdf(L+1);
    delta = 0.5*(y1-y3)/max(y1-2*y2+y3, eps); delta = max(min(delta,0.5),-0.5);
else
    delta = 0;
end
Lref = L + delta; fhat = Fs / Lref; pk = pkv;
end

function fhat = cand_comb_log_acf(r, Fs, tauMin, tauMax, K, wk, epslog)
rpos = max(r, 0);
S = -Inf(tauMax,1);
for L=tauMin:tauMax
    s=0; cnt=0;
    for k=1:K
        q = k*L; if q>tauMax, break; end
        s = s + wk(k)*log(epslog + rpos(q+1));
        cnt = cnt + 1;
    end
    if cnt>=2, S(L)=s; end
end
[~,L0] = max(S);
if L0>tauMin && L0<tauMax
    y1=S(L0-1); y2=S(L0); y3=S(L0+1);
    if isfinite(y1)&&isfinite(y2)&&isfinite(y3) && (y1-2*y2+y3)<0
        delta = 0.5*(y1-y3)/max(y1-2*y2+y3, eps);
        delta = max(min(delta,0.5),-0.5);
    else
        delta = 0;
    end
else
    delta = 0;
end
Lref = L0 + delta; fhat = Fs / Lref;
end

function fhat = cand_hps_from_psd(r, Fs, fmin, fmax, dfHz, K)
% PSD via ACF, puis Harmonic Product/Sum
[F,Pxx] = psd_from_acf(r, Fs);
fgrid = (fmin:dfHz:fmax).';
Spec = @(f) interp1(F, Pxx, f, 'linear', 0);
S = zeros(size(fgrid));
for i=1:numel(fgrid)
    f0 = fgrid(i); s=0; cnt=0;
    for k=1:K
        fk = k*f0; if fk > F(end), break; end
        s = s + Spec(fk);
        cnt = cnt + 1;
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
% Reprend la stratégie 2) (odd-even + ratio + spectral)
K = P.Comb_K; wk = weights_K(K, P.Comb_weight);
tauMax = numel(r)-1;  interpR = @(q) safe_interp((0:tauMax).', r, q);

% Scores peigne ACF
[sComb_T, Sodd_T, Seven_T] = comb_score(r, T_cand*Fs, K, wk, P.lambdaEven, P.UseLogComb);
[sComb_H, Sodd_H, Seven_H] = comb_score(r, (T_cand*Fs)/2, K, wk, P.lambdaEven, P.UseLogComb);

% Spectre & scores odd/even
[Fspec, Pxx] = psd_from_acf(r, Fs);
BW = max(P.SpecBW_Hz, Fspec(2)-Fspec(1)); 
Spec = @(f0) band_energy(Fspec, Pxx, f0, BW);
f0 = 1/T_cand; fH = 2/T_cand;
[sSpec_T, sSpec_H] = deal(odd_even_spec(Spec, f0, K, wk, P.lambdaEven), ...
                          odd_even_spec(Spec, fH, K, wk, P.lambdaEven));

% Ratio direct
RT = interpR(T_cand*Fs); RH = interpR((T_cand*Fs)/2);
sRatio_T = (RT - P.gammaHalf*RH); sRatio_H = (RH - P.gammaHalf*RT);
nr = max(abs([sRatio_T sRatio_H])) + eps; sRatio_T = sRatio_T/nr; sRatio_H = sRatio_H/nr;

S_T = P.wComb*sComb_T + P.wSpec*sSpec_T + P.wRratio*sRatio_T;
S_H = P.wComb*sComb_H + P.wSpec*sSpec_H + P.wRratio*sRatio_H;

if S_H > S_T
    f_corr = fH;   % utiliser T/2 -> double fréquence
else
    f_corr = f0;   % garder T
end
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

function c = odd_even_contrast_acf(r, Fs, T, K, wk, lambdaEven, useLog)
[score,~,~] = comb_score(r, T*Fs, K, wk, lambdaEven, useLog);
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

% ================== Kalman ==================
function fhat = kalman_f(f_meas, dt, sigma_fdot, sigma_meas, adapt, qvec)
A = [1 dt; 0 1];
Q = (sigma_fdot^2) * [dt^4/4 dt^3/2; dt^3/2 dt^2];
H = [1 0];
x = [max(f_meas(1),0); 0]; P = eye(2);
fhat = nan(size(f_meas));
for k=1:numel(f_meas)
    % prédiction
    x = A*x; P = A*P*A.' + Q;
    z = f_meas(k);
    if isfinite(z) && z>0
        if adapt
            q = max(0.05, min(1, qvec(k)));      % qualité [0..1]
            R = (sigma_meas^2) / (q^2);          % +qualité -> -bruit mesure
        else
            R = sigma_meas^2;
        end
        S = H*P*H.' + R; Kk = (P*H.')/S;
        x = x + Kk*(z - H*x); P = (eye(2)-Kk*H)*P;
    end
    fhat(k) = x(1);
end
end

% ================== Utilitaires ==================
function w = weights_K(K, mode)
switch lower(string(mode))
    case "equal", w = ones(1,K);
    otherwise,    w = 1./(1:K);
end
w = w / sum(w);
end

function [F, Pxx] = psd_from_acf(r, Fs)
rsym = [flipud(r(2:end)); r(:)];
N = numel(rsym);
w = 0.5 - 0.5*cos(2*pi*(0:N-1)'/(N-1));
S = fft(rsym .* w);
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
end

function [pks, locs] = findpeaks_safe(x, varargin)
try
    [pks, locs] = findpeaks(x, varargin{:});
catch
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
