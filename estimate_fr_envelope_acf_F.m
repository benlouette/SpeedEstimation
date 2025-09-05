function [fr_hat, OUT] = estimate_fr_envelope_acf_F(x, Fs, P)
% Solution F — ACF d'enveloppe (après non-linéarité)
% Étapes : bande résonante -> non-linéarité -> LPF -> décimation -> ACF ->
%          YIN/CMNDF (depuis ACF) -> anti-demi-tour (odd/even + R(T/2) vs R(T))

x = x(:);

% --------- Défaults ---------
DEF = struct('bandMode','auto','searchBand',[200 3000],'nBands',18,'bandManual',[800 2000], ...
    'nonlinearity','halfwave','gamma',0.7,'envFs',500, ...
    'fmin',0.5,'fmax',120,'maxLagSec',1.0,'yinThresh',0.15,'smoothMs',1.0, ...
    'Kharm',4,'lambdaEven',0.8,'gammaHalf',1.12,'plot',true);
fn = fieldnames(DEF);
for k=1:numel(fn), if ~isfield(P,fn{k}), P.(fn{k}) = DEF.(fn{k}); end, end

% --------- 1) Choix de la bande résonante ---------
switch lower(P.bandMode)
    case 'manual'
        f1f2 = P.bandManual;
    otherwise
        f1f2 = find_impulsive_band_kurt(x, Fs, P.searchBand, P.nBands);
end
y = bp_iir(x, Fs, f1f2(1), f1f2(2), 4);

% --------- 2) Non-linéarité -> enveloppe brute ---------
switch lower(P.nonlinearity)
    case 'halfwave'
        env0 = max(y, 0);               % rectification demi-onde
    case 'abs'
        env0 = abs(y);
    case 'power'
        g = max(P.gamma, 0.1);
        env0 = sign(y).*abs(y).^g;      % "soft-rect"
        env0 = max(env0, 0);            % garder partie positive
    case 'hilbert'
        env0 = abs(hilbert(y));
    otherwise
        env0 = max(y, 0);
end

% --------- 3) Lissage (LPF) & décimation ---------
M = max(1, floor(Fs / P.envFs));
Fs_env = Fs / M;
fc = min(0.4*Fs_env, 2*P.fmax);  % coupe à < Nyquist env et > fmax
env1 = lp_iir(env0, Fs, fc, 4);
env  = env1(1:M:end);

% --------- 4) ACF de l'enveloppe ---------
maxLag = min(round(P.maxLagSec*Fs_env), numel(env)-2);
acf_env = acf_fft_norm(env, maxLag);  % R_env(0)=1
lags = (0:maxLag).'; tau = lags/Fs_env;

% --------- 5) YIN/CMNDF (depuis ACF) ---------
tauMin = max(2, floor(Fs_env / P.fmax));
tauMax = min(maxLag-1, ceil(Fs_env / P.fmin));
if tauMin+2 >= tauMax
    fr_hat = 0; OUT = struct('reason','bad_range'); return;
end
% d(tau) = 2*(1 - r(tau)), CMNDF : d'(tau) = d(tau) * tau / sum_{j<=tau} d(j)
d = 2*(1 - acf_env(2:tauMax+1));
cum = cumsum(d); tIdx = (1:tauMax).';
dprime = d .* (tIdx ./ max(cum,eps));
% Lissage léger
w = max(1, round(P.smoothMs*1e-3*Fs_env));
dprime_s = movavg(dprime, w);

% Sélection YIN : première tau où d' < seuil
mask = true(size(dprime_s)); mask(1:tauMin-1) = false;
idx = find(dprime_s(mask) < P.yinThresh, 1, 'first');
if ~isempty(idx)
    Lc = find(mask,1,'first') + idx - 1;
    win = max(3, round(0.1*Lc));
    i1 = max(tauMin, Lc - win); i2 = min(tauMax, Lc + win);
    [~,rel] = min(dprime_s(i1:i2)); L = i1 + rel - 1;
else
    [~,L] = min(replace_inf(~mask, dprime_s, +Inf));
    if L<tauMin || L>tauMax, [~,L] = min(dprime_s(tauMin:tauMax)); L = L + tauMin - 1; end
end
% Interpolation parabolique
if L<2 || L>tauMax-1
    Lref = L;
else
    y1=dprime_s(L-1); y2=dprime_s(L); y3=dprime_s(L+1);
    delta = 0.5*(y1-y3)/max(y1-2*y2+y3, eps);
    delta = max(min(delta,0.5),-0.5);
    Lref = L + delta;
end
T_hat = Lref / Fs_env;

% --------- 6) Anti-demi-tour (sécurité) ---------
interpR = @(q) safe_interp((0:maxLag).', acf_env, q);
T  = T_hat * Fs_env;  % en samples "enveloppe"
R_T  = interpR(T);
R_T2 = interpR(T/2);
% peigne odd/even
S_odd=0; S_even=0;
for k=1:P.Kharm
    q = k*T; if q<=tauMax
        if mod(k,2)==1, S_odd=S_odd+interpR(q); else, S_even=S_even+interpR(q); end
    end
end
combContrast = (S_odd - P.lambdaEven*S_even)/max(S_odd + P.lambdaEven*S_even, eps);
isHalf = (R_T2 > P.gammaHalf*R_T) && (S_even > S_odd);
if isHalf, T_hat = T_hat/2; end

fr_hat = 1 / max(T_hat, eps);

% --------- 7) Sorties & graphes ---------
OUT = struct();
OUT.Fs_env = Fs_env;
OUT.band   = f1f2;
OUT.env    = env;
OUT.acf_env= acf_env;
OUT.nsdf_T = 1 - (interp_dprime(dprime_s, T_hat*Fs_env)/2); % approx
OUT.periodicity = interpR(T_hat*Fs_env);
OUT.combContrast= combContrast;
OUT.isHalf = isHalf;
OUT.T_hat = T_hat;

if P.plot
    t = (0:numel(x)-1)'/Fs;
    figure('Name','Signal & bande résonante');
    subplot(2,1,1); plot(t, x); grid on; xlabel('t (s)'); ylabel('x(t)'); title('Signal brut');
    subplot(2,1,2);
    y_show = bp_iir(x, Fs, f1f2(1), f1f2(2), 4);
    plot(t, y_show); grid on; xlabel('t (s)'); ylabel('y(t)'); 
    title(sprintf('Bande résonante [%.0f %.0f] Hz', f1f2(1), f1f2(2)));

    te = (0:numel(env)-1)'/Fs_env;
    figure('Name','Enveloppe');
    plot(te, env); grid on; xlabel('t (s)'); ylabel('enveloppe'); 
    title(sprintf('Non-linéarité: %s, Fs_{env}=%.0f Hz', P.nonlinearity, Fs_env));

    figure('Name','ACF enveloppe & YIN');
    subplot(2,1,1);
    plot(tau, acf_env); grid on; xlim([0 P.maxLagSec]);
    xline(T_hat,'--',sprintf('T=%.4f s  (f=%.2f Hz)', T_hat, fr_hat));
    xlabel('\tau (s)'); ylabel('R_{env}(\tau)'); title('ACF de l''enveloppe (norm.)');
    subplot(2,1,2);
    tt = (1:numel(dprime_s)).'/Fs_env;
    plot(tt, dprime_s); grid on; yline(P.yinThresh,'--','seuil');
    xline(T_hat,'--'); xlabel('\tau (s)'); ylabel('d''(\tau)'); 
    title('YIN/CMNDF sur ACF enveloppe');
end
end

% ================= Helpers =================
function f1f2 = find_impulsive_band_kurt(x, Fs, searchBand, nBands)
fmin = max(10, searchBand(1));
fmax = min(Fs/2-100, searchBand(2));
edges = logspace(log10(fmin), log10(fmax), nBands+1);
scores = zeros(nBands,1);
for k=1:nBands
    y = bp_iir(x, Fs, edges(k), edges(k+1), 4);
    scores(k) = kurtosis(y);
end
[~,i] = max(scores);
f1f2 = [edges(i) edges(i+1)];
end

function y = bp_iir(x, Fs, f1, f2, order)
W = sort([f1 f2])/(Fs/2);
W(1)=max(W(1),1e-6); W(2)=min(W(2),0.999);
[b,a]=butter(order, W, 'bandpass');
y=filtfilt(b,a,double(x));
end

function y = lp_iir(x, Fs, fc, order)
W = min(fc/(Fs/2), 0.999);
[b,a]=butter(order, W, 'low');
y=filtfilt(b,a,double(x));
end

function acf = acf_fft_norm(x, maxLag)
x = x(:) - mean(x);
L = numel(x);
nfft = 2^nextpow2(2*L-1);
X = fft(x, nfft);
S = X.*conj(X);
r = ifft(S, 'symmetric');
r = r(1:maxLag+1);
acf = r / max(r(1), eps);
end

function y = movavg(x,w)
if w<=1, y=x; return; end
k=ones(w,1)/w;
try, y=filtfilt(k,1,double(x)); catch, y=conv(double(x),k,'same'); end
end

function dp = interp_dprime(dprime, q)
N = numel(dprime);
if q<2 || q>N-1
    i0 = max(1,min(N-1,floor(q))); a = q - i0;
    dp = (1-a)*dprime(i0) + a*dprime(i0+1);
    return;
end
L = round(q);
y1=dprime(L-1); y2=dprime(L); y3=dprime(L+1);
delta = 0.5*(y1 - y3) / max(y1 - 2*y2 + y3, eps);
delta = max(min(delta,0.5),-0.5);
Lref = L + delta;
i0=floor(Lref); a=Lref-i0; dp=(1-a)*dprime(i0)+a*dprime(i0+1);
end

function v = safe_interp(x, y, q)
q = q(:); v = nan(size(q)); N = numel(x);
for i=1:numel(q)
    if q(i) < 1 || q(i) > N-2, v(i)=NaN;
    else
        i0=floor(q(i)); a=q(i)-i0; v(i)=(1-a)*y(i0+1)+a*y(i0+2);
    end
end
if isscalar(v), v=v(1); end
end

function y = replace_inf(mask, x, val)
y = x; y(mask) = val;
end

function s = tern(c,a,b)
if c, s=a; else, s=b; end
end
