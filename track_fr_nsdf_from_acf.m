function [tCenters, f_meas, f_filt, FRAMES] = track_fr_nsdf_from_acf(x, Fs, WinSec, HopSec, P, KAL)
% TRACK_FR_NSDF_FROM_ACF
%   Tracking f_r(t) par fenêtres glissantes :
%     - ACF par fenêtre (R(0) normalisée à 1)
%     - NSDF (approx. stationnaire) + McLeod peak picking (parabolique)
%     - Anti-demi-tour (R(T/2) vs R(T) + peigne impairs/pairs)
%     - Kalman 2x2 (état [f; fdot]) pour lisser
%
% In :
%   x       : signal temporel (brut ou déjà filtré/enveloppé)
%   Fs      : Hz
%   WinSec  : longueur fenêtre (s)
%   HopSec  : pas (s)
%   P       : params NSDF/anti-demi-tour (voir demo)
%   KAL     : params Kalman (dt, init, sigmas,...)
%
% Out:
%   tCenters: temps centre de fenêtre (s)
%   f_meas  : mesure par fenêtre (Hz)
%   f_filt  : estimation lissée (Kalman) (Hz)
%   FRAMES  : struct array diag (nsdf_peak, periodicity, isHalf, ...)

x = x(:);
N = numel(x);
win  = round(WinSec*Fs);
hop  = round(HopSec*Fs);
nFrm = 1 + floor((N - win)/hop);
tCenters = ((0:nFrm-1)*hop + (win/2))/Fs;

% Pré-allocation
f_meas = nan(nFrm,1);
FRAMES = repmat(struct('nsdf_peak',nan,'periodicity',nan,'isHalf',false,'T',nan,'f',nan), nFrm,1);

% Boucle fenêtres
for m = 1:nFrm
    i1 = (m-1)*hop + 1;
    i2 = i1 + win - 1;
    xw = x(i1:i2);

    % --- ACF par FFT (rapide), 0..maxLag ---
    maxLagSamp = min(round(P.MaxLagSec*Fs), win-2);
    acf = acf_frame_fft(xw, maxLagSamp);  % normalisée R(0)=1

    % --- Estimation NSDF/MPM (ACF-only) + anti-demi-tour ---
    [fr_hat, OUT] = acf_nsdf_estimator(acf, Fs, P);

    f_meas(m) = fr_hat;
    FRAMES(m).nsdf_peak   = OUT.nsdf_T;
    FRAMES(m).periodicity = OUT.periodicity;
    FRAMES(m).isHalf      = OUT.isHalf;
    FRAMES(m).T           = OUT.T_hat;
    FRAMES(m).f           = fr_hat;
end

% --- Kalman 2x2 sur f ---
f_filt = kalman_track(f_meas, KAL);

end

% ====== ACF (FFT-based, normalisée) ======
function acf = acf_frame_fft(xw, maxLag)
% ACF rapide (non biaisée approx.) via FFT sur fenêtre locale
xw = xw(:) - mean(xw);
L  = numel(xw);
nfft = 2^nextpow2(2*L-1);
X = fft(xw, nfft);
S = X .* conj(X);
r = ifft(S, 'symmetric');
r = r(1:L);                     % lags >= 0
% normalisation
acf = r / max(r(1), eps);
acf = acf(1:maxLag+1);
end

% ====== NSDF/MPM + anti-demi-tour depuis ACF ======
function [fr_hat, OUT] = acf_nsdf_estimator(acf, Fs, P)
acf = acf(:);
Nlag = numel(acf);
lags = (0:Nlag-1).';

% bornes
tauMin = max(2, floor(Fs/P.fmax));
tauMax = min(Nlag-2, ceil(Fs/P.fmin));
if tauMin+2 >= tauMax
    fr_hat = 0; OUT = struct('nsdf_T',0,'periodicity',0,'isHalf',false,'T_hat',NaN); return;
end

% NSDF approx. stationnaire : nsdf ~ r = acf (déjà R(0)=1)
nsdf = acf;
% lissage léger
w = max(1, round(P.SmoothMs*1e-3*Fs));
nsdf_s = movavg(nsdf, w);

% Cherche première crête > MinNSDF dans [tauMin..tauMax]
search = nsdf_s; search(1:tauMin-1) = -Inf; search(tauMax+1:end) = -Inf;
[pkVal, pkLoc] = findpeaks_safe(search, 'MinPeakHeight', P.MinNSDF);
if isempty(pkLoc)
    [pkVal, pkLoc] = max(search);   % fallback : max dans la bande
else
    pkVal = pkVal(1); pkLoc = pkLoc(1);
end

% Affinage parabolique
guard = max(2, round(P.PeakGuard*pkLoc));
i1 = max(2, pkLoc-guard); i2 = min(tauMax-1, pkLoc+guard);
[~, rel] = max(nsdf_s(i1:i2)); L = i1+rel-1;
if L<2 || L>tauMax-1, Lref = L;
else
    y1=nsdf_s(L-1); y2=nsdf_s(L); y3=nsdf_s(L+1);
    delta = 0.5*(y1-y3)/max(y1-2*y2+y3, eps);
    delta = max(min(delta,0.5),-0.5);
    Lref  = L + delta;
end
T_hat = Lref/Fs;

% Anti-demi-tour sur l'ACF
interpR = @(q) safe_interp((0:Nlag-1).', acf, q);
T  = T_hat*Fs;
R_T  = interpR(T);
R_T2 = interpR(T/2);

% peigne odd/even
S_odd=0; S_even=0;
for k=1:P.Kharm
    q = k*T; if q<=tauMax
        if mod(k,2)==1, S_odd=S_odd+interpR(q);
        else,           S_even=S_even+interpR(q);
        end
    end
end
combContrast = (S_odd - P.lambdaEven*S_even)/max(S_odd + P.lambdaEven*S_even, eps);

isHalf = (R_T2 > P.gammaHalf*R_T) && (S_even > S_odd);
if isHalf, T_hat = T_hat/2; end

fr_hat = 1/max(T_hat,eps);

OUT = struct();
OUT.nsdf_T = nsdf_s(max(1,round(T_hat*Fs)));
OUT.periodicity = interpR(T_hat*Fs);
OUT.isHalf = isHalf;
OUT.T_hat = T_hat;
OUT.combContrast = combContrast;
end

% ====== Kalman 2x2 (f, fdot) discret ======
function fhat = kalman_track(f_meas, K)
dt = K.dt;
A = [1 dt; 0 1];         % modèle vitesse constante
Q = (K.sigma_fdot^2) * [dt^4/4 dt^3/2; dt^3/2 dt^2];  % bruit de processus
% init
xhat = [K.f0; K.fdot0];
P    = eye(2);
fhat = nan(size(f_meas));
for k=1:numel(f_meas)
    % prediction
    xhat = A*xhat;
    P    = A*P*A.' + Q;
    z = f_meas(k);
    if isfinite(z) && z>0
        if K.adapt_from_nsdf && isfield(K,'nsdf_map') % (optionnel)
            R = K.sigma_meas^2;
        else
            R = K.sigma_meas^2;
        end
        H = [1 0]; S = H*P*H.' + R; Kk = (P*H.')/S;
        xhat = xhat + Kk*(z - H*xhat);
        P = (eye(2)-Kk*H)*P;
    end
    fhat(k) = xhat(1);
end
end

% ===== helpers =====
function y = movavg(x,w)
if w<=1, y=x; return; end
k = ones(w,1)/w;
try
    y = filtfilt(k,1,double(x));
catch
    y = conv(double(x), k, 'same');  % fallback si pas de filtfilt
end
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
