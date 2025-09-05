function [fr_hat, OUT] = estimate_fr_correntropy_1H(Fs, P, varargin)
% Estimation de la fréquence via "correntropy" (ACF non-linéaire).
%
% Deux modes d'entrée :
%   estimate_fr_correntropy_1H(Fs, P, 'x',  x)   % signal dispo -> vraie correntropy
%   estimate_fr_correntropy_1H(Fs, P, 'acf',acf) % ACF seule     -> proxy correntropy
%
% Sorties :
%   fr_hat : fréquence estimée (Hz)
%   OUT    : struct (V correntropy, diagnostics, etc.)

% --------- Parse inputs ----------
x = []; acf = [];
for k=1:2:numel(varargin)
    switch lower(varargin{k})
        case 'x',   x   = varargin{k+1}(:);
        case 'acf', acf = varargin{k+1}(:);
    end
end
haveX = ~isempty(x);
haveA = ~isempty(acf);

% --------- Defaults ----------
DEF = struct('MaxLagSec',1.0,'fmin',0.5,'fmax',120, ...
             'sigmaMode','mad','sigma',[],'sigmaProxy',0.25, ...
             'Threshold',0.15,'SmoothMs',1.0, ...
             'Kharm',4,'lambdaEven',0.8,'gammaHalf',1.12,'Plot',true);
fn = fieldnames(DEF); for i=1:numel(fn), if ~isfield(P,fn{i}), P.(fn{i})=DEF.(fn{i}); end, end

% --------- Bornes lag ----------
if haveX
    N = numel(x);
    maxLagAvail = N-2;
else
    maxLagAvail = numel(acf)-2;
end
tauMin = max(2, floor(Fs / P.fmax));
tauMax = min(maxLagAvail, ceil(Fs / P.fmin));
maxLagSamp = min(maxLagAvail, round(P.MaxLagSec*Fs));
tauMax = min(tauMax, maxLagSamp);
if tauMin+2 >= tauMax
    error('Plage de lag invalide: ajuste fmin/fmax/MaxLagSec au regard de la longueur.');
end

% --------- 1) Correntropy V(τ) ----------
% V(τ) = E[ exp( - (x[n]-x[n-τ])^2 / (2σ^2) ) ]
if haveX
    % Sigma : largeur du noyau (échelle de l'amplitude de x)
%     sig = select_sigma(x, P);
    %case 'silverman'
    %sig = 1.06*std(double(x))*numel(x)^(-1/5);
    %if sig<=0, sig = max(1e-3, std(double(x))+eps); end
    xm = median(x);
    sig = 1.4826*median(abs(x - xm)) + eps;
    sig = max(sig, 1e-6);
    V = compute_correntropy_signal(x, tauMax, sig);
else
    % ACF-only : proxy (stationnaire approx) => remplace le produit par kernel sur r(τ)
    r = acf / max(acf(1), eps);                 % normalise R(0)=1
    sigp = max(P.sigmaProxy, 1e-3);
    V = exp( -(1 - max(r(1:tauMax+1), -1)) / (sigp^2) );   % K(r) ~ exp(-(1-r)/σp^2)
end

% Normalise à V(0)=1 (sécurité)
V = V(:); V = V / max(V(1), eps);

% --------- 2) Sélection de période (CMNDF-like sur correntropy) ----------
% d(τ) = 1 - V(τ)  (analogue YIN via correntropy)
d = 1 - V(2:tauMax+1);
cum = cumsum(d); tIdx = (1:tauMax).';
dprime = d .* (tIdx ./ max(cum, eps));         % CMNDF
% lissage léger
w = max(1, round(P.SmoothMs*1e-3*Fs));
dprime_s = movavg(dprime, w);

% seuil YIN : première traversée sous Threshold
search = dprime_s;
search(1:tauMin-1) = +Inf;
idx = find(search < P.Threshold, 1, 'first');
if ~isempty(idx)
    Lc = idx; win = max(3, round(0.1*Lc));
    i1 = max(tauMin, Lc-win); i2 = min(tauMax, Lc+win);
    [~,rel] = min(dprime_s(i1:i2)); L = i1 + rel - 1;
else
    [~,L] = min(dprime_s);
    if L<tauMin || L>tauMax, [~,L] = min(dprime_s(tauMin:tauMax)); L = L + tauMin - 1; end
end
% interpolation parabolique
if L<2 || L>tauMax-1
    Lref = L;
else
    y1=dprime_s(L-1); y2=dprime_s(L); y3=dprime_s(L+1);
    delta = 0.5*(y1 - y3) / max(y1 - 2*y2 + y3, eps);
    delta = max(min(delta, 0.5), -0.5);
    Lref = L + delta;
end
T_hat = Lref / Fs;

% --------- 3) Anti demi-tour (sur V(τ)) ----------
interpV = @(q) safe_interp((0:tauMax).', V, q);
T  = T_hat * Fs;
V_T  = interpV(T);
V_T2 = interpV(T/2);

% Peigne impairs/pairs
S_odd=0; S_even=0;
for k=1:P.Kharm
    q = k*T; if q<=tauMax
        if mod(k,2)==1, S_odd  = S_odd  + interpV(q);
        else            S_even = S_even + interpV(q);
        end
    end
end
combContrast = (S_odd - P.lambdaEven*S_even)/max(S_odd + P.lambdaEven*S_even, eps);

isHalf = (V_T2 > P.gammaHalf*V_T) && (S_even > S_odd);
if isHalf
    T_hat = T_hat/2;
end

fr_hat = 1 / max(T_hat, eps);

% --------- 4) Sorties & graphes ----------
OUT = struct();
OUT.V = V;
OUT.T_hat = T_hat;
OUT.nsdf_like_T = 1 - interp_dprime(dprime_s, T_hat*Fs); % indicateur (petit d' ⇒ grande périodicité)
OUT.periodicity = interpV(T_hat*Fs);
OUT.combContrast= combContrast;
OUT.isHalf      = isHalf;
OUT.params      = P;

if P.Plot
    tau = (0:tauMax).'/Fs;
    figure('Name','Correntropy V(\tau)');
    plot(tau, V, 'b'); grid on; hold on;
    xline(T_hat, '--', sprintf('T=%.4f s (f=%.2f Hz)', T_hat, fr_hat), 'Color',[0.2 0.2 0.2]);
    xline(T_hat/2, ':', 'T/2'); xline(2*T_hat, ':', '2T');
    xlabel('\tau (s)'); ylabel('V(\tau)'); xlim([0 P.MaxLagSec]);
    title('Correntropy (normalisée)');

    figure('Name','CMNDF sur correntropy');
    tt = (1:numel(dprime_s)).'/Fs;
    plot(tt, dprime_s, 'b'); grid on; hold on;
    yline(P.Threshold, '--', 'Seuil'); xline(T_hat, '--');
    xlabel('\tau (s)'); ylabel('d''(\tau)'); xlim([0 P.MaxLagSec]);
    title('Sélection de période (YIN/CMNDF sur V)');

    figure('Name','Odd vs Even (sur V)');
    bar([S_odd, S_even]); grid on; set(gca,'XTickLabel',{'impairs','pairs'});
    ylabel('Somme V(kT)'); title(sprintf('Contraste=%.3f  %s', combContrast, tern(isHalf,'→ demi-tour corrigé','')));
end
end

% =================== Helpers ===================

function sig = select_sigma(x, P)
% Choix de σ pour le noyau gaussien (amplitude du signal)
switch lower(P.sigmaMode)
    case 'fixed'
        if isempty(P.sigma), error('sigmaMode=fixed mais P.sigma est vide.'); end
        sig = P.sigma;
    case 'silverman'
        sig = 1.06*std(double(x))*numel(x)^(-1/5);
        if sig<=0, sig = max(1e-3, std(double(x))+eps); end
    otherwise % 'mad' (robuste)
        xm = median(x);
        sig = 1.4826*median(abs(x - xm)) + eps;
end
% garde une valeur minimale non nulle
sig = max(sig, 1e-6);
end

function V = compute_correntropy_signal(x, tauMax, sigma)
% V(τ) = mean( exp( - (x[n]-x[n-τ])^2 / (2σ^2) ) )
x = double(x(:));
N = numel(x);
V = zeros(tauMax+1,1);
V(1) = 1;  % à τ=0, la correntropy vaut 1 (kernel à distance nulle)
for t = 1:tauMax
    M = N - t;
    if M<=0, break; end
    d = x(1+t:N) - x(1:N-t);
    V(1+t) = mean( exp( - (d.^2) / (2*sigma^2) ) );
end
end

function y = movavg(x,w)
if w<=1, y=x; return; end
k = ones(w,1)/w;
try, y = filtfilt(k,1,double(x));
catch, y = conv(double(x), k, 'same'); end
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

function dp = interp_dprime(dprime, q)
N = numel(dprime);
if q<2 || q>N-1
    i0 = max(1, min(N-1, floor(q))); a = q - i0;
    dp = (1-a)*dprime(i0) + a*dprime(i0+1);
    return;
end
L = round(q);
y1=dprime(L-1); y2=dprime(L); y3=dprime(L+1);
delta = 0.5*(y1 - y3) / max(y1 - 2*y2 + y3, eps);
delta = max(min(delta,0.5), -0.5);
Lref = L + delta;
i0 = floor(Lref); a = Lref - i0;
dp = (1-a)*dprime(i0) + a*dprime(i0+1);
end

function s = tern(c,a,b)
if c, s=a; else, s=b; end
end
