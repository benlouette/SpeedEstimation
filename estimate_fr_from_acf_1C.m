function [fr_hat, OUT] = estimate_fr_from_acf_1C(acf, Fs, P, x)
% Estimation de la fréquence de rotation via NSDF / McLeod (solution 1C).
% Fonctionne:
%   (A) En mode **exact** si le signal 'x' est fourni (NSDF McLeod classique)
%   (B) En mode **ACF-only** sinon, via une approximation stationnaire + correction
%       de bord optionnelle selon la normalisation ('biased'/'unbiased').
%
% In:
%   acf : vecteur ACF (lags >= 0, acf(1)=R(0) au sens de ta normalisation)
%   Fs  : fréquence d’échantillonnage (Hz)
%   P   : paramètres (voir demo)
%   x   : (optionnel) signal temporel original — permet NSDF exacte
%
% Out:
%   fr_hat : fréquence estimée (Hz)
%   OUT    : diagnostics (struct)

acf = acf(:);
Nlag = numel(acf);
lags = (0:Nlag-1).';
tau  = lags / Fs;

% ---------- paramètres par défaut ----------
def = struct('MaxLagSec',1.0,'fmin',0.5,'fmax',120,'MinNSDF',0.5, ...
             'SmoothMs',1.0,'PeakGuard',0.1,'Kharm',4,'lambdaEven',0.8, ...
             'gammaHalf',1.10,'ACFMode','auto','Nproxy',[],'Plot',true);
fn = fieldnames(def);
for k=1:numel(fn)
    if ~isfield(P, fn{k}), P.(fn{k}) = def.(fn{k}); end
end

% ---------- bornes de recherche ----------
tauMinSamp = max(2, floor(Fs / P.fmax)); % >=2 pour parabolique
tauMaxSamp = min(Nlag-2, ceil(Fs / P.fmin));
maxLagSamp = min(Nlag-1, round(P.MaxLagSec*Fs));
tauMaxSamp = min(tauMaxSamp, maxLagSamp);
if tauMinSamp+2 >= tauMaxSamp
    error('Plage de lag invalide : ajuste fmin/fmax/MaxLagSec.');
end

% ---------- NSDF (McLeod) ----------
haveX = (nargin>=4) && ~isempty(x);
if haveX
    % ===== Mode exact (McLeod) depuis x =====
    x = x(:);
    N  = numel(x);
    tauIdx = (1:tauMaxSamp).';
    nsdf = zeros(tauMaxSamp,1);
    % précompute énergies cumulées
    x2  = x.^2;
    E   = cumsum(x2);           % E[k] = sum_{n=1..k} x[n]^2
    E0  = E(end);               % énergie totale
    for i=1:numel(tauIdx)
        T = tauIdx(i);
        M = N - T;
        if M <= 0, break; end
        % num = 2 * sum_{n=T+1..N} x[n] x[n-T]
        num = 2 * sum( x(1+T:N) .* x(1:N-T) );
        % denom = sum_{n=1..N} x[n]^2 + sum_{n=T+1..N} x[n-T]^2 = E0 + E[N-T]
        denom = E0 + E(N-T);
        nsdf(i) = num / max(denom, eps);
    end
    % pad à la taille affichage
    nsdf = [0; nsdf];  %#ok<AGROW>
else
    % ===== Mode ACF-only (approximation stationnaire) =====
    % r = ACF normalisée à R(0)
    R0 = acf(1);
    r  = (R0~=0) * (acf / R0) + (R0==0) * acf;

    % Correction de bord optionnelle selon normalisation de l'ACF
    % - 'normalized' : déjà corrél. normalisée -> nsdf ≈ r
    % - 'biased'     : acf ≈ C(τ)/N      -> approx C(τ) ∝ acf, Eτ/E0 ≈ (N-τ)/N
    % - 'unbiased'   : acf ≈ C(τ)/(N-τ)  -> approx C(τ) ∝ acf*(Nproxy-τ)
    % Si P.Nproxy vide, on prend Nproxy = Nlag*2 (approx).
    mode = lower(string(P.ACFMode));
    if mode == "auto"
        mode = "normalized";  % choix sûr par défaut
    end
    nsdf = zeros(Nlag,1);
    switch mode
        case "normalized"
            % NSDF ≈ r (stationnaire)
            nsdf = r;
        case "biased"
            Np = P.Nproxy; if isempty(Np), Np = Nlag*2; end
            tauIdx = (0:Nlag-1).';
            % C(τ) ≈ acf * Np ; Eτ/E0 ≈ (Np-τ)/Np
            % nsdf ≈ 2*C(τ) / (E0 + Eτ) -> proportionnel à r(τ) * (Np/(2*Np - τ))
            w = Np ./ max(2*Np - tauIdx, 1); % fenêtre de correction
            nsdf = r .* w;
        case "unbiased"
            Np = P.Nproxy; if isempty(Np), Np = Nlag*2; end
            tauIdx = (0:Nlag-1).';
            % C(τ) ≈ acf*(Np - τ) ; Eτ/E0 ≈ (Np-τ)/Np
            % nsdf ≈ 2*C(τ) / (E0 + Eτ) -> ∝ r(τ) * (Np - τ) / (2*Np - τ)
            w = max(Np - tauIdx, 0) ./ max(2*Np - tauIdx, 1);
            nsdf = r .* w * (2);  % facteur 2 pour garder l'échelle ~[0..1]
        otherwise
            nsdf = r;  % fallback
    end
end

% ---------- Lissage léger ----------
w = max(1, round(P.SmoothMs*1e-3*Fs));
nsdf_s = movavg(nsdf, w);

% ---------- Recherche de la crête NSDF ----------
search = nsdf_s(1:tauMaxSamp+1);
search(1:tauMinSamp-1) = -Inf;              % exclure < tau_min
[pkVal, pkLoc] = findpeaks_safe(search, 'MinPeakHeight', P.MinNSDF);
if isempty(pkLoc)
    % pas de crête au-dessus du seuil -> prendre le max global dans la bande
    [pkVal, pkLoc] = max(search);
else
    % garder la première crête au-dessus du seuil (McLeod)
    pkVal = pkVal(1); pkLoc = pkLoc(1);
end

% ---------- Interpolation parabolique autour de la crête ----------
guard = max(2, round(P.PeakGuard * pkLoc));
i1 = max(2, pkLoc - guard);
i2 = min(tauMaxSamp-1, pkLoc + guard);
[~, rel] = max(nsdf_s(i1:i2));
L = i1 + rel - 1;
if L<2 || L>tauMaxSamp-1
    Lref = L;
else
    y1 = nsdf_s(L-1); y2 = nsdf_s(L); y3 = nsdf_s(L+1);
    delta = 0.5*(y1 - y3) / max(y1 - 2*y2 + y3, eps);
    delta = max(min(delta, 0.5), -0.5);
    Lref = L + delta;
end
T_hat = Lref / Fs;

% ---------- Anti-demi-tour (ACF) ----------
% Interp linéaire sur ACF normalisée pour les tests de parité
r_norm = acf / max(acf(1), eps);
interpR = @(q) safe_interp((0:Nlag-1).', r_norm, q);

T  = T_hat * Fs;   % en échantillons
R_T  = interpR(T); 
R_T2 = interpR(T/2);

% Peigne odd/even
K = P.Kharm;
S_odd = 0; S_even = 0;
for k=1:K
    q = k*T;
    if q <= tauMaxSamp
        if mod(k,2)==1, S_odd  = S_odd  + interpR(q);
        else            S_even = S_even + interpR(q);
        end
    end
end
combContrast = (S_odd - P.lambdaEven*S_even) / max(S_odd + P.lambdaEven*S_even, eps);

% Critère demi-tour
isHalf = (R_T2 > P.gammaHalf*R_T) && (S_even > S_odd);
if isHalf
    T_hat = T_hat/2;
end

fr_hat = 1 / max(T_hat, eps);

% ---------- sorties & tracés ----------
OUT = struct();
OUT.nsdf = nsdf_s; OUT.Lref = Lref; OUT.T_hat = T_hat;
OUT.nsdf_T = nsdf_s(max(1, round(T_hat*Fs)));
OUT.periodicity = R_T;
OUT.combContrast = combContrast; OUT.isHalf = isHalf;
OUT.params = P;

if P.Plot
    figure('Name','NSDF / McLeod');
    tt = (0:numel(nsdf_s)-1).'/Fs;
    plot(tt, nsdf_s, 'b'); grid on; hold on;
    yline(P.MinNSDF, '--', 'Seuil');
    xline(T_hat, '--', sprintf('T\\_hat=%.4f s  (f=%.2f Hz)', T_hat, fr_hat), 'Color',[0.2 0.2 0.2]);
    xlabel('\tau (s)'); ylabel('NSDF(\tau)'); xlim([0, P.MaxLagSec]);
    title('NSDF (lissée) & T estimé');

    figure('Name','ACF & périodes');
    plot(tau, r_norm, 'b'); grid on; hold on;
    xline(T_hat, '--', 'T'); xline(T_hat/2, ':', 'T/2');
    xlabel('\tau (s)'); ylabel('ACF (norm.)'); xlim([0, P.MaxLagSec]);
    title('ACF normalisée & marquage périodicités');

    figure('Name','Peigne odd/even sur ACF');
    tmarks = (1:K)*T_hat; yodd = nan(size(tmarks)); yeven = nan(size(tmarks));
    for k=1:K
        if mod(k,2)==1, yodd(k) = interpR(k*T_hat*Fs);
        else            yeven(k)= interpR(k*T_hat*Fs);
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

function [pks, locs] = findpeaks_safe(x, varargin)
% utilise findpeaks si dispo (Signal Toolbox), sinon fallback
try
    [pks, locs] = findpeaks(x, varargin{:});
catch
    p = inputParser;
    p.addParameter('MinPeakHeight', -Inf);
    p.parse(varargin{:});
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

function s = tern(c,a,b)
if c, s=a; else, s=b; end
end
