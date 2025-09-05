function [fr_hat, OUT] = estimate_fr_ar_yw_from_acf(acf, Fs, P)
% Estimation de la rotation via AR (Yule–Walker) à partir de l'ACF.
% Étapes:
%  1) Choix d'ordre AR(p) (auto par AICc ou fixé)
%  2) Levinson–Durbin sur r(0..p) -> polynôme AR A(z), erreur E
%  3) Spectre AR: Pxx(f) = E / |A(e^{-j2π f/Fs})|^2
%  4) Recherche de la fondamentale: Harmonic Sum sur Pxx(f)
%  5) Anti-demi-tour: impairs vs pairs + test R(T/2) vs R(T) (sur l'ACF)

acf = acf(:);
Nlag = numel(acf);
r0 = max(acf(1), eps);
r  = acf / r0;  % normalisation sûre

% -------- Params défaut --------
DEF = struct('Order','auto','pmax',60,'Neff',[], ...
             'fmin',0.5,'fmax',120,'dfHz',0.1,'Kharm',5,'WeightMode','1/k', ...
             'lambdaEven',0.8,'gammaHalf',1.12,'SpecNfft',16384,'Plot',true);
fn = fieldnames(DEF); for i=1:numel(fn), if ~isfield(P,fn{i}), P.(fn{i})=DEF.(fn{i}); end, end

% -------- Choix ordre p --------
if isnumeric(P.Order) && isscalar(P.Order)
    p = min(P.Order, Nlag-1);
else
    pmax = min(P.pmax, Nlag-1);
    if isempty(P.Neff)
        Neff = max(4*Nlag, 2000);  % proxy raisonnable si inconnu
    else
        Neff = P.Neff;
    end
    aic_best = Inf; p = 8;  % init
    for pp = 4:pmax
        [A,E] = levinson_safe(r(1:pp+1), pp);    % r(1)=r0 normalisé
        if ~isfinite(E) || E<=0, continue; end
        % AICc ~ Neff*log(E) + 2*pp + 2*pp*(pp+1)/(Neff-pp-1)
        aic = Neff*log(E) + 2*pp + 2*pp*(pp+1)/max(Neff-pp-1,1);
        if aic < aic_best
            aic_best = aic; p = pp; E_best = E; A_best = A; %#ok<NASGU>
        end
    end
end

% -------- Coefs AR (si pas déjà pris) --------
if exist('A_best','var')
    A = A_best; E = E_best;
else
    [A,E] = levinson_safe(r(1:p+1), p);
end
if ~isfinite(E) || E<=0
    warning('E (erreur de prédiction) non valide — ajustez p/ACF.'); E = max(E, 1e-6);
end

% -------- Spectre AR --------
nfft = max(2048, 2^nextpow2(P.SpecNfft));
w = 2*pi*(0:nfft/2)'/nfft;         % [0..pi]
f_spec = w * (Fs/(2*pi));          % [0..Fs/2]
H = freqz(1, A, w);                % 1/A(e^{jw})
Pxx = (E) ./ max(abs(H).^2, eps);  % PSD_AR (échelle relative)

% -------- Grille des fondamentales --------
fmin = max(P.fmin, 0.1);
fmax = min(P.fmax, f_spec(end)-1);
fgrid = (fmin:P.dfHz:fmax).';
K = P.Kharm;
wk = strcmpi(P.WeightMode,'equal')*ones(1,K) + ~strcmpi(P.WeightMode,'equal')*(1./(1:K));
wk = wk / sum(wk);

% Interpolation PSD_AR
Spec = @(fg) interp1(f_spec, Pxx, fg, 'linear', 'extrap');

% -------- Score Harmonic Sum --------
S = zeros(size(fgrid));
for i=1:numel(fgrid)
    f0 = fgrid(i);
    s = 0; cnt = 0;
    for k=1:K
        fk = k*f0;
        if fk > f_spec(end), break; end
        s = s + wk(k)*Spec(fk);
        cnt = cnt + 1;
    end
    if cnt>=2, S(i) = s; else, S(i) = -Inf; end
end

% Pic + raffinement parabolique
[~, idx] = max(S);
if idx>1 && idx<numel(S)
    y1=S(idx-1); y2=S(idx); y3=S(idx+1);
    delta = 0.5*(y1 - y3) / max(y1 - 2*y2 + y3, eps);
    delta = max(min(delta, 0.5), -0.5);
else
    delta = 0;
end
f_hat = fgrid(idx) + delta*P.dfHz;

% -------- Anti-demi-tour --------
% Odd vs even (sur PSD_AR)
Sodd = 0; Seven = 0;
for k=1:K
    fk = k*f_hat; if fk>f_spec(end), break; end
    if mod(k,2)==1, Sodd = Sodd + wk(k)*Spec(fk);
    else            Seven= Seven+ wk(k)*Spec(fk);
    end
end
combContrast = (Sodd - P.lambdaEven*Seven) / max(Sodd + P.lambdaEven*Seven, eps);

% Comparaison ACF: R(T/2) vs R(T)
T_hat_samp = Fs / max(f_hat, eps);
R_T  = interp_acf(r, T_hat_samp);
R_T2 = interp_acf(r, T_hat_samp/2);
isHalf = (R_T2 > P.gammaHalf*R_T) && (Seven > Sodd);
if isHalf
    f_hat = f_hat/2;
end

fr_hat = f_hat;

% -------- Sorties --------
OUT = struct();
OUT.p = p;
OUT.A = A; OUT.E = E;
OUT.f_spec = f_spec; OUT.Pxx = Pxx;
OUT.fgrid = fgrid; OUT.S = S;
OUT.R_T = R_T; OUT.R_T2 = R_T2;
OUT.combContrast = combContrast; OUT.isHalf = isHalf;
OUT.params = P;

end

% ========= Helpers =========
function [A,E] = levinson_safe(r, p)
% r: [r0 r1 ... rp], p: ordre
% Retourne A(z)=1 + a1 z^-1 + ... + ap z^-p  (convention MATLAB), E: erreur
try
    [A,E] = levinson(r, p);
catch
    % Fallback Levinson–Durbin
    A = zeros(p+1,1); A(1)=1; E = r(1);
    if E<=0, E = eps; end
    k = zeros(p,1);
    for m=1:p
        if m==1
            lambda = -r(2)/E;
        else
            lambda = -(r(m+1) + A(2:m).'*flipud(r(2:m))) / E;
        end
        k(m) = lambda;
        Aold = A(1:m);
        A(2:m) = A(2:m) + lambda*flipud(Aold(1:m-1));
        A(m+1) = lambda;
        E = E*(1 - lambda^2);
        if E<=eps, E=eps; end
    end
end
A = A(:);  % colonne
end

function val = interp_acf(r, q)
% Interpolation linéaire de l'ACF normalisée aux indices fractionnaires q (en échantillons)
N = numel(r);
if q < 1 || q > N-2, val = NaN; return; end
i0 = floor(q); a = q - i0;
val = (1-a)*r(i0+1) + a*r(i0+2);
end

function s = tern(c,a,b)
if c, s=a; else, s=b; end
end
