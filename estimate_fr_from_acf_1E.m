function [fr_hat, OUT] = estimate_fr_from_acf_1E(acf, Fs, P)
% Estimation de la frequence de rotation via Peigne-produit (log-domain)
%  ACF-only:
%    - r = ACF normalisee R(0)=1
%    - Score peigne S(T) = sum_k w_k * log(eps + max(0, r(kT)))
%    - Recherche T in [Fs/fmax .. Fs/fmin], max de S(T)
%    - Interpolation parabolique sur S(T) (sub-sample)
%    - Anti demi-tour: (R(T/2) > gamma*R(T) & pairs>impairs) => T <- T/2
%
% In:
%   acf : vecteur ACF (lags >= 0, acf(1) = R(0))
%   Fs  : frequence d'echantillonnage (Hz)
%   P   : struct params (voir demo)
%
% Out:
%   fr_hat : frequence estimee (Hz)
%   OUT    : diagnostics

acf = acf(:);
Nlag = numel(acf);
lags = (0:Nlag-1).';
tau  = lags / Fs;

% ---------- params par defaut ----------
def = struct('MaxLagSec',1.0,'fmin',0.5,'fmax',120,'Kharm',6, ...
             'WeightMode','1/k','Eps',1e-3,'SmoothMs',1.0,'PeakRefine',true, ...
             'lambdaEven',0.8,'gammaHalf',1.12,'Plot',true);
fn = fieldnames(def);
for k=1:numel(fn)
    if ~isfield(P, fn{k}), P.(fn{k}) = def.(fn{k}); end
end

% ---------- normalisation ACF ----------
R0 = acf(1);
r  = (R0~=0) * (acf / R0) + (R0==0) * acf;

% ---------- lissage leger (optionnel) ----------
w = max(1, round(P.SmoothMs*1e-3*Fs));
if w>1, r_eval = movavg(r, w); else, r_eval = r; end

% ---------- bornes en lag ----------
tauMin = max(2, floor(Fs / P.fmax)); % >=2 pour parabolique si refine
tauMax = min(Nlag-2, ceil(Fs / P.fmin));
maxLagSamp = min(Nlag-1, round(P.MaxLagSec*Fs));
tauMax = min(tauMax, maxLagSamp);
if tauMin+2 >= tauMax
    error('Plage invalide: ajuste fmin/fmax/MaxLagSec.');
end

% ---------- poids w_k ----------
K = P.Kharm;    % Nombre d’harmoniques considérées
switch lower(string(P.WeightMode))  % P.WeightMode = "equal" → tous les harmoniques pèsent autant.
    case "equal", wk = ones(1,K);
    otherwise,    wk = 1./(1:K);    % '1/k' par defaut Sinon (par défaut "1/k") → les harmoniques sont pondérés par 1/k.
        %Cela réduit progressivement l’importance des harmoniques lointains, plus sensibles au bruit et à la décroissance de l’ACF.
end
% Normalise les poids pour que leur somme fasse 1.
wk = wk / sum(wk);

% ---------- score peigne S(T) sur la grille entiere ----------
% On utilise interpolation lineaire aux positions fractionnaires k*T plus tard
% Pour l'eval discrete (entiers), on prend r_eval(k*L) aux indices entiers.
S = -inf(tauMax,1);% On prépare un vecteur S(L) (score pour une période candidate T=L échantillons)
%Initialisé à −Inf pour distinguer les L non évaluables (ou trop pauvres en harmoniques).
rpos = max(r_eval, 0);          % partie positive uniquement (AM -> pics positifs)
% On écrase les valeurs négatives de l’ACF lissée (r_eval) à 0.
% Intuition : pour une enveloppe AM (pics positifs attendus), on ne veut pas que des valeurs négatives “cassent” le produit (ou le log).
% Effet : on favorise les périodes où r(kL) est franchement positif.

logeps = @(x) log(P.Eps + x);   % Petite constante Eps pour éviter log(0) → stabilité numérique.

% for L = tauMin:tauMax
%     s = 0; count = 0;
%     for k=1:K
%         q = k*L;
%         if q>tauMax, break; end
%         s = s + wk(k) * logeps( rpos(q+1) );  % +1 car r(1)=lag0
%         count = count + 1;
%     end
%     if count >= 2     % au moins deux termes pour robustesse
%         S(L) = s;
%     end
% end
for L = tauMin:tauMax
    kmax = min(K, floor(tauMax / L));
    if kmax < 2, continue; end

    used   = 1:kmax;
    w_used = wk(used);
    w_used = w_used / sum(w_used);        % renormalisation locale

    idx  = used .* L + 1;                 % indices k*L (+1 car r(1)=lag0)
    vals = rpos(idx);                     % r(kL) >= 0

    % --- trois variantes équivalentes qui donnent un scalaire ---
    % 1) produit scalaire (robuste aux formes)
    s = dot(w_used(:), logeps(vals(:)));
    % 2) somme après alignement en colonne
    % s = sum( w_used(:) .* logeps(vals(:)) );
    % 3) w_used ligne × vals colonne -> scalaire
    % s = w_used(:).' * logeps(vals(:));

    S(L) = s;                             % <-- maintenant s est scalaire
end
% ---------- maximum grossier ----------
[~, L0] = max(S);
if ~isfinite(S(L0))
    % secours : restreindre K ou revenir a 1er pic ACF
    [~, L0] = max(r_eval(tauMin:tauMax)); L0 = L0 + tauMin - 1;
end

% ---------- raffinement parabolique sur S(L-1),S(L),S(L+1) ----------
if P.PeakRefine && L0>tauMin && L0<tauMax
    % evalue S a L0-1 et L0+1 en re-evaluant (deja calcule)
    y1 = S(L0-1); y2 = S(L0); y3 = S(L0+1);
    if isfinite(y1) && isfinite(y2) && isfinite(y3) && (y1 - 2*y2 + y3) < 0
        delta = 0.5*(y1 - y3) / max(y1 - 2*y2 + y3, eps);
        delta = max(min(delta, 0.5), -0.5);
    else
        delta = 0;
    end
    Lref = L0 + delta;
else
    Lref = L0;
end
T_hat = Lref / Fs;

% ---------- Anti-demi-tour ----------
interpR = @(q) safe_interp((0:Nlag-1).', r_eval, q);
T  = T_hat*Fs;                % en echantillons (fractionnaire)
R_T  = interpR(T);
R_T2 = interpR(T/2);

% peigne odd/even (log-domain aussi)
Sodd = 0; Seven = 0;
for k=1:K
    q = k*T;
    if q <= tauMax
        val = max(interpR(q), 0);
        if mod(k,2)==1, Sodd = Sodd + wk(k)*logeps(val);
        else            Seven= Seven+ wk(k)*logeps(val);
        end
    end
end
combContrast = (Sodd - P.lambdaEven*Seven) / max(abs(Sodd) + P.lambdaEven*abs(Seven), eps);

isHalf = (R_T2 > P.gammaHalf*R_T) && (Seven > Sodd);
if isHalf
    T_hat = T_hat/2;
end

fr_hat = 1 / max(T_hat, eps);

% ---------- sorties & plots ----------
OUT = struct();
OUT.S = S; OUT.tauGrid = (0:numel(S)-1).'/Fs; OUT.L0 = L0; OUT.Lref = Lref;
OUT.T_hat = T_hat; OUT.periodicity = R_T;
OUT.combContrast = combContrast; OUT.isHalf = isHalf; OUT.params = P;

if P.Plot
    % Courbe de score S(T)
    figure('Name','Peigne-produit (log) — Score S(T)');
    Taxis = (0:numel(S)-1).'/Fs;
    plot(Taxis, S, 'b'); grid on; hold on;
    xline(T_hat, '--', sprintf('T\\_hat=%.4f s  (f=%.2f Hz)', T_hat, fr_hat), 'Color',[0.2 0.2 0.2]);
    xlabel('T (s)'); ylabel('S(T) = \\sum w_k log(\\epsilon + R_+(kT))'); xlim([0, P.MaxLagSec]);
    title(sprintf('Score peigne-produit (K=%d, %s)', K, P.WeightMode));

    % ACF normalisee avec T et T/2
    figure('Name','ACF normalisee');
    plot(tau, r_eval, 'b'); grid on; hold on;
    xline(T_hat, '--', 'T'); xline(T_hat/2, ':', 'T/2');
    xlabel('\\tau (s)'); ylabel('ACF norm.'); xlim([0, P.MaxLagSec]);
    title('ACF normalisee & marquage T');

    % Odd vs Even (log-sum)
    figure('Name','Odd vs Even (log-sum)');
    bar([Sodd, Seven]); grid on; set(gca,'XTickLabel',{'odd','even'});
    ylabel('Somme log(\\epsilon + R_+)'); 
    title(sprintf('Contraste (odd-even)=%.3f  %s', combContrast, tern(isHalf,'\\rightarrow demi-tour corrige','')));
end
end

% ===== helpers =====

function y = movavg(x,w)
if w<=1, y=x; return; end
k = ones(w,1)/w;
try
    y = filtfilt(k,1,double(x));
catch
    y = conv(double(x), k, 'same');  % fallback si pas de Signal Toolbox
end
end

function v = safe_interp(x, y, q)
% interpolation lineaire aux indices fractionnaires q (en echantillons)
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
