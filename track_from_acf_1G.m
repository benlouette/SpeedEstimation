function [tFrames, f_meas, f_viterbi, f_kal, FR] = track_from_acf_1G(acf_frames, Fs, tFrames, P)
% TRACK_FROM_ACF_1G
%  Suivi de la période \tau(t) à partir d'ACF déjà calculées (colonnes).
%  - Mesure par frame: YIN/CMNDF sur l'ACF (solution 1B adaptée)
%  - Emission HMM: peigne log (solution 1E) -> coût = -score
%  - Viterbi sur grille \tau (bins en échantillons) avec transitions limitées ±D
%  - Kalman sur \tau (état [tau; tau_dot]) -> lissage
%
% In:
%   acf_frames : [Nlag x Nf] ou [Nlag x 1]
%   Fs         : Hz (du signal original utilisé pour générer l'ACF)
%   tFrames    : (1 x Nf) temps centre des frames (s) [optionnel]
%   P          : struct paramètres (voir demo)
%
% Out:
%   tFrames    : temps par frame (s)
%   f_meas     : mesure brute par frame (Hz)
%   f_viterbi  : chemin Viterbi (Hz)
%   f_kal      : Kalman filtré (Hz)
%   FR         : struct array diag par frame

acf_frames = double(acf_frames);
if isvector(acf_frames), acf_frames = acf_frames(:); end
[Nlag, Nf] = size(acf_frames);

% if nargin<3 || isempty(tFrames)
    tFrames = 0:(Nf-1);
% else
%     tFrames = tFrames(:).';
%     if numel(tFrames) ~= Nf
%         error('tFrames doit avoir Nf éléments.');
%     end
% end

% ----- Defaults -----
D = struct('fmin',0.5,'fmax',120,'maxLagSec',1.0,'yinThresh',0.15,'smoothMs',1.0, ...
           'Kharm',5,'WeightMode','1/k','Eps',1e-3, ...
           'gridStepSamp',1,'maxJumpBins',6,'sigmaTrans',2.0,'octavePenalty',6.0, ...
           'kalman_dt',[],'kal_sigma_taudot',5e-3,'kal_sigma_meas',1.5e-3, ...
           'plot',true);
fn = fieldnames(D);
for k=1:numel(fn), if ~isfield(P,fn{k}), P.(fn{k})=D.(fn{k}); end, end

% ----- Bornes lag utiles -----
tauMin = max(2, floor(Fs/P.fmax));
tauMax = min(Nlag-2, ceil(Fs/P.fmin));
maxLagSamp = min(Nlag-1, round(P.maxLagSec*Fs));
tauMax = min(tauMax, maxLagSamp);
if tauMin+2 >= tauMax
    error('Plage de lag invalide: ajuster fmin/fmax/maxLagSec vs Nlag.');
end

% ----- Grille HMM sur tau (en échantillons) -----
gridStep = max(1, round(P.gridStepSamp));
tauGrid  = (tauMin:gridStep:tauMax).';
S = numel(tauGrid);

% Poids peigne
K = P.Kharm;
wk = strcmpi(P.WeightMode,'equal') * ones(1,K) + ~strcmpi(P.WeightMode,'equal')*(1./(1:K));
wk = wk / sum(wk);

% Préallocs
f_meas    = nan(1,Nf);
tau_meas  = nan(1,Nf);
FR        = repmat(struct('nsdf_peak',nan,'periodicity',nan,'tau',nan), 1, Nf);

% ----- Mesure et émission HMM par frame -----
EM = nan(S,Nf);   % coût d'émission = -score peigne log
for j = 1:Nf
    r = acf_frames(:,j);
    % normalisation (sécurisée)
%     r = r / max(r(1), eps);
    % lissage léger
    w = max(1, round(P.smoothMs*1e-3*Fs));
    if w>1, r_eval = movavg(r, w); else, r_eval = r; end

    % --- Mesure brute: YIN/CMNDF-from-ACF (solution 1B condensée) ---
    d = 2*(1 - r_eval(2:tauMax+1));
    cum = cumsum(d); tIdx = (1:tauMax).';
    dprime = d .* (tIdx ./ max(cum,eps));
    % masque bande
    dprime(1:tauMin-1) = +Inf;
    % seuil YIN
    idx = find(dprime < P.yinThresh, 1, 'first');
    if ~isempty(idx)
        Lc = idx; win = max(3, round(0.1*Lc));
        i1 = max(tauMin, Lc-win); i2 = min(tauMax, Lc+win);
        [~,rel] = min(dprime(i1:i2)); L = i1+rel-1;
    else
        [~,L] = min(dprime);
        if L<tauMin || L>tauMax, [~,L]=min(dprime(tauMin:tauMax)); L=L+tauMin-1; end
    end
    % parabolique
    if L>=2 && L<=tauMax-1
        y1=dprime(L-1); y2=dprime(L); y3=dprime(L+1);
        delta = 0.5*(y1-y3)/max(y1-2*y2+y3, eps);
        delta = max(min(delta,0.5),-0.5);
        Lref = L + delta;
    else
        Lref = L;
    end
    tau_meas(j) = Lref;
    f_meas(j)   = Fs / Lref;
    % qualité
    nsdf_peak   = 1 - (interp_lin(dprime, Lref)/2);
    FR(j).nsdf_peak = max(0,min(1,nsdf_peak));
    FR(j).tau       = Lref;

    % --- Emission HMM: peigne log (solution 1E) ---
    rpos = max(r_eval, 0);
    for s = 1:S
        Lg = tauGrid(s);
        ssum = 0; count = 0;
        for k=1:K
            q = k*Lg;
            if q>tauMax, break; end
            ssum = ssum + wk(k) * log(P.Eps + rpos(q+1)); % r(1) = lag0
            count = count + 1;
        end
        if count>=2
            EM(s,j) = -ssum;   % coût = -score
        else
            EM(s,j) = +Inf;
        end
    end
end

% ----- Viterbi sur tauGrid -----
Dmax = max(1, round(P.maxJumpBins));    % portée de transition en bins
sig2 = P.sigmaTrans.^2;
octPen = P.octavePenalty;

% initialisation
phi = EM(:,1);
bp  = zeros(S,Nf);  % backpointers
for j=2:Nf
    prev = phi;
    cur  = EM(:,j);
    newphi = +Inf(S,1);
    newbp  = ones(S,1);
    for s=1:S
        p1 = max(1, s-Dmax); p2 = min(S, s+Dmax);
        idx = p1:p2;
        % coût transition quadratique
        dBins = (s - idx);
        Tcost = (dBins.^2) / (2*sig2);
        % pénalité octave si ~ double/half
        tau_s   = tauGrid(s);
        tau_idx = tauGrid(idx);
        ratio   = tau_s ./ tau_idx;
        is_oct  = (abs(ratio-2) < 0.08) | (abs(ratio-0.5) < 0.08);
        Tcost(is_oct) = Tcost(is_oct) + octPen;
        % choisir meilleur prédécesseur
        [val, kmin] = min(prev(idx) + Tcost);
        cand = val + cur(s);
        if cand < newphi(s)
            newphi(s) = cand;
            newbp(s)  = idx(kmin);
        end
    end
    phi = newphi; bp(:,j) = newbp;
end
% remonter le chemin
[~, sN] = min(phi);
path = zeros(1,Nf); path(Nf) = sN;
for j=Nf-1:-1:1
    path(j) = bp(path(j+1), j+1);
end
tau_vit = tauGrid(path).';
f_viterbi = Fs ./ tau_vit;

% ----- Kalman sur tau (état [tau; tau_dot]) -----
if isempty(P.kalman_dt)
    if numel(tFrames)>=2
        dt = median(diff(tFrames));
        if dt<=0, dt = 1; end
    else
        dt = 1;
    end
else
    dt = P.kalman_dt;
end
tau_meas_s = tau_meas / Fs;   % en secondes
tau_kal_s  = kalman_tau(tau_meas_s, dt, P.kal_sigma_taudot, P.kal_sigma_meas);
tau_kal    = tau_kal_s * Fs;
f_kal      = Fs ./ tau_kal;

% ----- Plot optionnel: carte d'émission -----
if P.plot && Nf>1
    figure('Name','Carte coût d''émission (HMM)');
    imagesc(tFrames, tauGrid/Fs, EM); axis xy; colorbar;
    hold on;
    plot(tFrames, tau_meas/Fs, 'w.-', 'DisplayName','mesure \tau'); 
    plot(tFrames, tau_vit /Fs, 'c-', 'LineWidth',1.5, 'DisplayName','Viterbi \tau');
    plot(tFrames, tau_kal /Fs, 'm-', 'LineWidth',2.0, 'DisplayName','Kalman \tau');
    xlabel('Temps (s)'); ylabel('\tau (s)'); legend('Location','best');
    title('HMM emissions (plus sombre = meilleur)');
end
end

% ===== Helpers =====
function y = movavg(x,w)
if w<=1, y=x; return; end
k = ones(w,1)/w;
try
    y = filtfilt(k,1,x);
catch
    y = conv(x,k,'same');
end
end

function v = interp_lin(x, q)
% interpole sur échelle d'échantillons (x: vecteur, q: réel)
N = numel(x);
if q<1, v=x(1); return; end
if q>N, v=x(end); return; end
i0 = floor(q);
if i0>=N, v=x(N); return; end
a  = q - i0;
v = (1-a)*x(i0) + a*x(i0+1);
end

function tau_hat = kalman_tau(z, dt, sigma_taudot, sigma_meas)
% Kalman 2x2 pour \tau (s) avec modèle vitesse constante
A = [1 dt; 0 1];
Q = (sigma_taudot^2) * [dt^4/4 dt^3/2; dt^3/2 dt^2];
H = [1 0];
R = sigma_meas^2;

x = [z(1); 0];  % init: tau=z0, taudot=0
P = eye(2);
tau_hat = nan(size(z));
for k=1:numel(z)
    % pred
    x = A*x;
    P = A*P*A.' + Q;
    if isfinite(z(k)) && z(k)>0
        S = H*P*H.' + R;
        K = (P*H.')/S;
        x = x + K*(z(k) - H*x);
        P = (eye(2)-K*H)*P;
    end
    tau_hat(k) = x(1);
end
end
