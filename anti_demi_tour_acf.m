function [f_corr, DEC] = anti_demi_tour_acf(acf, Fs, T_cand, P)
% Stratégie anti-demi-tour (ACF-only) : décide entre T, T/2 (et optionnellement 2T)
% à partir d'un candidat T_cand provenant de n'importe quel estimateur.
%
% In:
%   acf   : vecteur ACF (lags >= 0), acf(1) = R(0)
%   Fs    : Hz
%   T_cand: période candidate (s)
%   P     : struct paramètres (voir demo)
%
% Out:
%   f_corr : fréquence corrigée (Hz)
%   DEC    : diagnostics (scores, décisions, etc.)

% ----- défauts -----
D = struct('Kharm',5,'WeightMode','1/k','lambdaEven',0.8,'gammaHalf',1.12, ...
           'UseLogComb',true,'SpecBW_Hz',1.0,'wComb',0.45,'wSpec',0.40,'wRratio',0.15, ...
           'considerDouble',false,'Plot',true);
fn = fieldnames(D);
for k=1:numel(fn), if ~isfield(P,fn{k}), P.(fn{k})=D.(fn{k}); end, end

r = acf(:) / max(acf(1), eps);
Nlag = numel(r);
tauMax = Nlag-1;

% --- helpers ---
interpR = @(q) safe_interp((0:Nlag-1).', r, q);

% --- Energies peigne ACF (odd vs even) ---
function [scoreComb, Sodd, Seven] = comb_score(r, T, K, wk, lambdaEven, useLog)
    Sodd = 0; Seven = 0;
    for k = 1:K
        q = k*T;
        if q > tauMax, break; end
        val = interpR(q);
        if useLog, val = log(1e-3 + max(val,0)); else, val = max(val,0); end
        if mod(k,2)==1, Sodd  = Sodd  + wk(k)*val;
        else            Seven = Seven + wk(k)*val;
        end
    end
    scoreComb = (Sodd - lambdaEven*Seven) / max(abs(Sodd) + lambdaEven*abs(Seven), eps);
end

% --- Spectre via ACF (Wiener–Khinchin) ---
[Fspec, Pxx] = psd_from_acf(r, Fs);
BW = max(P.SpecBW_Hz, Fspec(2)-Fspec(1));   % largeur mini ~ 1 bin
Spec = @(f0) band_energy(Fspec, Pxx, f0, BW);

% Poids harmoniques
K = P.Kharm;
wk = strcmpi(P.WeightMode,'equal')*ones(1,K) + ~strcmpi(P.WeightMode,'equal')*(1./(1:K));
wk = wk / sum(wk);

% -------------------- Scores pour T et T/2 --------------------
T  = T_cand * Fs;     % en échantillons (fractionnaire)
T2 = T/2;

% (1) Peigne ACF (odd-even)
[sComb_T , Sodd_T , Seven_T ] = comb_score(r, T , K, wk, P.lambdaEven, P.UseLogComb);
[sComb_T2, Sodd_T2, Seven_T2] = comb_score(r, T2, K, wk, P.lambdaEven, P.UseLogComb);

% (2) Peigne spectral (odd-even)
f0  = Fs / max(T ,eps);
f02 = Fs / max(T2,eps);
Sodd_spec_T = 0; Seven_spec_T = 0;
Sodd_spec_T2= 0; Seven_spec_T2= 0;
for k=1:K
    if mod(k,2)==1
        Sodd_spec_T  = Sodd_spec_T  + wk(k)*Spec(k*f0);
        Sodd_spec_T2 = Sodd_spec_T2 + wk(k)*Spec(k*f02);
    else
        Seven_spec_T  = Seven_spec_T  + wk(k)*Spec(k*f0);
        Seven_spec_T2 = Seven_spec_T2 + wk(k)*Spec(k*f02);
    end
end
sSpec_T  = (Sodd_spec_T  - P.lambdaEven*Seven_spec_T ) / max(Sodd_spec_T  + P.lambdaEven*Seven_spec_T , eps);
sSpec_T2 = (Sodd_spec_T2 - P.lambdaEven*Seven_spec_T2) / max(Sodd_spec_T2 + P.lambdaEven*Seven_spec_T2, eps);

% (3) Ratio direct R(T) vs R(T/2)
RT  = interpR(T);
RT2 = interpR(T2);
sRatio_T  = (RT  - P.gammaHalf*RT2);   % >0 si T préférable
sRatio_T2 = (RT2 - P.gammaHalf*RT );   % >0 si T/2 préférable

% Normalisation douce des scores ratio ([-1,1] approx)
nr = max(abs([sRatio_T sRatio_T2])) + eps;
sRatio_T  = sRatio_T  / nr;
sRatio_T2 = sRatio_T2 / nr;

% Score global
S_T  = P.wComb*sComb_T  + P.wSpec*sSpec_T  + P.wRratio*sRatio_T;
S_T2 = P.wComb*sComb_T2 + P.wSpec*sSpec_T2 + P.wRratio*sRatio_T2;

decision = 'keep T';
T_best = T;
if S_T2 > S_T
    decision = 'use T/2';
    T_best = T2;
end

% (option) considérer 2T
S_T2x = -Inf;
if P.considerDouble
    T2x = 2*T;
    [sComb_T2x, ~, ~] = comb_score(r, T2x, K, wk, P.lambdaEven, P.UseLogComb);
    f2x = Fs / max(T2x,eps);
    Sodd_spec_T2x=0; Seven_spec_T2x=0;
    for k=1:K
        if mod(k,2)==1, Sodd_spec_T2x = Sodd_spec_T2x + wk(k)*Spec(k*f2x);
        else            Seven_spec_T2x= Seven_spec_T2x+ wk(k)*Spec(k*f2x);
        end
    end
    sSpec_T2x = (Sodd_spec_T2x - P.lambdaEven*Seven_spec_T2x) / ...
                max(Sodd_spec_T2x + P.lambdaEven*Seven_spec_T2x, eps);
    RT2x = interpR(T2x);
    sRatio_T2x = (RT2x - P.gammaHalf*RT) / (abs(RT2x - P.gammaHalf*RT) + eps);
    S_T2x = P.wComb*sComb_T2x + P.wSpec*sSpec_T2x + P.wRratio*sRatio_T2x;
    if S_T2x > max(S_T,S_T2)
        decision = 'use 2T';
        T_best = T2x;
    end
end

f_corr = Fs / max(T_best, eps);

% ----- Diagnostics & plots -----
DEC = struct();
DEC.sComb  = [sComb_T  sComb_T2  ];
DEC.sSpec  = [sSpec_T  sSpec_T2  ];
DEC.sRatio = [sRatio_T sRatio_T2 ];
DEC.S_all  = [S_T S_T2 S_T2x];
DEC.RT = RT; DEC.RT2 = RT2;
DEC.decisionStr = decision;
DEC.Kharm = P.Kharm; DEC.lambdaEven = P.lambdaEven;

if P.Plot
    tau = (0:Nlag-1).'/Fs;
    figure('Name','ACF & candidats');
    plot(tau, r); grid on; hold on;
    xline(T/Fs,  '--', sprintf('T  (%.2f Hz)', f0));
    xline(T2/Fs, ':',  sprintf('T/2(%.2f Hz)', f02));
    if P.considerDouble, xline(2*T/Fs, ':', '2T'); end
    xlabel('\tau (s)'); ylabel('ACF norm.'); xlim([0, min(1, (Nlag-1)/Fs)]);
    title(sprintf('Anti-demi-tour — décision: %s', decision));

    figure('Name','Scores (T vs T/2)');
    cats = categorical({'comb','spec','ratio'}); cats=reordercats(cats,{'comb','spec','ratio'});
    bar(cats, [ [sComb_T; sSpec_T; sRatio_T], [sComb_T2; sSpec_T2; sRatio_T2] ] );
    grid on; legend('T','T/2'); ylabel('score normalisé'); title('Votes des tests');
end
end

% ================= Helpers =================
function [F, Pxx] = psd_from_acf(r, Fs)
% PSD via Wiener–Khinchin (ACF -> FFT), ACF déjà normalisée
rsym = [flipud(r(2:end)); r(:)];
N = numel(rsym);
w = 0.5 - 0.5*cos(2*pi*(0:N-1)'/(N-1)); % Hann
S = fft(rsym .* w);
Pfull = real(S);
P = Pfull(1:floor(N/2)+1);
F = (0:floor(N/2))' * (Fs/N);
Pxx = max(P, 0);
end

function E = band_energy(F, Pxx, f0, bw)
if f0<=0 || f0>F(end), E=0; return; end
m = (F >= (f0-bw)) & (F <= (f0+bw));
if ~any(m), E=0; return; end
E = trapz(F(m), Pxx(m));
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
