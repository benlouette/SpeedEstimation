function [T_cand, q, info] = initial_candidate_T_from_acf(acf, Fs, fmin, fmax, P)
% Seed YIN/CMNDF "léger" depuis l'ACF (lags >= 0, acf(1)=R(0)).
% Renvoie:
%   T_cand : période candidate (s) ou NaN si invalide
%   q      : qualité ~ 1 - d'(T)/2 (dans [0,1])
%   info   : struct debug (indices, bornes)

if nargin<5, P = struct(); end
if ~isfield(P,'SmoothMs'),   P.SmoothMs = 1.0;  end   % lissage léger
if ~isfield(P,'MaxLagSec'),  P.MaxLagSec= 1.0;  end   % borne sécurité
if ~isfield(P,'UseHalfFix'), P.UseHalfFix = true; end % anti-demi-tour simple

r = acf(:);
if r(1) ~= 0, r = r / r(1); end

% (1) lissage léger de l’ACF (stabilise CMNDF)
w = max(1, round(P.SmoothMs*1e-3*Fs));
if w>1
    k = ones(w,1)/w;
    try r = filtfilt(k,1,double(r)); catch, r = conv(double(r),k,'same'); end
    r = r(:);
end

% (2) bornes de recherche en lag
Nlag = numel(r);
tauMin = max(2, floor(Fs / fmax));                         % >=2 pour parabole
tauMax = min([Nlag-2, ceil(Fs / fmin), round(P.MaxLagSec*Fs)]);
if tauMin+2 >= tauMax, T_cand = NaN; q = 0; info=struct('reason','bad_range'); return; end

% (3) CMNDF sur la bande
d = 2*(1 - r(2:tauMax+1));
cum = cumsum(d);
tIdx = (1:tauMax).';
dprime = d .* (tIdx ./ max(cum, eps));
dprime(1:tauMin-1) = +Inf;

% (4) minimum global + parabole (avec garde convexité)
[~,L] = min(dprime);
if L>=2 && L<=tauMax-1
    y1=dprime(L-1); y2=dprime(L); y3=dprime(L+1);
    denom = (y1 - 2*y2 + y3);
    if denom > 0                           % convexité ok
        delta = 0.5*(y1 - y3) / denom;
        delta = max(min(delta, 0.5), -0.5);
    else
        delta = 0;                         % pas de raffinement si non convexe
    end
else
    delta = 0;
end
Lref = L + delta;

% (5) qualité et anti-demi-tour (optionnel, simple)
q = max(0,min(1, 1 - interp1(1:numel(dprime), dprime, Lref, 'linear','extrap')/2));
if P.UseHalfFix
    L2 = 2*Lref;
    if L2 <= tauMax
        d1 = interp1(1:numel(dprime), dprime, Lref, 'linear','extrap');
        d2 = interp1(1:numel(dprime), dprime, L2,   'linear','extrap');
        if d2 < 0.85*d1                      % 2T clairement meilleur
            Lref = L2;
        end
    end
end

T_cand = Lref / Fs;
info = struct('tauMin',tauMin,'tauMax',tauMax,'L',L,'Lref',Lref,'q',q);
end