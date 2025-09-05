function [fr_acf, tau, R] = autocorr_envelope(env, Fs, frRange)
% ACF de l’enveloppe -> première période significative
env = env - mean(env);
R = xcorr(env, 'biased'); 
L = numel(env);
R = R(L:end); % tau >= 0
tau = (0:numel(R)-1).'/Fs;

% Cherche pic hors zéro delay
Tmin = 1/frRange(2); Tmax = 1/frRange(1);
mask = (tau>=Tmin) & (tau<=Tmax);
if ~any(mask), fr_acf=0; return; end
[~, idx] = max(R(mask));
iTau = find(mask,1,'first') + idx - 1;
fr_acf = 1 / tau(iTau);
end
