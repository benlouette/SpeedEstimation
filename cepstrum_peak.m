function [fr_cep, qvec, Cmag] = cepstrum_peak(env, Fs, frRange)
% Cepstre (log-FFT inverse) de l’enveloppe pour trouver la période
N  = 2^nextpow2(numel(env));
w  = hann(numel(env));
E  = fft(env(:).*w, N);
logS = log(abs(E)+eps);
c   = ifft(logS, 'symmetric');

% Quefrency
qvec = (0:numel(c)-1).'/Fs;

% Cherche pic dans [Tmin, Tmax]
Tmin = 1/frRange(2); Tmax = 1/frRange(1);
mask = (qvec>=Tmin) & (qvec<=Tmax);
Cmag = c;
if ~any(mask), fr_cep=0; return; end
[~, idx] = max(Cmag(mask));
iQ = find(mask,1,'first') + idx - 1;
fr_cep = 1/qvec(iQ);
end
