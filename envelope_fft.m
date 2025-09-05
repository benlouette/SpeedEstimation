function [F, A, fr_pk] = envelope_fft(env, Fs, frRange)
% FFT de l'enveloppe et pic principal dans [frRange]
N  = 2^nextpow2(numel(env));
w  = hann(numel(env));
E  = fft((env(:).*w), N);
A  = abs(E(1:floor(N/2)+1));
F  = (0:floor(N/2)).'/N*Fs;

% Chercher pic dans la plage basse fréquence
mask = (F>=frRange(1)) & (F<=frRange(2));
if ~any(mask)
    fr_pk = 0; return;
end
[Amax, idx] = max(A(mask)); %#ok<ASGLU>
fIdx = find(mask,1,'first') + idx - 1;
fr_pk = F(fIdx);
end
