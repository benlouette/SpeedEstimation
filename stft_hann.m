function [Z, fVec, tVec] = stft_hann(x, Fs, winLen, hop, nfft)
% STFT simple avec fenêtre de Hann (periodic), sans padding temporel entre trames
% In :
%   x      : signal (vecteur)
%   Fs     : fréquence d'échantillonnage [Hz]
%   winLen : longueur de fenêtre [échantillons]
%   hop    : pas entre trames [échantillons]
%   nfft   : taille FFT
% Out:
%   Z    : spectrogramme complexe (nFreq x nFrames)
%   fVec : fréquences [Hz] (nFreq x 1)
%   tVec : temps centre de chaque trame [s] (1 x nFrames)

x = x(:);
N = numel(x);

% Assure au moins 1 trame
if N < winLen
    x = [x; zeros(winLen - N, 1)];
    N = winLen;
end

w = hann(winLen, 'periodic');
nFrames = 1 + floor((N - winLen) / hop);
if nFrames <= 0
    nFrames = 1;
end

Z = zeros(nfft/2+1, nFrames);
tVec = zeros(1, nFrames);

for m = 1:nFrames
    i1 = (m-1)*hop + 1;
    i2 = i1 + winLen - 1;
    frame = x(i1:i2) .* w;
    X = fft(frame, nfft);
    Z(:, m) = X(1:nfft/2+1);
    tVec(m) = ((i1 + i2)/2 - 1) / Fs; % temps du centre de trame
end

fVec = (0:nfft/2).' * (Fs/nfft);
end
