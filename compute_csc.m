function [alphaVec_out, fVecValid, CSCabs, alphaScore] = compute_csc(Z, fVec, Fs, alphaVec_in, fBand)
% COMPUTE_CSC  Approximation de la Cyclic Spectral Coherence (CSC)
%   S_x(f; α) ≈ E_t[ Z(f+α/2, t) * conj(Z(f-α/2, t)) ]
%   CSC = |S_x| / sqrt( P(f+α/2) * P(f-α/2) ), moyenné sur t (trames STFT)
%
% In:
%   Z            : STFT complexe (nFreq x nFrames), voir stft_hann.m
%   fVec         : vecteur des fréquences (nFreq x 1), Hz
%   Fs           : fréquence d'échantillonnage (Hz) [non utilisé ici mais gardé pour cohérence]
%   alphaVec_in  : grille des fréquences cycliques α (Hz), e.g. min:step:max
%   fBand        : [f1 f2] bande (Hz) utilisée pour intégrer la CSC
%
% Out:
%   alphaVec_out : grille α (Hz) effectivement utilisée (identique à alphaVec_in)
%   fVecValid    : fréquences (Hz) où la CSC a été évaluée (dans la bande)
%   CSCabs       : |CSC(f, α)| de taille [numel(fVecValid) x numel(alphaVec_out)]
%   alphaScore   : score(α) = moyenne_f |CSC(f, α)| sur f ∈ fBand (NaN omis)

[nFreq, ~] = size(Z);
fVec = fVec(:);
df = fVec(2) - fVec(1);

% Puissance moyenne par fréquence (sur les trames)
Pbar = mean(abs(Z).^2, 2) + eps;

% Bande où intégrer la CSC
fMaskBand = (fVec >= fBand(1)) & (fVec <= fBand(2));
if ~any(fMaskBand), fMaskBand = true(size(fVec)); end
fBandVec = fVec(fMaskBand);

alphaVec = alphaVec_in(:).';
nAlpha = numel(alphaVec);

CSCabs = nan(numel(fBandVec), nAlpha);

for ia = 1:nAlpha
    a = alphaVec(ia);
    halfShiftBins = round((a/2) / df);   % décalage demi-alpha en bins
    if halfShiftBins < 1
        continue;
    end

    % Indices valides pour f ± a/2
    kmin = 1 + halfShiftBins;
    kmax = nFreq - halfShiftBins;
    if kmax <= kmin
        continue;
    end
    k = (kmin:kmax).';

    Zp = Z(k + halfShiftBins, :);   % f + a/2
    Zm = Z(k - halfShiftBins, :);   % f - a/2

    % Estimation S_x(f; α) = E_t[ Zp * conj(Zm) ]
    Sfa = mean(Zp .* conj(Zm), 2);

    % Normalisation
    denom = sqrt(Pbar(k + halfShiftBins) .* Pbar(k - halfShiftBins));
    CSC_line = Sfa ./ denom;

    % Interpole |CSC| sur la grille fBandVec
    CSCabs(:, ia) = interp1(fVec(k), abs(CSC_line), fBandVec, 'linear', NaN);
end

% Score(α) = moyenne sur f (en omettant les NaN)
alphaScore = mean(CSCabs, 1, 'omitnan');

alphaVec_out = alphaVec;
fVecValid    = fBandVec;
end
