function [fr_sb, bestScore] = sideband_grid_score(y, Fs, fCenter, winHz, frRange)
% Cherche un espacement de bandes latérales régulier autour d'une porteuse
% en maximisant un "peigne" de corrélation sur le spectre brut
N   = 2^nextpow2(numel(y));
w   = hann(numel(y));
Y   = fft(y(:).*w, N);
A   = abs(Y(1:floor(N/2)+1));
F   = (0:floor(N/2)).'/N*Fs;

% Fenêtre autour de la porteuse
maskC = (F>=fCenter-winHz) & (F<=fCenter+winHz);
Fc = F(maskC); Ac = A(maskC);
if ~any(maskC), fr_sb=0; bestScore=0; return; end

% grille d'espacements (fr) à tester
frGrid = linspace(frRange(1), frRange(2), 200);
scores = zeros(size(frGrid));

for i = 1:numel(frGrid)
    fr = frGrid(i);
    % positions de sidebands: fCenter +/- k*fr
    kmax = floor(winHz / fr);
    if kmax<1, scores(i)=0; continue; end
    fLines = fCenter + (-kmax:kmax)*fr;
    % somme des amplitudes près de chaque ligne
    score = 0;
    bwBin = max(1, round(1* numel(Fc)/(2*winHz))); % tolérance ~1 bin
    for fL = fLines
        [~, idx] = min(abs(Fc - fL));
        i1 = max(1, idx-bwBin); i2 = min(numel(Ac), idx+bwBin);
        score = score + max(Ac(i1:i2));
    end
    scores(i) = score / (2*kmax+1);
end

[bestScore, iBest] = max(scores);
fr_sb = frGrid(iBest);
end
