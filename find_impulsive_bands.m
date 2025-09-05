function [bands, scores] = find_impulsive_bands(x, Fs, searchBand, nBands)
% Découpe [fmin,fmax] en nBands sous-bandes et score par kurtosis
fmin = max(10, searchBand(1));
fmax = min(Fs/2-100, searchBand(2));
edges = logspace(log10(fmin), log10(fmax), nBands+1);
bands = [edges(1:end-1).' edges(2:end).'];
scores = zeros(nBands,1);

for k = 1:nBands
    f1 = bands(k,1); f2 = bands(k,2);
    y  = bandpass_iir(x, Fs, f1, f2, 4);
    scores(k) = kurtosis(y); % impulsivité (simple et efficace)
end
end
