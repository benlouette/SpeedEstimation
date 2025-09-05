function x = generate_bearing_like_signal(t, Fs, fr, f0, bw)
% Synthèse: résonance à f0 excitée par impacts modulés à fr
% bw ~ largeur de bande (Hz)
N   = numel(t);
x   = 0.02*randn(N,1);     % bruit

% Porteuse résonante (IIR peaking)
y = bandpass_iir(randn(N,1), Fs, max(10,f0-bw/2), f0+bw/2, 4);

% Enveloppe AM impulsive: trains d'impacts à fr
period = round(Fs/fr);
imp = zeros(N,1);
idx = 1:period:N;
imp(idx) = 1;
% rendre les impacts un peu "réalistes"
h = gausswin(round(0.004*Fs)); h = h/sum(h);
env = conv(imp, h, 'same');

% Modulation
x = x + 0.8 * (env .* y);

% Option: sidebands supplémentaires (harmoniques)
x = x + 0.3 * (env .* bandpass_iir(randn(N,1), Fs, max(10,f0-1.5*bw), f0+1.5*bw, 4));
end
