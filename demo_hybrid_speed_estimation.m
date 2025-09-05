% DEMO — Estimation tachless hybride (enveloppe + sidebands + CSC)
clear; close all; clc;

% ---------- Synthèse d’un signal de test ----------
Fs  = 51200;         % Hz
T   = 10;            % s
t   = (0:1/Fs:T-1/Fs).';
fr0 = 17.3;          % Hz, vitesse "vraie" pour validation
f0  = 6500;          % Hz, résonance
Q   = 25;            % ~ largeur de bande f0/Q
bw  = f0/Q;

% Utilise ta propre donnée en remplaçant la ligne suivante
x = generate_bearing_like_signal(t, Fs, fr0, f0, bw);

% ---------- Paramètres ----------
P = struct();
P.searchBand     = [500 20000];   % zone de recherche des résonances
P.nBands         = 24;            % nb sous-bandes pour la kurtosis
P.topK           = 2;             % nb de bandes retenues
P.minFr          = 0.2;           % Hz
P.maxFr          = 100;           % Hz
P.sidebandWinHz  = 600;           % fenêtre autour de la porteuse pour sidebands
P.plotting       = true;

% CSC (nécessaire au mode hybride)
P.alphaStep      = 0.2;           % pas en Hz pour la grille alpha
P.stft.winLen    = 4096;
P.stft.hop       = 1024;
P.stft.nfft      = 8192;

% ---------- Estimation ----------
[fr_hat, OUT] = estimate_speed_hybrid(x, Fs, P); %#ok<ASGLU>
fprintf('\nHybride: fr_hat = %.3f Hz (vrai = %.3f Hz)\n', fr_hat, fr0);
