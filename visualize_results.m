function visualize_results(x, Fs, P, D, fr_hat)
figure('Name','Time & Spectrum'); 
subplot(2,1,1); plot((0:numel(x)-1)/Fs, x); grid on; xlabel('Time (s)'); ylabel('x(t)');
title('Signal brut');
subplot(2,1,2); 
N=2^nextpow2(numel(x)); X=fft(x.*hann(numel(x)),N); 
F=(0:floor(N/2))/N*Fs; A=abs(X(1:floor(N/2)+1));
plot(F,20*log10(A+eps)); grid on; xlim([0 Fs/2]); xlabel('Hz'); ylabel('|X| dB'); 
title('Spectre brut');

figure('Name','Bands Kurtosis'); 
stairs([D.bands(:,1); D.bands(end,2)], [D.bandScores; D.bandScores(end)]); 
grid on; xlabel('Hz'); ylabel('Kurtosis'); title('Score d''impulsivité par bande');
hold on; for i=1:size(D.bandsKeep,1)
    xline(D.bandsKeep(i,1),'r--'); xline(D.bandsKeep(i,2),'r--');
end
legend('kurtosis par bande','bandes retenues');

% Détails par bande retenue
for i=1:numel(D.bandSummaries)
    S = D.bandSummaries{i};
    figure('Name',sprintf('Bande %d: [%.0f %.0f] Hz',i,S.f1,S.f2));
    subplot(3,1,1);
    plot(S.Fenv, S.Aenv); grid on; xlim([0 P.maxFr]); 
    xlabel('Hz'); ylabel('|FFT(env)|'); 
    title(sprintf('Enveloppe FFT — fr(env)=%.2f Hz', S.fr_env));
    xline(fr_hat,'r--','fr\_hat');

    subplot(3,1,2);
    plot(S.acLag, S.acCurve); grid on; 
    xlim([0 1/P.minFr]); xlabel('Tau (s)'); ylabel('ACF');
    title(sprintf('Autocorr enveloppe — fr(acf)=%.2f Hz', S.fr_acf));
    xline(1/fr_hat,'r--','T\_hat');

    subplot(3,1,3);
    plot(S.qvec, S.Cmag); grid on; 
    xlim([0 1/P.minFr]); xlabel('Quefrency (s)'); ylabel('Cepstrum');
    title(sprintf('Cepstre enveloppe — fr(cep)=%.2f Hz', S.fr_cep));
    xline(1/fr_hat,'r--','T\_hat');
end
end
