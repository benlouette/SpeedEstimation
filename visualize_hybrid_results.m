function visualize_hybrid_results(x, Fs, P, OUT, fr_hat)
% Vue d'ensemble: bandes retenues + indicateurs par bande

figure('Name','Bandes impulsives (kurtosis)');
stairs([OUT.bands(:,1); OUT.bands(end,2)], [OUT.bandScores; OUT.bandScores(end)]);
grid on; xlabel('Hz'); ylabel('Kurtosis'); title('Score impulsivité par bande');
hold on; for i=1:size(OUT.bandsKeep,1)
    xline(OUT.bandsKeep(i,1),'r--'); xline(OUT.bandsKeep(i,2),'r--');
end
legend('kurtosis','bandes retenues');

% Résumé par bande
for i=1:numel(OUT.Band)
    B = OUT.Band(i);
    figure('Name',sprintf('Bande %d — [%.0f %.0f] Hz',i,B.f1,B.f2));

    % Enveloppe FFT
    subplot(2,2,1);
    if ~isempty(B.Fenv)
        plot(B.Fenv, B.Aenv); grid on;
        xlim([0 P.maxFr]); xlabel('Hz'); ylabel('|FFT(env)|');
        title(sprintf('Enveloppe: fr(env)=%.2f Hz, SNR=%.1f dB',B.fr_env,B.snr_env));
        xline(fr_hat,'r--','fr\_hat');
    else
        text(0.1,0.5,'(pas de FFT enveloppe)','Units','normalized');
        axis off;
    end

    % ACF enveloppe
    subplot(2,2,2);
    if ~isempty(B.acLag)
        plot(B.acLag, B.acCurve); grid on;
        xlim([0 1/P.minFr]); xlabel('Tau (s)'); ylabel('ACF');
        title(sprintf('ACF: fr(acf)=%.2f Hz', B.fr_acf));
        xline(1/max(fr_hat,eps),'r--','T\_hat');
    else
        text(0.1,0.5,'(pas d''ACF)','Units','normalized'); axis off;
    end

    % Cepstre enveloppe
    subplot(2,2,3);
    if ~isempty(B.qvec)
        plot(B.qvec, B.Cmag); grid on;
        xlim([0 1/P.minFr]); xlabel('Quefrency (s)'); ylabel('Cepstrum');
        title(sprintf('Cepstre: fr(cep)=%.2f Hz', B.fr_cep));
        xline(1/max(fr_hat,eps),'r--','T\_hat');
    else
        text(0.1,0.5,'(pas de cepstre)','Units','normalized'); axis off;
    end

    % Récap des estimations numériques
    subplot(2,2,4);
    vals = [B.fr_env, B.fr_acf, B.fr_cep, B.fr_sb, B.fr_csc];
    names = {'env','acf','cep','sb','csc'};
    bar(vals); grid on; set(gca,'XTickLabel',names);
    ylabel('Hz'); title(sprintf('Estims bande %d (fr\\_hat=%.2f Hz)',i,fr_hat));
    yline(fr_hat,'r--','fr\_hat');
end

% Spectre brut pour repère
figure('Name','Signal brut & spectre');
subplot(2,1,1); plot((0:numel(x)-1)/Fs, x); grid on; xlabel('Time (s)'); ylabel('x(t)'); title('Signal brut');
subplot(2,1,2);
N=2^nextpow2(numel(x)); X=fft(x.*hann(numel(x)),N);
F=(0:floor(N/2))/N*Fs; A=abs(X(1:floor(N/2)+1));
plot(F,20*log10(A+eps)); grid on; xlim([0 Fs/2]); xlabel('Hz'); ylabel('|X| (dB)'); title('Spectre brut');
end
