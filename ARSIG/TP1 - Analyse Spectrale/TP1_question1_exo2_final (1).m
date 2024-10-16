% execrice 2
close all;
Nf = 4096;   N = 201;
axefreq = 0:Fe/Nf:(Nf-1)/Nf*Fe;
%fenêtre
figure()
fenetre_tf = fft(fen,Nf);
plot(axefreq,abs(fenetre_tf))
title(' TF du fênetre ')
xlabel('Fréquence Fe')
ylabel('|Fen(f)|')

%periodogramme signal
figure()
subplot 121

p_x = 1/N*abs(fft(x,Nf)).^2;
plot(axefreq,abs(p_x))
title(' Périodogramme du signal x ')
xlabel('Fréquence (cycles/jour)')
ylabel('|X(f)|')
subplot 122
plot(x)
title('Signal x ')
xlabel('Temps (jour)')
ylabel('vitesse radiale (km/s)')

%methode Clean
periodogramme = 1/N*(abs(fft(x, Nf))).^2;
residu = x;
k = 0;
while(k<8)
    TF_residu = fft(residu, Nf);
    [max_y, id_max] = max(periodogramme);
    frequence_max = axefreq(id_max);
    [amp_est,phi_est] = estim_amp_phase(residu,t,frequence_max);
    contribuition_estime = sin(2*pi*frequence_max*t+phi_est).*fen ;
    residu = residu - amp_est.*contribuition_estime;
    TF_residu_new = fft(residu, Nf);
    periodogramme = 1/N*(abs(TF_residu_new)).^2;

    %résidu avant et apres
    figure();
    subplot(131);
    plot(axefreq,abs(TF_residu));
    set(gca,'xlim',[0,Fe/2]);

    hold on
    plot(axefreq,abs(TF_residu_new), '-r');
    title('TF du résidu');
    legend('Résidu avant','Résidu aprés')
    xlabel('Fréquence (cycles par jour)'), ylabel('|résidu(f)|')
    hold off

    %Résidu domaine temporelle
    subplot(132);
    plot(residu);
    title('Résidu domaine temporelle')
    xlabel('Temps (jour)'), ylabel('|résidu(f)|')

    %composant detectée
    subplot(133);
    plot(axefreq,abs(fft(contribuition_estime,Nf)), '-r');
    set(gca,'xlim',[0,Fe/2]);
    title('TF composante detecté');
    xlabel('Fréquence (cycles par jour)'), ylabel('|résidu(f)|')
    k = k + 1;
end