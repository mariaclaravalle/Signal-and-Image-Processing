close all

Fs = 1;
N_options = {100, 100, 300, 100};

x_signals = {x1, x2, x3, x4};


for k = 1:4
    N = N_options{k};
    Nf = 4*2^nextpow2(N); %Nf>>N
    axefreq = 0:Fs/Nf:(Nf-1)/Nf*Fs;

    %Temps
    x = x_signals{k};
    figure()
    plot(x)
    title(['Signal x', num2str(k),' temporel'])
    xlabel('Nombre déchantillons N')
    ylabel(['x',num2str(k)])

     %Périodogramme standard
    TF_x = 1/N*abs(fft(x, Nf)).^2;
    figure()
    plot(axefreq, abs(TF_x))
    set(gca,'yscale','log') %Pour tracer dans l'echelle log
    set(gca,'xlim',[0,Fs/2]) %Restringe l'affichage à première motié
    title(['Périodogramme du signal x', num2str(k), ' standard'])
    xlabel('Fréquence lambda')
    ylabel(['log(x',num2str(k), ')'])

    %Périodogramme fenêtré 

    figure()
    w = hamming(N);
    TF_x_h = 1/N*abs(fft(x.*w, Nf)).^2;
    plot(axefreq, abs(TF_x_h))
    set(gca,'yscale','log')
    set(gca,'xlim',[0,Fs/2])
    title(['Périodogramme du signal x', num2str(k), ' fenêtre de Hamming'] )
    xlabel('Fréquence lambda')
    ylabel(['x',num2str(k)])

    %Périodogramme Welch
    
    Nwin = 60; %60
    Noverlap = 2; %2
    
    periodoWelch = pwelch(x,Nwin,Noverlap,Nf,Fs,'twosided');
    figure()
    plot(axefreq, abs(periodoWelch))
    set(gca,'yscale','log')
    set(gca,'xlim',[0,Fs/2])
    title(['Périodogramme du signal x', num2str(k), ' Welch'])
    xlabel('Fréquence lambda')
    ylabel(['log(x',num2str(k),')'])

    %Modélisation AR

    P = 7;
    [a,sigma2] = arcov(x,P);
    S_AR = sigma2./abs((fft(a,Nf))).^2;
    figure()
    plot(axefreq, S_AR);
    %set(gca,'yscale','log');
    set(gca,'xlim',[0,Fs/2])
    title(['Périodogramme du signal x',num2str(k), ' autorégressive (AR)'])
    xlabel('Fréquence lambda')
    ylabel(['x',num2str(k)])

    %Modélisation MUSIC

    S_MUSIC = pmusic(x,P,axefreq,Fs);
    figure()
    plot(axefreq, S_MUSIC);
    %set(gca,'yscale','log');
    set(gca,'xlim',[0,Fs/2])
    title(['Périodogramme du signal x',num2str(k), ' MUSIC'])
    xlabel('Fréquence lambda')
    ylabel(['x',num2str(k)])

end


%Pour sauvagader les figures comme .jpg--------------------------------

% % Encontra todos os handles de figuras abertas
% figHandles = findall(0, 'Type', 'figure');
% 
% % Loop para salvar cada figura como JPEG
% for i = 1:length(figHandles)
%     % Torna a figura atual ativa
%     figure(figHandles(i));
% 
%     % Define o nome do arquivo para cada figura
%     filename = ['figure_' num2str(i) '.jpg'];
% 
%     % Salva a figura como JPEG (usando o método 'print' para melhor qualidade)
%     print(figHandles(i), filename, '-djpeg', '-r300');  % '-r300' define a resolução (300 dpi)
% end