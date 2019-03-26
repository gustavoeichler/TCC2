%% Este � o programa que ser� utilizado no meu Trabalho de Conclus�o de Curso
%% Mantenha sempre organizado!
%% Siga boas pr�ticas de programa��o sempre!
%% Muita for�a de vontade, os erros v�o vir, mas no fim tudo vai dar certo.
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Gerando o sinal de teste                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%   Carregar par�metros do sinal de teste   %%%%%%%%%%%%%%%%%%

    f_i = 5;        %frequ�ncia inicial do sweep
    
    f_f = 20000;    %frequ�ncia final do sweep
    
    T =  20.0226;    %tempo de dura��o do sweep em segundos
    
    fs = 44100;     %frequ�ncia de amostragem


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%         Fun��o que gera o sinal           %%%%%%%%%%%%%%%%%% 

    [Sinal, t, L,fs] = Gerar_Sinal(f_i,f_f, T,fs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%        Fun��o que aplica um fade in e fade out no sinal           %%%%

    [Sinal_com_fade] = Fade(Sinal,fs/5,fs/5,t,fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%      Carregar o Sweep inverso   %%%%%%%%%%%%%%%%%%%%%%%%

    load('SweepInverso.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%      Carregar Sinal da sa�da do amplificador
[ Saida, f_ax, Saida_FFT, len] = Carregar_saida('Condenser 1sweep (cortado).wav', fs);

%%      Realizar a decomposi��o
[ Impulse_response ] = deconvolucao(Saida,invSweep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N= 3; %n�mero de impulse responses

hm = Separacao_Kernel(Impulse_response, N, L,fs);
kernel_1 = hm(:,1);
kernel_2 = hm(:,2);
kernel_3 = hm(:,3);


Hm = fft(hm);
figure;
plot(linspace(-fs/2,fs/2,length(Hm)),abs(fftshift(Hm)));


Senoide1 = audioread('Guitarra Solo(2).wav');
%Senoide1 = audioread('audiocheck.net_sin_440Hz_-3dBFS_5s.wav');
Senoide2 = audioread('audiocheck.net_sin_1000Hz_-3dBFS_5s.wav');

Senoide = Senoide1;
Senoide_fft = fft(Senoide);
figure;
plot(linspace(-fs/2,fs/2,length(Senoide_fft)),abs(fftshift(Senoide_fft)));
title('sinal de teste');

Senoide_2 = Senoide.*Senoide;
Senoide_3 = Senoide.*Senoide_2;


Saida_1 = conv(Senoide,kernel_1);
Saida_2 = conv(Senoide_2,kernel_2);
Saida_3 = conv(Senoide_3,kernel_3);




Saida_geral = Saida_1 + Saida_2 + Saida_3;
Saida_geral_FFT = fft(Saida_geral);
figure;
plot(linspace(-fs/2,fs/2,length(Saida_geral_FFT)),abs(fftshift(Saida_geral_FFT)));

Saida_fft = fft(Saida_1);
figure;
plot(linspace(-fs/2,fs/2,length(Saida_fft)),abs(fftshift(Saida_fft)));
title('Saida 1 Frequencia');
Saida_2_fft = fft(Saida_2); Saida_2_fft(1) = 0;
figure;
plot(linspace(-fs/2,fs/2,length(Saida_2_fft)),abs(fftshift(Saida_2_fft)));
title('Saida 2 Frequencia');
Saida_3_fft = fft(Saida_3);
figure;
plot(linspace(-fs/2,fs/2,length(Saida_3_fft)),abs(fftshift(Saida_3_fft)));
title('Saida 3 Frequencia');



figure;
plot(linspace((-fs/2),(fs/2), length(Saida_FFT)), abs(fftshift(Saida_FFT)));
title ('Resposta em frequencia da Saida gravada');
hold on;
plot(linspace(-fs/2,fs/2,length(Saida_geral_FFT)),abs(fftshift(Saida_geral_FFT)));





%%%%%%%%

A(N,N) = 0;
for n = 1:N
    for m = 1:N
        if ( (n>=m) && (mod(n+m,2)==0) )
            % Eq. (48) of [Novak et al. (2010)]
            A(n,m) = (((-1)^(2*n+(1-m)/2))/(2^(n-1)))*nchoosek(n,(n-m)/2) ;
        end
    end
end

% tranform the HHFR Hm to filters Gm
Gm = Hm/A;
Gm_1 = Gm(:,1);
Gm_2 = Gm(:,2);
Gm_3 = Gm(:,3);


figure;
plot(linspace(-fs/2,fs/2,length(Gm)),(abs(fftshift((Gm)))));
% hold on;
% plot(linspace(-fs/2,fs/2,length(Gm_2)),(abs(Gm_2)));
% hold on;
% plot(linspace(-fs/2,fs/2,length(Gm_3)),(abs(Gm_3)));
% 
% Saida = conv(Senoide,Gm_1,'same');
% Saida_2 = conv(Senoide_2,Gm_2,'same');
% Saida_3 = conv(Senoide_3,Gm_3,'same');
% Saida_4 = conv(Senoide_4,Gm_4,'same');
% Saida_5 = conv(Senoide_5,Gm_5,'same');
% Saida_6 = conv(Senoide_6,Gm_6,'same');
% 


gm_1 = ifft(Gm_1);
gm_2 = ifft(Gm_2);
gm_3 = ifft(Gm_3);

output1 = conv(Saida_1,gm_1,'valid');
output2 = conv(Saida_3,gm_3,'valid');
Output = output1 + output2;
FFT_OUT = fft(Output);
figure;
plot(linspace(-fs/2,fs/2,length(FFT_OUT)),(abs(fftshift((FFT_OUT)))));