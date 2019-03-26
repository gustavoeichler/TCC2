function [ Sinal,t, L, fs ] = Gerar_Sinal( f_i, f_f, T, fs )
%Função para gerar os sinais de teste
%   f_i = frequencia inicial do sinal
%   f_f = frequencia final do sinal
%   T = tempo do sinal em segundos
%   fs = frequencia de amostragem do sinal

t = (0:1/fs:T-1/fs);
L = T/log(f_f/f_i);

Sinal = sin(2*pi*f_i*L*(exp(t/L)));

FFT_Sinal = fft(Sinal);
figure;
subplot(2,1,1);
plot(t,Sinal);
title('Sinal exponencial no tempo');
grid on;
subplot(2,1,2);
plot(linspace((-fs/2),(fs/2),length(FFT_Sinal)), abs(fftshift(FFT_Sinal)));
title('Sinal exponencial na frequencia');
grid on;

end

