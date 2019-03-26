function [ Sinal_com_fade ] = Fade( x,n_in, n_out, t, fs )
%UNTITLED4 Summary of this function goes here
%   x = Sinal de entrada
%   n_in = número de samples que vão ser atenuados no início
%   n_out = número de samples que vão ser atenuados no final

x = x;
fade_in = (1-cos((0:n_in-1)/n_in*pi))/2;
index = 1:n_in;
x(index) = x(index).*fade_in;

% fade-out the input signal

fade_out = (1-cos((0:n_out-1)/n_out*pi))/2;
index = (1:n_out)-1;
x(end-index) = x(end-index).*fade_out;
Sinal_com_fade = x;


FFT_Sinal_com_fade = fft(Sinal_com_fade);
figure;
subplot(2,1,1);
plot(t,Sinal_com_fade);
title('Sinal exponencial com fade in e fade out no tempo');
grid on;
subplot(2,1,2);
plot(linspace((-fs/2),(fs/2),length(FFT_Sinal_com_fade)), abs(fftshift(FFT_Sinal_com_fade)));
title('Sinal exponencial com fade in e fade out na frequencia');
grid on;
end

