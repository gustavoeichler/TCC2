function [ Saida, f_ax, Saida_FFT, len] = Carregar_saida( audio,fs )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

Saida = audioread(audio);
%Saida = audio;

Saida_FFT = fft(Saida);
len = length(Saida_FFT);
f_ax = linspace(0,fs,len+1); 
f_ax(end) = [];

figure;
plot(linspace((-fs/2),(fs/2), length(Saida_FFT)), abs(fftshift(Saida_FFT)));
title ('Resposta em frequencia da Saida gravada');
grid on;


end

