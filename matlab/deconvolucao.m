function [ Impulse_response ] = deconvolucao( sinal_de_saida, sine_inverse)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
sinal_de_entrada = vec2mat(sine_inverse,1);
%sinal_de_entrada = sine_inverse;
IR = conv(sine_inverse, sinal_de_saida,'full');

IR_FFT = fft(IR);
IR_FFT(1) = 0;
Impulse_response=ifft(IR_FFT);

end

