clear
close all
clc;

%% Parametros para la identifiación:
Ts = 0.005;
conv_N = 1000;
Ts = Ts/conv_N;

%% Simulacion
pos = 1;
for i = 5:2:15
    vector = prbs15(100000, [1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1], i);
    vector = conv(vector, [1 zeros(1, conv_N)]);
    
    N = length(vector);
    % Especificaciones de Frecuencia
    Fs = 1/Ts;
    dF = Fs/length(vector);
    f = 2:dF:Fs/2;

    Fourier_vector = fftshift(fft(vector));

    mag = abs(Fourier_vector).^2/length(vector);
    mag = mag(N/2+2:N);
    
    subplot(2,3,pos)
    semilogx(f, 10*log10(mag))
    title('PRBS de orden '+string(i))
    xlabel('Frecuencia [Hz]')
    ylabel('dB')
    pos = pos + 1;
end
