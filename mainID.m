%% Limpieza de variables y consola
clear
close all
clc;

%% Parametros
Semilla_PRBS = [1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1];
M = 100000;                     % Numero de muestras
gamma = 5000;                   % Ventana hanning
Ts = 0.005;                     % Ts de enunciado
Ts_exp = 1/0.8;                 % Ts calculado experimental

%% Vectores para la identificación
tfinal = Ts*(M-1);              % Tiempo final
t = (0:Ts:tfinal)';             % Vector de tiempo
t_exp = (0:Ts_exp:tfinal)';     % Vector de tiempo sub muestreado
N_t = length(t);                % Largo vector tiempo
N_t_exp = length(t_exp);        % Largo vector tiempo sub muestreado
gamma_e = gamma*Ts/Ts_exp;      % gamma normalizado para exp

%% Vectores de entrada
% Vector ruidoso partir de PRBS (0 1) y pasarlo a NRZ
input = prbs(M, Semilla_PRBS)';
input = (input - 0.5)*2;

% Vector submuestreado
input_exp = input(1:Ts_exp/Ts:end);

% Entrada a modelo
u1     = [t, input];
u1_exp = [t_exp, input_exp];

%% Simulacion
% Correr modelo con entrada y experimental
[t1, ~, output] = sim('BlackBox', tfinal, [], u1);
[t2, ~, output_exp] = sim('BlackBox', tfinal, [], u1_exp);

% Quitamos Drift of parte DC.
output = detrend(output(1:N_t));

% Creamos salida submuestreada.
output_exp = detrend(output_exp(1:N_t_exp));

%% Correlaciones
% Correlaciones sobremuestreado
R_u  = xcorr(input,  input,  'biased', gamma);
R_yu = xcorr(output, input,  'biased', gamma);
R_y  = xcorr(output, output, 'biased', gamma);
N_u  = length(R_u);

% Correlaciones submuestreado
R_u_e  = xcorr(input_exp,  input_exp,  'biased', gamma_e);
R_yu_e = xcorr(output_exp, input_exp,  'biased', gamma_e);
R_y_e  = xcorr(output_exp, output_exp, 'biased', gamma_e);
N_u_e  = length(R_u_e);

%% Funcion de Enventanado
hann   = hann_ventana(N_u);     %Ventana normalizada a sobremuestreo
hann_e = hann_ventana(N_u_e);   %Ventana normalizada a submuestreo

%% Estimacion de Espectro
% Formula: fi_u^N (w) = Sum W_y(tau)*Ru*e^(-j*w*tau)
% Comentario: NO USAR FFT DIRECTAMENTE
w   = 2*pi*(0:2*gamma)/sqrt(N_u);
w_e = 2*pi*(0:2*gamma_e)/sqrt(N_u_e);

% Estimacion suave de entrada
fi_u   = exp(-1i*(w'*w))*R_u;
fi_u_e = exp(-1i*(w_e'*w_e))*R_u_e;

% Estimacion suave de salida.
fi_y   = exp(-1i*(w'*w))*(R_y.*hann);
fi_y_e = exp(-1i*(w_e'*w_e))*(R_y_e.*hann_e);

% Estimacion suave cruzada
fi_yu   = exp(-1i*(w'*w))*(R_yu.*hann);
fi_yu_e = exp(-1i*(w_e'*w_e))*(R_yu_e.*hann_e);

% Formula: G_n = fi^N_yu(w)/fi^N_u(w)
% Estimacion suave de respuesta en frecuencia
G_N   = fftshift(fi_yu./fi_u);
G_N_e = fftshift(fi_yu_e./fi_u_e);

% Espectro de perturbacion (Disturbance Spectrum)
fi_v   = fi_y - (fi_yu'*fi_yu)./fi_u;
fi_v_e = fi_y_e - (fi_yu_e'*fi_yu_e)./fi_u_e;

% Espectro de coherencia (Coherence Spectrum), valor entre 0 y 1
kappa   = sqrt((fi_yu'*fi_yu)/(fi_y'*fi_u));
kappa_e = sqrt((fi_yu_e'*fi_yu_e)/(fi_y_e'*fi_u_e));

%% Graficar respuesta en frecuencia
% Vectores de frecuencia sobre muestreada
Fs = [1/Ts; 1/Ts_exp];                      % Frecuencia de muestreo
dF = Fs.*[1/gamma; 1/gamma_e];              % Espaciado en frecuencia
f   = (1:gamma)*dF(1);                      % vector de frecuencia
f_e = (1:gamma_e)*dF(2);                    % vector de frecuencia exp
G_w   = G_N(gamma+1:2*gamma);               % Vector G_w como respuesta
G_w_e = G_N_e(gamma_e+1:2*gamma_e);         % Vector G_w_e como respuesta
mag   = 10*log10(abs(G_w));                 % dB
mag_e = 10*log10(abs(G_w_e));               % dB exp

%% Bode de sistema sobremuestreado
figure
subplot(2,1,1)
semilogx(f, mag)
xlim([dF(1) 10^(1.4)])
grid on
xlabel('frecuencia [Hz]')
ylabel('Magnitud [dB]')
title('Respuesta en magnitud')
subplot(2,1,2)
semilogx(f, phase(G_w))
xlim([dF(1) 10^(1.4)])
grid on
xlabel('frecuencia [Hz]')
ylabel('Fase \phi [radianes]')
title('Respuesta en fase')
sgtitle('Bode del sistema') 

%% Bode de sistema submuestreado
figure
subplot(2,1,1)
semilogx(f_e, mag_e)
xlim([dF(2) 10^(1.4)])
grid on
xlabel('frecuencia [Hz]')
ylabel('Magnitud [dB]')
title('Respuesta en magnitud')
subplot(2,1,2)
semilogx(f_e, phase(G_w_e))
xlim([dF(2) 10^(1.4)])
grid on
xlabel('frecuencia [Hz]')
ylabel('Fase \phi [radianes]')
title('Respuesta en fase')
sgtitle('Bode del sistema con frecuencia confiable')

%% Estimacion de Ruido presente
figure
f2   = (-Fs(1):dF(1):Fs(1));
f2_e = (-Fs(2):dF(2):Fs(2));
subplot(1,2,1)
plot(f2, abs(fftshift(fi_v)))
grid on
xlabel('frecuencia [Hz]')
ylabel('|\phi^N_v(\omega)|')
title('Espectro de perturbación')
subplot(1,2,2)
plot(f2_e, abs(fftshift(fi_v_e)))
grid on
xlabel('frecuencia [Hz]')
ylabel('|\phi^N_v(\omega)|')
title('Espectro de perturbación submuestreada')

