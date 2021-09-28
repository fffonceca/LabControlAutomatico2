%% Identificacion no paramerica con Blackman-Tukey
% G4: P. Ariztia, J. Contreras, C. Turrieta
clear 
clc

tfinal = 1.5;

% Tiempo de muestreo "continuo" para simulaciones iniciales "analogas"
Ts = 1e-4;
t_cont = (0:Ts:tfinal)';
n_cont = length(t_cont);



%% 1. Caracterizar respuesta en lazo abierto

% Respuesta con entrada 0 para obtener espectro del ruido
[t1, x_1, y_zero] = sim('BlackBox', tfinal, [], [t_cont, zeros(n_cont, 1)]);

% Escalones para ejemplificar respuesta en lazo abierto
[t2, x_2, y_step1] = sim('BlackBox', tfinal, [], [t_cont, 2*ones(n_cont, 1)]);
[t3, x_3, y_step2] = sim('BlackBox', tfinal, [], [t_cont, -2*ones(n_cont, 1)]);

% Graficos
subplot 221;plot(t1, y_zero);title(sprintf("Ruido: %0.3f mV (RMS)", 1e3*rms(y_zero)));
subplot 222;plot(t2, y_step1, t3, y_step2);title("Escalones")


%% 2. Determinar ancho de banda empiricamente
% Simulaciones con sinusoides para encontrar frecuencia de crossover y con
% ello el ancho de banda. Calculamos la amplitud tras cierto tiempo para
% evitar incluir la transiente inicial

%% 2.a Simular multiples sinusoides y medir sus amplitudes
f_sweep = logspace(0, 4, 40);
amp = zeros(1, length(f_sweep));
i = 0;
for f=f_sweep
    i = i+1;
    
    % Entrada
    u_sin = sin(2*pi*f*t_cont);
    
    % Simulacion
    [tf_1, x_4, yf] = sim('BlackBox', tfinal, [], [t_cont, u_sin]);
    
    % Calculo de amplitud 
    amp(i) = rms(detrend(yf(500:length(yf))))*sqrt(2);
end


%% 2.b Calcular ancho de banda segun la frecuencia de crossover

% Mostrar "Bode" obtenido
subplot 223;
figure; semilogx(f_sweep, 20*log10(amp))

% Consideramos el ancho de banda como la frecuencia de amplitud -3db
BW_i = interp1(20*log10(amp), f_sweep, -3);

% Para capturar altas frecuencias, muestreamos 4x el ancho de banda
Ts = 1/(BW_i*8); % 4x2 por Nyquist

line(-3, 'r');
line(BW_i, 'r');
xlabel("f");
ylabel("Amplitud (dB)");
title("Entradas sinusoidales")
grid on
text(BW_i+1, 10, sprintf("\\Leftarrow %0.3f \\rm{Hz}", BW_i));



%% 3. Simular sistema con entradas PRBS
% Usar idinput para generar PRBS con cierto periodo. Simular M periodos y 
% promediar las salidas para rechazar el ruido de mediciÃ³n

% n_samples debe ser suficientemente alto para que Ts*n_sample sea
% un tiempo "largo"
n_samples = 4095;%2047;
n_periods = 200;

tfinal = n_samples*n_periods*Ts;
t = linspace(0, tfinal, n_samples*n_periods)';

% Generar PRBS periodico y simular
PRBS = idinput([n_samples, 1, n_periods], 'prbs');
[tf, x, yf] =  sim('BlackBox', tfinal, [], [t, PRBS]);

% Truncar salidas ya que se añaden un par de muestras por alguna razon
tf = tf(1:n_periods*floor(length(tf)/n_periods));
yf = yf(1:n_periods*floor(length(yf)/n_periods));

% Quitar drift a yf
yf = detrend(yf);

% Reshape modulo el periodo y remover drift (suponiendo que sube mucho mas
% lento que un periodo)
tf_r = reshape(tf, [], n_periods);
yf_r = reshape(yf, [], n_periods);

% Promediar todos los periodos
y_final = mean(yf_r, 2);
y_std = std(yf_r, 0, 2);

N = length(y_final);


% Plotear salida y desviacion estandar
% https://www.mathworks.com/matlabcentral/answers/180829-shade-area-between-graphs#answer_169649
subplot 224;
tp = tf_r(:,1);
plot(tp, y_final); hold on;
c1 = y_final - y_std; c2 = y_final + y_std;
x2 = [tp', flipud(tp)']; inbetween = [c1', flipud(c2)'];
fill(x2, inbetween,'g');
title(sprintf("Respuesta promedio a PRBS (%i periodos)",n_periods));




%% 4. Calcular espectros suavizados
%% 4.a Variables auxiliares para los calculos
gamma = floor(N/2);
tau = -gamma:1:gamma;
Np = length(tau);
om = linspace(1,1/Ts,n_samples);
u = interp1(1:n_samples, PRBS(1:n_samples), 1:Np);
y = interp1(1:n_samples, y_final(1:n_samples), 1:Np);


% Señales shifteadas para calculos de correlación
Y = repmat(y, [Np, 1]); %dim1:tau, dim2:t
U = repmat(u, [Np, 1]);
Y_SHIFT = Y;
U_SHIFT = U;
for i=1:Np
    Y_SHIFT(i,:) = circshift(Y(i,:),i);
    U_SHIFT(i,:) = circshift(U(i,:),i);
end


%% 4.b Calcular correlaciones

% CALCULAR CORRELACIONES
Rsw_gen = @(w, s)(fftshift(sum(w.*s, 2)./Np));
Ryy = Rsw_gen(Y, Y_SHIFT);
Ruu = Rsw_gen(U, U_SHIFT);
Ryu = Rsw_gen(Y, U_SHIFT);

% foo_plot = @(f)(plot(tau, f, 'LineWidth', 2));
% subplot 311; foo_plot(Ryy); title("Ryy");
% subplot 312; foo_plot(Ruu); title("Ruu");
% subplot 313; foo_plot(Ryu); title("Ryu");


%% 4.c Ventanas
% En una version previa del codigo se hacia un for con todas las ventanas
% de a continuacion. Se decidio finalmente por la de hamming.

gen_win = @(f)(f(length(Y)));
%W1 = gen_win(@rectwin);
%W2 = gen_win(@gausswin);
W3 = gen_win(@hamming); % Decidí quedarnos con la hamming pq si nomas :v
%W4 = gen_win(@hann);

W = W3;

%% 4.d CALCULAR ESPECTROS
Phi_gen = @(R)((R'*exp(-1i*tau'*2*pi*(1:length(R))./length(R)))');
Phiyy = Phi_gen(W.*Ryy);
Phiuu = Phi_gen(W.*Ruu);
Phiyu = Phi_gen(W.*Ryu);

% Para validar:
% Phiyy = fft(Ryy.*W);
% Phiuu = fft(Ruu.*W);
% Phiyu = fft(Ryu.*W);
% foo_plot = @(f)(semilogx(om, abs(f), 'LineWidth', 2));
% figure
% subplot 311; foo_plot(Phiyy); hold on; title("Phiyy");
% subplot 312; foo_plot(Phiuu); title("Phiuu");
% subplot 313; foo_plot(Phiyu); title("Phiyu");


%% 4.e CALCULAR BODE SUAVISADO Y MEDIDAS DE DESEMPEÑO
GN = Phiyu./Phiuu;
Phiv = Phiyy - (abs(Phiyu).^2)./Phiuu;
kappa = sqrt(power(abs(Phiyu),2)./(Phiyy.*Phiuu));



%% 4.f PLOTEAR
figure
subplot 321; semilogx(om, 20*log10(abs(GN)));title("|GN|"); %ylim([1 1000]); 
subplot 323; semilogx(om, 20*log10(abs(Phiv))); title("|Phiv|"); %ylim([1 1e4]);
subplot 325; semilogx(om, abs(kappa));title("|K|"); ylim([0 1.1]);
subplot 322; semilogx(om, angle(GN));title("arg(GN)"); %ylim([1 1000]); 
subplot 324; semilogx(om, angle(Phiv)); title("arg(Phiv)"); %ylim([1 1e4]);
subplot 326; semilogx(om, angle(kappa));title("arg(K)"); %ylim([1e-2 1e-1]); 