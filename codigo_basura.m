%% Basura
%% Correlacion
y1 = y1(1:1:length(u1));
[Ruu_e, lagsU] = xcorr(u1(:,2), u1(:,2), 'biased', gamma);
[Ryu_e, lagsYU] = xcorr(u1(:,2), y1(:,1), 'biased', gamma);
[Ryy_e, lagsYY] = xcorr(y1(:,1), y1(:,1), 'biased', gamma);

%% DFT y G(e^jw)
fi_uu = fft(hann(2*gamma+1).*Ruu_e);
fi_yu = fft(hann(2*gamma+1).*Ryu_e);
fi_yy = fft(hann(2*gamma+1).*Ryy_e);

%% Especificaciones de Frecuencia
Fs = 1/Ts;
dF = Fs/length(Ruu_e);        % hertz
f = -Fs/2:dF:Fs/2-dF;           % hertz

%% Planta
G_n = fi_yu./fi_uu;
% Espectro de perturbacion
P_n = fi_yy - (fi_yu'*fi_yu)./fi_uu;
% Espectro de coherencia
K = sqrt((fi_yu'*fi_yu)/(fi_yy'*fi_uu));

%% Ploteo
figure 
grid on
X1 = G_n;
X2 = fftshift(P_n);
mag_1 = abs(X1);
mag2db(mag_1)
phase_1 = unwrap(angle(X1));
plot(1:length(mag_1),mag_1)
% plot(1:length(G_)

