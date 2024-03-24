clear all; close all; clc

[x, Fs] = audioread('ses.wav');
Ts = 1 / Fs;
t = (0:length(x) - 1) * Ts;

figure
plot(t, x)
grid on
xlabel('Zaman [sn]')
ylabel('Genlik')
title(['x(t) Ses İşareti, Fs: ', num2str(Fs)])

sound(x, Fs)
% < Ses çalmayı durdurmak için komutu seçip F9'a basın >

%% Frekans analizi
X = fftshift(abs(fft(x)));
F = linspace(-Fs/2, Fs/2, length(X)); % F vektörünün uzunluğunu X vektörüyle aynı yapın

figure
plot(F, X)
grid on
xlim([-Fs/2, Fs/2])
xlabel('Frekans [Hz]')
ylabel('Genlik [V]')

%% Basit filtre tasarimi
% Frekans domaininde filtreleme yapmak = merkez frekansı etrafında −Fc<F<Fc
% olan kısmın alınması, diğer kısımların sıfır yapılmasıdır. Bunun için F
% üzerinden bir kare işareti üretilebilir.

Fc = 1000;
H = zeros(size(F));
for i = 1:length(F)
    if abs(F(i)) < Fc
        H(i) = 1;
    end
end

figure
plot(F, X / max(X), F, H, 'm', 'linewidth', 2)
grid on
xlabel('Frekans [Hz]')
ylabel('Genlik [V]')
xlim([-Fs/2, Fs/2])
legend('FFT of x', 'Filter Response')

%%
X = abs(fft(x));
Hs = fftshift(H);
figure
plot(X / max(X))
hold on
plot(Hs, 'm', 'linewidth', 2)
grid on
xlabel('Frekans [Hz]')
ylabel('Genlik [V]') 


y = ifft(Hs .* fft(x));

yr = real(y);
Yr = fftshift(abs(fft(yr)));
X = fftshift(abs(fft(x)));
F = linspace(-Fs/2 , Fs/2 , length(X));

figure
subplot(2,1,1)
plot(F, X)
grid on
xlim([-Fs/2 , Fs/2])
xlabel('Frekans [Hz]')
ylabel('Genlik [V]')

subplot(2,1,2)
plot(F, Yr, 'r')
grid on
xlim([-Fs/2 , Fs/2])
xlabel('Frekans [Hz]')
ylabel('Genlik [V]')

sound(yr, Fs)

Fc = 5000;
H = zeros(numel(F),1);
H( abs(F) <= Fc ) = 1; % pratik kod yazimi
% Önemli Not: MATLAB indis olarak 0 kabul etmez
Hs = fftshift(H);
y = ifft(Hs .* fft(x));
yr = real(y);
Yr = fftshift(abs(fft(yr)));
X = fftshift(abs(fft(x)));
F = linspace(-Fs/2 , Fs/2 , numel(X));
figure
subplot(2,1,1)
plot(F, X/max(X), F, H, 'g-', 'linewidth', 2)
grid on
xlim([-Fs/2 , Fs/2])
xlabel('Frekans [Hz]')
ylabel('Genlik [V]')
legend('FFT of x', 'Filter Response')

subplot(2,1,2)
plot(F, Yr, 'r')
grid on
xlim([-Fs/2 , Fs/2])
xlabel('Frekans [Hz]')
ylabel('Genlik [V]')

yr = yr / max(abs(yr)); % normalizasyon işlemi çıkış [-1 1] arasında olsun istiyoruz
sound(yr, Fs)

Fc = 3e3;
c = cos(2*pi*Fc*t)'; % taşıyıcı işaret
m = yr .* c; % DSB modulasyon
s = m .* c; % Demodulasyon için modüle işareti lokal osilatör ile çarptık

M = fftshift(abs(fft(m)));
S = fftshift(abs(fft(s)));
F = linspace(-Fs/2 , Fs/2 , numel(M));

figure
subplot(3,1,1)
plot(F, Yr, 'b')
grid on
xlim([-Fs/2 , Fs/2])
xlabel('Frekans [Hz]')
ylabel('Genlik [V]')
title('Filtrelenmiş Sinyalin FFT')

subplot(3,1,2)
plot(F, M, 'r')
grid on
xlim([-Fs/2 , Fs/2])
xlabel('Frekans [Hz]')
ylabel('Genlik [V]')
title('Modüle Sinyalin FFT')

subplot(3,1,3)
plot(F, S, 'g-')
grid on
xlim([-Fs/2 , Fs/2])
xlabel('Frekans [Hz]')
ylabel('Genlik [V]')
title('Demodule Sinyalin FFT')

% İşaretleri dinleyelim...
sound(m, Fs) % Kapatmak için :: clear sound
sound(s, Fs) % Kapatmak için :: clear sound


figure
plot(F, M/max(M), 'r', F, H, 'g-', 'linewidth', 2)
grid on
xlim([-Fs/2 , Fs/2])
xlabel('Frekans [Hz]')
ylabel('Genlik [V]')
title('Modüle Sinyalin FFT ve Filtre Tepkisi')

m_recover = ifft(Hs .* fft(m));
m_recover = real(m_recover);
m_recover = m_recover / max(abs(m_recover));
sound(m_recover, Fs) % clear sound

%%
figure
plot(F, S/max(S), 'r', F, H, 'g-', 'linewidth', 2)
grid on
xlim([-Fs/2 , Fs/2])
xlabel('Frekans [Hz]')
ylabel('Genlik [V]')
title('Demodule Sinyalin FFT ve Filtre Tepkisi')

s_recover = ifft(Hs .* fft(s));
s_recover = real(s_recover);
s_recover = s_recover / max(abs(s_recover));
sound(s_recover, Fs) % clear sound

