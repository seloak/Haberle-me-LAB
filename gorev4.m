%%S_1a
clc; 
close all; 
clear all;

% Genliği 5 ve frekansı 10 olan cosinus giriş işaretini oluşturma
t = 0:0.01:1; % 0 ile 1 arasında 0.01 adım aralığında
A = 5; % Genlik
f = 10; % Frekans
cos_input = A * cos(2*pi*f*t);

% Giriş işaretinin zaman alanındaki grafiği
subplot(2,1,1);
plot(t, cos_input);
title('Zaman Alanında Giriş İşareti');
xlabel('Zaman');
ylabel('Genlik');

% Giriş işaretinin eksen düzenli genlik spektrumunu çizme
N = length(cos_input); % İşaretin uzunluğu
f_axis = (-N/2:N/2-1)*(1/N); % Frekans eksenini oluşturma
cos_input_fft = fftshift(fft(cos_input)); % FFT'yi alıp kaydırmak
mag_spectrum = abs(cos_input_fft); % Genlik spektrumunu hesaplama

% Genlik spektrumunun grafiği
subplot(2,1,2);
plot(f_axis, mag_spectrum);
title('Eksen Düzenli Genlik Spektrumu');
xlabel('Frekans');
ylabel('Genlik');

%%S_1b
clc; 
close all; 
clear all;

% 1A'da istenilen sinyalin oluşturulması
t = 0:0.01:1;
s = 5 * cos(10 * 2 * pi * t);

% Örneklenmiş sinyallerin oluşturulması
n1_values = [5, 10, 50];
sampled_signals = cell(1, length(n1_values));

for i = 1:length(n1_values)
    n1 = n1_values(i);
    sampled_signals{i} = 5 * cos(10 * 2 * pi * (0:1/(n1-1):1)); % Örneklenmiş sinyal oluşturulurken örneklenen nokta sayısına göre adım düzeltilmiştir
end

% Grafiklerin çizdirilmesi
figure;

% Analog sinyal grafiği
subplot(length(n1_values) + 1, 1, 1);
plot(t, s);
title('Analog Sinyal');
xlabel('Zaman');
ylabel('Genlik');

% Örneklenmiş sinyal grafiği
for i = 1:length(n1_values)
    subplot(length(n1_values) + 1, 1, i + 1);
    plot(0:1/(n1_values(i)-1):1, sampled_signals{i}, '-o'); % Örneklenmiş sinyal oluşturulurken örneklenen nokta sayısına göre adım düzeltilmiştir
    title(['Örneklenmiş Sinyal (n1 = ', num2str(n1_values(i)), ')']);
    xlabel('Zaman');
    ylabel('Genlik');
end

%%2_a


% Parametrelerin tanımlanması
A = 10; % Genlik
n = 20; % Bit sayısı
n1 = 10; % Örnek sayısı
L = 2^n; % Kuantalama seviyesi

% Sinyal oluşturma ve örnekleme
x = 0:2*pi/n1:4*pi; % n1 sayıda örnek alınır
s = A * sin(x); % Sinyalin genliği A'dır

% Kuantalama işlemi
Amax = A;
Amin = -A;
del = (Amax - Amin) / L; % Adım boyu
part = Amin:del:Amax; % Seviyeler
code = Amin - (del/2):del:Amax + (del/2); % Kuantalama değerleri
[ind, q] = quantiz(s, part, code); % Kuantalama süreci

% Kuantalı sinyalin çizdirilmesi
figure;
stem(q);
title('Kuantalı İşaret');
xlabel('Zaman');
ylabel('Genlik');

%%2_b
% Kodlama işlemi
code = de2bi(ind, n, 'left-msb'); % Decimali ikili sistemde kodlama
coded_signal = reshape(code', 1, []); % Kodlanmış sinyal vektörü oluşturma

% Kodlanmış sinyalin çizdirilmesi
figure;
stairs(coded_signal);
title('Kodlanmış İşaret');
xlabel('Zaman');
ylabel('Genlik');

%%2_c
% Demodülasyon işlemi
index = bi2de(code, 'left-msb'); % İkili sistemdeki kodları ondalık değerlere dönüştürme
demodulated_signal = del * index + Amin + (del / 2); % Kuantalama değerlerine geri alma

% Orijinal işaret ve demodule edilmiş işaretin çizdirilmesi
figure;
t = 0:0.01:12; % 0 ile 1 arasında 0.01 adım aralığında
A = 10; % Genlik
f = 0.16; % Frekans
cos_input = A * sin(2*pi*f*t);
plot(t, cos_input, 'g', 'LineWidth', 1.5); % Orijinal işaret yeşil renkte çizdirilir
hold on;
plot(x, demodulated_signal, 'r', 'LineWidth', 1.5); % Demodule edilmiş işaret kırmızı renkte çizdirilir
title('Demodule Edilmiş ve Orijinal İşaretler');
xlabel('Zaman');
ylabel('Genlik');
legend('Orijinal İşaret', 'Demodule Edilmiş İşaret');
grid on;
hold off;

%%3_a_b_c
clc;
close all;
clear all;

% Parametrelerin tanımlanması
A = 10; % İşaretin genliği
f = 0.2; % İşaretin frekansı
t = 0:0.01:12; % İşaretin zaman aralığı

% PCM giriş işareti
pcm_signal = A * sin(2*pi*f*t);

% Beyaz gürültü oluşturma
noise_power = 10^(0/10); % Gürültü gücü (0 dB)
noise = sqrt(noise_power) * randn(size(t)); % Beyaz gauss gürültüsü

% Giriş işareti ile gürültüyü toplama
pcm_with_noise_0dB = pcm_signal + noise;

% 20 dB gürültü gücü
noise_power_20dB = 10^(20/10); % Gürültü gücü (20 dB)
noise_20dB = sqrt(noise_power_20dB) * randn(size(t)); % Beyaz gauss gürültüsü

% Giriş işareti ile 20 dB gürültüyü toplama
pcm_with_noise_20dB = pcm_signal + noise_20dB;

% İki işareti aynı grafikte çizdirme
figure;
subplot(2,1,1);
plot(t, pcm_with_noise_0dB, 'g', 'LineWidth', 1.5); % 0 dB gürültülü PCM işareti
title('0 dB Gürültülü PCM İşareti');
xlabel('Zaman');
ylabel('Genlik');
grid on;

subplot(2,1,2);
plot(t, pcm_with_noise_20dB, 'r', 'LineWidth', 1.5); % 20 dB gürültülü PCM işareti
title('20 dB Gürültülü PCM İşareti');
xlabel('Zaman');
ylabel('Genlik');
grid on;

%%3_a
%% PCM Giriş İşareti Oluşturma ve Gürültü Eklenmesi
clc; 
close all; 
clear all;

% Parametreler
A = 5; % Genlik
f = 1000; % Frekans (Hz)
fs = 10*f; % Örnekleme frekansı
t = 0:1/fs:0.01; % Zaman aralığı

% PCM Giriş İşareti Oluşturma
pcm_input_signal = A * sin(2*pi*f*t);

% Beyaz Gauss Gürültüsü Oluşturma
noise_power = 10^(20/10); % Gürültü gücü (20 dB)
noise = sqrt(noise_power) * randn(size(pcm_input_signal)); % Beyaz Gauss gürültüsü

% Gürültüyü PCM giriş işaretine ekleyerek Gürültülü PCM İşareti Oluşturma
pcm_noisy_signal = pcm_input_signal + noise;

% PCM Giriş İşareti ve Gürültülü PCM İşaretinin Çizdirilmesi
figure;
plot(t, pcm_input_signal, 'b', 'LineWidth', 1.5); % PCM giriş işareti (mavi renkte)
hold on;
plot(t, pcm_noisy_signal, 'r', 'LineWidth', 1.5); % Gürültülü PCM işareti (kırmızı renkte)
title('PCM Giriş İşareti ve Gürültülü PCM İşareti');
xlabel('Zaman (s)');
ylabel('Genlik');
legend('PCM Giriş İşareti', 'Gürültülü PCM İşareti');
grid on;

%%3_b
%% Gürültülü PCM İşaretinin Demodülasyonu ve Karşılaştırılması
% Demodülasyon işlemi
index_pcm_with_noise = bi2de(code, 'left-msb'); % İkili sistemdeki kodları ondalık değerlere dönüştürme
demodulated_pcm_with_noise = del * index_pcm_with_noise + Amin + (del / 2); % Kuantalama değerlerine geri alma

% Giriş işareti ve demodüle edilmiş PCM işaretinin çizdirilmesi
figure;
plot(t, pcm_input_signal, 'b', 'LineWidth', 1.5); % PCM giriş işareti (mavi renkte)
hold on;
plot(t, demodulated_pcm_with_noise, 'r', 'LineWidth', 1.5); % Demodüle edilmiş PCM işareti (kırmızı renkte)
title('PCM Giriş İşareti ve Demodüle Edilmiş PCM İşareti');
xlabel('Zaman (s)');
ylabel('Genlik');
legend('PCM Giriş İşareti', 'Demodüle Edilmiş PCM İşareti');
grid on;




