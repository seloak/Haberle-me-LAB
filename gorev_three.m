%%S_1aCift Yanband Isareti (DSB) ve Yorumlanmasi
% Fs : Isareti MATLAB ortaminda analog gibi islemek icin
% kullanilan ornekleme frekansi (Fc'nin 100 kati)
% Fc : Tasiyicinin frekansi
close all,clear all,clc
Fs = 5000; Ts = 1/Fs;
Fc = 50;
A = 1;
t = -1:Ts:1;
input = A*exp(-5*t.^2); % bilgi isareti
carrier = A*cos(2*pi*Fc*t); % tasiyici
output = (input.*carrier);
figure,
subplot(311),plot(t,input),title('Bilgi İşareti'),grid on
xlabel('Zaman [sn]'),ylabel('Genlik [V]')
subplot(312),plot(t,carrier),title('Taşıyıcı İşareti'),grid on
xlabel('Zaman [sn]'),ylabel('Genlik [V]')
subplot(313),plot(t,output),title('Modüleli İşaret'),grid on
xlabel('Zaman [sn]'),ylabel('Genlik [V]')

%% 1b) Cift Yanband Isaretinin (DSB) Frekans Analizi
% Giriş işareti
input_fft = fftshift(abs(fft(input)));
f_input = linspace(-Fs/2, Fs/2, numel(input_fft));

% Taşıyıcı işareti
carrier_fft = fftshift(abs(fft(carrier)));
f_carrier = linspace(-Fs/2, Fs/2, numel(carrier_fft));

% Çıkış işareti
output_fft = fftshift(abs(fft(output)));
f_output = linspace(-Fs/2, Fs/2, numel(output_fft));

figure
subplot(3,1,1)
plot(f_input, input_fft)
title('Giriş İşaretinin Genlik Spektrumu')
xlabel('Frekans [Hz]')
ylabel('Genlik')

subplot(3,1,2)
plot(f_carrier, carrier_fft)
title('Taşıyıcı İşaretin Genlik Spektrumu')
xlabel('Frekans [Hz]')
ylabel('Genlik')

subplot(3,1,3)
plot(f_output, output_fft)
title('Modüleli İşaretin Genlik Spektrumu')
xlabel('Frekans [Hz]')
ylabel('Genlik')

%% 2a) Tasiyicili Genlik Modulasyonu (AM) ve Yorumlanmasi
close all, clear all, clc
Fs = 5000; Ts = 1/Fs;
Fc = 100; % taşıyıcı frekans - Hz
Fm = 10;
t = -1:Ts:1;
A = 1; ma = 1;
input = A*cos(2*pi*Fm*t); % bilgi işareti
carrier = cos(2*pi*Fc*t); % taşıyıcı
s_am = (1 + ma*cos(2*pi*Fm*t)).*carrier; % AM işareti

figure,
subplot(311), plot(t, input), title('Bilgi İşareti'), grid on
xlabel('Zaman [sn]'), ylabel('Genlik [V]')
subplot(312), plot(t, carrier), title('Taşıyıcı İşareti'), grid on
xlabel('Zaman [sn]'), ylabel('Genlik [V]')
subplot(313), plot(t, s_am), title('AM İşareti'), grid on
xlabel('Zaman [sn]'), ylabel('Genlik [V]'), xlim([-0.5 0.5])

%% 2b) Tasiyicili Genlik Modulasyonu (AM) ve Yorumlanmasi
close all, clear all, clc
Fs = 5000; Ts = 1/Fs;
Fc = 100; % taşıyıcı frekans - Hz
Fm = 10;
t = -1:Ts:1;
A = 1;
mas = [0.1, 0.5, 3];

figure;
for i = 1:length(mas)
    ma = mas(i);
    input = A*cos(2*pi*Fm*t); % bilgi işareti
    carrier = cos(2*pi*Fc*t); % taşıyıcı
    s_am = (1 + ma*cos(2*pi*Fm*t)).*carrier; % AM işareti
    
    subplot(length(mas), 1, i);
    plot(t, s_am);
    title(['AM İşareti, m_a = ', num2str(ma)]);
    xlabel('Zaman [sn]');
    ylabel('Genlik [V]');
    grid on;
    xlim([-0.5 0.5]);
end


%% 2C) Fourier Analizi
close all, clear all, clc
Fs = 5000; Ts = 1/Fs;
Fc = 100; % taşıyıcı frekans - Hz
Fm = 10;
t = -1:Ts:1;
A = 1;
ma = 1;

% Giriş işareti
input = A*cos(2*pi*Fm*t);

% Taşıyıcı işareti
carrier = cos(2*pi*Fc*t);

% Çıkış işareti (AM işareti)
s_am = (1 + ma*cos(2*pi*Fm*t)).*carrier;

% Genlik spektrumlarının hesaplanması
N = length(t); % Veri noktalarının sayısı
X_input = fftshift(abs(fft(input)))/N;
X_carrier = fftshift(abs(fft(carrier)))/N;
X_s_am = fftshift(abs(fft(s_am)))/N;

% Frekans vektörünün oluşturulması
f = linspace(-Fs/2, Fs/2, N);

% Grafiklerin çizdirilmesi
figure;

% Giriş işareti genlik spektrumu
subplot(3, 1, 1);
plot(f, X_input);
title('Giriş İşareti Genlik Spektrumu');
xlabel('Frekans [Hz]');
ylabel('Genlik');
xlim([-3*Fc, 3*Fc]);
grid on;

% Taşıyıcı işareti genlik spektrumu
subplot(3, 1, 2);
plot(f, X_carrier);
title('Taşıyıcı İşareti Genlik Spektrumu');
xlabel('Frekans [Hz]');
ylabel('Genlik');
xlim([-3*Fc, 3*Fc]);
grid on;

% Çıkış işareti (AM işareti) genlik spektrumu
subplot(3, 1, 3);
plot(f, X_s_am);
title('Çıkış İşareti (AM İşareti) Genlik Spektrumu');
xlabel('Frekans [Hz]');
ylabel('Genlik');
xlim([-3*Fc, 3*Fc]);
grid on;

%% 3a) Hilbert dönüşümü
close all, clear all, clc
Fs = 5000; Ts = 1/Fs;
Fc = 50; Fm = 10;
Am = 1; Ac = 5; 
t = -1:Ts:1;

% Mesaj işareti
m = Am*cos(2*pi*Fm*t);

% Taşıyıcı işareti
c = Ac*cos(2*pi*Fc*t);


mhat = hilbert(m); % Mesaj işaretinin Hilbert dönüşümü
usb = c.*exp(1j.*real(mhat)); % USB (Üst Yan Bant)
lsb = c.*exp(1j.*imag(mhat)); % LSB (Alt Yan Bant)

% USB ve LSB işaretlerinin genlik spektrumları
N = length(t);
USB_spectrum = fftshift(abs(fft(usb)))/N;
LSB_spectrum = fftshift(abs(fft(lsb)))/N;
f = linspace(-Fs/2, Fs/2, N); % Frekans vektörü oluşturulması

% Grafiklerin çizdirilmesi
figure;

% USB işaretinin genlik spektrumu
subplot(2, 1, 1);
plot(f, USB_spectrum);
title('USB İşaretinin Genlik Spektrumu');
xlabel('Frekans [Hz]');
ylabel('Genlik');
grid on;

% LSB işaretinin genlik spektrumu
subplot(2, 1, 2);
plot(f, LSB_spectrum);
title('LSB İşaretinin Genlik Spektrumu');
xlabel('Frekans [Hz]');
ylabel('Genlik');
grid on;


% 3B) DSB işaretin bir bandının bastırılması ile SSB modulasyonu
% 5.bölümde verilen frekans bölgesinde filtreleme benzeri bir filtreleme
% işlemi gerçekleştirmeniz beklenmektedir.

%% 3b) DSB isaretin bandlarindan birisinin bastirilmasi ile SSB modulasyonunu gerçekleştirme
% Verilen frekans bölgesinde filtreleme benzeri bir işlem gerçekleştirilir.

% Basit filtre tasarimi
close all, clear all, clc
Fc = 1000; % Kesim frekansı

H = zeros(numel(m),1); % Filtre tepkisi
for i = 1:numel(m)
    if abs(m(i)) < Fc
        H(i) = 1; % Kesim frekansı altında birim genlik
    end
end

% SSB modulasyonu
s_ssb_usb = fftshift(H) .* fft(usb); % USB (Upper Sideband)
s_ssb_lsb = fftshift(H) .* fft(lsb); % LSB (Lower Sideband)

% Genlik spektrumlarını çizdirme
figure;
subplot(211), plot(m, abs(s_ssb_usb)), title('USB Genlik Spektrumu'), xlabel('Frekans [Hz]'), ylabel('Genlik'), grid on;
subplot(212), plot(m, abs(s_ssb_lsb)), title('LSB Genlik Spektrumu'), xlabel('Frekans [Hz]'), ylabel('Genlik'), grid on;
