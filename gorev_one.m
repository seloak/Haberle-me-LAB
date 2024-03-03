% S_1a%
A = [42 15 40 20; 72 9 54 88; 0 19 42 3; 30 35 69 67];
a = A(2, 4); % İndeksleri düzeltildi
b = A(2, :);
c = A(end, :);
disp(['a = ', num2str(a)]);
disp(['b = ', num2str(b)]);
disp(['c = ', num2str(c)]);

% S_1b%
A = [42 15 40 20; 72 9 54 88; 0 19 42 3; 30 35 69 67];
A_transpose = A'; % Transpozunu al
result = A .* A_transpose; % Element çarpımı
disp(result);
% S_1c%
A = [42 15 40 20; 72 9 54 88; 0 19 42 3; 30 35 69 67];
A_transpose = A'; % Transpozunu al
result = A * A_transpose; % Matris çarpımı
disp(result);
%S_1d%
A = [42 15 40 20; 72 9 54 88; 0 19 42 3; 30 35 69 67];
result = A(2, :) .* A(:, 3);
disp(result);
%S_1e%
A = [42 15 40 20; 72 9 54 88; 0 19 42 3; 30 35 69 67];
result = A(1, :) * A(:, 4);
disp(result);
%S_1f%
A = [42 15 40 20; 72 9 54 88; 0 19 42 3; 30 35 69 67];
result1 = A(2:3, 2:3) * A(2:3, end);
disp(result1);
%S_1g%
A = [42 15 40 20; 72 9 54 88; 0 19 42 3; 30 35 69 67];
result = A * A(:, 2);
disp(result);
%S_1h%
A = [42 15 40 20; 72 9 54 88; 0 19 42 3; 30 35 69 67];
result = A(2, 2:end) .* A(2:end, 2);
disp(result);
%S_2a%

Fs = 100; % Örnekleme frekansı
t = (0:1/Fs:3); % Zaman aralığı ve adımlar
x1 = 2*cos(2*pi*t); % x1(t) işaretini hesapla

figure; % Yeni bir figür penceresi aç
subplot(3,1,1);
plot(t, x1, 'b', 'LineWidth', 1.5); % İşareti çiz
grid on; % Izgara çizgilerini aç
xlabel('Zaman [sn]'); % x ekseni etiketi
ylabel('Çıkış [V]'); % y ekseni etiketi
title('x_1(t) İşareti'); % Başlık
axis([0 3 -2 2]); % Eksen sınırlarını ayarla

%---------------------------------------

x2 = 4*square(2*pi*t); % x2(t) işaretini hesapla
subplot(3,1,2);
plot(t, x2, 'r', 'LineWidth', 1.5); % İşareti çiz
grid on; % Izgara çizgilerini aç
xlabel('Zaman [sn]'); % x ekseni etiketi
ylabel('Çıkış [V]'); % y ekseni etiketi
title('x_2(t) İşareti'); % Başlık
axis([0 3 -5 5]); % Eksen sınırlarını ayarla
%----------------------------
x3 = x2 + x1
subplot(3,1,3);
plot(t, x3, 'm', 'LineWidth', 1.5); % İşareti çiz
grid on; % Izgara çizgilerini aç
xlabel('Zaman [sn]'); % x ekseni etiketi
ylabel('Çıkış [V]'); % y ekseni etiketi
title('x_3(t)=x_2(t) + x_1(t) İşareti'); % Başlık
axis([0 3 -10 10]); % Eksen sınırlarını ayarla

%S_3a%
Fs = 100; % Örnekleme frekansı
t = -5:1/Fs:5; % Zaman aralığı ve adımlar
x = 3.2 * cos(2*pi*0.25*t) - 2.1 * square(2*pi*2*t) + 5.3 * sin(2*pi*0.5*t + pi/17);

% Gürültü ekleyin
snr = 5; % SNR değeri (dB)
y = awgn(x, snr, 'measured');

% Gürültüsüz işareti çizdirin
subplot(3,1,1);
plot(t, x, 'b', 'LineWidth', .5);
grid on;
xlabel('Zaman [sn]');
ylabel('Genlik [V]');
title('Gürültüsüz İşaret');

% Gürültülü işareti çizdirin
subplot(3,1,2);
plot(t, y, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Zaman [sn]');
ylabel('Genlik [V]');
title('Gürültülü İşaret');

% SNR değerini hesaplayın
snr_4B = snr_degeri(x, y);
disp(['4B için SNR değeri: ', num2str(snr_4B), ' dB']);

% Otokorelasyonları hesaplayın
Rx = xcorr(x, 'biased');
Ry = xcorr(y, 'biased');

% Grafikleri çizdirin
subplot(3,1,3);
plot(Rx, 'g', 'LineWidth', .5);
plot(Ry, 'm', 'LineWidth', .5);
grid on;
xlabel('Lag');
ylabel('Otokorelasyon');
title('x(t) ve y(t) İşaretlerinin Otokorelasyonları');
legend('Rx', 'Ry');




%S_3E%
% Define the given signal x(t)
Fs = 100; % Sampling frequency
t = -5:1/Fs:5; % Time vector
x_t = 3.2 * cos(2*pi*0.25*t) - 2.1 * square(2*pi*2*t) + 5.3 * sin(2*pi*0.5*t + pi/17);

% Generate Gaussian noise with mean 2 and standard deviation 1
n = 2 + 1 * randn(size(t));

% Add noise to the signal
z_t = x_t + n;

% Plot the noisy signal and the original signal
figure;
subplot(2,1,1);
plot(t, x_t, 'b');
title('Original Signal x(t)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, z_t, 'r');
title('Noisy Signal z(t)');
xlabel('Time (s)');
ylabel('Amplitude');

% Calculate SNR (Signal-to-Noise Ratio)
SNR = 10 * log10(sum(x_t.^2) / sum(n.^2));

disp(['SNR: ', num2str(SNR), ' dB']);


function snr_value = snr_degeri(x, y)
    % SNR hesaplamak için gerekli formül
    signal_power = sum(x.^2) / length(x);
    noise_power = sum((x - y).^2) / length(x);
    snr_value = 10 * log10(signal_power / noise_power);
end