x = 0:pi/1000:20*pi;
y = sin(x);
subplot(2,1,1);
plot(x, y);
xlabel('x = 0:2\pi');
ylabel('Genlik');
title('Sinüs İşareti, Bir Periyot Boyunca');

x = 0:pi/100:2*pi;
y1 = 2*cos(x); y2 = cos(x);
y3 = 0.5*cos(x);

subplot(2,1,2);
plot(x, y1, 'k--', ...
     x, y2, 'r--', ...
     x, y3, 'b--');
xlabel('0 \leq x \leq 2\pi')
ylabel('Kosinüs Fonksiyonları')
legend('2*cos(x)', 'cos(x)', '0.5*cos(x)')