clc; 
close all; 
clear all;

n=input('N−bit PCM sistemi icin n degerini girin: ');
n1=input('Bir periyottaki ornek sayisini girin: ');

L=2^n; % Kuantalama seviyesi

% Sinyal Uretimi ve Ornekleme
x=0:2*pi/n1:4*pi; % n1 sayida ornek secilmelidir
s=8*sin(x); % Sinyalin Genligi 8v'dir

subplot(5,1,1);
plot(x,s);
title('Analog Sinyal'); 
ylabel('Genlik−−−>'); 
xlabel('Zaman−−−>');

subplot(5,1,2);
stem(x,s); 
grid on;
title('Orneklenmis Sinyal'); 
ylabel('Genlik−−−>'); 
xlabel('Zaman−−−>');

% Kuantalama Islemi
Amax=8; 
Amin=-Amax;
del=(Amax-Amin)/L; % adim boyu
part=Amin:del:Amax; % Seviye, del farki ile Amin ve Amax arasindadir
code=Amin-(del/2):del:Amax+(del/2); % Kuantalama Degerleri
[ind,q]=quantiz(s,part,code); % Kuantalama Sureci
% ind indeks numarasini, q kuantalama degerlerini icerir
l1=length(ind); l2=length(q);
for i=1:l1
    if(ind(i)~=0) % Dizini ikili sistemde yapmak icin 0'dan N'ye baslar
        ind(i)=ind(i)-1;
    end
    i=i+1;
end
for i=1:l2
    if(q(i)==Amin-(del/2)) % Seviyeler arasi kuantalama yapilir
        q(i)=Amin+(del/2);
    end
end
subplot(5,1,3);
stem(q);grid on;
title('Kuantalanmis sinyal'); 
ylabel('Genlik−−−>'); 
xlabel('Zaman−−−>');

% Kodlama Islemi
code=de2bi(ind,'left-msb'); % Decimali ikili sisteme donusturme
k=1;
for i=1:l1
    for j=1:n
        coded(k)=code(i,j); % kod matrisini kodlanmis bir satir vektorune donusturme
        j=j+1;
        k=k+1;
    end
    i=i+1;
end
subplot(5,1,4); 
grid on;
stairs(coded); % Kodlanmis Sinyal
axis([0 100 -2 3]);
title('Kodlanmis Sinyal'); 
ylabel('Genlik−−−>'); 
xlabel('Zaman−−−>');


% Demodülasyon Islemi
qunt=reshape(coded,n,length(coded)/n);
index=bi2de(qunt','left-msb'); % Dizini ondalik degerlerine geri alma
q=del*index+Amin+(del/2); % Kuantalama degerlerine geri alma
subplot(5,1,5); 
grid on;
plot(q);
title('Demoduleli Sinyal'); 
ylabel('Genlik−−−>'); 
xlabel('Zaman−−−>');



