clear all;
close all;

%% Part 3
refinedFreqPS = 2*pi*xlsread('2_3_1 data.xlsx',1,'A10:A10009');
refinedPS = xlsread('2_3_1 data.xlsx',1,'B10:B10009');
refinedFreqMag = 2*pi*xlsread('2_3_1 data.xlsx',2,'A10:A10009');
refinedMag = xlsread('2_3_1 data.xlsx',2,'B10:B10009');
FreqPS = 2*pi*xlsread('2_3_1 data.xlsx',3,'A10:A10009');
PS = xlsread('2_3_1 data.xlsx',3,'B10:B10009');
FreqMag = 2*pi*xlsread('2_3_1 data.xlsx',4,'A10:A10009');
Mag = xlsread('2_3_1 data.xlsx',4,'B10:B10009');
%%

% 1st break frequency
x1data=FreqMag(80:400);
y1data=Mag(80:400);
x2data=FreqMag(5:10);
y2data=Mag(5:10);
p1=polyfit(log10(x1data),y1data,1);
wlog1=logspace(2,4,500);
mag1=p1(1)*log10(wlog1)+p1(2);
p2=polyfit(log10(x2data),y2data,1);
wlog2=logspace(2,4,500);
mag2=p2(1)*log10(wlog2)+p2(2);

%2nd break frequency
x3data=refinedFreqMag(30:60);
y3data=refinedMag(30:60);
x4data=refinedFreqMag(170:250);
y4data=refinedMag(170:250);
p3=polyfit(log10(x3data),y3data,1);
wlog3=logspace(4.5,5.6,500);
mag3=p3(1)*log10(wlog3)+p3(2);
p4=polyfit(log10(x4data),y4data,1);
wlog4=logspace(4.5,6,500);
mag4=p4(1)*log10(wlog4)+p4(2);

figure(1)
subplot(2,1,1);
semilogx(FreqMag,Mag);
hold on;
grid on;
semilogx(refinedFreqMag,refinedMag,'k');

xlim([200 500000]);
%ylim([-40 25]);
ylabel('Amplitude Ratio (dB)','FontSize',12);
title('Experimental Bode Response Plot','FontSize',14);
legend('Location','best','Large Time Division','Small Time Division');

subplot(2,1,2);
semilogx(FreqPS,PS);
hold on;
grid on;
semilogx(refinedFreqPS,refinedPS,'k');
xlim([200 500000]);
%ylim([-180 0]);
xlabel('Frequency (rad/s)','FontSize',12);
ylabel('Phase Shift (deg)','FontSize',12);
legend('Location','best','Large Time Division','Small Time Division');
%%
figure(2)

subplot(2,1,1);
semilogx(FreqMag,Mag);
hold on;
grid on;
semilogx(wlog1,mag1,'r')
semilogx(wlog2,mag2,'b')
Intersection1=find(abs(mag2-mag1)<=(0.02));
break1=wlog1(Intersection1);
%plot(break1,mag1(Intersection1),'o')
xlim([500 10000]);
%ylim([-40 25]);
ylabel('Amplitude Ratio (dB)','FontSize',12);
xlabel('Frequency (rad/s)','FontSize',12);
title('1st Break Frequency','FontSize',14);
%legend('Location','best','Simulated Bode Plot','Experimental Bode Plot');

figure(3)

subplot(2,1,1);
semilogx(refinedFreqMag,refinedMag,'k');
hold on;
grid on;
semilogx(wlog3,mag3,'r')
semilogx(wlog4,mag4,'b')
Intersection2=find(abs(mag4-mag3)<=(0.03));
break2=wlog3(Intersection2);
%plot(3.96e4,mag4(Intersection2),'o')
xlim([50000 400000]);
%ylim([-40 25]);
ylabel('Amplitude Ratio (dB)','FontSize',12);
xlabel('Frequency (rad/s)','FontSize',12);
title('2nd Break Frequency','FontSize',14);
%legend('Location','best','Simulated Bode Plot','Experimental Bode Plot');
%%
Rp=162; %ohms
Rs=408;
Rm=1e6;
Lp=Rp/(break1);
Lo=(Rm+2*Rs)/(2*pi*break2);

plateaudB=-22; %dB, approx
plateau=10^(plateaudB/20);
x=.1;
omega=2000*2*pi;
%sens=(sqrt(((2*Rm*Km*x*omega)/((Rm+Rs)*Rp))^2))/(sqrt((((2*Lp*Lo*omega^2)/(Rp*(Rm+2*Rs)))+1)^2)+sqrt((omega*((Lp/Rp)+((2*Lo)/(Rm+2*Rs))))^2))
sensTerm= sqrt((plateau*(sqrt(((-(2*Lp*Lo*omega^2)/(Rp*(Rm+2*Rs)))+1)^2)+sqrt((omega*((Lp/Rp)+((2*Lo)/(Rm+2*Rs))))^2)))^2);
Km=(sensTerm*(Rm+Rs)*Rp)/(2*Rm*x*omega);
Gain=((2*Rm*Km*x)/((Rm+2*Rs)*Rp));
%%
tau1=Lp/Rp;
tau2=(2*Lo)/(Rm+2*Rs);

figure(4)
sys=tf([Gain 0],[tau1*tau2 tau1+tau2 1]);

bode(sys)
% %xlim([100 100000]);
% %ylim([-40 25]);
% ylabel('Amplitude Ratio (dB)','FontSize',12);
% title('Experimental Bode Response Plot','FontSize',14);
% legend('Location','best','Large Time Division','Small Time Division');
% subplot(2,1,2);
% semilogx(FreqPS,PS);
% 
% %xlim([100 100000]);
% %ylim([-180 0]);
% xlabel('Frequency (Hz)','FontSize',12);
% ylabel('Phase Shift (deg)','FontSize',12);
% legend('Location','best','Large Time Division','Small Time Division');

