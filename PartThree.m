clear all;
close all;

%% Part 3
refinedFreqPS = xlsread('2_3_1 data.xlsx',1,'A10:A10009');
refinedPS = xlsread('2_3_1 data.xlsx',1,'B10:B10009');
refinedFreqMag = xlsread('2_3_1 data.xlsx',2,'A10:A10009');
refinedMag = xlsread('2_3_1 data.xlsx',2,'B10:B10009');
FreqPS = xlsread('2_3_1 data.xlsx',3,'A10:A10009');
PS = xlsread('2_3_1 data.xlsx',3,'B10:B10009');
FreqMag = xlsread('2_3_1 data.xlsx',4,'A10:A10009');
Mag = xlsread('2_3_1 data.xlsx',4,'B10:B10009');
%%

%p=polyfit(FreqMag(100:800),Mag(100:800),1);
xdata=FreqMag(100:800);
ydata=FreqMag(100:800);

p=polyfit(log10(xdata),ydata,1);
wlog=logspace(2,3);
mag=p(1)*log10(wlog)+p(2);

figure(1)

subplot(2,1,1);
semilogx(FreqMag,Mag);
hold on;
grid on;
semilogx(refinedFreqMag,refinedMag,'k');
semilogx(wlog,mag,'r')
%xlim([10 10000]);
%ylim([-40 25]);
ylabel('Amplitude Ratio (dB)','FontSize',12);
title('Experimental Bode Response Plot','FontSize',14);
%legend('Location','best','Simulated Bode Plot','Experimental Bode Plot');

subplot(2,1,2);
semilogx(FreqPS,PS);
hold on;
grid on;
semilogx(refinedFreqPS,refinedPS,'k');
%xlim([10 10000]);
%ylim([-180 0]);
xlabel('Frequency (Hz)','FontSize',12);
ylabel('Phase Shift (deg)','FontSize',12);
%legend('Location','best','Simulated Bode Plot','Experimental Bode Plot');



