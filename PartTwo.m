clear all;
close all;

%% Part 2.1 and Part 2.2

% loading files
down = importdata('Deflect down.lvm','\t',32);
downData = down.data;
up = importdata('Deflect up.lvm','\t',32);
upData = up.data;
null = importdata('Null.lvm','\t',32);
nullData = null.data;

downMod = importdata('demod and filter below null.lvm','\t',32);
downModData = downMod.data;
upMod = importdata('demod and filter above null.lvm','\t',32);
upModData = upMod.data;
nullMod = importdata('demod and filter null.lvm','\t',32);
nullModData = nullMod.data;

% deflect down breakdown
coilOneDownTime = downData(:,1);
coilOneDownVolt = downData(:,2);
coilTwoDownTime = downData(:,3);
coilTwoDownVolt = downData(:,4);
% mod:
downModTime = downModData(:,1);
downModVolt = downModData(:,2);
downFiltTime = downModData(:,3);
downFiltVolt = downModData(:,4);

% deflect up breakdown
coilOneUpTime = upData(:,1);
coilOneUpVolt = upData(:,2);
coilTwoUpTime = upData(:,3);
coilTwoUpVolt = upData(:,4);
% mod:
upModTime = upModData(:,1);
upModVolt = upModData(:,2);
upFiltTime = upModData(:,3);
upFiltVolt = upModData(:,4);

% null breakdown
coilOneNullTime = nullData(:,1);
coilOneNullVolt = nullData(:,2);
coilTwoNullTime = nullData(:,3);
coilTwoNullVolt = nullData(:,4);
% mod:
nullModTime = nullModData(:,1);
nullModVolt = nullModData(:,2);
nullFiltTime = nullModData(:,3);
nullFiltVolt = nullModData(:,4);

%% Part a), before mod

% plot of coil 1 and 2 in deflect down position and mod output
figure(1);
hold on;
plot(coilOneDownTime,coilOneDownVolt);
plot(coilTwoDownTime,coilTwoDownVolt,'--');
plot(downModTime,downModVolt,'-.');
xlabel('Time (s)','FontSize',12);
ylabel('Voltage (V)','FontSize',12);
legend('Location','best','Coil 1 before mod','Coil 2 before mod','Mod output');

% plot of coil 1 and 2 in deflect up position and mod output
figure(2);
hold on;
plot(coilOneUpTime,coilOneUpVolt);
plot(coilTwoUpTime,coilTwoUpVolt,'--');
plot(upModTime,upModVolt,'-.');
xlabel('Time (s)','FontSize',12);
ylabel('Voltage (V)','FontSize',12);
legend('Location','best','Coil 1 before mod','Coil 2 before mod','Mod Output');

% plot of coil 1 and 2 in null position and mod output
figure(3);
hold on;
plot(coilOneNullTime,coilOneNullVolt);
plot(coilTwoNullTime,coilTwoNullVolt,'--');
plot(nullModTime,nullModVolt,'-.');
xlabel('Time (s)','FontSize',12);
ylabel('Voltage (V)','FontSize',12);
legend('Location','best','Coil 1 before mod','Coil 2 before mod','Mod Output');

% plot of all mod outputs
figure(4);
hold on;
plot(downModTime,downModVolt);
plot(upModTime,upModVolt,'--');
plot(nullModTime,nullModVolt,'-.');
xlabel('Time (s)','FontSize',12);
ylabel('Modulated Voltage (V)','FontSize',12);
legend('Location','best','Deflect Down Mod Output','Deflect Up Mod Output','Null Mod Output');

% plot(downFiltTime,downFiltVolt,'--');
% plot(upFiltTime,upFiltVolt,'--');
% plot(nullFiltTime,nullFiltVolt,'--');