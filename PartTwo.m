clear all;
close all;

%% Part 2.1

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
downModTime = downModData(:,3);
downModVolt = downModData(:,4);
downFiltTime = downModData(:,1);
downFiltVolt = downModData(:,2);

% deflect up breakdown
coilOneUpTime = upData(:,1);
coilOneUpVolt = upData(:,2);
coilTwoUpTime = upData(:,3);
coilTwoUpVolt = upData(:,4);
% mod:
upModTime = upModData(:,3);
upModVolt = upModData(:,4);
upFiltTime = upModData(:,1);
upFiltVolt = upModData(:,2);

% null breakdown
coilOneNullTime = nullData(:,1);
coilOneNullVolt = nullData(:,2);
coilTwoNullTime = nullData(:,3);
coilTwoNullVolt = nullData(:,4);
% mod:
nullModTime = nullModData(:,3);
nullModVolt = nullModData(:,4);
nullFiltTime = nullModData(:,1);
nullFiltVolt = nullModData(:,2);

%% Part a), before mod

% plot of coil 1 and 2 in deflect down position and mod output
figure(1);
hold on;
plot(coilOneDownTime,coilOneDownVolt);
plot(coilTwoDownTime,coilTwoDownVolt,'--');
xlabel('Time (s)','FontSize',12);
ylabel('Voltage (V)','FontSize',12);
legend('Location','best','Coil 1 before mod','Coil 2 before mod');

% plot of coil 1 and 2 in deflect up position and mod output
figure(2);
hold on;
plot(coilOneUpTime,coilOneUpVolt);
plot(coilTwoUpTime,coilTwoUpVolt,'--');
xlabel('Time (s)','FontSize',12);
ylabel('Voltage (V)','FontSize',12);
legend('Location','best','Coil 1 before mod','Coil 2 before mod');

% plot of coil 1 and 2 in null position and mod output
figure(3);
hold on;
plot(coilOneNullTime,coilOneNullVolt);
plot(coilTwoNullTime,coilTwoNullVolt,'--');
xlabel('Time (s)','FontSize',12);
ylabel('Voltage (V)','FontSize',12);
legend('Location','best','Coil 1 before mod','Coil 2 before mod');

% plot of all mod outputs
figure(4);
hold on;
plot(downModTime,downModVolt);
plot(upModTime,upModVolt,'--');
plot(nullModTime,nullModVolt,'k-.');
xlabel('Time (s)','FontSize',12);
ylabel('Modulated Voltage (V)','FontSize',12);
legend('Location','best','Deflect Down Mod Output','Deflect Up Mod Output','Null Mod Output');

% plot of all filtered outputs
figure(5);
hold on;
plot(downFiltTime,downFiltVolt);
plot(upFiltTime,upFiltVolt,'--');
plot(nullFiltTime,nullFiltVolt,'k-.');
xlabel('Time (s)','FontSize',12);
ylabel('Modulated Voltage (V)','FontSize',12);
legend('Location','best','Deflect Down Filtered Output','Deflect Up Filtered Output','Null Filtered Output');

%% Part 2.2 a) - c)
% data from lab
mass_bucket = 7.65; % [kg]
total_masses = ([500,200,100,50,20,0,-mass_bucket] + mass_bucket)/1000*2.20462;
total_weight = total_masses*32.2; % [lbm -> lbf] * 32.2 ft/s^2
VoutWeight = [-643.76,-247.43,-140.55,-76.09,-28.59,-10.475,1.335]/1000; % [V]

deflections = [0,-0.025,-0.05,-0.075,-0.1,-0.125,-0.150]; % [in]
VoutDef = [-5.9836,-79.61,-153.93,-227.99,-302.31,-376.57,-450.23]/1000; % [V]

% linearity of the curves
% fit the curves
p1 = polyfit(total_weight,VoutWeight,1);
p2 = polyfit(deflections,VoutDef,1);

% make the fit lines
weight = total_weight;
Vw = weight*p1(1) + p1(2);

def = deflections;
Vdef = def*p2(1) + p2(2);

% R squared formula
Rsqw = 1 - sum((VoutWeight - Vw).^2)/sum((VoutWeight - mean(VoutWeight)).^2); % weight
RsqD = 1 - sum((VoutDef - Vdef).^2)/sum((VoutDef - mean(VoutDef)).^2); % deflection

% finding spring constant K
% match weight to displacement
def_equivalent = (Vw - p2(2))/p2(1);
% fit the matched displacement-weight data
p3 = polyfit(def_equivalent,weight,1);
K = p3(1);
def_eq = def_equivalent;
w = def_eq*p3(1) + p3(2); % make the best fit line

figure(6);
hold on;
plot(total_weight,VoutWeight,'X');
plot(weight,Vw,'--');
grid on;
xlabel('Weight (lb_f)','FontSize',12);
ylabel('Output Voltage (V)','FontSize',12);

figure(7);
hold on;
plot(deflections,VoutDef,'O');
plot(def,Vdef,'--');
grid on;
xlabel('Deflections (in)','FontSize',12);
ylabel('Output Voltage (V)','FontSize',12);

figure(8);
hold on;
plot(def_equivalent,weight,'x');
plot(def_eq,w,'--');
grid on;
xlabel('Displacements (in)','FontSize',12);
ylabel('Weight (lb_f)','FontSize',12);

%% Part 2.2 d)
% load data
vibrationData = importdata('2_2_4 Vibrating.lvm','\t',29);
vibrationT = vibrationData.data(5000:38000,1);
vibrationV = vibrationData.data(5000:38000,2);

% usual zeta and frequency calculations using log decrement
[maxV,~] = peakdet(vibrationV,0.1);
TV = diff(vibrationT(maxV(:,1)));
TV = TV(1);

n = 4; % using all the peaks
sigmaV = (1/(n-1))*log(maxV(1,2)/maxV(4,2));
zetaV = sigmaV/sqrt(4*pi^2 + sigmaV^2); % the damping ratio
wdV = 2*pi/TV; % The damped frequency in rad/s -> 2pi rad/cyc div s/cyc
wnV = wdV/sqrt(1 - zetaV^2); % undamped natural frequency

figure(9);
hold on;
plot(vibrationT,vibrationV);
plot(vibrationT(maxV(:,1)),vibrationV(maxV(:,1)),'O');
xlabel('Time (s)','FontSize',12);
ylabel('Voltage (V)','FontSize',12);

%% part e)

% natural frequency = sqrt(spring constant / effective mass)
Meff = -K/wnV^2;

%% part f)

% input - RI - 1/Cs I = 0; output = 1/Cs I
% I = CsOutput; input = (R+1/Cs)CsOutput
% input = RCsOutput + Output
% output/input = 1/(RCs+1)
R = 100e3;
C = 0.047e-6;

sys = tf(1,[R*C 1]);
bode(sys);