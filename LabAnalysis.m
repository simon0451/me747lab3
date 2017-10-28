clear all; close all;

%% Part 1

% loading data
lvtdata = importdata('LVT v(t).lvm','\t',29);
lvtT = lvtdata.data(:,1);
lvtV = lvtdata.data(:,2);

% cut out the extra parts
lvtT = lvtT(1:1500);
lvtV = lvtV(1:1500);

restVStd = std(lvtV(1:200));
zeroCrossInd = zeros(1,3);
% Find where the first free fall begin
iniSlopeStartInd = 0;
iniSlopeEndInd = 0;
for i = 1:length(lvtV)
    if (lvtV(i) < 0-restVStd)
        iniSlopeStartInd = i;
        
        % find where the first free fall ends
        for j = iniSlopeStartInd:length(lvtV)
            if (lvtV(j+1) - lvtV(j) > 0)
                iniSlopeEndInd = j;
                
                % Find the next three zero crossings
                n = 1;
                for k = iniSlopeEndInd:length(lvtV)
                    % switch signs due to data going up and down
                    if (mod(n,2) == 1 && lvtV(k) > 0)
                        zeroCrossInd(n) = k-1;
                        n = n + 1;
                    elseif (mod(n,2) == 0 && lvtV(k) < 0)
                        zeroCrossInd(n) = k-1;
                        n = n + 1;
                    elseif (n == 4)
                        break;
                    end
                end
                break;
            end
        end
        break;
    end
end

% 1 a) Plot LVT Output Voltage vs. Time
figure(1);
plot(lvtT,lvtV);
hold on;
plot(lvtT(zeroCrossInd),lvtV(zeroCrossInd),'X');
grid on;
xlim([-0.4 0.7]);
xlabel('Time (s)','FontSize',12);
ylabel('Voltage (V)','FontSize',12);
legend('Location','best','LVT Output','Zero Crossings');

%% finding the slopes
p = polyfit(lvtT(iniSlopeStartInd:iniSlopeEndInd),lvtV(iniSlopeStartInd:iniSlopeEndInd),1);

% see document for origin of formula
sens = p(1)/-386; % [V / in/s] LVT sensitivity
lvtVelocity = lvtV/sens;

% integrate the signal
lvtPosition = cumtrapz(lvtT,lvtVelocity);
% adjust to rest at 0
lvtPosition = lvtPosition - mean(lvtPosition(1200:end));

% 1 f)
figure(2);
grid on;
xlim([-0.5 0.75]);
yyaxis left;
plot(lvtT,lvtVelocity);
ylabel('Velocity (in/s)','FontSize',12);
ylim([-30 40]);
yyaxis right;
plot(lvtT,lvtPosition);
ylim([-0.8 1.4]);
xlabel('Time (s)','FontSize',12);
ylabel('Position (in)','FontSize',12);
legend('Location','best','LVT Velocity','Integrated Position');

%% 1 g) Find damping ratio, damped frequency, and spring constant
% Location the fourth zero crossing as start of damped response
for i = zeroCrossInd(end)+1:length(lvtV)
    if (lvtV(i) < 0)
        dampStart = i;
        break;
    end
end

% reassign data range
dampedT = lvtT(dampStart:end);
dampedV = lvtV(dampStart:end);
% find the peaks
[maxV,~] = peakdet(dampedV,0.1);
maxV = maxV(2:end,:); % remove the first peak

TV = diff(dampedT(maxV(:,1)));
TV = TV(1); % use the first damped period

n = 4; % use all 4 peaks

sigmaV = (1/(n-1))*log(maxV(1,2)/maxV(4,2));
zetaV = sigmaV/sqrt(4*pi^2 + sigmaV^2); % the damping ratio
wdV = 2*pi/TV; % The damped frequency in rad/s -> 2pi rad/cyc div s/cyc
wnV = wdV/sqrt(1 - zetaV^2); % undamped natural frequency
mass = 74.9/1000; % [kg]
mass_english = mass*2.20462;
k = wnV^2*mass; % 1/(wn^2) = m/k -> k/m = wn^2 -> k = m*wn^2
k_english = k*2.20462/32.2; % convert from kg/s^2 to lbf/ft

figure(3);
hold on;
plot(dampedT,dampedV);
plot(dampedT(maxV(:,1)),dampedV(maxV(:,1)),'O');
grid on;
xlabel('Time (s)','FontSize',12);
ylabel('Voltage (V)','FontSize',12);

%% 1 h)
B = 2*zetaV*k/wnV; % B/k = 2*zeta/wn, [N/m/s]
B_english = B/3.28084*0.224809; % converted to [lbf/ft/s]

% damping force: Favg = B*vi
[maxVel,minVel] = peakdet(lvtVelocity,1);

vi = minVel(1:3,2);
Favg = B_english*vi/12; % [lbf / ft/s] * [in/s] = lbf in/ft = lbf*12;

% delV = 18.53+29.5;%lvtVelocity(maxVel(2,1)) - lvtVelocity(minVel(1,1));
% delT = -0.03489+0.06897;%lvtT(maxVel(2,1))-lvtT(minVel(1,1));
% Favg2 = mass_english*delV/delT/12; % currently has units of lbm ft/s^2 = lbf

figure(4);
hold on;
plot(lvtT,lvtVelocity);
plot(lvtT(maxVel(:,1)),lvtVelocity(maxVel(:,1)),'O');
plot(lvtT(minVel(:,1)),lvtVelocity(minVel(:,1)),'O');