clear all; close all;

%% Part 1

% loading data
lvtdata = importdata('LVT v(t).lvm','\t',29);
lvtT = lvtdata.data(:,1);
lvtV = lvtdata.data(:,2);

% lvtP = cumtrapz(lvtT,lvtV);

% 1 a) Plot LVT Output Voltage vs. Time
figure(1);
% yyaxis left;
plot(lvtT,lvtV);
ylabel('Voltage (V)');
% yyaxis right;
% plot(lvtT,lvtP);
xlabel('Time (s)');
% ylabel('Integrated V (Vs)');

% Observation: Mass is free falling (linear slope) for the first 3 times
% i) first 200 data points are initial stuff
restVStd = std(lvtV(1:200));
% Find where the first free fall begin
startInd = 0;
for i = 1:length(lvtV)
    if (lvtV(i) < 0-restVStd)
        startInd = i;
        % find where the first free fall ends
        for j = startInd:length(lvtV)
            if (lvtV(j+1) - lvtV(j) > 0)
                endInd = j;
                break;
            end
        end
        break;
    end
end

% hold on;
% plot(lvtT(startInd),lvtV(startInd),'X',lvtT(endInd),lvtV(endInd),'O');