n = (1:100)';

% Bounds
Qmax = 25; % 450; % Qmax = 50;
Qmin = 0; % Qmin = 25;
SV = mean([11.55,15.00,12.70,10.83,11.92]) * 0.5;
Vbase = 0.075; % fraction of maximum

% Load flow data
data = readmatrix("QMed.txt");
tOrig = normalize(1:length(data),'range') * 99 + 1;
QOrig = data(:,1); QOrig = normalize(QOrig,'range') * (Qmax - Qmin) + Qmin;

% Uniformize data
Q = interp1(tOrig,QOrig,n);

% Integrate into Volume
V = normalize(cumtrapz(n,Q - mean(Q)) / length(n),'range'); 
V = V * SV + SV * Vbase;

% Plotting
figure;
plot(n,[Q,V])

% Save template
save('Qmed.mat','Q');
save('Vmed.mat','V');