n = (1:100)';

% Bounds
Qmax = 35; % 550; % Qmax = 50;
Qmin = 0; % Qmin = 0;
SV = mean([13.35,13.19,15.20,10.56,8.71]) * 0.5;
Vbase = 0.10; % fraction of maximum

% Load data
data = readmatrix("QLarge.txt");
tOrig = data(:,1); tOrig = normalize(tOrig,'range') * 99 + 1;
QOrig = data(:,2); QOrig = normalize(QOrig,'range') * (Qmax - Qmin) + Qmin;

% Uniformize data
Q = interp1(tOrig,QOrig,n);

% Integrate into Volume
V = normalize(cumtrapz(n,Q - mean(Q)) / length(n),'range'); 
V = V * SV + SV * Vbase;

% Plotting
figure;
plot(n,[Q,V]);

% Save template
save('Qlarge.mat','Q');
save('Vlarge.mat','V');