% fs = 50;
n = (1:100)';
% t = (n - 1) / fs;

% Bounds
Qmax = 35; % 550; % Qmax = 50;
Qmin = 0; % Qmin = 0;
Vmax = mean([13.35,13.19,15.20,10.56,8.71]);
Vbase = 0.7; % fraction of maximum

% Load data
data = readmatrix("QLarge.txt");
tOrig = data(:,1); tOrig = normalize(tOrig,'range') * 99 + 1;
QOrig = data(:,2); QOrig = normalize(QOrig,'range') * (Qmax - Qmin) + Qmin;

% Uniformize data
Q = interp1(tOrig,QOrig,n);

% Integrate into Volume
V = normalize(cumtrapz(n,Q - mean(Q)) / length(n),'range'); 
V = V * (1 - Vbase) * Vmax + Vmax * Vbase;

% Reobtain Q
Q = diff(V);
Q = [Q(1); Q] - min(Q);

% Plotting
figure;
plot(n,[Q,V]);

% Save template
save('Qlarge.mat','Q');
save('Vlarge.mat','V');