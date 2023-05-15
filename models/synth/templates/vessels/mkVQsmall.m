% fs = 50;
n = (1:100)';
% t = (n - 1) / fs;

% Bounds
Qmax = 5; % 250; % Qmax = 1;
Qmin = 0; % Qmin = 50;
Vmax = mean([3.47,3.68,4.50,7.56,7.84]);
Vbase = 0.7; % fraction of maximum

% Load flow data
data = readmatrix("QSmall.txt");
tOrig = normalize(1:length(data),'range') * 99 + 1;
QOrig = data(:,1); QOrig = normalize(QOrig,'range') * (Qmax - Qmin) + Qmin;

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
plot(n,[Q,V])

% Save template
save('Qsmall.mat','Q');
save('Vsmall.mat','V');