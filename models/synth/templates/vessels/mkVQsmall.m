n = (1:100)';

% Bounds
Qmax = 5; % 250; % Qmax = 1;
Qmin = 0; % Qmin = 50;
SV = mean([3.47,3.68,4.50,7.56,7.84]) * 0.5;
Vbase = 0.05; % fraction of maximum

% Load flow data
data = readmatrix("QSmall.txt");
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
save('Qsmall.mat','Q');
save('Vsmall.mat','V');