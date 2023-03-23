n = (1:100)';

% Bounds
Qmax = 550;
Qmin = 0;
Vmax = mean([13.35,13.19,15.20,10.56,8.71]);

% Load data
data = readmatrix("QLarge.txt");
tOrig = data(:,1); tOrig = normalize(tOrig,'range') * 99 + 1;
QOrig = data(:,2); QOrig = normalize(QOrig,'range') * (Qmax - Qmin) + Qmin;

% Uniformize data
Q = interp1(tOrig,QOrig,n);

% Integrate into Volume
V = cumtrapz(n,Q - mean(Q)) / length(n); 
V = V - min(V);

% Plotting
figure;
plot(n,[Q,V])

% Save template
save('Qlarge.mat','Q');
save('Vlarge.mat','V');