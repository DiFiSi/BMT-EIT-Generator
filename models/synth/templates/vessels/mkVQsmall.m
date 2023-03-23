n = (1:100)';

% Bounds
Qmax = 250;
Qmin = 50;
Vmax = mean([3.47,3.68,4.50,7.56,7.84]);

% Load flow data
data = readmatrix("QSmall.txt");
tOrig = normalize(1:length(data),'range') * 99 + 1;
QOrig = data(:,1); QOrig = normalize(QOrig,'range') * (Qmax - Qmin) + Qmin;

% Uniformize data
Q = interp1(tOrig,QOrig,n);

% Integrate into Volume
V = cumtrapz(n,Q - mean(Q)) / length(n); 
V = V - min(V);

% Plotting
figure;
plot(n,[Q,V])

% Save template
save('Qsmall.mat','Q');
save('Vsmall.mat','V');