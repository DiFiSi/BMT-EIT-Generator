n = (1:100)';

% Bounds
Qmax = 450;
Qmin = 25;
Vmax = mean([11.55,15.00,12.70,10.83,11.92]);

% Load flow data
data = readmatrix("QMed.txt");
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
save('Qmed.mat','Q');
save('Vmed.mat','V');