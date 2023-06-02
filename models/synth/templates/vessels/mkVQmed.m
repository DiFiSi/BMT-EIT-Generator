% fs = 50;
n = (1:100)';
% t = (n - 1) / fs;

% Bounds
% Qmax = 25; % 450; % Qmax = 50;
% Qmin = 0; % Qmin = 25;
Vmax = mean([11.55,15.00,12.70,10.83,11.92]);
Vbase = 0.7; % fraction of maximum

% Load flow data
data = readmatrix("QMed.txt");
tOrig = normalize(1:length(data),'range') * 99 + 1;
QOrig = data(:,1); QOrig = normalize(QOrig,'range'); % * (Qmax - Qmin) + Qmin;

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
yyaxis left;
plot(n,Q);
ylabel("Flow [cm3/s]")
yyaxis right;
plot(n,V);
ylabel("Volume [cm3]")

% Save template
save('Qmed.mat','Q');
save('Vmed.mat','V');