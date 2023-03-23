fs = 100; % target is 100 samples
fundFreq = 60; fundFreq = fundFreq / 60;
n = (1:1 / fundFreq * fs)';
t = (n - 1) / fs;

y = readmatrix("Vventricle.txt");
y = y(:,2);
tOrig = 1:length(y); tOrig = normalize(tOrig,'range');

% Normalize
y = (y - min(y)) / (max(y) - min(y)) * 2 - 1;

% Interpolate between t(1) and t(end)
y = interp1(tOrig,y,t,'spline');

% Plotting
figure;
plot(y)

% Save template
save('Vventricle.mat','y');