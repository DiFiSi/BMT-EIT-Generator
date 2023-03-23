% Define params
fs = 25; % target is 100 samples
Vt = 1;
tauFactor = 5;

fundFreq = 15; fundFreq = fundFreq / 60;
ie = 0.5;
tauIns = 1 / fundFreq * ie / tauFactor;
tauExp = 1 / fundFreq * (1 - ie) / tauFactor;

n = (1:1 / fundFreq * fs)';
nPeak = abs(n - tauFactor * tauIns * fs);
nPeak = find(nPeak == min(nPeak));
t = (n - 1) / fs;

% Define model
vInsFun = @(tau,t) (1 - exp(- t / tau));
vExpFun = @(tau,t,del) [ones(del,1);(exp(- t(1:del) / tau))];
          
% Build pulse
yIns = vInsFun(tauIns,t);% .* 
yExp = vExpFun(tauExp,t,nPeak);
y = Vt * yIns .* yExp;

% Normalize
y = (y - min(y)) / (max(y) - min(y)) * 2 - 1;

% Plotting
figure;
plot(t,y)

% Save template
save('VawMand.mat','y');