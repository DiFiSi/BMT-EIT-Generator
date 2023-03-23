% Define params
fs = 100;
fundFreq = 60; fundFreq = fundFreq / 60;

AHarms = [10;5;2;1;0.5];
phiHarms = [0; 0.4 * pi; 0.6 * pi; 0.8 * pi; pi];

n = (1:1 / fundFreq * fs)';
t = (n - 1) / fs;

% Define model
cosFun = @(A,phi,f,fs,n) A * cos(2 * pi * f / fs * n + phi);
          
% Build pulse harmonics
y = 0;
nHarms = length(AHarms);
for h = 1:nHarms
   A = AHarms(h);
   phi = phiHarms(h);
   f = fundFreq * h;
   
   y = y + cosFun(A,phi,f,fs,n); 
end

% Shift and normalize
y = (circshift(y,40) - min(y)) / (max(y) - min(y)) * 2 - 1;

% Plotting
figure;
plot(t,y);

% Save template
save('Vblood.mat','y');