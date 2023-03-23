% Load Hoogi's modulation template
load('templates.mat');
y = - T_resp_add_cell{6};

% Normalize
y = (y - min(y)) / (max(y) - min(y)) * 2 - 1;

% % Interpolate
% y = interp1(1:100,y,linspace(1,100,1000),'linear')';

% Plotting
figure;
plot(y)

% Save template
save('VawSpont.mat','y');