% Load Hoogi's modulation template
load('templates.mat');
y = T_resp_mod_cell{1};

% Normalize
y = (y - min(y)) / (max(y) - min(y)) * 2 - 1;

% Plotting
figure;
plot(y)

% Save template
save('Vmod.mat','y');