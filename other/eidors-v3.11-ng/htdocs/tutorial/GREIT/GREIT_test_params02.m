% Simulate obj $Id: GREIT_test_params02.m 4823 2015-03-29 15:17:16Z bgrychtol-ipa $

% Specify positions to simulate (only x-axis here)
r =  linspace(0,0.9,100);
xyzr = [r; zeros(1,100); ones(1,100);
     0.05*ones(1,100)];

[vh,vi] = simulate_movement(imgs, xyzr);

% Show GREIT images
opt.noise_figure = 0.5;
i_gr = mk_GREIT_model(fmdl,0.2,[],opt);
imgr = inv_solve(i_gr, vh, vi(:,1:5:100));
imgr.show_slices.img_cols = 5;
show_slices(imgr);

print_convert('GREIT_test_params02a.png','-density 60');
