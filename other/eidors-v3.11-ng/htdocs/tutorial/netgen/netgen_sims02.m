% Netgen simulation $Id: netgen_sims02.m 3791 2013-04-04 15:48:25Z aadler $

% Calculate stimulation pattern adjacent
stim = mk_stim_patterns(9,1,[0,1],[0,1],{},1);
img.fwd_model.stimulation = stim;

% Get all voltages so we can plot it
img.fwd_solve.get_all_meas = 1;

% Homogeneous model
img.elem_data(:) = 1;

vh = fwd_solve(img);

% Show Voltage for stim pattern #1
imgn = rmfield(img,'elem_data');
imgn.node_data = vh.volt(:,1);

h1= subplot(221);
show_fem(imgn);

% Show Voltage for stim pattern #2
imgn = rmfield(img,'elem_data');
imgn.node_data = vh.volt(:,2);

h2= subplot(222);
show_fem(imgn);

imgn.calc_colours.cb_shrink_move = [0.5,0.8,-.02];
common_colourbar([h1,h2],imgn);
print_convert netgen_sims02a.png '-density 100'
