% Explore Stimulation Patterns
% $Id: tutorial030b.m 5544 2017-06-16 16:30:20Z aadler $

% We have a 16 electrode EIT machine with adjacent drive
adjdrv= mk_stim_patterns(16,1, '{ad}','{ad}', ...
             {'no_meas_current','do_redundant'} );

% Arrange 16 electrodes in a zigzag
zigzag_mdl= fmdl; zigzag_mdl.stimulation= adjdrv;
zigzag_pat= [ 1:2:15;
             18:2:32]; 
zigzag_mdl.electrode= fmdl.electrode( zigzag_pat(:) );

subplot(131); show_stim_meas_pattern(zigzag_mdl,'ya')

% Arrange 16 electrodes as planar
planar_mdl= fmdl; planar_mdl.stimulation= adjdrv;
planar_pat= [ 1:2:15;
             17:2:32]'; 
planar_mdl.electrode= fmdl.electrode( planar_pat(:) );

subplot(132); show_stim_meas_pattern(planar_mdl,'ya')

% Arrange 16 electrodes as planar-opposite
pl_ops_mdl= fmdl; pl_ops_mdl.stimulation= adjdrv;
pl_ops_pat= [ rem( (0:7)*3, 8)*2+1;
              rem( (0:7)*3, 8)*2+17]'; 
pl_ops_mdl.electrode= fmdl.electrode( pl_ops_pat(:) );

subplot(133); show_stim_meas_pattern(pl_ops_mdl,'ya')

print_convert tutorial030b.png '-density 75';
