% Reconstruction Model $Id: basic_3d_05.m 3126 2012-06-08 16:17:56Z bgrychtol $
J = calc_jacobian( calc_jacobian_bkgnd( imdl) );
iRtR = inv(prior_noser( imdl ));
hp = 0.17;
iRN = hp^2 * speye(size(J,1));
RM = iRtR*J'/(J*iRtR*J' + iRN);
imdl.solve = @solve_use_matrix; 
imdl.solve_use_matrix.RM  = RM;
