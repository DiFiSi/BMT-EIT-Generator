function img= inv_solve_TV_pdipm( inv_model, data1, data2)
% INV_SOLVE_TV_PDIPM inverse solver for Andrea Borsic's
%   Total Variation solver for use with difference EIT
% img= inv_solve_TV_pdipm( inv_model, data1, data2)
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
% Parameters (see primaldual_tvrecon_lsearch for description)
%   alpha1 (imdl.inv_solve_TV_pdipm.alpha1)
%   alpha2 (set with imdl.hyperparameter.value)
%   beta   (imdl.inv_solve_TV_pdipm.beta)
%   want_dual_variable  (set to 1 if you want access to dual)
% Termination parameters
%  max_iters =  inv_model.parameters.max_iteration (default 10)
%      Max number of iterations before stopping
%  min change = inv_model.parameters.min_change   (default 0)
%      Min Change in objective fcn (norm(y-Jx)^2 + hp*TV(x)) before stopping

% (C) 2002-2009 Andrea Borsic and Andy Adler. License: GPL version 2 or version 3
% $Id: inv_solve_TV_pdipm.m 6247 2022-03-23 13:47:00Z aadler $


p= get_params(inv_model);

dva = calc_difference_data( data1, data2, inv_model.fwd_model);
% TEST CODE -> Put elsewhere
back_val = get_good_background(inv_model, data1);
inv_model.jacobian_bkgnd.value= back_val;

sol= [];
for i=1:size(dva,2)
   [soln,dual_x]=primaldual_tvrecon_lsearch(inv_model, dva(:,i), ...
        p.maxiter,p.alpha1,p.alpha2, p.tol, p.beta, p.min_change);

   if ~p.keepiters
      soln=soln(:,end);
   end

   sol=[sol, soln];
end

img.name= 'solved by inv_solve_TV_pdipm';
img.elem_data = sol;
img.fwd_model= inv_model.fwd_model;
try if inv_model.inv_solve_TV_pdipm.want_dual_variable
   img.dual_data = dual_x;
end; end

function p = get_params(inv_model);
   
   try   
       p.alpha1= inv_model.inv_solve_TV_pdipm.alpha1;
   catch
       p.alpha1= 1e-2;
   end

   try   
       p.beta= inv_model.inv_solve_TV_pdipm.beta;
   catch
       p.beta= 1e-4;
   end

   p.alpha2= calc_hyperparameter( inv_model);

   try   
       p.min_change = inv_model.parameters.min_change;
   catch
       p.min_change = 0;
   end

   try   
       p.maxiter = inv_model.parameters.max_iterations;
   catch
       p.maxiter= 10;
   end

   try   
       p.keepiters = inv_model.parameters.keep_iterations;
   catch
       p.keepiters= 0;
   end

   p.tol = 0; % TODO

function back_val = get_good_background(inv_mdl, data1);

   % Create homogeneous model
   IM= eidors_obj('image','');
   IM.fwd_model= inv_mdl.fwd_model;
   s= ones(size(IM.fwd_model.elems,1),1);
   IM.elem_data= s;

   vsim= fwd_solve( IM);
   back_val=abs( data1\vsim.meas );
   back_val=1;

