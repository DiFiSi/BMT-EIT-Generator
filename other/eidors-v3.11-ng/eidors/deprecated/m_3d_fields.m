function [v_f] = m_3d_fields(vtx,el_no,m_ind,E,tol,gnd_ind,v_f);
%function [v_f] = m_3d_fields(vtx,el_no,m_ind,E,tol,gnd_ind,v_f);
%
%This function calculates the measurement fields using preconditioned conjugate gradients.
%
%
%
%vtx     = The vertices
%el_no   = The total number of electrodes in the system
%m_ind   = The measurements matrix (indices of electrode pairs)
%E       = The full rank system matrix
%tol     = The tolerance in the forward solution
%gnd_ind = The ground index
%v_f     = The measurements fields

%FIXME: This code should call forward_solver, it can then decide
%what the best solver strategy is. Right now cgls is slower than \

% (C) Nick Polydorides GPL v2 or v3. $Id: m_3d_fields.m 6502 2022-12-30 14:29:12Z aadler $


warning('EIDORS:deprecated','M_3D_FIELDS is deprecated as of 06-Jun-2012. ');

[vr,vc] = size(vtx);

Is_supl = zeros(vr,size(m_ind,1));
%no of electrodes x no of measurements (now currents)!

MC = [];

for i=1:size(m_ind,1)

   m_n = zeros(el_no,1);

   m_n(m_ind(i,1)) = 1;
   m_n(m_ind(i,2)) = -1;

   MC = [MC,m_n];

end

I = [Is_supl;MC];
I(gnd_ind,:) = 0;

% stupidity to be matlab 6+7 compatible
if nargin < 7;  v_f = zeros(size(E,1),size(I,2)); end
if isempty(v_f);v_f = zeros(size(E,1),size(I,2)); end

maxiter=10000; % This should be high enough, but it may maybe this should
% depend on the number of measurements?


if isreal(E)==1

   ver = eidors_obj('interpreter_version');

   if ver.isoctave % OCtave doesn't have Cholinc yet (as of 2.9.13)
      M= [];
   else
      opts.droptol = tol/10;
      if ver.ver < 7.012
         M = cholinc(E,opts.droptol);
      else
         opts.droptol = 1e-6;
         opts.type = 'ict'; %otherwise droptol is ignored opts.type = 'nofill';

         %         M = ichol(E,opts);
         % ichol makes pcg even slower. It's better to use no pre-conditioner
         M = [];
      end
   end

   for y=1:size(MC,2)
      %Set this line to suit your approximation needs. ***************
      %for more details use help pcg on Matlab's command window.
      [v_f(:,y),flag,relres,iter,resvec] = pcg(E,I(:,y), ...
      tol*norm(I(:,y)),maxiter,M',M,v_f(:,y));
   end

else  %is real

   %Preconditioner
   % LUinc no longer available (matlab 6 or so)
%  [L,U] = luinc(E,tol/10);
%  [L,U] = ilu(E, struct('droptol',E/10));
   opts.droptol = tol/10 
   opts.type = 'nofill'; %otherwise droptol is ignored opts.type = 'nofill';
   [L,U] = ilu(E, opts);

   for y=1:size(MC,2)

      [v_f(:,y),flag,relres,iter,resvec] = bicgstab(E,I(:,y), ...
      tol*norm(I(:,y)),maxiter,L,U);

   end
end %is complex


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
