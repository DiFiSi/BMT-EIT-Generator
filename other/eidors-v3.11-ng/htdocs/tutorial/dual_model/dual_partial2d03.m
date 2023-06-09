% Dual Partial $Id: dual_partial2d03.m 4839 2015-03-30 07:44:50Z aadler $

% base model
imdl.RtR_prior = @prior_noser;
imdl.hyperparameter.value = .2;

% Reconstruction model - only fine reconstruction
frec_mdl= imdl;
frec_mdl.fwd_model= fmdl; % fine model

% Reconstruction model - dual model
drec_mdl = imdl;
drec_mdl.fwd_model= fmdl; % fine model
% coarse to fine mapping
nf_els= size(fmdl.elems,1); nc_els= size(cmdl.elems,1);
drec_mdl.fwd_model.coarse2fine= sparse(c2f_idx, 1:nc_els, 1, nf_els, nc_els);
