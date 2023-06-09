% Show lung images $Id: eidors_colours06.m 4855 2015-04-02 15:34:44Z aadler $
load montreal_data_1995
imdl = mk_common_model('c2t2',16);
imdl.fwd_model = mdl_normalize(imdl.fwd_model,1);
imdl.RtR_prior= @prior_gaussian_HPF;
imdl.hyperparameter.value = 0.45;
img = inv_solve(imdl, zc_resp(:,1), zc_resp(:,20));

img.calc_colours.ref_level= 0;
img.calc_colours.cb_shrink_move = [0.5,0.8,-.10];

clf; subplot(221);
show_fem(img,[1,1]);
axis equal; axis off; axis tight;

subplot(222);
img.calc_colours.ref_level =  0.1;
show_fem(img,[1,1]);
axis equal; axis off; axis tight;

print_convert eidors_colours06.png '-density 75'
