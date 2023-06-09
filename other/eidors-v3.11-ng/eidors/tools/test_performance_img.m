function [params_img] =  test_performance_img( imdls, fmdl );
% TEST_PERFORMANCE: test of difference reconstruction algorithms
%   in terms of the performance parameters defined in GREIT, image output
%
% Input:
%   imdls = cell array of inverse models to test on
%   fmdl  = fmdl of the shape to test on (default is 16 electrode tank).
%
%
% params_img =  test_performance( imdls );
%  params_img = {AR, PE, RES, SD, RNG, NOI}
%  is a cell array of 3D (Np x Np x N_models) arrays, each representing 
%  the given perfomance metrics as an image computed using regularly spaced
%  targets.
% 
%
% Example
%   i_bp = mk_common_gridmdl('backproj');
%   i_gp = mk_common_gridmdl('GREITc1');
%   test_performance( { i_gp, i_gr })

% (c) 2011 Andy Adler and Barlomiej Grychtol. Licenced under GPL v2 or v3
% $Id: test_performance_img.m 5112 2015-06-14 13:00:41Z aadler $

if ischar(imdls) && strcmp(imdls,'UNIT_TEST'); do_unit_test; return; end

if nargin==1
   fmdl = ng_mk_cyl_models([2,1,0.05],[16,1],[0.05]); 
   fmdl.stimulation = mk_stim_patterns(16,1,[0,1],[0,1],{},1);
end
switch fmdl.type
    case 'fwd_model'
        imgs= mk_image( fmdl, 1);
    case 'image'
        imgs = fmdl;
        fmdl = fmdl.fwd_model;
    otherwise
        error('Second parameter not understood');
end
% if ~iscell(imdls)
%     imdls = {imdls};
% end

if (0) % old code
Nsim = 100;
r =  linspace(0,0.9,Nsim);
ctr = [0,0, mean(fmdl.nodes(:,3))];  % Assume x,y centre is zero
maxx = max(abs(fmdl.nodes(:,1) - ctr(1)));
maxy = max(abs(fmdl.nodes(:,2) - ctr(2)));

th=  r*4321; % want object to jump around in radius
xyzr = [maxx*r.*cos(th); maxy*r.*sin(th); ctr(3)*ones(1,Nsim); 
        0.05/mean([maxx,maxy])*ones(1, Nsim)];
else 
%%
    N = 128;
     % this bit is copied from mdl_slice_mapper to make sure
     xmin = min(fmdl.nodes(:,1));    xmax = max(fmdl.nodes(:,1));
     xmean= mean([xmin,xmax]); xrange= xmax-xmin;

     ymin = min(fmdl.nodes(:,2));    ymax = max(fmdl.nodes(:,2));
     ymean= mean([ymin,ymax]); yrange= ymax-ymin;

     range= max([xrange, yrange]);
     xspace = linspace( xmean - range*0.5, xmean + range*0.5, N );
     yspace = linspace( ymean + range*0.5, ymean - range*0.5, N );
     % end of copied block
    bnd_nodes = unique(fmdl.boundary);
    min_bb = min(fmdl.nodes(bnd_nodes,:));
    max_bb = max(fmdl.nodes(bnd_nodes,:));
    h = mean([min_bb(3),max_bb(3)]);
    r = 0.025 * range;
%     xspace = linspace(min_bb(1),max_bb(1),N);
%     yspace = linspace(max_bb(2),min_bb(2),N);
    [X Y] = meshgrid(xspace,yspace);
    img = mk_image(fmdl,1);
    img.calc_colours.npoints = N;
    M = calc_slices(img,1);
    IN = M==1;
    n_px = nnz(IN);
    OUT = isnan(M);
    xyzr = [X(IN)'; Y(IN)'; h*ones(1,n_px);  r*ones(1,n_px); ];
%%
end
    
    
[vh,vi] = simulate_movement(imgs, xyzr);
IN6 = repmat(IN,[1 1 6]); IN6 = shiftdim(IN6,2);
for i= 1:length(imdls);
    if iscell(imdls)
        imgr = inv_solve(imdls{i}, vh, vi);
        pnoise = calc_noise_figure( imdls{i}, vh, vi );
    else
        imgr = inv_solve(imdls(i), vh, vi);
        pnoise = calc_noise_figure( imdls(i), vh, vi );
    end
   imgr.calc_colours.npoints = N;
   param_GR = eval_GREIT_fig_merit(imgr, xyzr);
   params{i} = [param_GR; pnoise];
   p_img = nan([6 size(IN)]);
   p_img(IN6) = params{i};
   params_img{i}=p_img;
end




return
clf; Nparams = 6;
for j=1:Nparams;
   p = [];
   for i= 1:length(params);
      p = [p, params{i}(j,:)'];
   end

%   subplot(5,1,j);
   hig = 0.95/Nparams;
   axes('position',[0.1,0.05 + hig*(Nparams-j),.88, hig]);
   plot(r, p);
   if j==1;     axis([0,0.9,0,2.1]);        ylabel('AR');
   elseif j==2; axis([0,0.9,-0.16,0.16]);   ylabel('PE');
   elseif j==3; axis([0,0.9,0,0.41]);       ylabel('RES');
   elseif j==4; axis([0,0.9,0,0.31]);       ylabel('SD');
   elseif j==5; axis([0,0.9,0,0.61]);       ylabel('RNG');
   elseif j==6; set(gca,'YScale','log');    ylabel('NF');
   end
   if j<Nparams; set(gca,'XTickLabel',[]);end
end









function do_unit_test;
% Reconstruct GREIT Images
imdl_gr = mk_common_gridmdl('GREITc1');

% Reconstruct backprojection Images
imdl_bp = mk_common_gridmdl('backproj');

% Reconstruct GN Images
imdl_gn = select_imdl( mk_common_model('d2c2', 16), {'Basic GN dif','Choose NF=0.5'});

test_performance( { imdl_gr, imdl_bp, imdl_gn } );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vh,vi,param]=test_performance_old(fmdl,imdl,r,N)

debug = false;

%% 1. Figure out the limits of fmdl
boundary = fmdl.boundary;
b_ver = unique(boundary);
boundary = fmdl.nodes(b_ver,:);
min_b = min(boundary)+repmat(0.5*r,1,3);
max_b = max(boundary)-repmat(0.5*r,1,3);
vol = get_elem_volume(fmdl);

%% 2. create N samples within the volume
rand('state', sum(100*clock));

n_elems = size(fmdl.elems,1);


c2f = eidors_obj('get-cache', {fmdl, N, r}, 'samples');
xyzr_pt = eidors_obj('get-cache', {fmdl, N, r}, 'points');

count = 0; %count points already created

if isempty(c2f) || isempty(xyzr_pt);
    disp('Generating sample points.');
    c2f = sparse(n_elems,N);
    while count < N
        n_pts = ceil(1.5 * (N - count) ); % more than needed to limit iterations
        vec = rand_pt ( repmat(min_b,n_pts,1), repmat(max_b,n_pts,1) );
        map = mk_c2f_circ_mapping(fmdl, [vec repmat(r,n_pts,1)]');
        frac = vol'*map/(4/3*pi*r^3);
        idx = find(frac > 0.9, N - count, 'first');
        
        xyzr_pt = [xyzr_pt, vec(idx,:)'];
        newcount = count + length(idx);
        c2f(:, (count+1):newcount ) = map(:,idx);
        count = newcount;
    end
    eidors_obj('set-cache', {fmdl, N, r}, 'samples', c2f);
    eidors_obj('set-cache', {fmdl, N, r}, 'points', xyzr_pt );
    disp('Done.');
end
 


% x=110;y=150;r=20;
% x=110;y=80;r=20;
% map = mk_c2f_circ_mapping(fmdl, [x,y,r]');
% img= mk_image(map);img.fwd_model = fmdl;
% phi= linspace(0,2*pi,20);xr=r*cos(phi)+x; yr=r*sin(phi)+y;
% show_fem(img,[0,1,1]); hold on; plot(xr,yr,'r');hold off
% vol = get_elem_volume(fmdl); vol'*map/(pi*r^2)

%% 3. reconstruct the images
% create a homegenous image
   elem_data = ones(size(fmdl.elems,1),1);
   img = mk_image(fmdl, elem_data);
   disp('Calculating Jacobian. This may take a (long) while.');
   J = calc_jacobian( img );
   disp('Done.');
   vh= fwd_solve(img);
   vh = vh.meas;
   vi= vh*ones(1,N) + J*c2f;

%% 4. calculate parameters


   
   img = inv_solve(imdl,vh,vi);
   img.calc_colours.npoints=32;
   imgs = calc_slices(img);

       param= [calc_noise_figure(imdl, vh, vi); ... % noise parameters
           GREIT_sim_params(  imgs, xyzr_pt)];            % image parameters
       
       r = sqrt(sum(xyzr_pt(1:2,:).^2));
       
       names = {'Noise','Amplitude','Position Error','Resolution', ...
           'Shape Deformation', 'Ringing'};
       for j=1:size(param,1)
           figure
           scatter(r, param(j,:));
           title(names{j});
       end

       
       
       
function x = rand_pt(x_min, x_max)
x = rand(size(x_min)) .* (x_max - x_min) + x_min;
