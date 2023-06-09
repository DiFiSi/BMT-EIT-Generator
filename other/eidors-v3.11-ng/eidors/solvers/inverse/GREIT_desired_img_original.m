function PSF= GREIT_desired_img_original(xyc, radius, opt)
%GREIT_DESIRED_IMG_ORIGINAL The original desired solution for GREIT
% PSF= GREIT_desired_img_original(xyc, radius, opt)
%   xyc     - array of point centers [2xNpoints]
%   radius  - the radius of the target on the desired image as a fraction
%             of the model radius (half the larger dimension in xy)
%   opt     - a struct with these mandatory fields:
%      .rec_model   a 2D model as generated by MK_GRID_MODEL, with some
%                   restrictions, see below. 
%
% The supplied opt.rec_model is only used to figure out the requested image
% size in pixels, and the coordinates of the model in space. If pixels are 
% removed from the model, a boolean opt.rec_model.inside must be specified
% indicating the removed pixels with false, as done by MK_PIXEL_SLICE.
% Nodes must not be removed. Only equal-size pixels are supported.
%
% NOTE that the amount of blur provided by this function decreases as the
% resolution of the image increases. 
%
% As of 2015-03-29, the default desired image function used by
% MK_GRID_MODEL is MK_DESIRED_IMAGE_SIGMOID.
%
% See also: CALC_GREIT_RM, MK_GREIT_MODEL, MK_PIXEL_SLICE, 
%           GREIT_DESIRED_IMG_SIGMOID

% (C) 2009-2015 Andy Adler and Bartlomiej Grychtol. 
% License: GPL v2 or v3
% $Id: GREIT_desired_img_original.m 4986 2015-05-11 20:09:28Z aadler $

   opt = parse_opt(opt.rec_model);
   
   copt.fstr = 'GREIT_desired_img_original';
   copt.cache_obj = {xyc, radius, opt};
   PSF = eidors_cache(@desired_soln,{xyc, radius, opt},copt);
end

function PSF = desired_soln(xyc,radius,opt)
   
   xsz = opt.imgsz(1); ysz = opt.imgsz(2);
   sz= xsz * ysz;
   xmin = opt.meshsz(1); xmax = opt.meshsz(2);
   ymin = opt.meshsz(3); ymax = opt.meshsz(4);
   % scale radius to half the greater dimension
   radius = radius * 0.5 * max(xmax-xmin, ymax-ymin);
   [x,y]= ndgrid(linspace(xmin,xmax,xsz), linspace(ymin,ymax,ysz));
   x_spc = (xmax-xmin)/(xsz-1) * 0.5;
   y_spc = (ymax-ymin)/(ysz-1) * 0.5;
   PSF = zeros(sz,size(xyc,2));
   for i=1:size(xyc,2);
      for dx = linspace(-x_spc, x_spc, 5)
         for dy = linspace(-y_spc, y_spc, 5)
            PSF(:,i) = PSF(:,i) +  1/25*( ...
               (dx+x(:)-xyc(1,i)).^2 + (dy+y(:)-xyc(2,i)).^2 ...
                        < radius^2 );
         end
      end
%     PSF(:,i) = PSF(:,i)/sum(PSF(:,i));
   end
   if ~isempty(opt.inside)
       PSF = PSF(opt.inside,:);
   end
   
end

function opt = parse_opt(mdl)
    opt.meshsz = [];    
    for i = 1:2
        opt.meshsz = [opt.meshsz min(mdl.nodes(:,i)) max(mdl.nodes(:,i))];
    end
    opt.imgsz(1) = numel(unique(mdl.nodes(:,1))) - 1;
    opt.imgsz(2) = numel(unique(mdl.nodes(:,2))) - 1;
%     opt.inside = [];
%     try
        opt.inside = mdl.inside;
%     end
end

