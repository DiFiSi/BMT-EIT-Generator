function [elem_data,ref_lev,max_scale] = scale_for_display( elem_data, pp)
% [elem_data,ref_lev,max_scale] = scale_for_display( elem_data, pp)
%
% PARAMETERS: elem_data
%  elem_data: data for fem elements or image pixels
%
% PARAMETERS: ref_lev
%  ref_lev:   if param provided, use it,
%               otherwise use the global value
%             Can be numeric or 'auto' or 'use_global' 
%
% PARAMETERS: clim
%    clim - colour limit. Colours more different from ref_level are cropped.
%         - if not specified or scale==[] => no limit
%
% OUTPUT: 
%    ref_lev, max_scale - the centre and max of the colour scale
%    elem_data - data scaled in the range [-1 .. 1]
%
% $Id: scale_for_display.m 6148 2021-10-05 11:46:08Z bgrychtol $
% (C) 2006 Andy Adler. Licensed under GPL v2

%FIXME - set to use the colours in the img.calc_colours fields

   global eidors_colours;
   clim = [];
   if nargin <=1
      ref_lev = eidors_colours.ref_level;
      component = 'real';
   elseif ischar(pp)  && strcmp(pp, 'use_global' );
      ref_lev = eidors_colours.ref_level;
      component = 'real';
   else
      ref_lev = pp.ref_level;
      clim    = pp.clim;
      component = pp.component;
   end

   if ~isnumeric(ref_lev)
      if ~strcmp(ref_lev, 'auto')
          error('ref_level parameter must be "auto" or numeric');
      end
      s_ed= elem_data(:);
      s_ed(isnan(s_ed)) = [];
      s_ed= sort(s_ed);
      e= length(s_ed);
      if e==0;
         error('Can''t display. All values NaN. Is raw data 0?')
      end
      % ensure symmetric rejection of data for small data sets
      % This means that 1-2*.35 = .3 of the data will be used
      %  to take the mean. It should make for better centring
      %  of reconstructions - AA may2015
      rej_vals = floor(.35*e);
      ref_lev = mean(s_ed( (rej_vals+1):(end-rej_vals) ));
   end

   elem_data = elem_data - ref_lev;

   switch component;
      case 'real'; elem_data = real(elem_data);
      case 'imag'; elem_data = imag(elem_data);
      otherwise;   error('specified component not real or imag');
   end
   
   max_scale = max(abs(elem_data(:))) + eps;

   % Crop output to the colour limit
   if ~isempty(clim)
      elem_data( elem_data> clim)=  clim;
      elem_data( elem_data<-clim)= -clim;
      max_scale = clim;
   end

   if isfield(eidors_colours,'image_field')
      elem_data(eidors_colours.image_field_idx) = eidors_colours.image_field_val;
   end

   elem_data = elem_data/max_scale;
