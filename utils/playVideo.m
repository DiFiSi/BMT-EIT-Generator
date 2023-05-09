function h = playVideo(vid, fr, colormap)
    % vid - a 3D matrix (2D over time)
    % fr - framerate
    % colormap - colormap string, e.g. 'autumn' or 'grayscale' or 'jet'
    % rangeMin, rangeMax - lower and higher limits of pixel values in
    % matrix
    
    % make sure data is double
    vid = double(vid);
    
    if ~exist('fr','var')
       fr = 30;
    end
    
    if ~exist('colormap','var')
       colormap = jet(256);
    end
    
    h = implay(vid,fr);  
    
    h.Visual.ColorMap.UserRange = 1;
    h.Visual.ColorMap.UserRangeMin = min(vid(:),[],'omitnan');
    h.Visual.ColorMap.UserRangeMax = max(vid(:),[],'omitnan');
    h.Visual.ColorMap.Map = colormap;
end