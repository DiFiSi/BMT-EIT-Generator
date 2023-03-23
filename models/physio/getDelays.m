function [lags,lagsN,dels] = getDelays(fmdl,pwv,origin,staticLlung,staticRlung,fs)
    [x,y,z] = getCoords(fmdl);
    gridX = mean(x,2);
    gridY = mean(y,2);
    gridZ = mean(z,2);

    gridX(~(staticLlung | staticRlung)) = nan;
    gridY(~(staticLlung | staticRlung)) = nan;
    gridZ(~(staticLlung | staticRlung)) = nan;

    dists = sqrt((origin(1) - gridX) .^ 2 + (origin(2) - gridY) .^ 2 + (origin(3) - gridZ) .^ 2);
    dels = dists / pwv; dels = dels - min(dels);
    
    lags = round(dels * fs);
    
    [lagsN,lagBounds] = histcounts(lags);
    lagBounds = lagBounds(1:end-1) + diff(lagBounds) / 2;

    lagsN = interp1(lagBounds,lagsN,lags,'spline');
    
%     lags(isnan(lags)) = 0;
%     dels(isnan(dels)) = 0;
%     lagsN(isnan(lagsN)) = 0;
end