function [dists,nLungElems,lagsCard,lagsResp] = getDelays(fmdl,v,origin,llungIdxs,rlungIdxs,fs)
    pwv = v(1); % pulse wave velocity
    vwv = v(2); % ventilatory wave velocity

    [x,y,z] = getCoords(fmdl);
    gridX = mean(x,2);
    gridY = mean(y,2);
    gridZ = mean(z,2);

    gridX(~(llungIdxs | rlungIdxs)) = nan;
    gridY(~(llungIdxs | rlungIdxs)) = nan;
    gridZ(~(llungIdxs | rlungIdxs)) = nan;

    dists = sqrt((origin(1) - gridX) .^ 2 + (origin(2) - gridY) .^ 2 + (origin(3) - gridZ) .^ 2);
    
    delsCard = dists / pwv; delsCard = delsCard - min(delsCard);
    lagsCard = round(delsCard * fs);
    
    delsResp = dists / vwv; delsResp = delsResp - min(delsResp);
    lagsResp = round(delsResp * fs);
    
    lagsCard(isnan(lagsCard)) = 0;
    lagsResp(isnan(lagsResp)) = 0;
    dists(isnan(dists)) = 0;
    
    [uniqueDists,~,uniqueIdxs] = unique(dists);
    uniqueCounts = accumarray(uniqueIdxs,1);
    uniqueCounts(uniqueDists == 0) = 0;
    
    nLungElems = uniqueCounts(uniqueIdxs);
end