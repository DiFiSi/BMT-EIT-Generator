function [rhoLargeFull,rhoMedFull,rhoSmallFull] = getVesselWeightsNorm(dists,nLungElems,nVessels)
    emptyIdxs = nLungElems == 0 | isnan(nLungElems) | isnan(dists);
    fullDists = dists;
    fullDists(emptyIdxs) = [];
    fullNLungElems = nLungElems;
    fullNLungElems(emptyIdxs) = [];
    wX = (fullDists - min(fullDists)) ./ (max(fullDists) - min(fullDists));

    % Large
    rhoLarge = normpdf(wX,0,0.25);
    rhoLarge = rhoLarge ./ sum(rhoLarge,'omitnan');
    rhoLarge = nVessels(1) * rhoLarge;

    % Medium
    rhoMed = normpdf(wX,0.5,0.15);
    rhoMed = rhoMed ./ sum(rhoMed,'omitnan');
    rhoMed = nVessels(2) * rhoMed;

    % Small
    rhoSmall = normpdf(wX,1,0.25);
    rhoSmall = rhoSmall ./ sum(rhoSmall,'omitnan');
    rhoSmall = nVessels(3) * rhoSmall;
    
%     figure;
%     plot(wX,[rhoLarge,rhoMed,rhoSmall],'.','MarkerSize',10);
%     xlabel("Normalized distance from bifurcation [1]");
%     ylabel("Number of vessels per distance [1]");
%     legend(["Large","Medium","Small"]);
    
    % Get N/Nelems
    rhoLarge = rhoLarge ./ fullNLungElems;
    rhoMed = rhoMed ./ fullNLungElems;
    rhoSmall = rhoSmall ./ fullNLungElems;
    
    rhoLargeFull = zeros(size(dists));
    rhoLargeFull(~emptyIdxs) = rhoLarge;
    rhoMedFull = zeros(size(dists));
    rhoMedFull(~emptyIdxs) = rhoMed;
    rhoSmallFull = zeros(size(dists));
    rhoSmallFull(~emptyIdxs) = rhoSmall;
    
%     figure;
%     plot(wX,[rhoLarge,rhoMed,rhoSmall],'.','MarkerSize',10);
%     xlabel("Normalized distance from bifurcation [1]");
%     ylabel("Number of vessels per distance / per element [1]");
%     legend(["Large","Medium","Small"]);
end