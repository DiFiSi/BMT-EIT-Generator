function [nLargeElemFull,rhoLargeFull,nMedElemFull,rhoMedFull,nSmallElemFull,rhoSmallFull] = getVesselWeightsNorm(dists,nLungElems,nVessels)
    emptyIdxs = nLungElems == 0 | isnan(nLungElems) | isnan(dists);
    fullDists = dists;
    fullDists(emptyIdxs) = [];
    fullNLungElems = nLungElems;
    fullNLungElems(emptyIdxs) = [];
    wX = (fullDists - min(fullDists)) ./ (max(fullDists) - min(fullDists));

    % Large
    rhoLarge = normpdf(wX,0,0.25);

    % Medium
    rhoMed = normpdf(wX,0.5,0.15);

    % Small
    rhoSmall = normpdf(wX,1,0.25);
    
    % Relative density of vessels
    rhoSum = sum([rhoLarge,rhoMed,rhoSmall],2);
    rhoLarge = rhoLarge ./ rhoSum;
    rhoMed = rhoMed ./ rhoSum;
    rhoSmall = rhoSmall ./ rhoSum;
    
    rhoLargeFull = zeros(size(dists));
    rhoLargeFull(~emptyIdxs) = rhoLarge;
    rhoMedFull = zeros(size(dists));
    rhoMedFull(~emptyIdxs) = rhoMed;
    rhoSmallFull = zeros(size(dists));
    rhoSmallFull(~emptyIdxs) = rhoSmall;
    
%     figure;
%     plot(wX,[rhoLarge,rhoMed,rhoSmall],'.','MarkerSize',10);
%     xlabel("Normalized distance from bifurcation [1]");
%     ylabel("Number of vessels per distance [1]");
%     legend(["Large","Medium","Small"]);
    
    % Get N/Nelems
    nLarge = nVessels(1) * (rhoLarge ./ sum(rhoLarge,'omitnan'));
    nMed = nVessels(2) * (rhoMed ./ sum(rhoMed,'omitnan'));
    nSmall = nVessels(3) * (rhoSmall ./ sum(rhoSmall,'omitnan'));
    
    nLargeElem = nLarge ./ fullNLungElems;
    nMedElem = nMed ./ fullNLungElems;
    nSmallElem = nSmall ./ fullNLungElems;
    
    nLargeElemFull = zeros(size(dists));
    nLargeElemFull(~emptyIdxs) = nLargeElem;
    nMedElemFull = zeros(size(dists));
    nMedElemFull(~emptyIdxs) = nMedElem;
    nSmallElemFull = zeros(size(dists));
    nSmallElemFull(~emptyIdxs) = nSmallElem;
    
%     figure;
%     plot(wX,[rhoLarge,rhoMed,rhoSmall],'.','MarkerSize',10);
%     xlabel("Normalized distance from bifurcation [1]");
%     ylabel("Number of vessels per distance / per element [1]");
%     legend(["Large","Medium","Small"]);
end