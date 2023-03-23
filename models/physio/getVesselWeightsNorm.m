function [wLarge,wMed,wSmall,lags,lagsN] = getVesselWeightsNorm(lags,lagsN,n)
    wX = (lags - min(lags)) ./ (max(lags) - min(lags));

    % Large
    wLarge = normpdf(wX,0,0.25);
    wLarge = n(1) * wLarge ./ sum(wLarge,'omitnan');

    % Medium
    wMed = normpdf(wX,0.5,0.15);
    wMed = n(2) * wMed ./ sum(wMed,'omitnan');

    % Small
    wSmall = normpdf(wX,1,0.25);
    wSmall = n(3) * wSmall ./ sum(wSmall,'omitnan');
    
    % Get N/Nelems
    wLarge = wLarge ./ lagsN;
    wMed = wMed ./ lagsN;
    wSmall = wSmall ./ lagsN;
    
    lags(isnan(lags)) = 0;
    lagsN(isnan(lagsN)) = 0;
    wLarge(isnan(wLarge)) = 0;
    wMed(isnan(wMed)) = 0;
    wSmall(isnan(wSmall)) = 0;
end