function [wLarge,wMed,wSmall,lags,lagsN] = getVesselWeightsExp(lags,lagsN,n)
    % Spatial presence of vessel classes scaled by number of vessels
    tau = 10;
    wX = (lags - min(lags)) ./ (max(lags) - min(lags));

    % Large
    wLargeFun = @(d) exp(-tau * d);
    wLarge = wLargeFun(wX);
    wLarge = n(1) * wLarge ./ sum(wLarge,'omitnan');

    % Medium
    wMed1Fun = @(d) exp((tau * 2) * (d - 0.5));
    wMed2Fun = @(d) exp(-(tau * 2) * (d - 0.5));
    wMed = zeros(size(wX));
    wMed(wX < 0.5) = wMed1Fun(wX(wX < 0.5));...
    wMed(wX >= 0.5) = wMed2Fun(wX(wX >= 0.5));...
    wMed = n(2) * wMed ./ sum(wMed,'omitnan');

    % Small
    wSmallFun = @(d) exp(tau * (d - 1));
    wSmall = wSmallFun(wX);
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

