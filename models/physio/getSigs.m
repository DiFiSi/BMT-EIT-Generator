function [sigCard,sigResp] = getSigs(phi,phiCard,phiResp,tempCard,tempResp)
    sigCard = interp1(phi,tempCard,phiCard,'spline');
    sigResp = interp1(phi,tempResp,phiResp,'spline');
end

