function [sigCard0,sigCard90,sigResp] = getSigs(phi,phiCard,phiResp,tempCard0,tempCard90,tempResp)
    phiCard0 = phiCard;
    phiCard90 = phiCard;
    phiCard90(phiCard0 >= 0) = phiCard0(phiCard0 >= 0) - pi;
    phiCard90(phiCard0 < 0) = phiCard0(phiCard0 < 0) + pi;

    sigCard0 = interp1(phi,tempCard0,phiCard0,'spline');
    sigCard90 = interp1(phi,tempCard90,phiCard90,'spline');
    sigResp = interp1(phi,tempResp,phiResp,'spline');
end