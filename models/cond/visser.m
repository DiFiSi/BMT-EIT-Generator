function deltaSigma = visser(v,R,HCT)
    %VISSER Summary of this function goes here
    %   v - average spatial velocity over time
    %   R - radius over time
    %   HCT - hematocrit
    
    deltaSigma = -0.45 .* HCT .* (1 - exp(-0.26 .* (v ./ R) .^ 0.39));
end

