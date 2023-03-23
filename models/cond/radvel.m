function [r,v] = radvel(vol,flow,L)
    r = sqrt(vol ./ (L * pi));
    v = flow ./ (pi * r .^ 2);
end

