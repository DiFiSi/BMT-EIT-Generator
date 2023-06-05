function V = mkTempVentr()
    fs = 100;
    fundFreq = 60; fundFreq = fundFreq / 60;
    n = (1:1 / fundFreq * fs)';
    t = (n - 1) / fs;

    V = readmatrix("Vventricle.txt");
    V = V(:,2);
    tOrig = 1:length(V); tOrig = normalize(tOrig,'range');

    % Normalize
    V = (V - min(V)) / (max(V) - min(V)) * 2 - 1;

    % Interpolate between t(1) and t(end)
    V = -interp1(tOrig,V,t,'spline');

    % Plotting
    figure;
    plot(V);
    ylabel("Volume [cm3]");
end

