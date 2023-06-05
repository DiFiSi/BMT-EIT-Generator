function [V,Q] = mkTempVess(Vmax,Vbase,tp)
    n = (1:100)';
    
    switch(tp)
        case "large"
            data = readmatrix("QLarge.txt");
        case "med"
            data = readmatrix("QMed.txt");
        case "small"
            data = readmatrix("QSmall.txt");
    end

    tOrig = normalize(1:size(data,1),'range') * 99 + 1;
    QOrig = data(:,1); QOrig = normalize(QOrig,'range');

    % Uniformize data
    Q = interp1(tOrig,QOrig,n);

    % Integrate into Volume
    V = normalize(cumtrapz(n,Q - mean(Q)) / length(n),'range'); 
    V = V * (1 - Vbase) * Vmax + Vmax * Vbase;

    % Reobtain Q
    Q = diff(V);
    Q = [Q(1); Q] - min(Q);

    % Plotting
    figure;
    yyaxis left;
    plot(n,Q);
    ylabel("Flow [cm3/s]")
    yyaxis right;
    plot(n,V);
    ylabel("Volume [cm3]")
end

