function sigDel = delaySigs(i, sig, dels)
    sigDel = zeros(length(dels),size(sig,2));
    for d = 1:length(dels)
        sigDel(d,:) = sig(i - dels(d),:);
    end
end