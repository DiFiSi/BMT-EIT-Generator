function [myocIdxs,chbrIdxs,llungIdxs,rlungIdxs,bkgIdxs] = getOrganIdxs(imgd,centMyoc,centChbr,centLlung,centRlung,rMyoc,rChbr,rLung,rLungMax,geoFun,upp,bot)
    selectMyoc = inline(geoFun(centMyoc,rMyoc),'x','y','z');
    selectChbr = inline(geoFun(centChbr,rChbr),'x','y','z');
    
    rlungDelta = (rLungMax - rLung') / 2;
    centL = centLlung + rlungDelta;
    centR = centRlung + rlungDelta;
    selectLlung = inline(geoFun(centL,rLung),'x','y','z');
    selectRlung = inline(geoFun(centR,rLung),'x','y','z');
    
    % Get indices of all organs and background
    myocIdxs = elem_select(imgd.fwd_model, selectMyoc) & upp;
    chbrIdxs = elem_select(imgd.fwd_model, selectChbr) & upp;
    llungIdxs = (elem_select(imgd.fwd_model, selectLlung) & bot) & ~myocIdxs;
    rlungIdxs = (elem_select(imgd.fwd_model, selectRlung) & bot) & ~myocIdxs;
    bkgIdxs = ~(myocIdxs | chbrIdxs | llungIdxs | rlungIdxs);
end