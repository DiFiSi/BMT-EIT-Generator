function [phiCard,fHR,phiResp,fHF,phi,t] = synthesizer(fHR,fHF,fLF,kHF,kLF,sdHF,sdLF,tSim,fs)
    % synthesizer_demo.m          19/06/2017
    % written by C. Hoog Antink
    % Copyright: Chair for Medical Information Technology
    %            Pauwelsstr. 20
    %            52074 Aachen
    %
    % This code demonstrates the synthesizer described in
    %
    % C. Hoog Antink, S. Leonhardt, and M. Walter: “A Synthesizer Framework
    % for Multimodal Cardiorespiratory Signals”. Biomedical Engineering
    % and Physics Express, June 2017
    % See in particular Figure 14.
    % We kindly ask you to cite the original publication when using this code.

    t = (0:1 / fs:tSim)';
    tFinal  = tSim * fs;

    [bBaro, aBaro] = butter(4, max(fLF) / 60 / (fs / 2));
    [bRR, aRR]     = butter(4, max(fHF) / 60 / (fs / 2));

    options = odeset('RelTol',1e-5,'AbsTol',1e-5);

    rng(1)
    randParamLF = filtfilt(bBaro,aBaro,randn(tFinal+1,1));
    randParamHF = filtfilt(bRR,aRR,randn(tFinal+1,1));

    if length(fHR) > 1 && length(fHF) > 1
        % Interpolate frequency arrays to desired simulation time
        fHR = interp1(linspace(0,tFinal+1,length(fHR)),fHR,(1:(tFinal+1)),'pchip');
        fHF = interp1(linspace(0,tFinal+1,length(fHF)),fHF,(1:(tFinal+1)),'pchip');

        % Remove transients introduced by interpolation
        [bTransient, aTransient] = butter(4, 0.4 / (fs / 2), 'low');
        fHR = filtfilt(bTransient,aTransient,fHR);
        fHF = filtfilt(bTransient,aTransient,fHF);
    else
        fHR = repmat(fHR,tFinal+1,1);
        fHF = repmat(fHF,tFinal+1,1);
    end

    [T,Y] = ode45(@(t,y)rigid(t,y,fs,randParamLF,randParamHF,fHR,fHF,fLF,sdLF,sdHF,kLF,kHF),[0 tFinal],[1 0 1 0 1 0],options);

    delta = 1;

    interpEqui = (0:delta:T(end))';

    sig2aEqui = interp1(T,Y(:,3),interpEqui);
    sig2bEqui = interp1(T,Y(:,4),interpEqui);
    sig3aEqui = interp1(T,Y(:,5),interpEqui);
    sig3bEqui = interp1(T,Y(:,6),interpEqui);

    phiResp = atan2(sig2bEqui,sig2aEqui);
    phiCard = atan2(sig3bEqui,sig3aEqui);
    phi = ((1:100) - 1) / 100 * 2 * pi - pi;
end

function dy = rigid(t,y,fs,randParam1,randParam2,fHR,fHF,fLF,sdLF,sdHF,kLF,kHF)
    m1   = randParam1(round(t)+1);
    m2   = randParam2(round(t)+1);
    fHR = fHR(round(t)+1);
    fHF = fHF(round(t)+1);

    dy = zeros(2,1);
    dy(1) = -2*pi*(fLF / 60 / fs + sdLF*m1)*y(2); % mayer oscillator 
    dy(2) =  2*pi*(fLF / 60 / fs + sdLF*m1)*y(1); % mayer oscillator

    dy(3) = -2*pi*(fHF / 60 / fs + sdHF*m2)*y(4); % lung oscillator
    dy(4) =  2*pi*(fHF / 60 / fs + sdHF*m2)*y(3); % lung oscillator

    dy(5) = -2*pi*(fHR / 60 / fs + kLF*y(1)+kHF*y(3))*y(6); % heart oscillator
    dy(6) =  2*pi*(fHR / 60 / fs + kLF*y(1)+kHF*y(3))*y(5); % heart oscillator
end