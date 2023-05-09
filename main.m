%% Prepare workspace, load models and load EIDORS
cd(fileparts(matlab.desktop.editor.getActiveFilename));
addpath(genpath(".\other"));
addpath(genpath(".\utils"));
addpath(genpath(".\models"));
resultsFolder = ".\results\";

run .\other\eidors-v3.11-ng\eidors\eidors_startup.m

%% Load FMDL and IMDL models
load('.\models\fem\ellip_fine.mat');
load('.\models\fem\ellip_fine_imdl.mat');

%% Load templates
% Template folder
tempFolder = "models\synth\templates\";

% Ventricular volume
vVentrTemp = load(fullfile(tempFolder,"heart\Vventricle.mat"));
vVentrTemp = vVentrTemp.y;

% Blood flow and volume
% Large
data = load(fullfile(tempFolder,"\vessels\Qlarge.mat"));
qLargeTemp = data.Q;
data = load(fullfile(tempFolder,"\vessels\Vlarge.mat"));
vLargeTemp = data.V;

% Medium
data = load(fullfile(tempFolder,"\vessels\Qmed.mat"));
qMedTemp = data.Q;
data = load(fullfile(tempFolder,"\vessels\Vmed.mat"));
vMedTemp = data.V;

% Small
data = load(fullfile(tempFolder,"\vessels\Qsmall.mat"));
qSmallTemp = data.Q;
data = load(fullfile(tempFolder,"\vessels\Vsmall.mat"));
vSmallTemp = data.V;

% Air volume
spontVentTemp = load(fullfile(tempFolder,"\ventilation\VawMand.mat"));
spontVentTemp = spontVentTemp.y;
mandVentTemp = load(fullfile(tempFolder,"\ventilation\VawSpont.mat"));
mandVentTemp = mandVentTemp.y;

%% Physiological constants
% System
fs     = 50;            % Simulated Sampling Frequency [Hz]
volElem = 0.0284; % cm^3
nElems = size(fmdl.elems,1);

% Heart
Vmax = 360;
Vs = 140;
Vedv = 200;

% Pulmonary Ventilation (mL)
Vfrc = 3000;
Vt = 500;
vwv = 1000;       % (instantaneous) cm/s

% Pulmonary perfusion
Vpul = 350;
pulFrac = 0.14;         % 14% of pumped blood into lungs
htc = 0.43;             % fraction
pwv = 2.00 * 100;       % cm/s

% Conductivities (S m-1 @ 100 kHz)
compOffset = 0.1; % offset of conductivity for computational stability
                  % (does not change tdEIT results)
sigmaBlood = 0.703;
sigmaHeart = 0.215;
sigmaLungDef = 0.270;
sigmaLungInf = 0.110;
sigmaFat = 0.0244;
sigmaMuscle = 0.362;
sigmaBone = 0.0211;
sigmaBkg = mean([sigmaFat,sigmaMuscle,sigmaBone]);

% Cardiorespiratory coupling
fHR    = [60; 60; 100; 150; 80; 50; 50];            % Mean Heart Rate                      [BPM]
fHF    = [15; 20; 50;  60;  40; 30; 20];            % Mean Respiratory Frequency           [BPM]
fLF    = 6;             % Mean Frequency of Mayer Oscillations [BPM]

kHF    = 0.001;         % Modulation Constant of Respiratory Oscillator
kLF    = 0.0005;        % Modulation Constant of Mayer Oscillator

sdHF   = 0.008;         % Standard Deviation of Respiratory 'Jitter' Component
sdLF   = 0.008;         % Standard Deviation of Mayer 'Jitter' Component

tSim   = 30;            % Simulation Time              [s]

%% Anatomical constants (cm)
% Heart
centChbr = [1.5,4,20]';     % heart chamber center
centMyoc = [1.5,4,20]';     % myocardium center
ratioMyoc = [0.7297,1.75]'; % ratio of ellipsoidal axes of myocardium shape
ratioChbr = [0.7297,1.75]'; % ratio of ellipsoidal axes of heart chamber shape

% Lungs
centBifur = [0,0,15]';  % bifrucation center
centLlung = [7,0,5]';   % left lung center
centRlung = [-7,0,5]';  % right lung center
ratioLung = [1.5,3.5]'; % ratio of ellipsoidal axes of lung shapes

Llarge = 0.1 * (12.35 + 18.07 + 25.97 + 35.69 + 25.30) / 5; % abg. length of large vessels
Lmed = 0.1 * (1.08 + 1.92 + 2.81 + 3.73 + 6.58) / 5;        % avg. length of medium vessels
Lsmall = 0.1 * (0.22 + 0.26 + 0.36 + 0.45 + 0.68) / 5;      % avg. length of small vessels
nLarge = 1 * (2 + 5 + 16 + 26 + 84);      % number of large vessels in 2 lungs
nMed = 1 * (142 + 182 + 227 + 315 + 315); % number of medium vessels in 2 lungs
nSmall = 1 * (180 + 50 + 48 + 81 + 113);  % number of small vessels in 2 lungs

%% Define geometric shape and volume functions
% Shape
sphFun = @(cent,r) strcat('(x-', num2str(cent(1)),').^2+(y-', num2str(cent(2)),').^2+(z-', num2str(cent(3)),').^ 2<', num2str(r),'.^2');
cylFun = @(cent,r) strcat('(x-', num2str(cent(1)),').^2+(y-', num2str(cent(2)),').^2<', num2str(r),'.^2');
eliFun = @(cent,r) strcat('(x-', num2str(cent(1)),').^2/',num2str(r(1)),'.^2+', ...
                          '(y-', num2str(cent(2)),').^2/',num2str(r(2)),'.^2+', ...
                          '(z-', num2str(cent(3)),').^2/',num2str(r(3)),'.^2<1');
                      
uppFun = @(bound) strcat('z<=',num2str(bound));
botFun = @(bound) strcat('z>=',num2str(bound));

% Volume
cylVol = @(r) pi * r(1) ^ 2 * r(2);
eliVol = @(r) 4 / 3 * pi * r(1) * r(2) * r(3);
sphVol = @(r) 4 / 3 * pi * r(1) ^ 3;

% Radius a, based on V, and remaining shape parameters
cylR = @(vol,ratio) (vol / (pi * ratio(1))) .^ (1 / 3);
eliR = @(vol,ratio) (3 * vol / (4 * pi * ratio(1) * ratio(2))) .^ (1 / 3);
sphR = @(vol,ratio) (3 * vol / (4 * pi)) .^ (1 / 3); 

%% Define HR and RR
% Couple frequencies
[phiCard,HR,phiResp,RR,phi,t] = synthesizer(fHR,fHF,fLF,kHF,kLF,sdHF,sdLF,tSim,fs);

% Interpolate templates according to phase signals
% sigCard0 (in-phase):
%   1 - vVentr
% sigCard90 (out-of-phase):
%   1 - vLarge
%   2 - qLarge
%   3 - vMed
%   4 - qMed
%   5 - vSmall
%   6 - qSmall
% sigResp:
%   1 - spontVent
%   2 - mandVent

[sigCard0,sigCard90,sigResp] = getSigs(phi,phiCard,phiResp,vVentrTemp,...
                                                           [vLargeTemp,qLargeTemp,vMedTemp,qMedTemp,vSmallTemp,qSmallTemp],...
                                                           [spontVentTemp,mandVentTemp]);

nSamples = length(sigCard0);                                                    
n = (1:nSamples)';
                                                    
VVentr = sigCard0;

VLarge = sigCard90(:,1);
QLarge = sigCard90(:,2);
QLarge(QLarge <= 0) = 0;

VMed = sigCard90(:,3);
QMed = sigCard90(:,4);
QMed(QMed <= 0) = 0;

VSmall = sigCard90(:,5);
QSmall = sigCard90(:,6);
QSmall(QSmall <= 0) = 0;

VVent = sigResp(:,1); % 1 - spontaenous; 2 - mandatory 

%% Define sinusoidal functions for movement
Vmyoc = Vmax - 1 / 2 * (Vs * VVentr + Vs);
rmyoc = eliR(Vmyoc,ratioMyoc);
rMyoc = rmyoc .* [1, ratioMyoc'];

Vchbr = Vedv - 1 / 2 * (Vs * VVentr + Vs);
rchbr = eliR(Vchbr,ratioChbr);
rChbr = rchbr .* [1, ratioChbr'];

Vvent = Vfrc + 1 / 2 * (Vt * VVent + Vt);
Vlung = Vvent + 0.25 * VLarge; % 0.25 to mimic dampening of the arterial pulsation
rlung = eliR(Vlung,ratioLung);
rLung = rlung .* [1, ratioLung'];
rLungMax = max(rLung)';

%% Define sinusoidal functions for conductivity
sigmaLungVent = (sigmaLungDef + sigmaLungInf) / 2 - (sigmaLungDef - sigmaLungInf) / 2 * normalize(Vlung,'range');
sigmaVent = (sigmaLungDef + sigmaLungInf) / 2 - (sigmaLungDef - sigmaLungInf) / 2 * normalize(Vvent,'range');

%% Initialize EIDORS model
imgd = mk_image(fmdl,1);
[stim,msel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
imgd.fwd_model.stimulation = stim;

% Hard bounds z-axis
uppBound = inline(uppFun(centChbr(3)),'x','y','z');
upp = logical(elem_select(fmdl, uppBound));
botBound = inline(botFun(centLlung(3)),'x','y','z');
bot = logical(elem_select(fmdl, botBound)); 

%% Debug vars
refHeart = [1.5,4,18];
refLungEarly = [7,0,12.5];
refLungLate = [10,0,12.5];

[xFMDL,yFMDL,zFMDL] = getCoords(fmdl);
gridX = mean(xFMDL,2);
gridY = mean(yFMDL,2);
gridZ = mean(zFMDL,2);

[~,heartIdx] = min(sqrt((refHeart(1) - gridX) .^ 2 + (refHeart(2) - gridY) .^ 2 + (refHeart(3) - gridZ) .^ 2));
[~,lungEarlyIdx] = min(sqrt((refLungEarly(1) - gridX) .^ 2 + (refLungEarly(2) - gridY) .^ 2 + (refLungEarly(3) - gridZ) .^ 2));
[~,lungLateIdx] = min(sqrt((refLungLate(1) - gridX) .^ 2 + (refLungLate(2) - gridY) .^ 2 + (refLungLate(3) - gridZ) .^ 2));

VbloodSig = [];
VairSig = [];
sigmaBloodSig = [];
sigmaLungVentSig = [];
fiSig = [];
NelemsSig = [];
sigmaLungSig = [];
sigmaLargeSig = [];
sigmaMedSig = [];
sigmaSmallSig = [];
rhoLargeSig = [];
rhoMedSig = [];
rhoSmallSig = [];
velLargeSig = [];
velMedSig = [];
velSmallSig = [];
radLargeSig = [];
radMedSig = [];
radSmallSig = [];

% Build signal with spontaenous ventilation
tSize = length(t);
delOffset = round((20 / min(pwv,vwv)) * fs);

vMIX = nan(208,tSize - delOffset);
vCRS = nan(208,tSize - delOffset);
vVRS = nan(208,tSize - delOffset);
for i = delOffset + 1:tSize
    tNow = t(i);
    
    selectMyoc = inline(eliFun(centMyoc,rMyoc(i,:)),'x','y','z');
    selectChbr = inline(eliFun(centChbr,rChbr(i,:)),'x','y','z');
    
    rlungDelta = (rLungMax - rLung(i,:)') / 2;
    centL = centLlung + rlungDelta;
    centR = centRlung + rlungDelta;
    selectLlung = inline(eliFun(centL,rLung(i,:)),'x','y','z');
    selectRlung = inline(eliFun(centR,rLung(i,:)),'x','y','z');
    
    % Get indices of all organs and background
    myocIdxs = elem_select(imgd.fwd_model, selectMyoc) & upp;
    chbrIdxs = elem_select(imgd.fwd_model, selectChbr) & upp;
    llungIdxs = (elem_select(imgd.fwd_model, selectLlung) & bot) & ~myocIdxs;
    rlungIdxs = (elem_select(imgd.fwd_model, selectRlung) & bot) & ~myocIdxs;
    bkgIdxs = ~(myocIdxs | chbrIdxs | llungIdxs | rlungIdxs);
    
    % Get delays and vessel weights
    [dists,nLungElems,lagsCard,lagsResp] = getDelays(fmdl,[pwv,vwv],centBifur,llungIdxs,rlungIdxs,fs);
    [rhoLarge,rhoMed,rhoSmall] = getVesselWeightsNorm(dists,nLungElems,[nLarge;nMed;nSmall]);
    
    % Calculate current volume and flow states of vessels
    delSigs = delaySigs(i,[VLarge,QLarge,VMed,QMed,VSmall,QSmall],lagsCard);
    VLargeNow = delSigs(:,1);
    QLargeNow = delSigs(:,2);
    VMedNow = delSigs(:,3);
    QMedNow = delSigs(:,4);
    VSmallNow = delSigs(:,5);
    QSmallNow = delSigs(:,6);
    
    % Get delayed ventilation signal
    delSigs = delaySigs(i,[Vlung,sigmaLungVent],lagsResp);
    VlungNow = delSigs(:,1);
    sigmaLungVentNow = delSigs(:,2);
    
    % Calculate current conductivity state of blood in vessels
    [radLarge,velLarge] = radvel(VLargeNow,QLargeNow,Llarge);
    sigmaLarge = sigmaBlood + visser(velLarge,radLarge,htc);
    [radMed,velMed] = radvel(VMedNow,QMedNow,Lmed);
    sigmaMed = sigmaBlood + visser(velMed,radMed,htc); 
    [radSmall,velSmall] = radvel(VSmallNow,QSmallNow,Lsmall);
    sigmaSmall = sigmaBlood + visser(velSmall,radSmall,htc); 
    sigmaBloodNow = rhoLarge    .* sigmaLarge   + ...
                    rhoMed      .* sigmaMed     + ...
                    rhoSmall    .* sigmaSmall;           
    
    % Calculate fraction of blood and air volume
    VbloodNow = rhoLarge   .* VLargeNow    + ...
                rhoMed     .* VMedNow      + ...
                rhoSmall   .* VSmallNow;
    VairNow = 2 * VlungNow / mean([sum(llungIdxs),sum(rlungIdxs)]);
    fi = VbloodNow ./ (VbloodNow + VairNow);
    
    % Use mixing formula to calculate conductivity of FEM element
    sigmaLungNow = maxgar(sigmaLungVentNow,sigmaBloodNow,fi);
    
    VbloodSig = [VbloodSig; [VbloodNow(heartIdx),VbloodNow(lungEarlyIdx),VbloodNow(lungLateIdx)]];
    VairSig = [VairSig; [VairNow(heartIdx),VairNow(lungEarlyIdx),VairNow(lungLateIdx)]];
    fiSig = [fiSig; [fi(heartIdx),fi(lungEarlyIdx),fi(lungLateIdx)]];
    NelemsSig = [NelemsSig; [sum(llungIdxs),sum(rlungIdxs),sum(myocIdxs),sum(chbrIdxs)]];
    
    sigmaBloodSig = [sigmaBloodSig; [sigmaBloodNow(heartIdx),sigmaBloodNow(lungEarlyIdx),sigmaBloodNow(lungLateIdx)]];
    sigmaLungSig = [sigmaLungSig; [sigmaLungNow(heartIdx),sigmaLungNow(lungEarlyIdx),sigmaLungNow(lungLateIdx)]];
    sigmaLungVentSig = [sigmaLungVentSig; [sigmaLungVent(i),sigmaLungVent(i),sigmaLungVent(i)]];

    sigmaLargeSig = [sigmaLargeSig; [sigmaLarge(heartIdx),sigmaLarge(lungEarlyIdx),sigmaLarge(lungLateIdx)]];
    sigmaMedSig = [sigmaMedSig; [sigmaMed(heartIdx),sigmaMed(lungEarlyIdx),sigmaMed(lungLateIdx)]];
    sigmaSmallSig = [sigmaSmallSig; [sigmaSmall(heartIdx),sigmaSmall(lungEarlyIdx),sigmaSmall(lungLateIdx)]];
    
    rhoLargeSig = [rhoLargeSig; [rhoLarge(heartIdx),rhoLarge(lungEarlyIdx),rhoLarge(lungLateIdx)]];
    rhoMedSig = [rhoMedSig; [rhoMed(heartIdx),rhoMed(lungEarlyIdx),rhoMed(lungLateIdx)]];
    rhoSmallSig = [rhoSmallSig; [rhoSmall(heartIdx),rhoSmall(lungEarlyIdx),rhoSmall(lungLateIdx)]];
    
    velLargeSig = [velLargeSig; [velLarge(heartIdx),velLarge(lungEarlyIdx),velLarge(lungLateIdx)]];
    velMedSig = [velMedSig; [velMed(heartIdx),velMed(lungEarlyIdx),velMed(lungLateIdx)]];
    velSmallSig = [velSmallSig; [velSmall(heartIdx),velSmall(lungEarlyIdx),velSmall(lungLateIdx)]];
    
    radLargeSig = [radLargeSig; [radLarge(heartIdx),radLarge(lungEarlyIdx),radLarge(lungLateIdx)]];
    radMedSig = [radMedSig; [radMed(heartIdx),radMed(lungEarlyIdx),radMed(lungLateIdx)]];
    radSmallSig = [radSmallSig; [radSmall(heartIdx),radSmall(lungEarlyIdx),radSmall(lungLateIdx)]];
    
%     % MIX
%     imgd.elem_data = compOffset     + ...
%                      + sigmaBkg       * bkgIdxs...
%                      + sigmaLungNow   .* llungIdxs...
%                      + sigmaLungNow   .* rlungIdxs...
%                      + sigmaBlood     * chbrIdxs...
%                      + sigmaMuscle    * myocIdxs;
%                           
%     % Get voltages
%     vi = fwd_solve(imgd);
%     vMIX(:,i - delOffset) = vi.meas;
%     
%     % VRS
%     sigmaLungNow = sigmaVent(i);
%     imgd.elem_data = compOffset     + ...
%                      + sigmaBkg       * (bkgIdxs | chbrIdxs | myocIdxs)...
%                      + sigmaLungNow   .* llungIdxs...
%                      + sigmaLungNow   .* rlungIdxs;
%                
%     % Get voltages
%     vi = fwd_solve(imgd);
%     vVRS(:,i - delOffset) = vi.meas;
%     
%     % CRS
%     % Assume no lung movement
%     rlungDelta = 0;
%     centL = centLlung + rlungDelta;
%     centR = centRlung + rlungDelta;
%     selectLlung = inline(eliFun(centL,rLungMax),'x','y','z');
%     selectRlung = inline(eliFun(centR,rLungMax),'x','y','z');
%     
%     % Get new indices of all organs and background
%     llungIdxs = (elem_select(imgd.fwd_model, selectLlung) & bot) & ~myocIdxs;
%     rlungIdxs = (elem_select(imgd.fwd_model, selectRlung) & bot) & ~myocIdxs;
%     bkgIdxs = ~(llungIdxs | rlungIdxs); % myocIdxs | chbrIdxs | 
%     
%     % No air volume in mixing formula
%     sigmaLungNow = maxgar(sigmaLungVent(i),sigmaBloodNow,1);
%     
%     imgd.elem_data = compOffset     + ...
%                    + sigmaBkg       * bkgIdxs...
%                    + sigmaLungNow   .* llungIdxs...
%                    + sigmaLungNow   .* rlungIdxs ...
%                    + sigmaBlood     * chbrIdxs...
%                    + sigmaMuscle    * myocIdxs;
%                
%     % Get voltages
%     vi = fwd_solve(imgd);
%     vCRS(:,i - delOffset) = vi.meas;
    
    disp(["Sample " + num2str(i) + "/" + num2str(tSize)]);
end

%% Relevant debug graphs
figure;
plot([VairSig(:,3),VbloodSig(:,2)]);
hold on; plot(Vmyoc);
legend(["Air","Blood","Myocardium"])
title("Air and blood volumes in single lung element");

figure;
plot([VairSig(:,3),rhoMedSig(:,2)]);
legend(["VAir","$\rho_{med}$"])
title("Air volume and number of medium sized lung vessels");

figure;
plot([rhoLargeSig(:,2),rhoMedSig(:,2),rhoSmallSig(:,2)]);
legend(["$\rho_{large}$","$\rho_{med}$","$\rho_{small}$"])
title("Number of differently-sized lung vessels in element");

figure;
plot([sigmaMedSig(:,2) - sigmaBlood,velMedSig(:,2),radSmallSig(:,2)]);
legend(["$\Delta\sigma_{med}$","$v_{med}$","$r_{med}$"])
title("Visser cond. change with blood velocity and radius");

figure;
plot(sigmaLungSig(:,2)); hold on;
plot(VairSig(:,3));
plot(sigmaBloodSig(:,2));
legend(["$\Delta\sigma_{lung}$","$\Delta\sigma_{vent}$","$\Delta\sigma_{blood}$"])
title("Visser cond. change with blood velocity and radius");

%% Reconstruct and save
% Remove NaNs
vMIXclean = rmmissing(real(vMIX),2);
vVRSclean = rmmissing(real(vVRS),2);
vCRSclean = rmmissing(real(vCRS),2);

% Reconstruct MIX
vh = mean(vMIXclean,2);
imgr = inv_solve(imdl, vh, vMIXclean);
imgMIX = calc_slices(imgr);

% Reconstruct VRS
vh = mean(vVRSclean,2);
imgr = inv_solve(imdl, vh, vVRSclean);
imgVRS = calc_slices(imgr);

% Reconstruct CRS
vh = mean(vCRSclean,2);
imgr = inv_solve(imdl, vh, vCRSclean);
imgCRS = calc_slices(imgr);

% Save
save(fullfile(resultsFolder,"abschlussbericht.mat"));

%% Test signals
heartPix = [24,35];
lungPix = [30,49];

% MIX
heartSig = squeeze(imgMIX(heartPix(1),heartPix(2),:));
lungSig = squeeze(imgMIX(lungPix(1),lungPix(2),:));
figure; plot([heartSig,lungSig])
title("MIX")

% VRS
heartSig = squeeze(imgVRS(heartPix(1),heartPix(2),:));
lungSig = squeeze(imgVRS(lungPix(1),lungPix(2),:));
figure; plot([heartSig,lungSig])
title("VRS")

% CRS
heartSig = squeeze(imgCRS(heartPix(1),heartPix(2),:));
lungSig = squeeze(imgCRS(lungPix(1),lungPix(2),:));
figure; plot([heartSig,lungSig])
title("CRS")
