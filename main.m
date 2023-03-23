%% EIDORS
run 'C:\Users\dfpsi\Documents\MedIT\EITGeneratorFull\EITGenerator\other\eidors-v3.11-ng\eidors\startup.m';

%% Load models
load('models\fem\ellip_fine.mat');
load('models\fem\ellip_fine_imdl.mat');

%% Load templates
% Template folder
tempFolder = "models\synth\templates\";

% Ventricular volume
vVentr = load(fullfile(tempFolder,"heart\Vventricle.mat"));
vVentr = vVentr.y;

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
nElems = 829360;

% Heart
Vmax = 360;
Vs = 140;
Vedv = 240;

% Pulmonary Ventilation (mL)
Vfrc = 3000;
Vt = 500;

% Pulmonary perfusion
Vpul = 350;
pulFrac = 0.14;         % 14% of pumped blood into lungs
htc = 0.43;             % fraction
pwv = 2.00 * 100;       % cm/s

% Conductivities (S m-1 @ 100 kHz)
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
centChbr = [2,8,20]';
centMyoc = [2,8,20]';
ratioMyoc = [0.7297,3.1622]';
ratioChbr = [0.7297,3.1622]';

% Lungs
centBifur = [0,0,15]';
centLlung = [10,0,5]';
centRlung = [-10,0,5]';
ratioLung = [1.5,3.5]';
Llarge = 0.1 * (12.35 + 18.07 + 25.97 + 35.69 + 25.30) / 5;
Lmed = 0.1 * (1.08 + 1.92 + 2.81 + 3.73 + 6.58) / 5;
Lsmall = 0.1 * (0.22 + 0.26 + 0.36 + 0.45 + 0.68) / 5;
nLarge = 2 + 5 + 16 + 26 + 84;
nMed = 142 + 182 + 227 + 315 + 315;
nSmall = 180 + 50 + 48 + 81 + 113;

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
% sigCard:
%   1 - vVentr
%   2 - vLarge
%   3 - vMed
%   4 - vSmall
%   5 - vLarge
%   6 - vMed
%   7 - vSmall
% sigResp:
%   1 - spontVent
%   2 - mandVent
[sigCard,sigResp] = getSigs(phi,phiCard,phiResp,[vVentr,vLargeTemp,vMedTemp,vSmallTemp,...
                                                        qLargeTemp,qMedTemp,qSmallTemp],...
                                                        [spontVentTemp,mandVentTemp]);

vVentr = sigCard(:,1);
                                                    
% vLarge = sigCard(:,2);
% vMed = sigCard(:,3);
% vSmall = sigCard(:,4);

VLarge = sigCard(:,2);
VMed = sigCard(:,3);
VSmall = sigCard(:,4);
QLarge = sigCard(:,5);
QMed = sigCard(:,6);
QSmall = sigCard(:,7);

spontVent = sigResp(:,1);
mandVent = sigResp(:,2);

%% Define sinusoidal functions for movement
Vmyoc = Vmax - 1 / 2 * (Vs * vVentr + Vs);
rmyoc = eliR(Vmyoc,ratioMyoc);
rMyoc = rmyoc .* [1, ratioMyoc'];

Vchbr = Vedv - 1 / 2 * (Vs * vVentr + Vs);
rchbr = eliR(Vchbr,ratioChbr);
rChbr = rchbr .* [1, ratioChbr'];

Vvent = Vfrc + 1 / 2 * (Vt * mandVent + Vt);
Vlung = VLarge + Vvent;

rlung = eliR(Vlung,ratioLung);
rLung = rlung .* [1, ratioLung'];
rLungMax = max(rLung)';

%% Define sinusoidal functions for conductivity
sigmaLungVent = (sigmaLungDef + sigmaLungInf) / 2 - (sigmaLungDef - sigmaLungInf) / 2 * mandVent;

%% Initialize EIDORS model
% imgd - dynamic model
imgd = mk_image(fmdl,1);
[stim,msel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'},1);
imgd.fwd_model.stimulation = stim;

% Hard bounds z-axis
uppBound = inline(uppFun(centChbr(3)),'x','y','z');
upp = elem_select(fmdl, uppBound);
botBound = inline(botFun(centLlung(3)),'x','y','z');
bot = elem_select(fmdl, botBound); 

%% Build single cycle
tSize = length(t);
delOffset = round(0.200 * fs); % limit of 200 ms of delay

vMIX = nan(208,tSize - delOffset);
vCRS = nan(208,tSize - delOffset);
vVRS = nan(208,tSize - delOffset);
for i = delOffset + 1:100
    tNow = t(i);
    
    selectMyoc = inline(eliFun(centMyoc,rMyoc(i,:)),'x','y','z');
    selectChbr = inline(eliFun(centChbr,rChbr(i,:)),'x','y','z');
    
    rlungDelta = (rLungMax - rLung(i,:)') / 2;
    centL = centLlung + rlungDelta;
    centR = centRlung + rlungDelta;
    selectLlung = inline(eliFun(centL,rLung(i,:)),'x','y','z');
    selectRlung = inline(eliFun(centR,rLung(i,:)),'x','y','z');
    
    % Get indices of lungs
    myocIdxs = elem_select(imgd.fwd_model, selectMyoc) .* upp;
    chbrIdxs = elem_select(imgd.fwd_model, selectChbr) .* upp;
    llungIdxs = elem_select(imgd.fwd_model, selectLlung) .* bot;
    rlungIdxs = elem_select(imgd.fwd_model, selectRlung) .* bot;
    
    % Get delays and vessel weights
    [lags,lagsN,dels] = getDelays(fmdl,pwv,centBifur,llungIdxs,rlungIdxs,fs);
    [wLarge,wMed,wSmall,lags,lagsN] = getVesselWeightsNorm(lags,lagsN,[nLarge;nMed;nSmall]);
    
    % Calculate current volume and flow states of vessels
    delSigs = delaySigs(i, [VLarge,QLarge,VMed,QMed,VSmall,QSmall],lags);
    VLargeNow = delSigs(:,1);
    QLargeNow = delSigs(:,2);
    VMedNow = delSigs(:,3);
    QMedNow = delSigs(:,4);
    VSmallNow = delSigs(:,5);
    QSmallNow = delSigs(:,6);
    
    % Calculate current conductivity state of blood in vessels
    [rad,vel] = radvel(VLargeNow,QLargeNow,Llarge);
    sigmaLarge = sigmaBlood + visser(vel,rad,htc);
    [rad,vel] = radvel(VMedNow,QMedNow,Lmed);
    sigmaMed = sigmaBlood + visser(vel,rad,htc);
    [rad,vel] = radvel(VSmallNow,QSmallNow,Lsmall);
    sigmaSmall = sigmaBlood + visser(vel,rad,htc);
    sigmaBloodNow = wLarge .* sigmaLarge + ...
                    wMed .* sigmaMed + ...
                    wSmall .* sigmaSmall;
                
    % Calculate fraction of blood and air volume            
    Vblood = wLarge .* VLargeNow + ...
             wMed .* VMedNow + ...
             wSmall .* VSmallNow;
    Vair = Vvent(i) / sum(llungIdxs | rlungIdxs);
    fi = Vblood ./ (Vblood + Vair);
    
    % Use mixing formula to calculate conductivity of FEM element
    sigmaLungNow = maxgar(sigmaLungVent(i),sigmaBloodNow,fi);
    
    % MIX
    imgd.elem_data = sigmaBkg + sigmaLungNow .* llungIdxs...
                              + sigmaLungNow .* rlungIdxs...
                              + sigmaBlood   * chbrIdxs...
                              + sigmaMuscle  * myocIdxs;
                          
    % Get voltages
    vi = fwd_solve(imgd);
    vMIX(:,i - delOffset) = vi.meas;
    
    % VRS
    sigmaLungNow = sigmaLungVent(i);
    imgd.elem_data = sigmaBkg + sigmaLungNow .* llungIdxs...
                              + sigmaLungNow .* rlungIdxs;
    % Get voltages
    vi = fwd_solve(imgd);
    vVRS(:,i - delOffset) = vi.meas;
    
    % CRS
    sigmaLungNow = sigmaBloodNow .* Vblood;
    imgd.elem_data = sigmaBkg + sigmaLungNow .* llungIdxs...
                              + sigmaLungNow .* rlungIdxs...
                              + sigmaBlood   * chbrIdxs...
                              + sigmaMuscle  * myocIdxs;
    % Get voltages
    vi = fwd_solve(imgd);
    vCRS(:,i - delOffset) = vi.meas;
    
    disp(["Sample " + num2str(i) + "/" + num2str(tSize)]);
end

%%
% Reconstruct MIX
vh = mean(vMIX,2);
imgr = inv_solve(imdl, vh, vMIX);
imgMIX = calc_slices(imgr);

% Reconstruct VRS
vh = mean(vVRS,2);
imgr = inv_solve(imdl, vh, vVRS);
imgVRS = calc_slices(imgr);

% Reconstruct CRS
vh = mean(vCRS,2);
imgr = inv_solve(imdl, vh, vCRS);
imgCRS = calc_slices(imgr);
