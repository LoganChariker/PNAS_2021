%ON+OFF
kernels = { [2 1   1 33 0] ...
            [1 2   1 33 0] ...
            };
dONs = [0];
grDirs = [1 -1];

tfsLevitt      = [1  2  4  6  8  10  12  14  16  18  20  25  30  35  ];
numPeriodsPerTF= [10 20 40 60 80 100 120 160 180 200 220 250 300 300 ]*3;

ampGainPrms  = zeros(length(kernels),length(tfsLevitt),length(grDirs),length(dONs));
phasePrms    = zeros(length(kernels),length(tfsLevitt),length(grDirs),length(dONs));
ampGainOffs  = zeros(1,length(tfsLevitt));
phaseOffs    = zeros(1,length(tfsLevitt));

for dirInd=1:length(grDirs)
for dONInd=1:length(dONs)
for kernelInd=1:length(kernels)
for tfInd=1:length(tfsLevitt)

tf = tfsLevitt(tfInd) * 2*pi /1000;

% implement different types of LGN
disp('Differentiating ON and OFF cell parameters');

%get grating parameters

sf=2.5;

%Timesteps
dt=.1;
T=1000;
ts = 0:dt:T;

%Kernels
scaleFac = 1;
OFFKernels = { [1 1 1  33  0] }; % [vert-stretch, neg-area-vert-stretch, neg-area-horiz-stretch, shift]
%ONKernels = { [posAdj negAdj negHorAdj ONShift] };

ONKernel=kernels{kernelInd};
OFFKernel=OFFKernels{1};
[gsOFF,gs]=getONOFFKernels(ts,ONKernel,OFFKernel);
dir=grDirs(dirInd);
d=0;
dON=dONs(dONInd);
[~, phi, B]    = LGNResponse(ts, gs, tf, sf, dir, dON, 'ON');
[~, ~, BAt4Hz] = LGNResponse(ts, gsOFF, 4*2*pi/1000, sf, dir, d, 'ON');
[~, phiOFF, BOFF] = LGNResponse(ts, gsOFF, tf, sf, dir, d, 'ON');
I_0=1/10;
contr=1.1;
lgnRedFac=1.5;
msPerSec=1000;
ampGainPrms(kernelInd,tfInd,dirInd,dONInd) = msPerSec * I_0 * lgnRedFac * contr * LGNSpFreqDep(sf,1.0) * .9 * (B/BAt4Hz);
phasePrms(kernelInd,tfInd,dirInd,dONInd) = phi;
ampGainOffs(tfInd) = msPerSec * I_0 * lgnRedFac * contr * LGNSpFreqDep(sf,1.0) * .9 * (BOFF/BAt4Hz);
phaseOffs(tfInd) = phiOFF;

end
end
end
end

disp('kernel prms set. now simulating...')

%%

figure;
nrows=length(kernels);
ncols=length(dONs)+1;

for kernelInd=1:length(kernels)
    row=kernelInd;
    col=1;
    subplot(nrows,ncols,(row-1)*ncols+col);
    %plot kernels
    ONKernel=kernels{kernelInd};
    OFFKernel=OFFKernels{1};
    [gsOFF,gsON]=getONOFFKernels(ts,ONKernel,OFFKernel);
    
    plot(ts,gsOFF,'k--');
    hold on;
    plot(ts,gsON,'k');
    xlim([0 100]);
    grid on;
    grid minor;
    if row==1
        legend('OFF','ON');
    end
    for dONInd=1:length(dONs)
        row=kernelInd;
        col=dONInd+1;
        subplot(nrows,ncols,(row-1)*ncols+col);
        
        dirIndL=1;
        dirIndR=2;
        ampsONL=ampGainPrms(kernelInd,:,dirIndL,dONInd);
        phasesONL=phasePrms(kernelInd,:,dirIndL,dONInd);
        ampsONR=ampGainPrms(kernelInd,:,dirIndR,dONInd);
        phasesONR=phasePrms(kernelInd,:,dirIndR,dONInd);
        cpxYOffs = [ampGainOffs'.*cos(phaseOffs'), ampGainOffs'.*sin(phaseOffs')];
        cpxYONRs = [ampsONR'.*cos(phasesONR'), ampsONR'.*sin(phasesONR')];
        cpxYONLs = [ampsONL'.*cos(phasesONL'), ampsONL'.*sin(phasesONL')];
        
        tfToShow=2;
        ampON=ampsONL(tfsLevitt==tfToShow);
        phaseON=phasesONL(tfsLevitt==tfToShow);
        ampOFF=ampGainOffs(tfsLevitt==tfToShow);
        phaseOFF=phaseOffs(tfsLevitt==tfToShow);
        
        tsGrating=0:2000;
        ysON =ampON*sin(2*pi*tfToShow/1000*tsGrating+phaseON);
        ysOFF=ampOFF*sin(2*pi*tfToShow/1000*tsGrating+phaseOFF);
        plot(tsGrating,ysOFF,'k--');
        hold on;
        plot(tsGrating,ysON,'k');
        title(['\phi_{dashed}=' num2str(round(phaseOFF*180/pi)) ...
               '   \phi_{solid}=' num2str(round(phaseON*180/pi))]);
        periodMs=1000/tfToShow;
        xlim([0 periodMs])
        %ylim([0 4])
        grid on;
        grid minor;

    end
end

%% fns

function [gsOFF,gsON]=getONOFFKernels(ts,ONKernel,OFFKernel)
%compute and implement OFFKernels
tau0OFF=2.3*OFFKernel(4)/22;
tau1OFF=4.5*OFFKernel(4)/22;
shift=OFFKernel(5);
stretchNegHoriz=OFFKernel(3);
stretchVert=OFFKernel(1);
stretchNegVert=OFFKernel(2)/stretchVert;
gsOFF = stretchVert*LGNKernel(ts,tau0OFF,tau1OFF,shift,stretchNegVert,stretchNegHoriz);
tau0=2.3*ONKernel(4)/22;
tau1=4.5*ONKernel(4)/22;
shift=ONKernel(5);
stretchNegHoriz=ONKernel(3);
stretchVert=ONKernel(1);
stretchNegVert=ONKernel(2)/stretchVert;
gsON = stretchVert*LGNKernel(ts,tau0,tau1,shift,stretchNegVert,stretchNegHoriz);
end