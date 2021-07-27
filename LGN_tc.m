%ON+OFF
kernels = { [2 1 1 35] ...
            };

tfs = 0:.1:30;

ampGainPrms  = zeros(length(kernels),length(tfs));

for kernelInd=1:length(kernels)
for tfInd=1:length(tfs)

tf = tfs(tfInd) * 2*pi /1000;

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
posAdj = kernels{kernelInd}(1);
negAdj = kernels{kernelInd}(2);
negHorAdj = kernels{kernelInd}(3);
OFFKernels = { [1   1  1  0] }; % [vert-stretch, neg-area-vert-stretch, neg-area-horiz-stretch, shift]
ONKernels = { [posAdj negAdj negHorAdj 0] };

%compute and implement OFFKernels
tau0=2.3*kernels{kernelInd}(4)/22;
tau1=4.5*kernels{kernelInd}(4)/22;

shift=ONKernels{1}(4);
stretchNegHoriz=ONKernels{1}(3);
stretchVert=ONKernels{1}(1);
stretchNegVert=ONKernels{1}(2)/stretchVert;
gs = stretchVert*LGNKernel(ts,tau0,tau1,shift,stretchNegVert,stretchNegHoriz);
dir=1;
d=0;
dON=0;
[~, phi, B]    = LGNResponse(ts, gs, tf, sf, dir, dON, 'ON');
ampGainPrms(kernelInd,tfInd) = B;
end
end

disp('kernel prms set. now simulating...')

%%

figure;
nrows=length(kernels);
ncols=1;

for kernelInd=1:length(kernels)
    row=kernelInd;
    amps=ampGainPrms(kernelInd,:);
    col=1;
    subplot(nrows,ncols,(row-1)*ncols+col);
    plot(tfs,amps);
end

