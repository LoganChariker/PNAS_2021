%ON+OFF
kernels = { [1.1 1 1 33 0] ...
            };

tfsLevitt      = [1];
numPeriodsPerTF= 1;

kernelInd=1;
tfInd=1;

tf = tfsLevitt(tfInd) * 2*pi /1000;

% implement different types of LGN
disp('Differentiating ON and OFF cell parameters');

%get grating parameters

sf=2.5;

%Timesteps
dt=.1;
T=1000;
ts = 0:dt:T;

figure;
subplot(1,3,1);
plot(ts, sin(tf*ts),'k');

%Kernels
scaleFac = 1;

ONKernel=kernels{kernelInd};
[~,gs]=getONOFFKernels(ts,ONKernel,[1 1 1 1 1]);
subplot(1,3,2);
plot(ts, gs);

dir=1;
d=0;
dON=0;
[~, phi, B]    = LGNResponse(ts, gs, tf, sf, dir, dON, 'ON');
subplot(1,3,3);
plot(ts, B*sin(tf*ts+phi));

disp('kernel prms set. now simulating...')


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