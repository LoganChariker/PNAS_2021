%ON+OFF
% kernels = { [1.0  0.5  1.0  30] ...
%             [1.0  0.5  1.0/0.5  30] ...
%             [5/6*1.0  5/6*0.5  1.0/0.5  5/6*30] ...            
%             [22/30*1.0  22/30*0.5  1.0/0.5  22/30*30] ...
%             [0.5*1.0  0.5*0.5  1.0/0.5  0.5*30] ...
%             };
kernels = { [1.0  1.0  1.0  00  0     -1] ...
            [1.0  1.0  1.0  00  0.1   1] ...
            };
        
gs = cell(1,length(kernels));
As = cell(1,length(kernels));
spTempKs = cell(1,length(kernels));

%Timesteps
dt=.1;
T=1000;
ts = 0:dt:T;        

for kernelInd=1:length(kernels);
%Kernels
scaleFac = 1;
posAdj = kernels{kernelInd}(1);
negAdj = kernels{kernelInd}(2);
negHorAdj = kernels{kernelInd}(3);
ONKernels = { [1*posAdj 1*negAdj 1*negHorAdj  10] };

%compute and implement OFFKernels
crossover=36;
tau0=2.3*crossover/22;
tau1=4.5*crossover/22;
shift=kernels{kernelInd}(4);
stretchNegHoriz=kernels{kernelInd}(3);
stretchVert=kernels{kernelInd}(1);
stretchNegVert=kernels{kernelInd}(2)/stretchVert;
gs{kernelInd} = stretchVert*LGNKernel(ts,tau0,tau1,shift,stretchNegVert,stretchNegHoriz);

sa = 0.066 * 2.0 * 1.2;
sb = 0.093 * 2.0 * 1.2;
a=1;
b=0.74;

A=@(x) a/(pi*sa^2)*exp(-(x/sa).^2)-b/(pi*sb^2)*exp(-(x/sb).^2);

d=kernels{kernelInd}(5);
% subplot(1,3,2);
xs=-1:.01:1;
ys=A(xs-d)*kernels{kernelInd}(6);
As{kernelInd}=ys;
% plot(xs,ys);

%s=subplot(1,3,3);

spTempKs{kernelInd} = (ys' * gs{kernelInd})';
end

spTempKsSummed = zeros(size(spTempKs{1}));
for i=1:length(spTempKs)
    spTempKsSummed = spTempKsSummed + spTempKs{i};
end

figure;
imagesc([xs(1) xs(end)],[ts(1) ts(end)],spTempKsSummed);
%ylim([size(im,1)-1000 size(im,1)])
set(gca,'ydir','normal');
ylim([0 100])


ylabel('t (ms)')
xlabel('x (deg)')
title('(1,1) OFF at d=0 + (2,1) ON unnormalized, 10 ms delay, at d=0.15')

cax=caxis;
%cax=[-max(cax) max(cax)];
cmsz=size(colormap,1);
cmsz=320;

cs=linspace(cax(1),cax(2),cmsz);
bys=1-max(0,cs)/cax(2);
rys=1-min(0,cs)/cax(1);
gys=1-(max(0,cs)/cax(2)+min(0,cs)/cax(1));
colormap([rys(:),gys(:),bys(:)]);
%%

revCorFreq = 1/30;
revCorPeriod = 1/revCorFreq;
lightMap = zeros(size(spTempKsSummed));
lightMap( floor(mod(ts,revCorPeriod*4)/revCorPeriod) == 0, : ) = -1;
lightMap( floor(mod(ts,revCorPeriod*4)/revCorPeriod) == 1, : ) = 0;
lightMap( floor(mod(ts,revCorPeriod*4)/revCorPeriod) == 2, : ) = 1;
lightMap( floor(mod(ts,revCorPeriod*4)/revCorPeriod) == 3, : ) = 0;

revCor = zeros(size(lightMap));

for xInd=1:size(lightMap,2);
    ysLight = lightMap(:,xInd);
    ysK = spTempKsSummed(1:150*10,xInd);
    
    convOut = conv(ysLight,ysK);    
    revCor(:,xInd) = convOut(1:length(ysLight));
end

%%
figure;
imagesc([xs(1) xs(end)],[ts(1) ts(end)],revCor)

set(gca,'ydir','normal');
ylim([0 100])

cax=caxis;
cmsz=320;
cs=linspace(cax(1),cax(2),cmsz);
bys=1-max(0,cs)/cax(2);
rys=1-min(0,cs)/cax(1);
gys=1-(max(0,cs)/cax(2)+min(0,cs)/cax(1));
colormap([rys(:),gys(:),bys(:)]);

figure;imagesc([xs(1) xs(end)],[ts(1) ts(end)],lightMap);set(gca,'ydir','normal')