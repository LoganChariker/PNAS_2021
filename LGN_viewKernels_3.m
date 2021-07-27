%ON+OFF
% kernels = { [1.0  0.5  1.0  30] ...
%             [1.0  0.5  1.0/0.5  30] ...
%             [5/6*1.0  5/6*0.5  1.0/0.5  5/6*30] ...            
%             [22/30*1.0  22/30*0.5  1.0/0.5  22/30*30] ...
%             [0.5*1.0  0.5*0.5  1.0/0.5  0.5*30] ...
%             };
kernels = { [1.0    1.0  1.0  00  0     -1] ...
            [1.1    0.3  1.0  10  0.12  1] ...
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
crossover=33;
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

tf = 4/1000;
sf = 4;
gratingX  = exp(1i*2*pi*sf*xs);
gratingTL = exp(1i*2*pi*tf*ts);
gratingTR = exp(-1i*2*pi*tf*ts);
gratingL = real( gratingTL'*gratingX );
gratingR = real( gratingTR'*gratingX );

kernToConv = spTempKsSummed(1:1500,end:-1:1);
convOutL = conv2( gratingL, kernToConv );
convOutL = convOutL(:,floor(length(xs)/2) : size(convOutL,2)-floor(length(xs)/2));
convOutR = conv2( gratingR, kernToConv );
convOutR = convOutR(:,floor(length(xs)/2) : size(convOutR,2)-floor(length(xs)/2));

%%
figure('position',[231   261   757   596]);
subplot(2,2,1);
imagesc([xs(1) xs(end)],[ts(1) ts(end)],convOutL/1e5)

set(gca,'ydir','normal');
ylim([0 100])
setRBColor();
ylim([500 750])
title('Response to right moving grating');
ylabel('t (ms)')

subplot(2,2,2);
imagesc([xs(1) xs(end)],[ts(1) ts(end)],gratingL);set(gca,'ydir','normal')
ylim([500 750])
title('Right moving grating');

subplot(2,2,3);
imagesc([xs(1) xs(end)],[ts(1) ts(end)],convOutR/1e5)

set(gca,'ydir','normal');
ylim([500 750])
setRBColor();
ylabel('t (ms)')
xlabel('kernel loc (deg)');
title('Response to left moving grating');

subplot(2,2,4);
imagesc([xs(1) xs(end)],[ts(1) ts(end)],gratingR);set(gca,'ydir','normal')
ylim([500 750])
xlabel('kernel loc (deg)');
title('Response to left moving grating');

%% fns

function setRBColor()
cax=caxis;
cmsz=320;
cs=linspace(cax(1),cax(2),cmsz);
bys=1-max(0,cs)/cax(2);
rys=1-min(0,cs)/cax(1);
gys=1-(max(0,cs)/cax(2)+min(0,cs)/cax(1));
colormap([rys(:),gys(:),bys(:)]);
end