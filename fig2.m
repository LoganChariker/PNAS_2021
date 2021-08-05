% Plot temporal and spatial kernels for fig 2

% Temporal kernel form
K = @(t,tau0,tau1) t.^6 /tau0^7 .* exp( -t/tau0 ) - ...
                   t.^6 /tau1^7 .* exp( -t/tau1 );

% Temporal kernel parameters
tau0=3.66;
tau1=7.16;

% Plot temporal kernel
figure;
subplot(1,2,1);
ts = 0:.01:100;
ks = K(ts,tau0,tau1);
plot(ts,ks,'k');

% Spatial kernel form
A = @(x,alpha,beta,sigmaa,sigmab) ...
      alpha/(pi*sigmaa^2) * exp(-x.^2/sigmaa^2) - ...
      beta /(pi*sigmab^2) * exp(-x.^2/sigmab^2);
  
% Spatial kernel parameters
alpha=1.0;
beta=0.74;
sigmaa=0.0894;
sigmab=0.1259;

% Plot spatial kernel
subplot(1,2,2);
xs = -.3:.01:.3;
as = A(xs,alpha,beta,sigmaa,sigmab);
plot(xs,as,'k');
