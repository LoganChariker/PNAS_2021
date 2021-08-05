% Plot the TF response panels in fig 5

% LGN parameters common to all panels
OFF_loc=0;
ON_delay=10;
ON_kernel_ab = [1.6 0.7];

% Plot parameters
TFs = [2, 4, 10, 16, 25];
ON_locs = [0.2, 0.15, 0.10, 0.05];

% Produce response panels
figure;
for TF_ind=1:length(TFs)
for ON_loc_ind=1:length(ON_locs)
    % Parameters unique to plot
    gr_tf = TFs(TF_ind);
    ON_loc = ON_locs(ON_loc_ind);

    % Plot placement
    row = ON_loc_ind;
    col = TF_ind;
    nrows = length(ON_locs);
    ncols = length(TFs);
    subplot(nrows,ncols,(row-1)*ncols+col);

    % Compute SF responses
    gr_sfs = 1 : .1 : 10;
    summed_resps_right = zeros(size(gr_sfs));
    summed_resps_left  = zeros(size(gr_sfs));
    for sf_ind=1:length(gr_sfs)
        gr_sf=gr_sfs(sf_ind);

        gr_dir = +1;
        summed_resps_right(sf_ind) = get_summed_response_amplitude( OFF_loc, ON_loc, ON_delay, gr_sf, gr_tf, gr_dir, ON_kernel_ab );

        gr_dir = -1;
        summed_resps_left(sf_ind)  = get_summed_response_amplitude( OFF_loc, ON_loc, ON_delay, gr_sf, gr_tf, gr_dir, ON_kernel_ab );
    end
    
    % Plot SF responses
    plot(gr_sfs,summed_resps_right,'k')
    hold on;
    plot(gr_sfs,summed_resps_left,'k--')
    set(gca,'xtick',2.^(0:8),...
            'xscale','log');
        
    if row==4 && col==1
        xlabel('SF (c/d)');
        ylabel('Current (sec^{-1})');
    end
end
end
%%

disp(get_SF_factor(2.5))

%% functions

function ks = get_tk_time_series( duration, dt, t_delay, ONOFF_type, kernel_ab )
    % ab notation
    a = kernel_ab(1);
    b = kernel_ab(2);
    
    % Temporal kernel form
    K = @(t,tau0,tau1) t.^6 /tau0^7 .* exp( -t/tau0 ) - ...
                       t.^6 /tau1^7 .* exp( -t/tau1 );

    % Temporal kernel parameters
    tau0=3.66;
    tau1=7.16;
    
    ts = 0 : dt : duration;
    if isequal(ONOFF_type, 'ON')
        ks = K(ts-t_delay,tau0,tau1);
        pos_part = ks>0;
        ks( pos_part ) = a * ks( pos_part );
        neg_part = ks<0;
        ks( neg_part ) = b * ks( neg_part );
    else
        ks = -K(ts-t_delay,tau0,tau1);
        neg_part = ks>0;
        ks( neg_part ) = a * ks( neg_part );
        pos_part = ks<0;
        ks( pos_part ) = b * ks( pos_part );
    end
    ks( ts < t_delay ) = 0;
end

function amp = get_summed_response_amplitude( OFF_loc, ON_loc, ON_delay, gr_sf, gr_tf, gr_dir, ON_kernel_ab )
    % Get initial phase of light intensity map at OFF and ON cell locations
    period = 1 / gr_sf;
    OFF_init_phase = 2*pi * OFF_loc / (gr_dir*period);
    ON_init_phase  = 2*pi * ON_loc  / (gr_dir*period);
    
    % Get temporal kernel time series
    duration_kern = 200;
    dt = 0.1;
    t_delay = 0;
    OFF_kernel_ab = [1 1];
    ks_OFF = get_tk_time_series( duration_kern, dt, t_delay, 'OFF', OFF_kernel_ab );
    t_delay = ON_delay;
    ks_ON  = get_tk_time_series( duration_kern, dt, t_delay, 'ON', ON_kernel_ab );
    
    % Compute the amplitude of summed responses of the OFF,ON cells
    C=get_SF_factor(gr_sf);
    ts = 0 : dt : duration_kern;
    cx_rep_OFF = exp(-1i*OFF_init_phase) * sum( ks_OFF .* exp(-2*pi*1i*(gr_tf/1000)*ts) ) * dt;
    cx_rep_ON  = exp(-1i*ON_init_phase)  * sum( ks_ON  .* exp(-2*pi*1i*(gr_tf/1000)*ts) ) * dt;
    amp = C*abs( cx_rep_OFF + cx_rep_ON );
end

function C = get_SF_factor( sf )
    % Spatial kernel form
    A = @(x,y,alpha,beta,sigmaa,sigmab) ...
          alpha/(pi*sigmaa^2) * exp(-(x.^2+y.^2)/sigmaa^2) - ...
          beta /(pi*sigmab^2) * exp(-(x.^2+y.^2)/sigmab^2);

    % Spatial kernel parameters
    alpha=1.0;
    beta=0.74;
    sigmaa=0.0894;
    sigmab=0.1259;
    
    % Store SF kernel in 2d
    dx = 0.01;
    xs = -.4:dx:.4;
    ys = xs;
    [Xs,Ys] = meshgrid(xs,ys);
    As = A(Xs,Ys,alpha,beta,sigmaa,sigmab);
    
    % Integrate against light intensity
    Ls = exp(2*pi*1i*sf*Xs);
    C = abs( sum( As(:) .* Ls(:) ) * dx^2 );
end
    