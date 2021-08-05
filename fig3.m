% Plot the response panels in fig 3A,B

% Time series 
duration=1000;
dt=.1;
ts=0:dt:duration;

% Grating parameters and LGN placement
gr_sf=2.5;
gr_tf=10;
gr_dir=1;
OFF_loc=0;
ON_loc=0.1;

% ON delay and grating direction in each of the 4 response plots of A,B
ON_delays = [0,  0,  10, 10];
gr_dirs   = [+1, -1, +1, -1];

% Produce each of the response plots of A,B
figure;
for i=1:4
    subplot(2,2,i);

    ONOFF_type='OFF';
    t_delay=0;
    loc_x=OFF_loc;
    gr_dir=gr_dirs(i);
    ys_OFF=get_LGN_response(duration, dt, t_delay, loc_x, gr_sf, gr_tf, gr_dir, ONOFF_type );
    plot(ts,ys_OFF,'color',[0 0 0],'linewidth',1)

    hold on;
    ONOFF_type='ON';
    t_delay=ON_delays(i);
    loc_x=ON_loc;
    gr_dir=gr_dirs(i);
    ys_ON=get_LGN_response(duration, dt, t_delay, loc_x, gr_sf, gr_tf, gr_dir, ONOFF_type );
    plot(ts,ys_ON,'color',[1 1 1],'linewidth',1)
    
    plot(ts,ys_OFF+ys_ON,'r','linewidth',2)
    
    xlim([100 300]);
    set(gca,'color',.75*[1 1 1])
    xlabel('Time (ms)')
    ylabel('Current (sec^{-1})')
end

%% Plot panel D
gr_tfs = 1:32;
summed_resps_right = zeros(size(gr_tfs));
summed_resps_left  = zeros(size(gr_tfs));
for tf_ind=1:length(gr_tfs)
    gr_tf=gr_tfs(tf_ind);
    
    OFF_loc = 0;
    ON_loc = 0.1;
    ON_delay = 10;
    gr_sf = 2.5;
    gr_dir = +1;
    summed_resps_right(tf_ind) = get_summed_response_amplitude( OFF_loc, ON_loc, ON_delay, gr_sf, gr_tf, gr_dir );
    
    gr_dir = -1;
    summed_resps_left(tf_ind)  = get_summed_response_amplitude( OFF_loc, ON_loc, ON_delay, gr_sf, gr_tf, gr_dir );
end

figure;
plot(gr_tfs,summed_resps_right,'k-');
hold on;
plot(gr_tfs,summed_resps_left,'k--');
set(gca,'xtick',2.^(0:5),...
        'xscale','log')
xlabel('TF (Hz)')
ylabel('Current (sec^{-1})');
xlim([1 32])
        

%% functions

function ks = get_tk_time_series( duration, dt, t_delay, ONOFF_type )
    % Temporal kernel form
    K = @(t,tau0,tau1) t.^6 /tau0^7 .* exp( -t/tau0 ) - ...
                       t.^6 /tau1^7 .* exp( -t/tau1 );

    % Temporal kernel parameters
    tau0=3.66;
    tau1=7.16;
    
    ts = 0 : dt : duration;
    if isequal(ONOFF_type, 'ON')
        ks = K(ts-t_delay,tau0,tau1);
    else
        ks = -K(ts-t_delay,tau0,tau1);
    end
    ks( ts < t_delay ) = 0;
end

function ys = get_LGN_response( duration, dt, t_delay, loc_x, gr_sf, gr_tf, gr_dir, ONOFF_type )
    ts = 0 : dt : duration;
    
    % Get light intensity time series at loc_x
    period = 1 / gr_sf;
    init_phase = 2*pi * loc_x / (gr_dir*period);
    Ls = -sin( -2*pi*(gr_tf/1000)*ts + init_phase );
    
    % Get temporal kernel time series
    duration_kern = 200;
    ks = get_tk_time_series( duration_kern, dt, t_delay, ONOFF_type );
    
    % Convolve light intensity with kernel
    C = .34;
    ys_uncut = C*conv( Ls, ks ) * dt;
    ys = ys_uncut(1:length(ts));
end

function amp = get_summed_response_amplitude( OFF_loc, ON_loc, ON_delay, gr_sf, gr_tf, gr_dir )
    % Get initial phase of light intensity map at OFF and ON cell locations
    period = 1 / gr_sf;
    OFF_init_phase = 2*pi * OFF_loc / (gr_dir*period);
    ON_init_phase  = 2*pi * ON_loc  / (gr_dir*period);
    
    % Get temporal kernel time series
    duration_kern = 200;
    dt = 0.1;
    t_delay = 0;
    ks_OFF = get_tk_time_series( duration_kern, dt, t_delay, 'OFF' );
    t_delay = ON_delay;
    ks_ON  = get_tk_time_series( duration_kern, dt, t_delay, 'ON' );
    
    % Compute the amplitude of summed responses of the OFF,ON cells
    C=.34;
    ts = 0 : dt : duration_kern;
    cx_rep_OFF = exp(-1i*OFF_init_phase) * sum( ks_OFF .* exp(-2*pi*1i*(gr_tf/1000)*ts) ) * dt;
    cx_rep_ON  = exp(-1i*ON_init_phase)  * sum( ks_ON  .* exp(-2*pi*1i*(gr_tf/1000)*ts) ) * dt;
    amp = C*abs( cx_rep_OFF + cx_rep_ON );
end