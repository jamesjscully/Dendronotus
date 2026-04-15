% Compare three neuromodulation models: Control, Conductance, and Alpha modulation
% Measures inter-burst periods and creates superimposed scatterplots

clear all; close all;

tmax=200; % max time in seconds
gelec=0.002;

%-----excitation from Si1----
alpha1=0.01;
beta1= 0.002;
scale1=alpha1/(alpha1+beta1);

%----Neuromodulation parameters---------------
alpham_g = 0.005;       % Forward rate for conductance modulation
betam_g = 0.0001;       % Backward rate for conductance modulation
scalem_g = (alpham_g - betam_g) / alpham_g;
gm = 3;                 % Conductance modulation factor

alpham_a = 0.005;       % Forward rate for alpha modulation
betam_a = 0.0001;       % Backward rate for alpha modulation
scalem_a = (alpham_a - betam_a) / alpham_a;
gm_a = 3;               % Alpha modulation factor

%----excitation from Si3s---------------
alphax = 0.011;
betax =  0.001;
scalex=(alphax-betax)/alphax;

%Inhibition from Si2;
alphai = .015;
betai = .0005;
scalei=(alphai-betai)/alphai;

%----- inhibition between Si3s---------------
alpha3=0.01;
beta3= 0.002;
scale3= alpha3/(alpha3+beta3);

%-----inhibition between Si2s----
alpha2=0.01;
beta2= 0.01;
scale2= alpha2/(alpha2+beta2);

%----------------------
cutoff=16.0; % cutoff frequency in Hz
noise=.0;

%Parameters4
Ca_shift0 = -10;
Ca_shift4 = -25; Ca_shift1 = -25; Ca_shift3 = Ca_shift4; Ca_shift2 = Ca_shift1;
x_shift =  -4;

Iapp=0.;
t1=1*1000; % this time is in msec
t2=35*1000; % this time is in msec
      
% Integration parameters    
tf=tmax*1000;   % max time in msec   
step=0.1;       % time step in msec
tsamp=2;        % sampling interval to save data in msec

% H-current      
gh = 0.0005;   Vhh = -53;

cut1=cutoff/1000; % cutoff frequency converted to msec

% Prepare for three model runs
model_names = {'Control (No Modulation)', 'Conductance Modulation', 'Alpha Modulation'};
colors = {[0.3 0.3 0.3], [0.8 0.2 0.2], [0.2 0.2 0.8]}; % Gray, Red, Blue

% Storage for results from all models
all_IBP_Si3 = cell(3,1);  % Inter-burst periods for Si3
all_IBP_Si2 = cell(3,1);  % Inter-burst periods for Si2
all_burst_num_Si3 = cell(3,1); % Burst numbers for Si3
all_burst_num_Si2 = cell(3,1); % Burst numbers for Si2
all_si1_freq_Si3 = cell(3,1); % Si1 frequency at Si3 bursts
all_si1_freq_Si2 = cell(3,1); % Si1 frequency at Si2 bursts
all_implied_freq_Si3 = cell(3,1); % Implied burst frequency for Si3
all_implied_freq_Si2 = cell(3,1); % Implied burst frequency for Si2

fprintf('Running three model comparisons...\n');

%% Run all three models
for model = 1:3
    fprintf('\n=== Running Model %d: %s ===\n', model, model_names{model});
    
    % Set model-specific parameters
    if model == 1  % Control - no modulation
        g0 = 0.001;
        g41 = 0.001;
        g32 = g41;
        g14 = 0.005;
        g23 = g14;
        g34 = 0.005;
        g43 = g34;
        g21 = 0.01;
        g12 = g21;
        use_g_mod = false;
        use_a_mod = false;
    elseif model == 2  % Conductance modulation
        g0 = 0.002;
        g41 = 0.04;
        g32 = g41;
        g14 = 0.005;
        g23 = g14;
        g34 = 0.005;
        g43 = g34;
        g21 = 0.01;
        g12 = g21;
        use_g_mod = true;
        use_a_mod = false;
    else  % model == 3, Alpha modulation
        g0 = 0.001;
        g41 = 0.001;
        g32 = g41;
        g14 = 0.005;
        g23 = g14;
        g34 = 0.005;
        g43 = g34;
        g21 = 0.01;
        g12 = g21;
        use_g_mod = false;
        use_a_mod = true;
    end
    
    % Initial values
    V0 = -44; V1 = -44; V2 = -44; V3 = -44; V4 = -44;
    if model == 1
        x0 = 0;
    else
        x0 = 0.9;
    end
    x1 = .6; x2 = .6; x3 = .5; x4 = .6;
    Ca0 = .3; Ca1 = 1.; Ca2 = 1.; Ca3 = 1.1; Ca4 = 1.;
    h0 = 0; h1 = 0; h2 = 0; h3 = 0; h4 = 0;
    n0 = 0; n1 = 0; n2 = 0; n3 = 0; n4 = 0;
    y1 = 0; y2 = 0;
    s0 = 0; s1 = 0; s2 = 0; s3 = 0; s4 = 0;
    s12 = 0; s21 = 0; s34 = 0.1; s43 = 0.1;
    
    if model == 2  % Conductance modulation
        sm_g = .99;
        sm_a = 0;
    elseif model == 3  % Alpha modulation
        sm_g = 0;
        sm_a = .8;
    else  % Control
        sm_g = 0;
        sm_a = 0;
    end
    
    tic
    [rnd1,nt] = bandlimnoise(cut1,step,tf);
    [rnd2,nt] = bandlimnoise(cut1,step,tf);
    rnd1=rnd1*noise; rnd2=rnd2*noise;
    nt1=round(tmax/tsamp)+1; is=round(tsamp/step);
    
    time=zeros(nt1,1);
    vv4=time;vv1=time;vv3=time;vv2=time;vv0=time;
    
    % Heaviside functions
    heav4=zeros(nt,1); heav1=heav4;
    j1=round(t1/step)+1; j2=round(t2/step)+1;
    heav4(j1:nt)=1; heav1(j2:nt)=1; heav=heav4-heav1;
    clear heav4 heav1;
    
    heav3=zeros(nt,1); heav2=heav3;
    heav3(j1:nt)=1; heav2(j2:nt)=1;
    clear heav3 heav2;
    
    j=0;
    for i=1:nt
        tt=(i-1)*step;
        
        V0new= V0 +step*(Iapp*heaviside(tt-t1)*heaviside(t2+1250-tt)+4*((0.1*(50-(127*V0/105+8265/105))/(exp((50 - (127*V0/105 ...
            +8265/105))/10) - 1))/((0.1*(50 - (127*V0/105 + 8265/105))/(exp((50 - (127*V0/105 + 8265/105))/10) - 1))+...
            (4*exp((25 - (127*V0/105 + 8265/105))/18))))^3*h0*(30 - V0) + 0.3*n0^4*(-75 - V0)+0.01*x0*(30-V0) ...
            +0.03*Ca0/(.5 + Ca0)*(-75 - V0)+0.003*(-40 - V0));
        
        scalex2 = 1;
        
        % Apply modulation based on model type
        if use_g_mod
            g41_eff = g41*(1+gm*sm_g/scalem_g);
            g32_eff = g32*(1+gm*sm_g/scalem_g);
        else
            g41_eff = g41;
            g32_eff = g32;
        end
        
        V1new= V1 +step*(4*((0.1*(50-(127*V1/105+8265/105))/(exp((50 - (127*V1/105 ...
            +8265/105))/10) - 1))/((0.1*(50 - (127*V1/105 + 8265/105))/(exp((50 - (127*V1/105 + 8265/105))/10) - 1))+...
            (4*exp((25 - (127*V1/105 + 8265/105))/18))))^3*h1*(30 - V1) + 0.3*n1^4*(-75 - V1)+0.01*x1*(30-V1) ...
            +0.03*Ca1/(.5 + Ca1)*(-75 - V1)+0.003*(-40 - V1)...
            -g41_eff*(V1-30)*s4/scalex2 ...
            -g12*(V1+80)*s21/scale2 ...
            +gelec*(V4-V1));
        
        V2new= V2 +step*(4*((0.1*(50-(127*V2/105+8265/105))/(exp((50 - (127*V2/105 ...
            +8265/105))/10) - 1))/((0.1*(50 - (127*V2/105 + 8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))+...
            (4*exp((25 - (127*V2/105 + 8265/105))/18))))^3*h2*(30 - V2) + 0.3*n2^4*(-75 - V2)+0.01*x2*(30-V2) ...
            +0.03*Ca2/(.5 + Ca2)*(-75 - V2)+0.003*(-40 - V2)...
            -g32_eff*(V2-30)*s3/scalex2 ...
            -g12*(V2+80)*s12/scale2 ...
            +gelec*(V3-V2));
        
        V3new =V3 +step*(4*((0.1*(50-(127*V3/105+8265/105))/(exp((50 - (127*V3/105 ...
            +8265/105))/10) - 1))/((0.1*(50 - (127*V3/105 + 8265/105))/(exp((50 - (127*V3/105 + 8265/105))/10) - 1))+...
            (4*exp((25 - (127*V3/105 + 8265/105))/18))))^3*h3*(30 - V3) + 0.3*n3^4*(-75 - V3)+0.01*x3*(30-V3) ...
            +0.03*Ca3/(.5 + Ca3)*(-75 - V3)+0.003*(-40 - V3) ...
            -g23*(V3+80)*s2/scalei ...
            -g34*(V3+80)*s43/scale3 ...
            +gelec*(V2-V3) ...
            -g0*(V3-30)*s0/scale1);
        
        V4new =V4 +step*(4*((0.1*(50-(127*V4/105+8265/105))/(exp((50 - (127*V4/105 ...
            +8265/105))/10) - 1))/((0.1*(50 - (127*V4/105 + 8265/105))/(exp((50 - (127*V4/105 + 8265/105))/10) - 1))+...
            (4*exp((25 - (127*V4/105 + 8265/105))/18))))^3*h4*(30 - V4) + 0.3*n4^4*(-75 - V4)+0.01*x4*(30-V4) ...
            +0.03*Ca4/(.5 + Ca4)*(-75 - V4)+0.003*(-40 - V4) ...
            -g14*(V4+80)*s1/scalei ...
            -g34*(V4+80)*s34/scale3 ...
            +gelec*(V1-V4) ...
            -g0*(V4-30)*s0/scale1);
        
        Ca0new=Ca0+step*(0.000012*(0.0085*x0*(140-V0+Ca_shift0)-Ca0));
        Ca4new=Ca4+step*(0.0003*(0.0085*x4*(140-V4+Ca_shift4)-Ca4));
        Ca1new=Ca1+step*(0.0003*(0.0085*x1*(140-V1+Ca_shift1)-Ca1));
        Ca3new=Ca3+step*(0.0003*(0.0085*x3*(140-V3+Ca_shift3)-Ca3));
        Ca2new=Ca2+step*(0.0003*(0.0085*x2*(140-V2+Ca_shift2)-Ca2));
        
        x0new =x0+step*(((1/(exp(0.15*(-V0-50+x_shift))+1))-x0)/100);
        x4new =x4+step*(((1/(exp(0.15*(-V4-50+x_shift))+1))-x4)/100);
        x1new =x1+step*(((1/(exp(0.15*(-V1-50+x_shift))+1))-x1)/100);
        x3new =x3+step*(((1/(exp(0.15*(-V3-50+x_shift))+1))-x3)/100);
        x2new =x2+step*(((1/(exp(0.15*(-V2-50+x_shift))+1))-x2)/100);
        
        h0new =h0+step*(((1-h0)*(0.07*exp((25 - (127*V0/105 + 8265/105))/20))-h0*(1.0/(1 + exp((55 - (127*V0/105 + 8265/105))/10))))/12.5);
        h4new =h4+step*(((1-h4)*(0.07*exp((25 - (127*V4/105 + 8265/105))/20))-h4*(1.0/(1 + exp((55 - (127*V4/105 + 8265/105))/10))))/12.5);
        h1new =h1+step*(((1-h1)*(0.07*exp((25 - (127*V1/105 + 8265/105))/20))-h1*(1.0/(1 + exp((55 - (127*V1/105 + 8265/105))/10))))/12.5);
        h3new =h3+step*(((1-h3)*(0.07*exp((25 - (127*V3/105 + 8265/105))/20))-h3*(1.0/(1 + exp((55 - (127*V3/105 + 8265/105))/10))))/12.5);
        h2new =h2+step*(((1-h2)*(0.07*exp((25 - (127*V2/105 + 8265/105))/20))-h2*(1.0/(1 + exp((55 - (127*V2/105 + 8265/105))/10))))/12.5);
        
        n0new =n0+step*(((1-n0)*(0.01*(55 - (127*V0/105 + 8265/105))/(exp((55 - (127*V0/105 + 8265/105))/10) - 1))-n0*(0.125*exp((45 - (127*V0/105 + 8265/105))/80)))/12.5);
        n4new =n4+step*(((1-n4)*(0.01*(55 - (127*V4/105 + 8265/105))/(exp((55 - (127*V4/105 + 8265/105))/10) - 1))-n4*(0.125*exp((45 - (127*V4/105 + 8265/105))/80)))/12.5);
        n1new =n1+step*(((1-n1)*(0.01*(55 - (127*V1/105 + 8265/105))/(exp((55 - (127*V1/105 + 8265/105))/10) - 1))-n1*(0.125*exp((45 - (127*V1/105 + 8265/105))/80)))/12.5);
        n3new =n3+step*(((1-n3)*(0.01*(55 - (127*V3/105 + 8265/105))/(exp((55 - (127*V3/105 + 8265/105))/10) - 1))-n3*(0.125*exp((45 - (127*V3/105 + 8265/105))/80)))/12.5);
        n2new =n2+step*(((1-n2)*(0.01*(55 - (127*V2/105 + 8265/105))/(exp((55 - (127*V2/105 + 8265/105))/10) - 1))-n2*(0.125*exp((45 - (127*V2/105 + 8265/105))/80)))/12.5);
        
        % Modulation variables
        s0new =s0+step*(alpha1*(1-s0)/(1+exp(-20*(V0+20)))-beta1*s0);
        
        if use_g_mod
            sm_g_new = sm_g+step*(alpham_g*sm_g*(1-sm_g)/(1+exp(-20*(V0+20)))-betam_g*(sm_g-0.0001));
        else
            sm_g_new = sm_g;
        end
        
        if use_a_mod
            sm_a_new = sm_a+step*(alpham_a*sm_a*(1-sm_a)/(1+exp(-20*(V0+20)))-betam_a*(sm_a-0.0001));
            alphax_eff = alphax*(1+gm_a*sm_a/scalem_a);
        else
            sm_a_new = sm_a;
            alphax_eff = alphax;
        end
        
        % Synaptic variables
        s4new =s4+step*(alphax_eff*s4*(1-s4)/(1+exp(-20*(V4+20)))-betax*(s4-0.0001));
        s3new =s3+step*(alphax_eff*s3*(1-s3)/(1+exp(-20*(V3+20)))-betax*(s3-0.0001));
        s1new =s1+step*(alphai*s1*(1-s1)/(1+exp(-20*(V1+20)))-betai*(s1-0.0001));
        s2new =s2+step*(alphai*s2*(1-s2)/(1+exp(-20*(V2+20)))-betai*(s2-0.0001));
        s12new =s12+step*(alpha2*(1-s12)/(1+exp(-20*(V1+20)))-beta2*s12);
        s21new =s21+step*(alpha2*(1-s21)/(1+exp(-20*(V2+20)))-beta2*s21);
        s34new =s34+step*(alpha3*(1-s34)/(1+exp(-20*(V3+20)))-beta3*s34);
        s43new =s43+step*(alpha3*(1-s43)/(1+exp(-20*(V4+20)))-beta3*s43);
        
        y1new =y1+step*(.5*((1/(1+exp(10*(V1-Vhh))))-y1)/(7.1+10.4/(1+exp((V1+68)/2.2))));
        y2new =y2+step*(.5*((1/(1+exp(10*(V2-Vhh))))-y2)/(7.1+10.4/(1+exp((V2+68)/2.2))));
        
        % Update all variables
        V0 = V0new; V1 = V1new; V2 = V2new; V3 = V3new; V4 = V4new;
        x0 = x0new; x1 = x1new; x2 = x2new; x3 = x3new; x4 = x4new;
        Ca0 = Ca0new; Ca1 = Ca1new; Ca2 = Ca2new; Ca3 = Ca3new; Ca4 = Ca4new;
        h0 = h0new; h1 = h1new; h2 = h2new; h3 = h3new; h4 = h4new;
        n0 = n0new; n1 = n1new; n2 = n2new; n3 = n3new; n4 = n4new;
        y1 = y1new; y2 = y2new;
        s0 = s0new; sm_g = sm_g_new; sm_a = sm_a_new;
        s1 = s1new; s2 = s2new; s3 = s3new; s4 = s4new;
        s12 = s12new; s21 = s21new; s34 = s34new; s43 = s43new;
        
        % Save data at sampling intervals
        if mod(i,is) ==0
            j=j+1;
            time(j)=tt;
            vv0(j)=V0; vv1(j)=V1; vv2(j)=V2; vv3(j)=V3; vv4(j)=V4;
        end
    end
    
    elapsed = toc;
    fprintf('Simulation completed in %.2f seconds\n', elapsed);
    
    % Convert time to seconds
    time=time/1000;
    
    %% Detect bursts using zero crossing logic
    % Threshold voltage for zero crossing
    V_threshold = -40;  % mV
    
    % Detect zero crossings for Si3R (V4)
    above_thresh_4 = vv4 > V_threshold;
    crossings_4 = diff(above_thresh_4);
    up_crossings_4_idx = find(crossings_4 == 1);  % Transitions from below to above threshold
    down_crossings_4_idx = find(crossings_4 == -1);  % Transitions from above to below threshold
    
    % Detect zero crossings for Si3L (V3)
    above_thresh_3 = vv3 > V_threshold;
    crossings_3 = diff(above_thresh_3);
    up_crossings_3_idx = find(crossings_3 == 1);
    down_crossings_3_idx = find(crossings_3 == -1);
    
    % Detect zero crossings for Si2L (V1)
    above_thresh_1 = vv1 > V_threshold;
    crossings_1 = diff(above_thresh_1);
    up_crossings_1_idx = find(crossings_1 == 1);
    down_crossings_1_idx = find(crossings_1 == -1);
    
    % Detect zero crossings for Si2R (V2)
    above_thresh_2 = vv2 > V_threshold;
    crossings_2 = diff(above_thresh_2);
    up_crossings_2_idx = find(crossings_2 == 1);
    down_crossings_2_idx = find(crossings_2 == -1);
    
    % Function to identify burst starts (up crossings after >1000ms down-cycle)
    min_down_time = 1.0;  % minimum down-cycle duration in seconds
    
    % Process Si3R (V4) bursts
    burst_starts_4 = [];
    if ~isempty(up_crossings_4_idx)
        for k = 1:length(up_crossings_4_idx)
            up_idx = up_crossings_4_idx(k);
            % Find the most recent down crossing before this up crossing
            prev_down_idx = down_crossings_4_idx(down_crossings_4_idx < up_idx);
            if ~isempty(prev_down_idx)
                down_idx = prev_down_idx(end);
                down_time = time(up_idx) - time(down_idx);
                if down_time >= min_down_time
                    burst_starts_4 = [burst_starts_4; time(up_idx)];
                end
            elseif time(up_idx) >= min_down_time  % First burst, check if enough time from start
                burst_starts_4 = [burst_starts_4; time(up_idx)];
            end
        end
    end
    
    % Process Si3L (V3) bursts
    burst_starts_3 = [];
    if ~isempty(up_crossings_3_idx)
        for k = 1:length(up_crossings_3_idx)
            up_idx = up_crossings_3_idx(k);
            prev_down_idx = down_crossings_3_idx(down_crossings_3_idx < up_idx);
            if ~isempty(prev_down_idx)
                down_idx = prev_down_idx(end);
                down_time = time(up_idx) - time(down_idx);
                if down_time >= min_down_time
                    burst_starts_3 = [burst_starts_3; time(up_idx)];
                end
            elseif time(up_idx) >= min_down_time
                burst_starts_3 = [burst_starts_3; time(up_idx)];
            end
        end
    end
    
    % Process Si2L (V1) bursts
    burst_starts_1 = [];
    if ~isempty(up_crossings_1_idx)
        for k = 1:length(up_crossings_1_idx)
            up_idx = up_crossings_1_idx(k);
            prev_down_idx = down_crossings_1_idx(down_crossings_1_idx < up_idx);
            if ~isempty(prev_down_idx)
                down_idx = prev_down_idx(end);
                down_time = time(up_idx) - time(down_idx);
                if down_time >= min_down_time
                    burst_starts_1 = [burst_starts_1; time(up_idx)];
                end
            elseif time(up_idx) >= min_down_time
                burst_starts_1 = [burst_starts_1; time(up_idx)];
            end
        end
    end
    
    % Process Si2R (V2) bursts
    burst_starts_2 = [];
    if ~isempty(up_crossings_2_idx)
        for k = 1:length(up_crossings_2_idx)
            up_idx = up_crossings_2_idx(k);
            prev_down_idx = down_crossings_2_idx(down_crossings_2_idx < up_idx);
            if ~isempty(prev_down_idx)
                down_idx = prev_down_idx(end);
                down_time = time(up_idx) - time(down_idx);
                if down_time >= min_down_time
                    burst_starts_2 = [burst_starts_2; time(up_idx)];
                end
            elseif time(up_idx) >= min_down_time
                burst_starts_2 = [burst_starts_2; time(up_idx)];
            end
        end
    end
    
    % Calculate inter-burst periods for Si3 (using Si3R)
    if length(burst_starts_4) > 1
        IBP_Si3 = diff(burst_starts_4);
        burst_num_Si3 = (2:length(burst_starts_4))';
    else
        IBP_Si3 = [];
        burst_num_Si3 = [];
    end
    
    % Calculate inter-burst periods for Si2 (using Si2L)
    if length(burst_starts_1) > 1
        IBP_Si2 = diff(burst_starts_1);
        burst_num_Si2 = (2:length(burst_starts_1))';
    else
        IBP_Si2 = [];
        burst_num_Si2 = [];
    end
    
    % Store results
    all_IBP_Si3{model} = IBP_Si3;
    all_IBP_Si2{model} = IBP_Si2;
    all_burst_num_Si3{model} = burst_num_Si3;
    all_burst_num_Si2{model} = burst_num_Si2;
    
    %% Detect Si1 spikes and calculate instantaneous frequency
    % Detect spikes in Si1 (V0) - use higher threshold since these are action potentials
    V_spike_threshold = -20;  % mV threshold for spike detection
    
    above_spike = vv0 > V_spike_threshold;
    spike_crossings = diff(above_spike);
    spike_times_idx = find(spike_crossings == 1);
    spike_times = time(spike_times_idx);
    
    % Calculate instantaneous Si1 spike frequency (Hz)
    if length(spike_times) > 1
        spike_intervals = diff(spike_times);  % in seconds
        spike_freq = 1 ./ spike_intervals;  % Convert to Hz
        spike_freq_times = spike_times(2:end);  % Time points for frequency measurements
    else
        spike_freq = [];
        spike_freq_times = [];
    end
    
    %% Calculate Si1 frequency at each burst start
    % For Si3 bursts
    if ~isempty(burst_starts_4) && ~isempty(spike_freq)
        si1_freq_at_si3_bursts = zeros(length(burst_starts_4), 1);
        for k = 1:length(burst_starts_4)
            % Find the Si1 frequency measurement closest in time to this burst
            [~, closest_idx] = min(abs(spike_freq_times - burst_starts_4(k)));
            if ~isempty(closest_idx)
                si1_freq_at_si3_bursts(k) = spike_freq(closest_idx);
            else
                si1_freq_at_si3_bursts(k) = NaN;
            end
        end
    else
        si1_freq_at_si3_bursts = [];
    end
    
    % For Si2 bursts
    if ~isempty(burst_starts_1) && ~isempty(spike_freq)
        si1_freq_at_si2_bursts = zeros(length(burst_starts_1), 1);
        for k = 1:length(burst_starts_1)
            [~, closest_idx] = min(abs(spike_freq_times - burst_starts_1(k)));
            if ~isempty(closest_idx)
                si1_freq_at_si2_bursts(k) = spike_freq(closest_idx);
            else
                si1_freq_at_si2_bursts(k) = NaN;
            end
        end
    else
        si1_freq_at_si2_bursts = [];
    end
    
    % Calculate implied burst frequencies from IBPs
    if ~isempty(IBP_Si3)
        implied_freq_Si3 = 1 ./ IBP_Si3;  % Convert IBP to frequency (Hz)
        % Get Si1 frequency at these burst times (skip first burst since IBP starts at burst 2)
        si1_freq_for_Si3_IBP = si1_freq_at_si3_bursts(2:end);
    else
        implied_freq_Si3 = [];
        si1_freq_for_Si3_IBP = [];
    end
    
    if ~isempty(IBP_Si2)
        implied_freq_Si2 = 1 ./ IBP_Si2;
        si1_freq_for_Si2_IBP = si1_freq_at_si2_bursts(2:end);
    else
        implied_freq_Si2 = [];
        si1_freq_for_Si2_IBP = [];
    end
    
    % Store Si1 frequency data
    all_si1_freq_Si3{model} = si1_freq_for_Si3_IBP;
    all_si1_freq_Si2{model} = si1_freq_for_Si2_IBP;
    all_implied_freq_Si3{model} = implied_freq_Si3;
    all_implied_freq_Si2{model} = implied_freq_Si2;
    
    fprintf('Model %d: Found %d bursts in Si3R, %d bursts in Si2L\n', ...
        model, length(burst_starts_4), length(burst_starts_1));
    fprintf('  Si1: %d spikes detected\n', length(spike_times));
    if ~isempty(spike_freq)
        fprintf('  Si1 mean frequency: %.2f Hz (std: %.2f Hz)\n', mean(spike_freq), std(spike_freq));
    end
    if ~isempty(IBP_Si3)
        fprintf('  Si3 mean IBP: %.2f s (std: %.2f s)\n', mean(IBP_Si3), std(IBP_Si3));
        fprintf('  Si3 implied mean frequency: %.3f Hz\n', mean(implied_freq_Si3));
    end
    if ~isempty(IBP_Si2)
        fprintf('  Si2 mean IBP: %.2f s (std: %.2f s)\n', mean(IBP_Si2), std(IBP_Si2));
        fprintf('  Si2 implied mean frequency: %.3f Hz\n', mean(implied_freq_Si2));
    end
end

%% Create superimposed scatterplots
fprintf('\nCreating comparison plots...\n');

% Check current directory
fprintf('Current directory: %s\n', pwd);

% Set renderer to painters for better compatibility
set(0, 'DefaultFigureRenderer', 'painters');

figure('Position', [100, 100, 800, 600]);
fprintf('Figure created with handle: %d\n', gcf);

% Combine Si2 and Si3 data for denser visualization
hold on;
for model = 1:3
    % Combine IBP data from both Si2 and Si3
    combined_IBP = [];
    combined_si1_freq = [];
    
    % Add Si3 data
    if ~isempty(all_IBP_Si3{model})
        combined_IBP = [combined_IBP; all_IBP_Si3{model}];
        combined_si1_freq = [combined_si1_freq; all_si1_freq_Si3{model}];
    end
    
    % Add Si2 data  
    if ~isempty(all_IBP_Si2{model})
        combined_IBP = [combined_IBP; all_IBP_Si2{model}];
        combined_si1_freq = [combined_si1_freq; all_si1_freq_Si2{model}];
    end
    
    % Convert Si1 frequency to period (1/frequency)
    if ~isempty(combined_si1_freq)
        combined_si1_period = 1 ./ combined_si1_freq;  % Convert Hz to seconds
        
        scatter(combined_si1_period, combined_IBP, 80, ...
            'MarkerFaceColor', colors{model}, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceAlpha', 0.6, ...
            'LineWidth', 0.5);
    end
end

xlabel('Si1 Spike Period (s)', 'FontSize', 14);
ylabel('Inter-Burst Period (s)', 'FontSize', 14);
title('Inter-Burst Period vs Si1 Spike Period', 'FontSize', 16, 'FontWeight', 'bold');
legend(model_names, 'Location', 'best', 'FontSize', 12);
grid on;
box on;
set(gca, 'FontSize', 12);

% Save figure
print(gcf, 'neuromod_comparison_IBP', '-dpng', '-r300');
fprintf('Figure saved as neuromod_comparison_IBP.png\n');

% Additional diagnostics
fig_file = fullfile(pwd, 'neuromod_comparison_IBP.png');
if exist(fig_file, 'file')
    fprintf('✓ PNG file confirmed at: %s\n', fig_file);
    file_info = dir(fig_file);
    fprintf('  File size: %d bytes\n', file_info.bytes);
else
    fprintf('✗ PNG file not found at expected location: %s\n', fig_file);
end

% Force figure to front and ensure it's visible
figure(gcf);
drawnow;
fprintf('Figure forced to front. Check if MATLAB figure window is visible.\n');

fprintf('\nComparison complete!\n');
