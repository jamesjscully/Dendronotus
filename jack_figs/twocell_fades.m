clear all 

tmax=50; % max time in seconds

%--------- green, blue: ts ---------------
alpha4 = 0.0175;
beta4 =  0.001;
scale4 = alpha4/(alpha4+beta4);

alpha1 = 0.01;
beta1 =  0.0012;
scale1 = alpha1/(alpha1+beta1);

g41 = .17;  % Synaptic coupling from neuron 4 to neuron 1
g14 = .1;   % Synaptic coupling from neuron 1 to neuron 4

gelec = 0.0001;  % Electrical coupling

%----------------------
cutoff = 16.0; % cutoff frequency in Hz
noise = .00;   % Noise amplitude
%----------------------

% Calcium and x-variable parameters
Ca_shift4 = -50;
Ca_shift1 = -20;
x_shift = -4;

Iapp = 0.0;  % Applied current

t1 = 00*1000; % this time is in msec
t2 = 00*1000; % this time is in msec
      
% H-current parameters     
gh = 0.000;  % h-current conductance
Vhh = -55;   % h-current half-activation

% Integration parameters 
tf = tmax*1000;     % max time in msec   
step = 0.1;         % time step in msec
tsamp = 2;          % sampling interval to save data in msec

% Setup noise
cut1 = cutoff/1000; % cutoff frequency converted to msec

% Initial values
V4 = -44; V1 = -44; Ca4 = .6; Ca1 = 1.3; h4 = 0; h1 = 0; n4 = 0; n1 = 0; 
x4 = .8; x1 = .92; s4 = 0; s1 = 0; 

tic
[rnd1, nt] = bandlimnoise(cut1, step, tf);
[rnd2, nt] = bandlimnoise(cut1, step, tf);
rnd1 = rnd1*noise; rnd2 = rnd2*noise;

nt1 = round(tmax/tsamp)+1;
is = round(tsamp/step);

time = zeros(nt1,1); vv4 = time; vv1 = time; ss4 = time; ss1 = time; Caa4 = time;
Caa1 = time; xx4 = time; xx1 = time;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% CALCULATE heaviside functions in advance! This speeds up calculation!
heaV4 = zeros(nt,1); heaV1 = heaV4;
j1 = round(t1/step)+1; j2 = round(t2/step)+1;
heaV4(j1:nt) = 1; heaV1(j2:nt) = 1; heav = heaV4-heaV1; 
clear heaV4 heaV1;

% Main simulation loop
j = 0;
for i = 1:nt
    tt = (i-1)*step;
    
    % Neuron 4 voltage dynamics
    V4 = V4 + step*(4*((0.1*(50-(127*V4/105+8265/105))/(exp((50 - (127*V4/105 ...
        +8265/105))/10) - 1))/((0.1*(50 - (127*V4/105 + 8265/105))/(exp((50 - (127*V4/105 + 8265/105))/10) - 1))+...
        (4*exp((25 - (127*V4/105 + 8265/105))/18))))^3*h4*(30 - V4) + 0.3*n4^4*(-75 - V4)+0.01*x4*(30-V4) ...
        +0.03*Ca4/(.5 + Ca4)*(-75 - V4)+0.003*(-40 - V4) -g14*(V4+80)*s1/scale4+gelec*(V1-V4)+rnd1(i));
    
    % Neuron 1 voltage dynamics
    V1 = V1 + step*(4*((0.1*(50-(127*V1/105+8265/105))/(exp((50 - (127*V1/105 ...
        +8265/105))/10) - 1))/((0.1*(50 - (127*V1/105 + 8265/105))/(exp((50 - (127*V1/105 + 8265/105))/10) - 1))+...
        (4*exp((25 - (127*V1/105 + 8265/105))/18))))^3*h1*(30 - V1) + 0.3*n1^4*(-75 - V1)+0.01*x1*(30-V1) ...
        +0.03*Ca1/(.5 + Ca1)*(-75 - V1)+0.003*(-40 - V1) -g41*(V1-40)*s4/scale4+gelec*(V4-V1)+rnd2(i));
    
    % Calcium dynamics
    Ca4 = Ca4 + step*(0.0003*(0.0085*x4*(140-V4+Ca_shift4)-Ca4));
    Ca1 = Ca1 + step*(0.0003*(0.0085*x1*(140-V1+Ca_shift1)-Ca1));
    
    % x-variable dynamics (potassium gating)
    x4 = x4 + step*(((1/(exp(0.15*(-V4-50+x_shift))+1))-x4)/100);
    x1 = x1 + step*(((1/(exp(0.15*(-V1-50+x_shift))+1))-x1)/100);
    
    % h-variable dynamics (sodium inactivation)
    h4 = h4 + step*(((1-h4)*(0.07*exp((25 - (127*V4/105 + 8265/105))/20))-h4*(1.0/(1 + exp((55 - (127*V4/105 + 8265/105))/10))))/12.5);
    h1 = h1 + step*(((1-h1)*(0.07*exp((25 - (127*V1/105 + 8265/105))/20))-h1*(1.0/(1 + exp((55 - (127*V1/105 + 8265/105))/10))))/12.5);
    
    % n-variable dynamics (potassium activation)
    n4 = n4 + step*(((1-n4)*(0.01*(55 - (127*V4/105 + 8265/105))/(exp((55 - (127*V4/105 + 8265/105))/10) - 1))-n4*(0.125*exp((45 - (127*V4/105 + 8265/105))/80)))/12.5);
    n1 = n1 + step*(((1-n1)*(0.01*(55 - (127*V1/105 + 8265/105))/(exp((55 - (127*V1/105 + 8265/105))/10) - 1))-n1*(0.125*exp((45 - (127*V1/105 + 8265/105))/80)))/12.5);
    
    % Synaptic dynamics
    s4 = s4 + step*(alpha4*s4*(1-s4)/(1+exp(-20*(V4+20)))-beta4*(s4-0.0001));
    s1 = s1 + step*(alpha1*s1*(1-s1)/(1+exp(-20*(V1+20)))-beta1*(s1-0.0001));

    % Save data at specified sampling interval
    if mod(i,is) == 0
        j = j + 1;
        time(j) = tt; 
        vv4(j) = V4; 
        vv1(j) = V1; 
        ss4(j) = s4;
        ss1(j) = s1; 
        Caa4(j) = Ca4;
        Caa1(j) = Ca1; 
        xx4(j) = x4; 
        xx1(j) = x1; 
        heavy(j) = heav(i);
    end
end  
toc    

% Convert time from msec to sec
time = time/1000;

figure(1)
clf
set(gcf, 'Renderer', 'painters')            % Use the Painters renderer for smoother lines
set(0, 'DefaultLineLineSmoothing', 'on')    % Turn on line smoothing
set(gcf, 'Position', [100, 100, 1600, 300]) % Larger figure size
set(gcf, 'PaperPositionMode', 'auto')
hold on
% Define offset and scaling factors for staggered plot
offset = 130;          % Vertical offset between traces
scale_factor_syn = 150; % Scaling factor for synaptic variables

% Plot neuron 1 (top trace) - voltage and synaptic variable
plot(time, vv1 + offset, 'Color', [0, 0.5, 0], 'LineWidth', 1.5)
plot(time, ss1/scale1*scale_factor_syn + offset - 85, 'Color', [0, 0.5, 0], 'LineWidth', 0.8)

% Plot neuron 4 (bottom trace) - voltage and synaptic variable
plot(time, vv4, 'Color', [0, 0, 0.6], 'LineWidth', 1.5)
plot(time, ss4/scale4*scale_factor_syn*2 - 85, 'Color', [0, 0, 0.6], 'LineWidth', 0.8)

% Set axis limits and labels
xlim([0 tmax*1.05])               % Extended to show labels
ylim([-95 offset+60])             % Adjusted to fit the traces
xlabel('Time (s)', 'Fontsize', 14)
ylabel('Voltage (mV)', 'Fontsize', 14)
box on

set(gca, 'YTick', []) % Remove y-axis ticks
set(gca, 'XTick', []) % Remove x-axis ticks

% Add L-shaped scale indicator in bottom right corner
% Define scale bar position and size
x_start = tmax * 0.95;              % X position for start of scale
y_start = -70;                      % Y position for start of scale
time_length = 5;                    % Time scale in seconds
voltage_height = 50;                % Voltage scale in mV


% Add trace labels
text(tmax*0.01, offset+42, 'Inhibitory', 'FontSize', 14, 'Color', [0, 0.5, 0])
text(tmax*0.01, 43, 'Excitatory', 'FontSize', 14, 'Color', [0, 0, 0.6])

%saveas(gcf, 'twocell_fades.jpg');
