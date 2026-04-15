%% Four-Neuron Network with Improved Plotting
% This script simulates a four-neuron network with calcium dynamics
% and creates publication-quality staggered plots in the requested style
clear all 
tmax=50; % max time in seconds

gelec=0.003;
% normal burst 4-5 sec 
% gelec=0 makes it shorter 3 sec
%gelec=0.00;
% gelec=0.006 makes burts longer less stable and 0.008 breaks all down
%gelec=0.006

%STEP 1
g41=0;
g41=.007; g32=g41;

%STEP 2
g14=.00;
g14=.008; g23=g14;

%STEP3 
g34=0.00;
g34=0.015; g43=g34;

% STEP4
g21=0.00; 
g21=0.035; g12=g21;

%--exc--41--Si4/3--> Si1/2-------------
% fast 3s burst
alpha4 = 0.08;  beta4 =  0.005;  scale4=alpha4/(alpha4+beta4);
alpha3 = alpha4; beta3 =  beta4; scale3=alpha3/(alpha3+beta3);


%- inh---14--Si1/2-->Si4/3-------------
% fast burts 5 sec 
alpha1 = 0.015; beta1 =  0.001; scale1=alpha1/(alpha1+beta1);
alpha1 = 0.010; beta1 =  0.001; scale1=alpha1/(alpha1+beta1);
alpha2 = alpha1; beta2 = beta1 ; scale2=alpha2/(alpha2+beta2);

%--inh---43 Si3/4---between------------
% burst ~ 3 sec 
alpha34=0.005; beta34=0.005;   scale34=alpha34/(alpha34+beta34); % inh 34 and 43 same

%--inh-between--12 Si1/2---------
alpha12=0.005; beta12= 0.005; scale12=alpha12/(alpha12+beta12);  % inh 12 and 21 same

%----------------------
cutoff=16.0; % cutoff frequency in Hz
noise=.0002   ;
%----------------------

%Parameters4
Ca_shift4 = -90; Ca_shift3 = Ca_shift4;
Ca_shift1 = -20; Ca_shift2 = Ca_shift1;
x_shift =  -4;

Iapp=0.0;

t1=00*1000; % this time is in msec
t2=00*1000; % this time is in msec
      
% H-current      
    gh    = 0.0005;
    Vhh   = -55;

% Integration    
tf = tmax*1000;   % max time in msec   
step = 0.1;       % time step in msec
tsamp = 2;        % sampling interval to save data in msec

cut1 = cutoff/1000; % cutoff frequency converted to msec

% Initial values

cut1=cutoff/1000; % cutoff frequency converted to msec
% Initial values
V4= -44; V1= -44; Ca4=.1; Ca1=1.3; h4 =0; h1 =0; n4 =0; n1 =0; 
x4 =.8; x1 =.92; s4 =0; s1 =0; s12=0; s43=0;

V3= -44; V2= -44; Ca3=.11; Ca2=1.3; h3 =0; h2 =0; n3 =0; n2 =0; 
x3 =.8; x2 =.92; s3 =0; s2 =0; s21=0; s34=0; y1=0; y2=0; y3=0; y4=0;

tic
[rnd1, nt] = bandlimnoise(cut1, step, tf);
[rnd2, nt] = bandlimnoise(cut1, step, tf);
rnd1 = rnd1*noise; rnd2 = rnd2*noise;

nt1 = round(tmax/tsamp)+1;
is = round(tsamp/step);

time = zeros(nt1,1); vv4 = time; vv1 = time; ss4 = time; ss1 = time; Caa4 = time;
Caa1 = time; xx4 = time; xx1 = time; vv3 = time; vv2 = time; ss3 = time; ss2 = time; Caa3 = time;
Caa2 = time; xx3 = time; xx2 = time; ss12 = time; ss21 = time; ss34 = time; ss43 = time;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% CALCULATE heaviside functions in advance! This speeds up calculation!

heav4 = zeros(nt,1); heav1 = heav4;
j1 = round(t1/step)+1; j2 = round(t2/step)+1;
heav4(j1:nt) = 1; heav1(j2:nt) = 1; heav = heav4-heav1; 
clear heav4 heav1;

heav3 = zeros(nt,1); heav2 = heav3;
j1 = round(t1/step)+1; j2 = round(t2/step)+1;
heav3(j1:nt) = 1; heav2(j2:nt) = 1; heav = heav3-heav2; 
clear heav3 heav2;

% Main simulation loop
j = 0;
for i = 1:nt
    tt = (i-1)*step;

    V1 = V1 + step*(Iapp*heaviside(tt-t1)*heaviside(t2-tt)+4*((0.1*(50-(127*V1/105+8265/105))/(exp((50 - (127*V1/105 ...
        +8265/105))/10) - 1))/((0.1*(50 - (127*V1/105 + 8265/105))/(exp((50 - (127*V1/105 + 8265/105))/10) - 1))+...
        (4*exp((25 - (127*V1/105 + 8265/105))/18))))^3*h1*(30 - V1) + 0.3*n1^4*(-75 - V1)+0.01*x1*(30-V1) ...
        +0.03*Ca1/(.5 + Ca1)*(-75 - V1)+0.003*(-40 - V1)+gh*((1/(1+exp(-(V1+73)/7.8)))^3)*y1*(-V2+120) ....
        -g41*(V1-30)*s4/scale4  -g12*(V1+80)*s21/scale12  +gelec*(V4-V1)+rnd2(i));

    V2 = V2 + step*(Iapp*heaviside(tt-t1)*heaviside(t2+50-tt)+4*((0.1*(50-(127*V2/105+8265/105))/(exp((50 - (127*V2/105 ...
        +8265/105))/10) - 1))/((0.1*(50 - (127*V2/105 + 8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))+...
        (4*exp((25 - (127*V2/105 + 8265/105))/18))))^3*h2*(30 - V2) + 0.3*n2^4*(-75 - V2)+0.01*x2*(30-V2) ...
        +0.03*Ca2/(.5 + Ca2)*(-75 - V2)+0.003*(-40 - V2)  ...
        -g32*(V2-30)*s3/scale3 -g12*(V2+80)*s12/scale12 +gelec*(V3-V2)+rnd2(i)); 
     
    V3 = V3 + step*(4*((0.1*(50-(127*V3/105+8265/105))/(exp((50 - (127*V3/105 ...
        +8265/105))/10) - 1))/((0.1*(50 - (127*V3/105 + 8265/105))/(exp((50 - (127*V3/105 + 8265/105))/10) - 1))+...
        (4*exp((25 - (127*V3/105 + 8265/105))/18))))^3*h3*(30 - V3) + 0.3*n3^4*(-75 - V3)+0.01*x3*(30-V3) ...
        +0.03*Ca3/(.5 + Ca3)*(-75 - V3)+0.003*(-40 - V3) +gh*((1/(1+exp(-(V3+73)/7.8)))^3)*y3*(-V3+120) ...
        -g23*(V3+80)*s2/scale2 -g34*(V3+80)*s43/scale34  +gelec*(V2-V3)+rnd1(i));

    V4 = V4 + step*(4*((0.1*(50-(127*V4/105+8265/105))/(exp((50 - (127*V4/105 ...
        +8265/105))/10) - 1))/((0.1*(50 - (127*V4/105 + 8265/105))/(exp((50 - (127*V4/105 + 8265/105))/10) - 1))+...
        (4*exp((25 - (127*V4/105 + 8265/105))/18))))^3*h4*(30 - V4) + 0.3*n4^4*(-75 - V4)+0.01*x4*(30-V4) ...
        +0.03*Ca4/(.5 + Ca4)*(-75 - V4)+0.003*(-40 - V4) +gh*((1/(1+exp(-(V4+73)/7.8)))^3)*y4*(-V4+120) ... 
        -g14*(V4+80)*s1/scale1 -g34*(V4+80)*s34/scale34 +gelec*(V1-V4)+rnd1(i));

    Ca4 = Ca4 + step*(0.0003*(0.0085*x4*(140-V4+Ca_shift4)-Ca4));
    Ca1 = Ca1 + step*(0.0003*(0.0085*x1*(140-V1+Ca_shift1)-Ca1));
    Ca3 = Ca3 + step*(0.0003*(0.0085*x3*(140-V3+Ca_shift3)-Ca3));
    Ca2 = Ca2 + step*(0.0003*(0.0085*x2*(140-V2+Ca_shift2)-Ca2));
    x4 = x4 + step*(((1/(exp(0.15*(-V4-50+x_shift))+1))-x4)/100);
    x1 = x1 + step*(((1/(exp(0.15*(-V1-50+x_shift))+1))-x1)/100);
    x3 = x3 + step*(((1/(exp(0.15*(-V3-50+x_shift))+1))-x3)/100);
    x2 = x2 + step*(((1/(exp(0.15*(-V2-50+x_shift))+1))-x2)/100);
    h4 = h4 + step*(((1-h4)*(0.07*exp((25 - (127*V4/105 + 8265/105))/20))-h4*(1.0/(1 + exp((55 - (127*V4/105 + 8265/105))/10))))/12.5);
    h1 = h1 + step*(((1-h1)*(0.07*exp((25 - (127*V1/105 + 8265/105))/20))-h1*(1.0/(1 + exp((55 - (127*V1/105 + 8265/105))/10))))/12.5);
    h3 = h3 + step*(((1-h3)*(0.07*exp((25 - (127*V3/105 + 8265/105))/20))-h3*(1.0/(1 + exp((55 - (127*V3/105 + 8265/105))/10))))/12.5);
    h2 = h2 + step*(((1-h2)*(0.07*exp((25 - (127*V2/105 + 8265/105))/20))-h2*(1.0/(1 + exp((55 - (127*V2/105 + 8265/105))/10))))/12.5);
    n4 = n4 + step*(((1-n4)*(0.01*(55 - (127*V4/105 + 8265/105))/(exp((55 - (127*V4/105 + 8265/105))/10) - 1))-n4*(0.125*exp((45 - (127*V4/105 + 8265/105))/80)))/12.5);
    n1 = n1 + step*(((1-n1)*(0.01*(55 - (127*V1/105 + 8265/105))/(exp((55 - (127*V1/105 + 8265/105))/10) - 1))-n1*(0.125*exp((45 - (127*V1/105 + 8265/105))/80)))/12.5);
    n3 = n3 + step*(((1-n3)*(0.01*(55 - (127*V3/105 + 8265/105))/(exp((55 - (127*V3/105 + 8265/105))/10) - 1))-n3*(0.125*exp((45 - (127*V3/105 + 8265/105))/80)))/12.5);
    n2 = n2 + step*(((1-n2)*(0.01*(55 - (127*V2/105 + 8265/105))/(exp((55 - (127*V2/105 + 8265/105))/10) - 1))-n2*(0.125*exp((45 - (127*V2/105 + 8265/105))/80)))/12.5);
    s4 = s4 + step*(alpha4*s4*(1-s4)/(1+exp(-20*(V4+20)))-beta4*(s4-0.0001));
    s1 = s1 + step*(alpha1*s1*(1-s1)/(1+exp(-20*(V1+20)))-beta1*(s1-0.0001));
    s3 = s3 + step*(alpha3*s3*(1-s3)/(1+exp(-20*(V3+20)))-beta3*(s3-0.0001));
    s2 = s2 + step*(alpha2*s2*(1-s2)/(1+exp(-20*(V2+20)))-beta2*(s2-0.0001));
    s12 = s12 + step*(alpha12*(1-s12)/(1+exp(-20*(V1+20)))-beta12*s12);
    s21 = s21 + step*(alpha12*(1-s21)/(1+exp(-20*(V2+20)))-beta12*s21);
    s34 = s34 + step*(alpha34*(1-s34)/(1+exp(-20*(V3+20)))-beta34*s34);
    s43 = s43 + step*(alpha34*(1-s43)/(1+exp(-20*(V4+20)))-beta34*s43);
    y1 = y1 + step*(.5*((1/(1+exp(10*(V1-Vhh))))-y1)/(7.1+10.4/(1+exp((V1+68)/2.2))));
    y2 = y2 + step*(.5*((1/(1+exp(10*(V2-Vhh))))-y2)/(7.1+10.4/(1+exp((V2+68)/2.2))));
    y3 = y3 + step*(.5*((1/(1+exp(10*(V3-Vhh))))-y3)/(7.1+10.4/(1+exp((V3+68)/2.2))));
    y4 = y4 + step*(.5*((1/(1+exp(10*(V4-Vhh))))-y4)/(7.1+10.4/(1+exp((V4+68)/2.2))));


    % Save data at specified sampling interval
    if mod(i,is) == 0
        j = j + 1;
        time(j) = tt; 
        vv4(j) = V4; vv1(j) = V1; ss4(j) = s4; ss1(j) = s1; Caa4(j) = Ca4;
        Caa1(j) = Ca1; xx4(j) = x4; xx1(j) = x1; heavy(j) = heav(i);
        vv3(j) = V3; vv2(j) = V2; ss3(j) = s3; ss2(j) = s2; Caa3(j) = Ca3;
        Caa2(j) = Ca2; xx3(j) = x3; xx2(j) = x2; ss12(j) = s12; ss21(j) = s21; 
        ss34(j) = s34; ss43(j) = s43;
    end
end  
toc    

% Convert time from msec to sec
time = time/1000;

%% Create enhanced figures with improved styling

%% Second plot: Full 4-cell staggered view
fig=figure(2);
clf
set(gcf, 'Renderer', 'painters') ;           % Use the Painters renderer for smoother lines
set(0, 'DefaultLineLineSmoothing', 'on')    % Turn on line smoothing
set(gcf, 'Position',  [100, 100, 800, 600]); % Larger figure size for 4 cells
set(gcf, 'PaperPositionMode', 'auto');
hold on

% Define offset and scaling factors for staggered plot
offset_step = 135;  % Vertical offset between traces
synscale = 30;     % Base scaling factor for synaptic variables

% Start with highest offset for top trace (Si2L)
offset = 520;       % Starting from the top

% Plot Si2L (top trace)
plot(time, vv1 + offset, 'Color', [0, 0, 0.6], 'LineWidth', 1.5)
plot(time, ss1/scale1*synscale + offset - 85, 'Color', [0, 0, 0.6], 'LineWidth', 1.5)
plot(time, ss12/scale12*synscale + offset - 85, 'Color', [0.3, 0.3, 0.3], 'LineWidth', 1)

% Decrease offset for Si3R
offset = offset - offset_step;

% Plot Si3R trace
plot(time, vv4 + offset, 'Color', [0.9, 0, 0], 'LineWidth', 1.5)
plot(time, ss4/scale4*synscale + offset - 85, 'Color', [0.9, 0, 0], 'LineWidth', 1.5)
plot(time, ss43/scale34*synscale + offset - 85, 'Color', [0.3, 0.3, 0.3], 'LineWidth', 1)

% Decrease offset for Si2R
offset = offset - offset_step;

% Plot Si2R trace
plot(time, vv2 + offset, 'Color', [0, 0, 0.99], 'LineWidth', 1.5)
plot(time, ss2/scale2*synscale + offset - 85, 'Color', [0, 0, 0.99], 'LineWidth', 1.5)
plot(time, ss21/scale12*synscale + offset - 85, 'Color', [0.3, 0.3, 0.3], 'LineWidth', 1)

% Decrease offset for Si3L
offset = offset - offset_step;

% Plot Si3L trace
plot(time, vv3 + offset, 'Color', [0.6, 0, 0], 'LineWidth', 1.5)
plot(time, ss3/scale3*synscale + offset - 85, 'Color', [0.6, 0, 0], 'LineWidth', 1.5)
plot(time, ss34/scale34*synscale + offset - 85, 'Color', [0.3 0.3, 0.3], 'LineWidth', 1)



% Remove axis ticks
set(gca, 'YTick', [])   % Remove y-axis ticks
set(gca, 'XTick', [])   % Remove x-axis ticks

% % % Add L‐shaped scale bar
x_start   = 12;
y_start   = 225;
yh_start   = y_start+5;
time_len  = 2;
volt_h    = 50;
plot([x_start, x_start], [yh_start, yh_start + volt_h], 'k', 'LineWidth', 2);
plot([x_start, x_start - time_len], [y_start, y_start], 'k', 'LineWidth', 2);
text(x_start - time_len/2, y_start - 15, [num2str(time_len) ' s'], 'HorizontalAlignment', 'center', 'FontSize', 18);
text(x_start + .2, yh_start + volt_h/2, [num2str(volt_h) ' mV'], 'VerticalAlignment', 'middle', 'FontSize', 18);

% Add trace labels
% text(tmax*0.01, 520, 'Si2L', 'FontSize', 14, 'Color', [0, 0, 0.6])
% text(tmax*0.01, 390, 'Si2R', 'FontSize', 14, 'Color', [0, 0, 0.99])
% text(tmax*0.01, 260, 'Si3L', 'FontSize', 14, 'Color', [0.6, 0, 0])
% text(tmax*0.01, 130, 'Si3R', 'FontSize', 14, 'Color', [0.9, 0, 0])
% Save the figure
%saveas(gcf, 'assembly_line.jpg');

% Set axis limits and labels
xlim([7.9 tmax-30.7])     % Extended to show labels
ylim([20 560])         % Adjusted to fit all traces
xlabel('Time [s]', 'Fontsize', 23)
ylabel('Voltage [mV]', 'Fontsize', 23)
box on

% End of stacked‐trace plotting
%exportgraphics(fig, 'norm_swim_5sec_LR.jpg', 'Resolution', 75,'ContentType', 'image', 'BackgroundColor', 'white')

%% This function generates band-limited noise
function [out, n] = bandlimnoise(cutoff, dt, tmax)
    n = round(tmax/dt)+1;
    noise = randn(1, n);
    
    % Create filter
    filt = zeros(1, n);
    if(mod(n, 2) == 0)
        filt(1:n*cutoff) = 1;
        filt(n-n*cutoff+1:n) = 1;
    else
        filt(1:ceil(n*cutoff)) = 1;
        filt(n-floor(n*cutoff)+1:n) = 1;
    end
    
    % Filter the noise in frequency domain
    noise_fft = fft(noise);
    out = real(ifft(noise_fft.*filt));
    out = out';
end