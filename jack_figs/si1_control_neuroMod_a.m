
clear all 

tmax=200; % max time in seconds

gelec=0.002;

%-----excitation from Si1----

alpha1=0.01;
beta1= 0.002;
scale1=alpha1/(alpha1+beta1);
g0=0.001;

%----excitation from Si3s---------------
alpham=0.005;
betam= 0.0001;
scalem=(alpham-betam)/alpham;
gm = 3;

alphax = 0.011;
betax =  0.001;
scalex=(alphax-betax)/alphax;

g41=.001;
g32=g41;

%Inhibition from Si2;
g14=.005;
g23=g14;
alphai = .015;
betai = .0005;
scalei=(alphai-betai)/alphai;

%----- inhibition between Si3s---------------
alpha3=0.01;
beta3= 0.002;
scale3= alpha3/(alpha3+beta3);

g34=0.005;
g43=g34;

%-----inhibition between Si2s----
alpha2=0.01;
beta2= 0.01;
scale2= alpha2/(alpha2+beta2);

g21=0.01;
g12=g21;

%----------------------
cutoff=16.0; % cutoff frequency in Hz
noise=.0;
%----------------------

%Parameters4

Ca_shift0 = -10;

%Parameters4
Ca_shift4 = -25; Ca_shift1 = -25; Ca_shift3 = Ca_shift4; Ca_shift2 = Ca_shift1;
x_shift =  -4;

Iapp=0.;
t1=1*1000; % this time is in msec
t2=35*1000; % this time is in msec
      
% Intergration    
tf=tmax*1000;   % max time in msec   
step=0.1;       % time step in msec
tsamp=2; % sampling interval to save data in msec

% H-current      
    gh    = 0.0005;   Vhh   = -53;

cut1=cutoff/1000; % cutoff frequency converted to msec

% Initial values
V0 = -44; V1 = -44; V2 = -44; V3 = -44; V4 = -44; 
x0 = .9; x1 = .6; x2 = .6; x3 = .5; x4 = .6;
Ca0 = .3; Ca1 = 1.; Ca2 = 1.; Ca3 = 1.1; Ca4 = 1.; 
h0 = 0; h1 = 0; h2 = 0; h3 = 0; h4 = 0; n0 = 0; n1 = 0; n2 = 0; n3 = 0; n4 = 0; y1 = 0; y2 = 0; s0 = 0;s1 = 0;s2 = 0;s3 = 0;
s4 = 0;s12 = 0;s21 = 0;s34 = 0.1;s43 = 0.1;sm = .8;

tic
[rnd1,nt] = bandlimnoise (cut1,step,tf); [rnd2,nt] = bandlimnoise (cut1,step,tf); rnd1=rnd1*noise; rnd2=rnd2*noise;
nt1=round(tmax/tsamp)+1; is=round(tsamp/step);

time=zeros(nt1,1);vv4=time;vv1=time;ss4=time;ss1=time;Caa4=time; Caa1=time;xx4=time;xx1=time;
vv3=time;vv2=time;ss3=time;ss2=time;Caa3=time; Caa2=time;xx3=time;xx2=time; ss12=time; ss21=time; ss34=time; ss43=time;
Caa0=time;xx0=time;vv0=time;ss0=time; 

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% CALCULATE heaviside functions in advance! This speeds up calculation!

heav4=zeros(nt,1); heav1=heav4; j1=round(t1/step)+1; j2=round(t2/step)+1; heav4(j1:nt)=1; heav1(j2:nt)=1; heav=heav4-heav1; 
clear heav4 heav1;

heav3=zeros(nt,1); heav2=heav3; j1=round(t1/step)+1; j2=round(t2/step)+1; heav3(j1:nt)=1; heav2(j2:nt)=1; heav=heav3-heav2; 
clear heav3 heav2;

j=0;
for i=1:nt
 tt=(i-1)*step;

V0new= V0 +step*(Iapp*heaviside(tt-t1)*heaviside(t2+1250-tt)+4*((0.1*(50-(127*V0/105+8265/105))/(exp((50 - (127*V0/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V0/105 + 8265/105))/(exp((50 - (127*V0/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V0/105 + 8265/105))/18))))^3*h0*(30 - V0) + 0.3*n0^4*(-75 - V0)+0.01*x0*(30-V0) ...
    +0.03*Ca0/(.5 + Ca0)*(-75 - V0)+0.003*(-40 - V0));

scalex2 = 1;% (alphax*(1 + gm*sm) -betax)/(alphax*(1 + gm*sm));

V1new= V1 +step*(4*((0.1*(50-(127*V1/105+8265/105))/(exp((50 - (127*V1/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V1/105 + 8265/105))/(exp((50 - (127*V1/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V1/105 + 8265/105))/18))))^3*h1*(30 - V1) + 0.3*n1^4*(-75 - V1)+0.01*x1*(30-V1) ...
    +0.03*Ca1/(.5 + Ca1)*(-75 - V1)+0.003*(-40 - V1)...%+gh*((1/(1+exp(-(V1+73)/7.8)))^3)*y1*(-V2+120) ....
    -g41*(V1-30)*s4/scalex2 ... % excitation from si3
    -g12*(V1+80)*s21/scale2  ... % reciprocal inhibition
    +gelec*(V4-V1));

V2new= V2 +step*(4*((0.1*(50-(127*V2/105+8265/105))/(exp((50 - (127*V2/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V2/105 + 8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V2/105 + 8265/105))/18))))^3*h2*(30 - V2) + 0.3*n2^4*(-75 - V2)+0.01*x2*(30-V2) ...
    +0.03*Ca2/(.5 + Ca2)*(-75 - V2)+0.003*(-40 - V2)... %+gh*((1/(1+exp(-(V2+73)/7.8)))^3)*y2*(-V2+120) ...
    -g32*(V2-30)*s3/scalex2 ... % excitation from si3
    -g12*(V2+80)*s12/scale2 ... % reciprocal inhibition
    +gelec*(V3-V2)); 

V3new =V3 +step*(4*((0.1*(50-(127*V3/105+8265/105))/(exp((50 - (127*V3/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V3/105 + 8265/105))/(exp((50 - (127*V3/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V3/105 + 8265/105))/18))))^3*h3*(30 - V3) + 0.3*n3^4*(-75 - V3)+0.01*x3*(30-V3) ...
    +0.03*Ca3/(.5 + Ca3)*(-75 - V3)+0.003*(-40 - V3)  ...
    -g23*(V3+80)*s2/scalei ... # slow inhibition
    -g34*(V3+80)*s43/scale3 ...
    +gelec*(V2-V3) ...
    -g0*(V3-30)*s0/scale1); % excitation from si1

 V4new =V4 +step*(4*((0.1*(50-(127*V4/105+8265/105))/(exp((50 - (127*V4/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V4/105 + 8265/105))/(exp((50 - (127*V4/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V4/105 + 8265/105))/18))))^3*h4*(30 - V4) + 0.3*n4^4*(-75 - V4)+0.01*x4*(30-V4) ...
    +0.03*Ca4/(.5 + Ca4)*(-75 - V4)+0.003*(-40 - V4)  ... 
    -g14*(V4+80)*s1/scalei ...
    -g34*(V4+80)*s34/scale3 ...
    +gelec*(V1-V4)   ...
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

% from si1
s0new =s0+step*(alpha1*(1-s0)/(1+exp(-20*(V0+20)))-beta1*s0); 
smnew = sm+step*(alpham*sm*(1-sm)/(1+exp(-20*(V0+20)))-betam*(sm-0.0001)); % from si1 for modulation

%excitation

s4new =s4+step*(alphax*(1+gm*sm)/scalem*s4*(1-s4)/(1+exp(-20*(V4+20)))-betax*(s4-0.0001));
s3new =s3+step*(alphax*(1+gm*sm)/scalem*s3*(1-s3)/(1+exp(-20*(V3+20)))-betax*(s3-0.0001));

%slow inhibition
s1new =s1+step*(alphai*s1*(1-s1)/(1+exp(-20*(V1+20)))-betai*(s1-0.0001));
s2new =s2+step*(alphai*s2*(1-s2)/(1+exp(-20*(V2+20)))-betai*(s2-0.0001));

%reciprocal inhibition
s12new =s12+step*(alpha2*(1-s12)/(1+exp(-20*(V1+20)))-beta2*s12);
s21new =s21+step*(alpha2*(1-s21)/(1+exp(-20*(V2+20)))-beta2*s21);
s34new =s34+step*(alpha3*(1-s34)/(1+exp(-20*(V3+20)))-beta3*s34);
s43new =s43+step*(alpha3*(1-s43)/(1+exp(-20*(V4+20)))-beta3*s43);

y1new =y1+step*(.5*((1/(1+exp(10*(V1-Vhh))))-y1)/(7.1+10.4/(1+exp((V1+68)/2.2))));
y2new =y2+step*(.5*((1/(1+exp(10*(V2-Vhh))))-y2)/(7.1+10.4/(1+exp((V2+68)/2.2))));

V0 = V0new; V1 = V1new; V2 = V2new; V3 = V3new; V4 = V4new; x0 = x0new; x1 = x1new; x2 = x2new; x3 = x3new; x4 = x4new; Ca0 = Ca0new; Ca1 = Ca1new; Ca2 = Ca2new;
Ca3 = Ca3new; Ca4 = Ca4new; h0 = h0new; h1 = h1new; h2 = h2new; h3 = h3new; h4 = h4new; n0 = n0new; n1 = n1new; n2 = n2new; n3 = n3new; n4 = n4new; y1 = y1new;
y2 = y2new; s0 = s0new; sm = smnew; s1 = s1new; s2 = s2new; s3 = s3new; s4 = s4new; s12 = s12new; s21 = s21new; s34 = s34new; s43 = s43new;


% Do not need to save every point!
if mod(i,is) ==0
    j=j+1;
    time(j)=tt; vv4(j)=V4; vv1(j)=V1; ss4(j)=s4; ss1(j)=s1; Caa4(j)=Ca4;
    Caa1(j)=Ca1; xx4(j)=x4; xx1(j)=x1; heavy(j)=heav(i);
    vv3(j)=V3; vv2(j)=V2; ss3(j)=s3; ss2(j)=s2; Caa3(j)=Ca3;
    Caa2(j)=Ca2; xx3(j)=x3; xx2(j)=x2; ss12(j)=s12; ss21(j)=s21; ss34(j)=s34; ss43(j)=s43;
    vv0(j)=V0; ss0(j)=s0; ssm(j)= sm;
end 
end  
toc    

%  
%  [vp1,tp1] = findpeaks(vV4,time,'MinPeakHeight',-20);
%  fr1=length(tp1)/(tp1(end)-tp1(1))*1000
%  
%  % convert back to sec
%  tp1=tp1/1000;
 time=time/1000;
%%
 %% Second plot: Full 5-cell staggered view
fig=figure(3);
clf
set(gcf, 'Renderer', 'painters') ;           % Use the Painters renderer for smoother lines
set(0, 'DefaultLineLineSmoothing', 'on')    % Turn on line smoothing
set(gcf, 'Position',  [10, 300, 1500, 650]); % Larger figure size for 4 cells
set(gcf, 'PaperPositionMode', 'auto');
hold on

%
% Define offset and scaling factors for staggered plot
offset_step = 130;  % Vertical offset between traces
synscale = 30;     % Base scaling factor for synaptic variables

% Start with highest offset for top trace (Si2L)
offset = 700;       % Starting from the top

% Plot Si2L (top trace)
plot(time, vv1 + offset, 'Color', [0, 0, 0.6], 'LineWidth', 1.5)
%plot(time, ss1/scale1*synscale + offset - 85, 'Color', [0, 0, 0.6], 'LineWidth', 1.5)
plot(time, ss1/scale2*synscale + offset - 80, 'Color', [0.0, 0.0, 0.6], 'LineWidth', 1)

% Decrease offset for Si2R
offset = offset - offset_step;

% Plot Si2R trace
plot(time, vv2 + offset, 'Color', [0, 0, 0.99], 'LineWidth', 1.5)
hold on
plot(time, ss2/scale2*synscale + offset - 85, 'Color', [0, 0, 0.99], 'LineWidth', 1.1)
%plot(time, ss21/scale12*synscale + offset - 85, 'Color', [0., 0.0, 0.99], 'LineWidth', 1)
hold on

% Decrease offset for Si3L
offset = offset - offset_step-20;

% Plot Si3L trace
plot(time, vv3 + offset, 'Color', [0.6, 0, 0], 'LineWidth', 1.5)
hold on
%plot(time, ss3/scale3*synscale + offset - 85, 'Color', [0.6, 0, 0], 'LineWidth', 1.1)
plot(time, ss34/scale3*synscale + offset - 85, 'Color', [0.6 0.0, 0.0], 'LineWidth', 1)

% Decrease offset for Si3R
offset = offset - offset_step-20;

% Plot Si3R trace
plot(time, vv4 + offset, 'Color', [0.9, 0, 0], 'LineWidth', 1.5);
hold on;
plot(time, ss43/scale3*synscale + offset - 85, 'Color', [0.9, 0, 0], 'LineWidth', 1.5)
%plot(time, ss43/scale34*synscale + offset - 85, 'Color', [0.9, 0., 0.0], 'LineWidth', 1)

% Decrease offset for Si1
offset = offset - offset_step-10;
% Plot Si1R trace
plot(time, vv0 + offset, 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.5)
plot(time, ss0/scale1*synscale + offset - 80, 'Color', [0.4, .4, .4], 'LineWidth', 1.5)

% Remove axis ticks
set(gca, 'YTick', [])   % Remove y-axis ticks
set(gca, 'XTick', [])   % Remove x-axis ticks

% % Add L‐shaped scale bar
x_start   = tmax-18.5;
y_start   = 110;
yh_start   = y_start-0;
time_len  = 5;
volt_h    = 50;
plot([x_start, x_start], [yh_start, yh_start + volt_h], 'k', 'LineWidth', 2);
plot([x_start, x_start - time_len], [y_start, y_start], 'k', 'LineWidth', 2);
text(x_start - time_len/2, y_start + 15, [num2str(time_len) ' s'], 'HorizontalAlignment', 'center', 'FontSize', 14);
text(x_start + .9, yh_start + volt_h/2, [num2str(volt_h) ' mV'], 'VerticalAlignment', 'middle', 'FontSize', 14);


% Set axis limits and labels
xlim([50 tmax-10])     % Extended to show labels
ylim([40 750])         % Adjusted to fit all traces
xlabel('Time [s]', 'Fontsize', 23)
ylabel('Voltage [mV]', 'Fontsize', 23)
box on

% End of stacked‐trace plotting
%period 10=12 sec
%exportgraphics(fig, 's1_Fig3B_No_neuroM.jpg', 'Resolution', 75,'ContentType', 'image', 'BackgroundColor', 'white')