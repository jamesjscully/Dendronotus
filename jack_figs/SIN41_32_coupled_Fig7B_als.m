
clear all 

tmax=100; % max time in seconds

gelec=0.0005;

%----excitation from Si3s---------------
alpha4 = .3;
beta4 =  0.06;
scale4=alpha4/(alpha4+beta4);

% inhibition from Si2s 
alpha1 = 0.012;
beta1 =  0.001;
scale1=alpha1/(alpha1+beta1);

g41=.015;
g14=.024;

g32=g41;
g23=g14;
alpha3 = alpha4;
beta3 =  beta4;
scale3=alpha3/(alpha3+beta3);
alpha2 = alpha1;
beta2 = beta1 ;
scale2=alpha2/(alpha2+beta2);


%----- inhibition between Si3s---------------
alpha34=0.01;
beta34= 0.002;
scale34=alpha34/(alpha34+beta34);

g34=0.005;
g43=g34;


%-----inhibition between Si2s----
alpha12=0.025;
beta12= 0.01;
scale12=alpha12/(alpha12+beta12);

g21=0.01;
g12=g21;

%----------------------
cutoff=16.0; % cutoff frequency in Hz
noise=.01   ;
%----------------------

%Parameters4
Ca_shift4 = 5;
Ca_shift3 = Ca_shift4+5;
Ca_shift1 = -10;

Ca_shift2 = Ca_shift1+5;
x_shift =  -4;

Iapp=.15;

t1=15*1000; % this time is in msec
t2=85*1000; % this time is in msec
      
% Intergration    
tf=tmax*1000;   % max time in msec   
step=0.1;       % time step in msec
tsamp=2; % sampling interval to save data in msec



cut1=cutoff/1000; % cutoff frequency converted to msec
% Initial values
V4= -40; V1= -40; Ca4=1.2; Ca1=1.3; h4 =0; h1 =0; n4 =0; n1 =0; 
x4 =.8; x1 =.92; s4 =0; s1 =0; s12=0; s43=0;

V3= -40; V2= -40; Ca3=1.2; Ca2=1.4; h3 =0; h2 =0; n3 =0; n2 =0; 
x3 =.8; x2 =.92; s3 =0; s2 =0; s21=0; s34=0;

tic
[rnd1,nt] = bandlimnoise (cut1,step,tf);
[rnd2,nt] = bandlimnoise (cut1,step,tf);
rnd1=rnd1*noise; rnd2=rnd2*noise;

nt1=round(tmax/tsamp)+1;
is=round(tsamp/step);

time=zeros(nt1,1);vv4=time;vv1=time;ss4=time;ss1=time;Caa4=time;
Caa1=time;xx4=time;xx1=time;
vv3=time;vv2=time;ss3=time;ss2=time;Caa3=time;
Caa2=time;xx3=time;xx2=time; ss12=time; ss21=time; ss34=time; ss43=time;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% CALCULATE heaviside functions in advance! This speeds up calculation!

heav4=zeros(nt,1); heav1=heav4;
j1=round(t1/step)+1; j2=round(t2/step)+1;
heav4(j1:nt)=1; heav1(j2:nt)=1; heav=heav4-heav1; 
clear heav4 heav1;

heav3=zeros(nt,1); heav2=heav3;
j1=round(t1/step)+1; j2=round(t2/step)+1;
heav3(j1:nt)=1; heav2(j2:nt)=1; heav=heav3-heav2; 
clear heav3 heav2;

j=0;
for i=1:nt
 tt=(i-1)*step;

V1= V1 +step*(0.*Iapp*heaviside(tt-t1)*heaviside(t2+250-tt)+4*((0.1*(50-(127*V1/105+8265/105))/(exp((50 - (127*V1/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V1/105 + 8265/105))/(exp((50 - (127*V1/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V1/105 + 8265/105))/18))))^3*h1*(30 - V1) + 0.3*n1^4*(-75 - V1)+0.01*x1*(30-V1) ...
    +0.03*Ca1/(.5 + Ca1)*(-75 - V1)+0.003*(-40 - V1) ....
    -g41*(V1-30)*s4/scale4  -g12*(V1+80)*s21/scale12  +gelec*(V4-V1)+rnd2(i));

V2= V2 +step*(0*Iapp*heaviside(tt-t1)*heaviside(t2-tt)+4*((0.1*(50-(127*V2/105+8265/105))/(exp((50 - (127*V2/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V2/105 + 8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V2/105 + 8265/105))/18))))^3*h2*(30 - V2) + 0.3*n2^4*(-75 - V2)+0.01*x2*(30-V2) ...
    +0.03*Ca2/(.5 + Ca2)*(-75 - V2)+0.003*(-40 - V2) ...
    -g32*(V2-30)*s3/scale3 -g12*(V2+80)*s12/scale12 +gelec*(V3-V2)+rnd2(i));

V3 =V3 +step*(0.95*Iapp*heaviside(tt-t1)*heaviside(t2-tt)+4*((0.1*(50-(127*V3/105+8265/105))/(exp((50 - (127*V3/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V3/105 + 8265/105))/(exp((50 - (127*V3/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V3/105 + 8265/105))/18))))^3*h3*(30 - V3) + 0.3*n3^4*(-75 - V3)+0.01*x3*(30-V3) ...
    +0.03*Ca3/(.5 + Ca3)*(-75 - V3)+0.003*(-40 - V3) -g23*(V3+80)*s2/scale2 -g34*(V3+80)*s43/scale34  +gelec*(V2-V3)+rnd1(i));


V4 =V4 +step*(Iapp*heaviside(tt+500-t1)*heaviside(t2+20-tt)+4*((0.1*(50-(127*V4/105+8265/105))/(exp((50 - (127*V4/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V4/105 + 8265/105))/(exp((50 - (127*V4/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V4/105 + 8265/105))/18))))^3*h4*(30 - V4) + 0.3*n4^4*(-75 - V4)+0.01*x4*(30-V4) ...
    +0.03*Ca4/(.5 + Ca4)*(-75 - V4)+0.003*(-40 - V4) -g14*(V4+80)*s1/scale1 -g34*(V4+80)*s34/scale34 +gelec*(V1-V4)+rnd1(i));



Ca4=Ca4+step*(0.0003*(0.0085*x4*(140-V4+Ca_shift4)-Ca4));
Ca1=Ca1+step*(0.0003*(0.0085*x1*(140-V1+Ca_shift1)-Ca1));
Ca3=Ca3+step*(0.0003*(0.0085*x3*(140-V3+Ca_shift3)-Ca3));
Ca2=Ca2+step*(0.0003*(0.0085*x2*(140-V2+Ca_shift2)-Ca2));

x4 =x4+step*(((1/(exp(0.15*(-V4-50+x_shift))+1))-x4)/100);
x1 =x1+step*(((1/(exp(0.15*(-V1-50+x_shift))+1))-x1)/100);
x3 =x3+step*(((1/(exp(0.15*(-V3-50+x_shift))+1))-x3)/100);
x2 =x2+step*(((1/(exp(0.15*(-V2-50+x_shift))+1))-x2)/100);

h4 =h4+step*(((1-h4)*(0.07*exp((25 - (127*V4/105 + 8265/105))/20))-h4*(1.0/(1 + exp((55 - (127*V4/105 + 8265/105))/10))))/12.5);
h1 =h1+step*(((1-h1)*(0.07*exp((25 - (127*V1/105 + 8265/105))/20))-h1*(1.0/(1 + exp((55 - (127*V1/105 + 8265/105))/10))))/12.5);
h3 =h3+step*(((1-h3)*(0.07*exp((25 - (127*V3/105 + 8265/105))/20))-h3*(1.0/(1 + exp((55 - (127*V3/105 + 8265/105))/10))))/12.5);
h2 =h2+step*(((1-h2)*(0.07*exp((25 - (127*V2/105 + 8265/105))/20))-h2*(1.0/(1 + exp((55 - (127*V2/105 + 8265/105))/10))))/12.5);

n4 =n4+step*(((1-n4)*(0.01*(55 - (127*V4/105 + 8265/105))/(exp((55 - (127*V4/105 + 8265/105))/10) - 1))-n4*(0.125*exp((45 - (127*V4/105 + 8265/105))/80)))/12.5);
n1 =n1+step*(((1-n1)*(0.01*(55 - (127*V1/105 + 8265/105))/(exp((55 - (127*V1/105 + 8265/105))/10) - 1))-n1*(0.125*exp((45 - (127*V1/105 + 8265/105))/80)))/12.5);
n3 =n3+step*(((1-n3)*(0.01*(55 - (127*V3/105 + 8265/105))/(exp((55 - (127*V3/105 + 8265/105))/10) - 1))-n3*(0.125*exp((45 - (127*V3/105 + 8265/105))/80)))/12.5);
n2 =n2+step*(((1-n2)*(0.01*(55 - (127*V2/105 + 8265/105))/(exp((55 - (127*V2/105 + 8265/105))/10) - 1))-n2*(0.125*exp((45 - (127*V2/105 + 8265/105))/80)))/12.5);

s4 =s4+step*(alpha4*s4*(1-s4)/(1+exp(-20*(V4+20)))-beta4*(s4-0.0001));
s1 =s1+step*(alpha1*s1*(1-s1)/(1+exp(-20*(V1+20)))-beta1*(s1-0.0001));
s3 =s3+step*(alpha3*s3*(1-s3)/(1+exp(-20*(V3+20)))-beta3*(s3-0.0001));
s2 =s2+step*(alpha2*s2*(1-s2)/(1+exp(-20*(V2+20)))-beta2*(s2-0.0001));

s12 =s12+step*(alpha12*(1-s12)/(1+exp(-20*(V1+20)))-beta12*s12);
s21 =s21+step*(alpha12*(1-s21)/(1+exp(-20*(V2+20)))-beta12*s21);
s34 =s34+step*(alpha34*(1-s34)/(1+exp(-20*(V3+20)))-beta34*s34);
s43 =s43+step*(alpha34*(1-s43)/(1+exp(-20*(V4+20)))-beta34*s43);


% Do not need to save every point!
if mod(i,is) ==0
    j=j+1;
    time(j)=tt; vv4(j)=V4; vv1(j)=V1; ss4(j)=s4; ss1(j)=s1; Caa4(j)=Ca4;
    Caa1(j)=Ca1; xx4(j)=x4; xx1(j)=x1; heavy(j)=heav(i);
    vv3(j)=V3; vv2(j)=V2; ss3(j)=s3; ss2(j)=s2; Caa3(j)=Ca3;
    Caa2(j)=Ca2; xx3(j)=x3; xx2(j)=x2; ss12(j)=s12; ss21(j)=s21; ss34(j)=s34; ss43(j)=s43;
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
 %     
fig=figure(7);
clf;
set(gcf, 'Renderer', 'painters');              % smoother lines
set(gcf, 'Position', [100, 100, 1000, 500]);  % larger window
hold on;

% Define offsets, scale factors, and colors for each trace
offsets      = [520, 390, 260, 130];
scale_factors= [100, 100,  40,  40];
colors       = {[0,0,0.6], [0,0,0.99], [0.6,0,0], [0.9,0,0]};
syn_traces   = {vv1, vv2, vv3, vv4};
ss_traces    = {ss1, ss2, ss34, ss43};
scales       = [scale1, scale2, scale34, scale34];

% Plot each synaptic voltage and corresponding activation
for k = 1:4
    off = offsets(k);
    plot(time, syn_traces{k} + off, 'Color', colors{k}, 'LineWidth', 1.5);
end


% % Add L‐shaped scale bar
% x_start   = tmax -39;
% y_start   = 220;
% yh_start   = y_start+10;
% time_len  = 2;
% volt_h    = 50;
% plot([x_start, x_start], [yh_start, yh_start + volt_h], 'k', 'LineWidth', 2);
% plot([x_start, x_start - time_len], [y_start, y_start], 'k', 'LineWidth', 2);
% text(x_start - time_len/2, y_start - 15, [num2str(time_len) ' s'], 'HorizontalAlignment', 'center', 'FontSize', 14);
% text(x_start + .2, yh_start + volt_h/2, [num2str(volt_h) ' mV'], 'VerticalAlignment', 'middle', 'FontSize', 14);

% End of stacked‐trace plotting
xc=0.2;
a=annotation('arrow', [xc xc], [190/560 210/560]); 
a.Color = 'red';
a.LineWidth = 4; 
a=annotation('arrow', [xc xc], [80/560 100/560]); 
a.Color = 'red';
a.LineWidth = 4; 


xc=0.813;
a=annotation('arrow', [xc xc], [210/560 190/560]); 
a.Color = 'red';
a.LineWidth = 4; 
a=annotation('arrow', [xc xc], [100/560 80/560]); 
a.Color = 'red';
a.LineWidth = 4; 


% Axes limits and labels
xlim([6, tmax-5]);
ylim([25, 560]);
xlabel('Time [s]', 'FontSize', 23);
ylabel('Voltage [mV]', 'FontSize', 23);
box on;

% Remove axis ticks
set(gca, 'YTick', []) % Remove y-axis ticks
set(gca, 'XTick', []) % Remove x-axis ticks


% End of stacked‐trace plotting
%saveas(gcf, '5f.jpg');
%exportgraphics(fig, '7b_als.jpg', 'Resolution', 75,'ContentType', 'image', 'BackgroundColor', 'white')