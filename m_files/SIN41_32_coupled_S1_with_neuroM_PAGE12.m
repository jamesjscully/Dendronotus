
clear all 

tmax=100; % max time in seconds

gelec=0.002;

%----excitation from Si3s---------------
alpha4 = 0.05;
beta4 =  0.008;
scale4=alpha4/(alpha4+beta4);

alpha1 = 0.012;
beta1 =  0.001;
scale1=alpha1/(alpha1+beta1);

g41=.01;
g32=g41;

%Inhibition from Si2;
g14=.02;

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

g34=0.01;
g43=g34;


%-----inhibition between Si2s----
alpha12=0.01;
beta12= 0.01;
scale12=alpha12/(alpha12+beta12);

g21=0.012;
g12=g21;

%----------------------
cutoff=16.0; % cutoff frequency in Hz
noise=.0  ;
%----------------------


%-----excitation from Si1----

alpha0=0.006;
beta0= 0.00018;
scale0=alpha0/(alpha0+beta0);
g0=0.006;

%----------------------
cutoff=16.0; % cutoff frequency in Hz
noise=.0  ;
%----------------------

%Parameters4

Ca_shift0 = -30;

%Parameters4
Ca_shift4 = -25;
Ca_shift1 = -24;
Ca_shift3 = Ca_shift4+5;
Ca_shift2 = Ca_shift1+8;
x_shift =  -4;

% H-current      
    gh    = 0.0005;
    Vhh   = -53;

Iapp=0.25;    
t1=15*1000; % this time is in msec
t2=55*1000; % this time is in msec
      
% Intergration    
tf=tmax*1000;   % max time in msec   
step=0.1;       % time step in msec
tsamp=2; % sampling interval to save data in msec

      
% Intergration    
tf=tmax*1000;   % max time in msec   
step=0.1;       % time step in msec
tsamp=2; % sampling interval to save data in msec


cut1=cutoff/1000; % cutoff frequency converted to msec
% Initial values
V4= -44; V1= -44; Ca4=1.1; Ca1=1.3; h4 =0; h1 =0; n4 =0; n1 =0; y1=0;
x4 =.99; x1 =.92; s4 =0; s1 =0; s12=0; s43=0;

V3= -44; V2= -44; Ca3=1.2; Ca2=1.3; h3 =0; h2 =0; n3 =0; n2 =0; y2=0;
x3 =.89; x2 =.92; s3 =0; s2 =0; s21=0; s34=0;
V0= -44; Ca0=1.3; h0 =0; n0 =0; s0=0; x0=0;


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
Caa0=time;xx0=time;vv0=time;ss0=time;

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


V0= V0 +step*(Iapp*heaviside(tt-t1)*heaviside(t2+1250-tt)+4*((0.1*(50-(127*V0/105+8265/105))/(exp((50 - (127*V0/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V0/105 + 8265/105))/(exp((50 - (127*V0/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V0/105 + 8265/105))/18))))^3*h0*(30 - V0) + 0.3*n0^4*(-75 - V0)+0.01*x0*(30-V0) ...
    +0.03*Ca0/(.5 + Ca0)*(-75 - V0)+0.003*(-40 - V0));

V1= V1 +step*(4*((0.1*(50-(127*V1/105+8265/105))/(exp((50 - (127*V1/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V1/105 + 8265/105))/(exp((50 - (127*V1/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V1/105 + 8265/105))/18))))^3*h1*(30 - V1) + 0.3*n1^4*(-75 - V1)+0.01*x1*(30-V1) ...
    +0.03*Ca1/(.5 + Ca1)*(-75 - V1)+0.003*(-40 - V1)+gh*((1/(1+exp(-(V1+73)/7.8)))^3)*y1*(-V2+120) ....
    -g41*(V1-30)*s4/scale4*(1+s0/scale0)  -g12*(V1+80)*s21/scale12  +gelec*(V4-V1)+rnd2(i));

V2= V2 +step*(4*((0.1*(50-(127*V2/105+8265/105))/(exp((50 - (127*V2/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V2/105 + 8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V2/105 + 8265/105))/18))))^3*h2*(30 - V2) + 0.3*n2^4*(-75 - V2)+0.01*x2*(30-V2) ...
    +0.03*Ca2/(.5 + Ca2)*(-75 - V2)+0.003*(-40 - V2) +gh*((1/(1+exp(-(V2+73)/7.8)))^3)*y2*(-V2+120) ...
    -g32*(V2-30)*s3/scale3*(1+s0/scale0) -g12*(V2+80)*s12/scale12 +gelec*(V3-V2)+rnd2(i)); 
 

V3 =V3 +step*(4*((0.1*(50-(127*V3/105+8265/105))/(exp((50 - (127*V3/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V3/105 + 8265/105))/(exp((50 - (127*V3/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V3/105 + 8265/105))/18))))^3*h3*(30 - V3) + 0.3*n3^4*(-75 - V3)+0.01*x3*(30-V3) ...
    +0.03*Ca3/(.5 + Ca3)*(-75 - V3)+0.003*(-40 - V3)  ...
    -g23*(V3+80)*s2/scale2*(1+0*s0/scale0) -g34*(V3+80)*s43/scale34  +gelec*(V2-V3)+rnd1(i)  -g0*(V3-30)*s0/scale0) ;

 V4 =V4 +step*(4*((0.1*(50-(127*V4/105+8265/105))/(exp((50 - (127*V4/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V4/105 + 8265/105))/(exp((50 - (127*V4/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V4/105 + 8265/105))/18))))^3*h4*(30 - V4) + 0.3*n4^4*(-75 - V4)+0.01*x4*(30-V4) ...
    +0.03*Ca4/(.5 + Ca4)*(-75 - V4)+0.003*(-40 - V4)  ... 
    -g14*(V4+80)*s1/scale1*(1+0*s0/scale0) -g34*(V4+80)*s34/scale34 +gelec*(V1-V4)+rnd1(i)   -g0*(V4-30)*s0/scale0);


Ca0=Ca0+step*(0.000012*(0.0085*x0*(140-V0+Ca_shift0)-Ca0)); 
Ca4=Ca4+step*(0.0003*(0.0085*x4*(140-V4+Ca_shift4)-Ca4));
Ca1=Ca1+step*(0.0003*(0.0085*x1*(140-V1+Ca_shift1)-Ca1));
Ca3=Ca3+step*(0.0003*(0.0085*x3*(140-V3+Ca_shift3)-Ca3));
Ca2=Ca2+step*(0.0003*(0.0085*x2*(140-V2+Ca_shift2)-Ca2));

x0 =x0+step*(((1/(exp(0.15*(-V0-50+x_shift))+1))-x0)/100);
x4 =x4+step*(((1/(exp(0.15*(-V4-50+x_shift))+1))-x4)/100);
x1 =x1+step*(((1/(exp(0.15*(-V1-50+x_shift))+1))-x1)/100);
x3 =x3+step*(((1/(exp(0.15*(-V3-50+x_shift))+1))-x3)/100);
x2 =x2+step*(((1/(exp(0.15*(-V2-50+x_shift))+1))-x2)/100);

h0 =h0+step*(((1-h0)*(0.07*exp((25 - (127*V0/105 + 8265/105))/20))-h0*(1.0/(1 + exp((55 - (127*V0/105 + 8265/105))/10))))/12.5);
h4 =h4+step*(((1-h4)*(0.07*exp((25 - (127*V4/105 + 8265/105))/20))-h4*(1.0/(1 + exp((55 - (127*V4/105 + 8265/105))/10))))/12.5);
h1 =h1+step*(((1-h1)*(0.07*exp((25 - (127*V1/105 + 8265/105))/20))-h1*(1.0/(1 + exp((55 - (127*V1/105 + 8265/105))/10))))/12.5);
h3 =h3+step*(((1-h3)*(0.07*exp((25 - (127*V3/105 + 8265/105))/20))-h3*(1.0/(1 + exp((55 - (127*V3/105 + 8265/105))/10))))/12.5);
h2 =h2+step*(((1-h2)*(0.07*exp((25 - (127*V2/105 + 8265/105))/20))-h2*(1.0/(1 + exp((55 - (127*V2/105 + 8265/105))/10))))/12.5);

n0 =n0+step*(((1-n0)*(0.01*(55 - (127*V0/105 + 8265/105))/(exp((55 - (127*V0/105 + 8265/105))/10) - 1))-n0*(0.125*exp((45 - (127*V0/105 + 8265/105))/80)))/12.5);
n4 =n4+step*(((1-n4)*(0.01*(55 - (127*V4/105 + 8265/105))/(exp((55 - (127*V4/105 + 8265/105))/10) - 1))-n4*(0.125*exp((45 - (127*V4/105 + 8265/105))/80)))/12.5);
n1 =n1+step*(((1-n1)*(0.01*(55 - (127*V1/105 + 8265/105))/(exp((55 - (127*V1/105 + 8265/105))/10) - 1))-n1*(0.125*exp((45 - (127*V1/105 + 8265/105))/80)))/12.5);
n3 =n3+step*(((1-n3)*(0.01*(55 - (127*V3/105 + 8265/105))/(exp((55 - (127*V3/105 + 8265/105))/10) - 1))-n3*(0.125*exp((45 - (127*V3/105 + 8265/105))/80)))/12.5);
n2 =n2+step*(((1-n2)*(0.01*(55 - (127*V2/105 + 8265/105))/(exp((55 - (127*V2/105 + 8265/105))/10) - 1))-n2*(0.125*exp((45 - (127*V2/105 + 8265/105))/80)))/12.5);

s0 =s0+step*(alpha0*s0*(1-s0)/(1+exp(-20*(V0+20)))-beta0*(s0-0.0001));
s4 =s4+step*(alpha4*s4*(1-s4)/(1+exp(-20*(V4+20)))-beta4*(s4-0.0001));
s1 =s1+step*(alpha1*s1*(1-s1)/(1+exp(-20*(V1+20)))-beta1*(s1-0.0001));
s3 =s3+step*(alpha3*s3*(1-s3)/(1+exp(-20*(V3+20)))-beta3*(s3-0.0001));
s2 =s2+step*(alpha2*s2*(1-s2)/(1+exp(-20*(V2+20)))-beta2*(s2-0.0001));

s12 =s12+step*(alpha12*(1-s12)/(1+exp(-20*(V1+20)))-beta12*s12);
s21 =s21+step*(alpha12*(1-s21)/(1+exp(-20*(V2+20)))-beta12*s21);
s34 =s34+step*(alpha34*(1-s34)/(1+exp(-20*(V3+20)))-beta34*s34);
s43 =s43+step*(alpha34*(1-s43)/(1+exp(-20*(V4+20)))-beta34*s43);

y1 =y1+step*(.5*((1/(1+exp(10*(V1-Vhh))))-y1)/(7.1+10.4/(1+exp((V1+68)/2.2))));
y2 =y2+step*(.5*((1/(1+exp(10*(V2-Vhh))))-y2)/(7.1+10.4/(1+exp((V2+68)/2.2))));

% Do not need to save every point!
if mod(i,is) ==0
    j=j+1;
    time(j)=tt; vv4(j)=V4; vv1(j)=V1; ss4(j)=s4; ss1(j)=s1; Caa4(j)=Ca4;
    Caa1(j)=Ca1; xx4(j)=x4; xx1(j)=x1; heavy(j)=heav(i);
    vv3(j)=V3; vv2(j)=V2; ss3(j)=s3; ss2(j)=s2; Caa3(j)=Ca3;
    Caa2(j)=Ca2; xx3(j)=x3; xx2(j)=x2; ss12(j)=s12; ss21(j)=s21; ss34(j)=s34; ss43(j)=s43;
    vv0(j)=V0; ss0(j)=s0;
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
 figure(4)
 clf
 subplot(6,1,1)
 plot(time,vv1,'Color',[0 .0 0.6]','LineWidth',1.5)
 hold on
 xlim([0 tmax]) 
 ylim([-75 50])
 xlabel('Time [sec]','Fontsize', 16),ylabel('Si2L','Fontsize', 16)


 subplot(6,1,2)
 plot(time,vv2,'Color',[0 .0 0.99]','LineWidth',1.5)
 hold on
 xlim([0 tmax]) 
 ylim([-75 50])
 xlabel('Time [sec]','Fontsize', 16),ylabel('Si2R','Fontsize', 16)


 subplot(6,1,3)
 plot(time,vv3,'Color',[0.6 0  0]','LineWidth',1.5)
 hold on
 xlim([0 tmax]) 
 ylim([-75 50])
xlabel('Time [sec]','Fontsize', 16),ylabel('Si3L','Fontsize', 16)


subplot(6,1,4)
 plot(time,vv4,'Color',[0.9 0  0]','LineWidth',1.5)
 hold on
 xlim([0 tmax]) 
 ylim([-75 50])
xlabel('Time [sec]','Fontsize', 16),ylabel('Si3R','Fontsize', 16)
%

subplot(6,1,6)
 plot(time,vv0,'Color',[0.5 0.5 0.5]','LineWidth',1.5)
 hold on
 xlim([0 tmax]) 
 ylim([-75 50])
xlabel('Time [sec]','Fontsize', 16),ylabel('Si1','Fontsize', 16)
%


%fftt = 28;
%txpo = text(-4.4,25,'B','Fontsize',fftt,'Color','black','FontName','Arial','FontWeight','bold')


 subplot(6,1,5)
 plot(time,ss4/scale4,'Color',[0.6 0.0 0.6],'LineWidth',1.5)
 hold on
 plot(time,ss1/scale1,'Color',[0.0 0.  0.6],'LineWidth',1.5)
 hold on
 plot(time,ss12/scale12,'Color',[0.0 0.5 0.0],'LineWidth',0.5)
 hold on
 plot(time,ss43/scale34,'Color',[0.1 0.1 0.1],'LineWidth',0.5)
 hold on
 plot(time,ss0/scale0,'Color',[0.7 0.7 0.7],'LineWidth',1.5)
 hold on
 ylim([-0.2  .9])
 xlim([0 tmax])
 box on
 xlabel('Time [sec]','Fontsize', 16),ylabel('rates S4,S1','Fontsize', 16)




 % subplot(6,1,6)
 % plot(time,ss3/scale3,'Color',[0.6 0.0 0.0],'LineWidth',1.5)
 % hold on
 % plot(time,ss2/scale2,'Color',[0.0 0. 0.6],'LineWidth',1.5)
 % hold on
 % plot(time,ss21/scale12,'Color',[0.0 0.5 0.0],'LineWidth',0.5)
 % hold on
 % plot(time,ss34/scale34,'Color',[0.1 0.1 0.1],'LineWidth',0.5)
 % hold on
 % ylim([-0.1  .9])
 % xlim([0 tmax])
 % box on
 % xlabel('Time [sec]','Fontsize', 16),ylabel('rates S4,S1','Fontsize', 16)


% subplot(6,1,5)
%  plot(time,Caa4,'Color',[0.6 0.0 0.0],'Linewidth',1.5)
%  hold on
%  plot(time,Caa1,'Color',[0.0 0. 0.6],'Linewidth',1.5)
%  hold on
%  ylim([0.1 1.4])
%  xlim([0 tmax])
%  xlabel('Time [sec]','Fontsize', 16),ylabel('[Ca 4/1]','Fontsize', 16)
% 
% subplot(6,1,6)
%  plot(time,Caa3,'Color',[0.99 0.0 0.0],'Linewidth',1.5)
%  hold on
%  plot(time,Caa1,'Color',[0.0 0. 0.99],'Linewidth',1.5)
%  hold on
%  ylim([0.1 1.4])
%  xlim([0 tmax])
%  xlabel('Time [sec]','Fontsize', 16),ylabel('[Ca 4/1]','Fontsize', 16)



  %print(gcf,'-djpeg','-r600' ,'fredo3b.jpeg');

  toc 
   figure(5)
 clf
 st=1;
 plot(Caa4(st:end),xx4(st:end),'Color',[0.6 0. 0.],'LineWidth',2)
  hold on
  plot(Caa1(st:end),xx1(st:end),'Color',[0.0 0.0 0.6],'LineWidth',2)
  hold on
  plot(Caa3(st:end),xx3(st:end),'Color',[0.99 0. 0.],'LineWidth',2)
  hold on
  plot(Caa2(st:end),xx2(st:end),'Color',[0.0 0.0 0.99],'LineWidth',2)
  hold on
  %title('Ca vs x', 'Fontsize', 11);
  

 xlabel('[Ca]-variable','Fontsize', 16),
  ylabel('x-variable','Fontsize', 16)  
  %title('Ca vs Ca-activaterd possasion gating varaible', 'Fontsize', 11);
  
  %SNIC part 
%--------------------- SNIC whole system gh = 0 -------------------------------
data1 = load('LP_whole_bk_gh0.mat','x');
cod_back = data1.x;
data1 = load('LP_whole_fwd_gh0.mat','x');
cod_fwd = data1.x;
x_cod_bkgh0 = cod_back(end-3,:);
Ca_cod_bkgh0 = cod_back(end-2,:);
x_cod_fdgh0 =  cod_fwd(end-3,1:150);
Ca_cod_fdgh0 = cod_fwd(end-2,1:150);

%------------------- plot SNIC ---------------------
plot(Ca_cod_bkgh0,x_cod_bkgh0-0.02,'-.','Color',[.85 .85 .85 ],'LineWidth',4)
hold on
plot(Ca_cod_fdgh0,x_cod_fdgh0-0.02,'-.','Color',[.85 .85 .85 ],'LineWidth',4)
hold on

axis([0.1 1.4 0.01 1.0])

fftt = 28;
txpo = text(0.1,.98,'А','Fontsize',fftt,'Color','black','FontName','Arial','FontWeight','bold')

fftt = 15;
txpo = text(0.22,.15,'[Ca]^\prime=0','Fontsize',fftt,'Color','red','FontName','Arial')

fftt = 16;
txpo = text(1.26,.17,'x^\prime=0','Fontsize',fftt,'Color','black','FontName','Arial')

fftt = 16;
txpo = text(.22,.44,'SNIC','Fontsize',fftt,'Color',[0.6 0.6 0.6],'FontName','Arial')



pos=length(xx4);
  plot(Caa4(pos),xx4(pos),'.','Color',[0. 0.0 1],'MarkerSize',50)
  hold on
  plot(Caa1(pos),xx1(pos),'.','Color',[0. 0.99 0.],'MarkerSize',50)
  hold on
plot(Caa4(pos),xx4(pos),'.','Color',[0.5 0.7 1],'MarkerSize',40)
  hold on
  plot(Caa1(pos),xx1(pos),'.','Color',[0. 0.7 0.],'MarkerSize',40)
  hold on
   plot(Caa4(pos),xx4(pos),'.','Color',[0.1 0.1 .5],'MarkerSize',20)
  hold on
  plot(Caa1(pos),xx1(pos),'.','Color',[0. 0.5 0.],'MarkerSize',20)
  hold on 
  

dca=-50;
dx=-4;
% 
%x4=0.1:0.01:.999;
%V4=dx-50-log(1./x4-1)/0.15;
 V4=-68:.1:-30;
x4=1./(1+exp(-0.15*(V4+50-dx)));
% 
gsyn=0.000;
Iapp=+0.00;

VI=30;
VK=-75;
Vs=127*V4/(VI-VK)-(115*VK+VI*12)/(VI-VK); 
amm=0.1*(50-Vs)./(exp((50-Vs)/10)-1);
bmm=4*exp((25-Vs)/18);                                                                                                                                         
ah=0.07*exp((25-Vs)/20);   
bh=1.0./(1+exp((55-Vs)/10));          
an=0.01*(55-Vs)./(exp((55-Vs)/10)-1);
bn=0.125*exp((45-Vs)/80);  
minf=amm./(amm+bmm);  
hinf=ah./(ah+bh); 
ninf=an./(an+bn);  
tauh=12.5./(ah+bh);  
taun=12.5./(an+bn);                                                                                                                                             
taux=100;                                                                                                                                                                                                                                                                                                                                                                                                     
h4=hinf;                                                                                                                                             
n4=ninf;  
Vhh=-50;
y1=1./(1+exp(10*(V4-Vhh)));   

%Jack 
% r = 127*V4/105 + 8265/105;
% amm=0.1*(50-r)./(exp((50-r)./10)-1);
% bmm=4*exp((25-r)./18);
% minfty = amm./(amm+bmm);

% TTX current
gT=0.01;
% little change from -0.15 as above
xinf=1./(1+exp(-0.15*(V4+50-dx)));  
Ittx=gT*xinf.*(V4-VI);

% leak current 
gL=0.003; VL=-40; 
Ileak=gL*(V4-VL);

%Calcium current 
rho=0.0003; gKCa=0.03; VCa=140;
Ca4=0.0085*x4.*(VCa-V4+dca);
ICa=gKCa*Ca4.*(V4-VK)./(0.5+Ca4);

%very small conntributions 
% Potassium currebt
gK=0.3; 
IK=gK*n4.^4.*(V4-VK);
%Sodium current
gI=4;
INa=-gI*minf.^3.*h4.*(30-V4);

% h-current
gh=0.0000;
Ihh=gh./((1.+exp(-(V4+63.)/7.8))).^3.*y1.*(V4-70);


%all currents  
Q1 =-gsyn*(V4-20) -Iapp -INa  -IK -Ittx -Ileak -Ihh;
 
% solve for Q1 = -0.03*Ca4./(Ca4+0.5).*(V4+75)
Ca1=0.5*Q1./(0.03*(V4-VK)-Q1);

%x22=Ca4./(0.0085*(140-V4+dca));
plot (Ca1,x4,'Color',[0 0.99 0.0],'LineWidth',1)
hold on

% nullcline Ca'=0
Q11=-gsyn*(V4-20)-Iapp+gI*minf.^3.*h4.*(VI-V4) + gK*n4.^4.*(VK-V4) +0.03*Ca4.*(V4+75)./(0.5+Ca4) -gL*(V4-VL) -gh*(((1./(1.+exp(-(V4+63.)/7.8))).^3).*y1.*(V4-70));
x41=Q11./(gT*(VI-V4));
Ca4=0.0085*x41.*(140-V4+dca);
plot (Ca4,x41,'r')
hold on



dca=-55;
dx=-4;
% 
%x4=0.1:0.01:.999;
%V4=dx-50-log(1./x4-1)/0.15;
 V4=-68:1:-26;
x4=1./(1+exp(-0.15*(V4+50-dx)));
% 
gsyn=0.009;
Iapp=+0.00;

VI=30;
VK=-75;
Vs=127*V4/(VI-VK)-(115*VK+VI*12)/(VI-VK); 
amm=0.1*(50-Vs)./(exp((50-Vs)/10)-1);
bmm=4*exp((25-Vs)/18);                                                                                                                                         
ah=0.07*exp((25-Vs)/20);   
bh=1.0./(1+exp((55-Vs)/10));          
an=0.01*(55-Vs)./(exp((55-Vs)/10)-1);
bn=0.125*exp((45-Vs)/80);  
minf=amm./(amm+bmm);  
hinf=ah./(ah+bh); 
ninf=an./(an+bn);  
tauh=12.5./(ah+bh);  
taun=12.5./(an+bn);                                                                                                                                             
taux=100;                                                                                                                                                                                                                                                                                                                                                                                                     
h4=hinf;                                                                                                                                             
n4=ninf;  
Vhh=-50;
y1=1./(1+exp(10*(V4-Vhh)));   

%Jack 
% r = 127*V4/105 + 8265/105;
% amm=0.1*(50-r)./(exp((50-r)./10)-1);
% bmm=4*exp((25-r)./18);
% minfty = amm./(amm+bmm);

% TTX current
gT=0.01;
% little change from -0.15 as above
xinf=1./(1+exp(-0.15*(V4+50-dx)));  
Ittx=gT*xinf.*(V4-VI);

% leak current 
gL=0.003; VL=-40; 
Ileak=gL*(V4-VL);

%Calcium current 
rho=0.0003; gKCa=0.03; VCa=140;
Ca4=0.0085*x4.*(VCa-V4+dca);
ICa=gKCa*Ca4.*(V4-VK)./(0.5+Ca4);

%very small conntributions 
% Potassium currebt
gK=0.3; 
IK=gK*n4.^4.*(V4-VK);
%Sodium current
gI=4;
INa=-gI*minf.^3.*h4.*(30-V4);

% h-current
gh=0.0000;
Ihh=gh./((1.+exp(-(V4+63.)/7.8))).^3.*y1.*(V4-70);


%all currents  
Q1 =-gsyn*(V4+70) -Iapp -INa  -IK -Ittx -Ileak -Ihh;
 
% solve for Q1 = -0.03*Ca4./(Ca4+0.5).*(V4+75)
Ca1=0.5*Q1./(0.03*(V4-VK)-Q1);

%x22=Ca4./(0.0085*(140-V4+dca));
plot (Ca1,x4,'Color',[0 0.1 0.99],'LineWidth',1)
hold on

% nullcline Ca'=0
Q11=-gsyn*(V4+70)-Iapp+gI*minf.^3.*h4.*(VI-V4) + gK*n4.^4.*(VK-V4) +0.03*Ca4.*(V4+75)./(0.5+Ca4) -gL*(V4-VL) -gh*(((1./(1.+exp(-(V4+63.)/7.8))).^3).*y1.*(V4-70));
x41=Q11./(gT*(VI-V4));
Ca4=0.0085*x41.*(140-V4+dca);
plot (Ca4,x41,'r')
hold on

%print(gcf,'-djpeg','-r600' ,'fredo3a.jpeg');
