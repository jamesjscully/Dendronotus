
clear all 

tmax=50; % max time in seconds


%--------- green, blue: ts ---------------
alpha4 = 0.0175;
beta4 =  0.001;
scale4=alpha4/(alpha4+beta4)

alpha1 = 0.01;
beta1 =  0.0012;
scale1=alpha1/(alpha1+beta1)

g41=.17;
g14=.1;

gelec=0.0001;


%----------------------
cutoff=16.0; % cutoff frequency in Hz
noise=.00   ;
%----------------------

%Parameters4
Ca_shift4 = -50;
Ca_shift1 = -20;

x_shift =  -4;

Iapp=0.0;

t1=00*1000; % this time is in msec
t2=00*1000; % this time is in msec
      
% H-current      
    gh    = 0.000;
    Vhh   = -55;

% Intergration    
tf=tmax*1000;   % max time in msec   
step=0.1;       % time step in msec
tsamp=2; % sampling interval to save data in msec



cut1=cutoff/1000; % cutoff frequency converted to msec
% Initial values
V4= -44; V1= -44; Ca4=.6; Ca1=1.3; h4 =0; h1 =0; n4 =0; n1 =0; 
x4 =.8; x1 =.92; s4 =0; s1 =0; 

tic
[rnd1,nt] = bandlimnoise (cut1,step,tf);
[rnd2,nt] = bandlimnoise (cut1,step,tf);
rnd1=rnd1*noise; rnd2=rnd2*noise;

nt1=round(tmax/tsamp)+1;
is=round(tsamp/step);

time=zeros(nt1,1);vv4=time;vv1=time;ss4=time;ss1=time;Caa4=time;
Caa1=time;xx4=time;xx1=time;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% CALCULATE heaviside functions in advance! This speeds up calculation!

heaV4=zeros(nt,1); heaV1=heaV4;
j1=round(t1/step)+1; j2=round(t2/step)+1;
heaV4(j1:nt)=1; heaV1(j2:nt)=1; heav=heaV4-heaV1; 
clear heaV4 heaV1;

j=0;
for i=1:nt
 tt=(i-1)*step;
V4 =V4 +step*(4*((0.1*(50-(127*V4/105+8265/105))/(exp((50 - (127*V4/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V4/105 + 8265/105))/(exp((50 - (127*V4/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V4/105 + 8265/105))/18))))^3*h4*(30 - V4) + 0.3*n4^4*(-75 - V4)+0.01*x4*(30-V4) ...
    +0.03*Ca4/(.5 + Ca4)*(-75 - V4)+0.003*(-40 - V4) -g14*(V4+80)*s1/scale4+gelec*(V1-V4)+rnd1(i));
V1= V1 +step*(4*((0.1*(50-(127*V1/105+8265/105))/(exp((50 - (127*V1/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V1/105 + 8265/105))/(exp((50 - (127*V1/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V1/105 + 8265/105))/18))))^3*h1*(30 - V1) + 0.3*n1^4*(-75 - V1)+0.01*x1*(30-V1) ...
    +0.03*Ca1/(.5 + Ca1)*(-75 - V1)+0.003*(-40 - V1) -g41*(V1-40)*s4/scale4+gelec*(V4-V1)+rnd2(i));
Ca4=Ca4+step*(0.0003*(0.0085*x4*(140-V4+Ca_shift4)-Ca4));
Ca1=Ca1+step*(0.0003*(0.0085*x1*(140-V1+Ca_shift1)-Ca1));
x4 =x4+step*(((1/(exp(0.15*(-V4-50+x_shift))+1))-x4)/100);
x1 =x1+step*(((1/(exp(0.15*(-V1-50+x_shift))+1))-x1)/100);
h4 =h4+step*(((1-h4)*(0.07*exp((25 - (127*V4/105 + 8265/105))/20))-h4*(1.0/(1 + exp((55 - (127*V4/105 + 8265/105))/10))))/12.5);
h1 =h1+step*(((1-h1)*(0.07*exp((25 - (127*V1/105 + 8265/105))/20))-h1*(1.0/(1 + exp((55 - (127*V1/105 + 8265/105))/10))))/12.5);
n4 =n4+step*(((1-n4)*(0.01*(55 - (127*V4/105 + 8265/105))/(exp((55 - (127*V4/105 + 8265/105))/10) - 1))-n4*(0.125*exp((45 - (127*V4/105 + 8265/105))/80)))/12.5);
n1 =n1+step*(((1-n1)*(0.01*(55 - (127*V1/105 + 8265/105))/(exp((55 - (127*V1/105 + 8265/105))/10) - 1))-n1*(0.125*exp((45 - (127*V1/105 + 8265/105))/80)))/12.5);
s4 =s4+step*(alpha4*s4*(1-s4)/(1+exp(-20*(V4+20)))-beta4*(s4-0.0001));
s1 =s1+step*(alpha1*s1*(1-s1)/(1+exp(-20*(V1+20)))-beta1*(s1-0.0001));


% Do not need to save every point!
if mod(i,is) ==0
    j=j+1;
    time(j)=tt; vv4(j)=V4; vv1(j)=V1; ss4(j)=s4;
    ss1(j)=s1; Caa4(j)=Ca4;
    Caa1(j)=Ca1; xx4(j)=x4; xx1(j)=x1; heavy(j)=heav(i);
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
 subplot(4,1,1)

 plot(time,vv4,'Color',[0 0  .6]','LineWidth',1.5)
 hold on
%  plot(time,heavy*10,'Color',[1 0  1]','LineWidth',1.5)
%  hold on
 xlim([0 tmax]) 
 ylim([-80 40])
  fftt = 28;
txpo = text(-4.4,25,'B','Fontsize',fftt,'Color','black','FontName','Arial','FontWeight','bold')
   
   xlabel('Time [sec]','Fontsize', 16),ylabel('Voltage','Fontsize', 16)

 subplot(4,1,2)
 plot(time,vv1,'Color',[0 .5 0]','LineWidth',1.5)
 hold on
  xlabel('Time [sec]','Fontsize', 16),ylabel('Voltage','Fontsize', 16)
 xlim([0 tmax]) 
 ylim([-65 40])
 
  %xlim([350 time(end)])

  subplot(4,1,3)
 % plot(time,10*mm1.*ss4,'blue')
  hold on
  plot(time,ss4/scale4,'Color',[0.1 0.1 0.6],'LineWidth',1.5)
  plot(time,ss1/scale1,'Color',[0.1 0.6 0.1],'LineWidth',1.5)
  hold on
  ylim([-0.02  .1])
  xlim([0 tmax])
  box on
   xlabel('Time [sec]','Fontsize', 16),ylabel('Release rate','Fontsize', 16)
  
subplot(4,1,4)
  plot(time,Caa4,'Color',[0. 0.0 0.8],'Linewidth',1.5)
  hold on
  plot(time,Caa1,'Color',[0.0 0.7 0.0],'Linewidth',1.5)
  hold on
  %plot(time,xx4,'Color',[0 .1 1])
  hold on
  %plot(time,xx2,'Color',[0 1 .1])
  hold on
  ylim([0.3 1.4])
  xlim([0 tmax])
 xlabel('Time [sec]','Fontsize', 16),ylabel('[Ca]','Fontsize', 16)
  %title('Ca and Ca-activaterd possasion gating varaible in time', 'Fontsize', 11);
  
  %print(gcf,'-djpeg','-r600' ,'fredo3b.jpeg');

  toc 
   figure(5)
 clf
 st=1;
 plot(Caa4(st:end),xx4(st:end),'Color',[0. 0. 0.6],'LineWidth',2)
  hold on
  plot(Caa1(st:end),xx1(st:end),'Color',[0.1 0.6 0.1],'LineWidth',2)
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

axis([0.2 1.4 0.01 1.0])

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
  

dca=-20;
dx=-4;
% 
%x4=0.1:0.01:.999;
%V4=dx-50-log(1./x4-1)/0.15;
 V4=-68:.1:-30;
x4=1./(1+exp(-0.15*(V4+50-dx)));
% 
gsyn=0.00025;
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
