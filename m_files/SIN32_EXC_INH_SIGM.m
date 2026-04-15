clear all 

tmax=50; % max time in seconds

%--------- green, blue: ts ---------------
alpha3 = 0.0175;
beta3 =  0.001;
scale3=alpha3/(alpha3+beta3)

alpha2 = 0.01;
beta2 =  0.0012;
scale2=alpha2/(alpha2+beta2)

g32=.17;
g23=.1;
gelec=0.0001;


%----------------------
cutoff=16.0; % cutoff frequency in Hz
noise=.00   ;
%----------------------

%Parameters3
Ca_shift3 = -50;
Ca_shift2 = -20;

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
V3= -44; V2= -44; Ca3=.6; Ca2=1.3; h3 =0; h2 =0; n3 =0; n2 =0; 
x3 =.8; x2 =.92; s3 =0; s2 =0; 

tic
[rnd1,nt] = bandlimnoise (cut1,step,tf);
[rnd2,nt] = bandlimnoise (cut1,step,tf);
rnd1=rnd1*noise; rnd2=rnd2*noise;

nt1=round(tmax/tsamp)+1;
is=round(tsamp/step);

time=zeros(nt1,1);vv3=time;vv2=time;ss3=time;ss2=time;Caa3=time;
Caa2=time;xx3=time;xx2=time;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% CALCULATE heaviside functions in advance! This speeds up calculation!

heaV3=zeros(nt,1); heav2=heaV3;
j1=round(t1/step)+1; j2=round(t2/step)+1;
heaV3(j1:nt)=1; heav2(j2:nt)=1; heav=heaV3-heav2; 
clear heaV3 heav2;

j=0;
for i=1:nt
 tt=(i-1)*step;
V3 =V3 +step*(4*((0.1*(50-(127*V3/105+8265/105))/(exp((50 - (127*V3/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V3/105 + 8265/105))/(exp((50 - (127*V3/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V3/105 + 8265/105))/18))))^3*h3*(30 - V3) + 0.3*n3^4*(-75 - V3)+0.01*x3*(30-V3) ...
    +0.03*Ca3/(.5 + Ca3)*(-75 - V3)+0.003*(-40 - V3) -g23*(V3+80)*s2/scale2+gelec*(V2-V3)+rnd1(i));
V2= V2 +step*(4*((0.1*(50-(127*V2/105+8265/105))/(exp((50 - (127*V2/105 ...
    +8265/105))/10) - 1))/((0.1*(50 - (127*V2/105 + 8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V2/105 + 8265/105))/18))))^3*h2*(30 - V2) + 0.3*n2^4*(-75 - V2)+0.01*x2*(30-V2) ...
    +0.03*Ca2/(.5 + Ca2)*(-75 - V2)+0.003*(-40 - V2) -g32*(V2-40)*s3/scale3+gelec*(V3-V2)+rnd2(i));
Ca3=Ca3+step*(0.0003*(0.0085*x3*(140-V3+Ca_shift3)-Ca3));
Ca2=Ca2+step*(0.0003*(0.0085*x2*(140-V2+Ca_shift2)-Ca2));
x3 =x3+step*(((1/(exp(0.15*(-V3-50+x_shift))+1))-x3)/100);
x2 =x2+step*(((1/(exp(0.15*(-V2-50+x_shift))+1))-x2)/100);
h3 =h3+step*(((1-h3)*(0.07*exp((25 - (127*V3/105 + 8265/105))/20))-h3*(1.0/(1 + exp((55 - (127*V3/105 + 8265/105))/10))))/12.5);
h2 =h2+step*(((1-h2)*(0.07*exp((25 - (127*V2/105 + 8265/105))/20))-h2*(1.0/(1 + exp((55 - (127*V2/105 + 8265/105))/10))))/12.5);
n3 =n3+step*(((1-n3)*(0.01*(55 - (127*V3/105 + 8265/105))/(exp((55 - (127*V3/105 + 8265/105))/10) - 1))-n3*(0.125*exp((45 - (127*V3/105 + 8265/105))/80)))/12.5);
n2 =n2+step*(((1-n2)*(0.01*(55 - (127*V2/105 + 8265/105))/(exp((55 - (127*V2/105 + 8265/105))/10) - 1))-n2*(0.125*exp((45 - (127*V2/105 + 8265/105))/80)))/12.5);
s3 =s3+step*(alpha3*s3*(1-s3)/(1+exp(-20*(V3+20)))-beta3*(s3-0.0001));
s2 =s2+step*(alpha2*s2*(1-s2)/(1+exp(-20*(V2+20)))-beta2*(s2-0.0001));


% Do not need to save every point!
if mod(i,is) ==0
    j=j+1;
    time(j)=tt; vv3(j)=V3; vv2(j)=V2; ss3(j)=s3;
    ss2(j)=s2; Caa3(j)=Ca3;
    Caa2(j)=Ca2; xx3(j)=x3; xx2(j)=x2; heavy(j)=heav(i);
end
end  
toc    

%  
%  [vp1,tp1] = findpeaks(vV3,time,'MinPeakHeight',-20);
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

 plot(time,vv3,'Color',[0 0  .6]','LineWidth',1.5)
 hold on
%  plot(time,heavy*10,'Color',[1 0  1]','LineWidth',1.5)
%  hold on
 xlim([0 tmax]) 
 ylim([-80 40])
  fftt = 28;
txpo = text(-4.4,25,'B','Fontsize',fftt,'Color','black','FontName','Arial','FontWeight','bold')
   
   xlabel('Time [sec]','Fontsize', 16),ylabel('Voltage','Fontsize', 16)

 subplot(4,1,2)
 plot(time,vv2,'Color',[0 .5 0]','LineWidth',1.5)
 hold on
  xlabel('Time [sec]','Fontsize', 16),ylabel('Voltage','Fontsize', 16)
 xlim([0 tmax]) 
 ylim([-65 40])
 
  %xlim([350 time(end)])

  subplot(4,1,3)
 % plot(time,10*mm1.*ss3,'blue')
  hold on
  plot(time,ss3/scale3,'Color',[0.1 0.1 0.6],'LineWidth',1.5)
  plot(time,ss2/scale2,'Color',[0.1 0.6 0.1],'LineWidth',1.5)
  hold on
  ylim([-0.02  .1])
  xlim([0 tmax])
  box on
   xlabel('Time [sec]','Fontsize', 16),ylabel('Release rate','Fontsize', 16)
  
subplot(4,1,4)
  plot(time,Caa3,'Color',[0. 0.0 0.8],'Linewidth',1.5)
  hold on
  plot(time,Caa2,'Color',[0.0 0.7 0.0],'Linewidth',1.5)
  hold on
  %plot(time,xx3,'Color',[0 .1 1])
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
 plot(Caa3(st:end),xx3(st:end),'Color',[0. 0. 0.6],'LineWidth',2)
  hold on
  plot(Caa2(st:end),xx2(st:end),'Color',[0.1 0.6 0.1],'LineWidth',2)
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



pos=length(xx3);
  plot(Caa3(pos),xx3(pos),'.','Color',[0. 0.0 1],'MarkerSize',50)
  hold on
  plot(Caa2(pos),xx2(pos),'.','Color',[0. 0.99 0.],'MarkerSize',50)
  hold on
plot(Caa3(pos),xx3(pos),'.','Color',[0.5 0.7 1],'MarkerSize',40)
  hold on
  plot(Caa2(pos),xx2(pos),'.','Color',[0. 0.7 0.],'MarkerSize',40)
  hold on
   plot(Caa3(pos),xx3(pos),'.','Color',[0.1 0.1 .5],'MarkerSize',20)
  hold on
  plot(Caa2(pos),xx2(pos),'.','Color',[0. 0.5 0.],'MarkerSize',20)
  hold on 
  

dca=-20;
dx=-4;
% 
%x3=0.1:0.01:.999;
%V3=dx-50-log(1./x3-1)/0.15;
 V3=-68:.1:-30;
x3=1./(1+exp(-0.15*(V3+50-dx)));
% 
gsyn=0.00025;
Iapp=+0.00;

VI=30;
VK=-75;
Vs=127*V3/(VI-VK)-(115*VK+VI*12)/(VI-VK); 
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
h3=hinf;                                                                                                                                             
n3=ninf;  
Vhh=-50;
y1=1./(1+exp(10*(V3-Vhh)));   

%Jack 
% r = 127*V3/105 + 8265/105;
% amm=0.1*(50-r)./(exp((50-r)./10)-1);
% bmm=4*exp((25-r)./18);
% minfty = amm./(amm+bmm);

% TTX current
gT=0.01;
% little change from -0.15 as above
xinf=1./(1+exp(-0.15*(V3+50-dx)));  
Ittx=gT*xinf.*(V3-VI);

% leak current 
gL=0.003; VL=-40; 
Ileak=gL*(V3-VL);

%Calcium current 
rho=0.0003; gKCa=0.03; VCa=140;
Ca3=0.0085*x3.*(VCa-V3+dca);
ICa=gKCa*Ca3.*(V3-VK)./(0.5+Ca3);

%very small conntributions 
% Potassium currebt
gK=0.3; 
IK=gK*n3.^4.*(V3-VK);
%Sodium current
gI=4;
INa=-gI*minf.^3.*h3.*(30-V3);

% h-current
gh=0.0000;
Ihh=gh./((1.+exp(-(V3+63.)/7.8))).^3.*y1.*(V3-70);


%all currents  
Q1 =-gsyn*(V3-20) -Iapp -INa  -IK -Ittx -Ileak -Ihh;
 
% solve for Q1 = -0.03*Ca3./(Ca3+0.5).*(V3+75)
Ca2=0.5*Q1./(0.03*(V3-VK)-Q1);

%x22=Ca3./(0.0085*(140-V3+dca));
plot (Ca2,x3,'Color',[0 0.99 0.0],'LineWidth',1)
hold on

% nullcline Ca'=0
Q11=-gsyn*(V3-20)-Iapp+gI*minf.^3.*h3.*(VI-V3) + gK*n3.^4.*(VK-V3) +0.03*Ca3.*(V3+75)./(0.5+Ca3) -gL*(V3-VL) -gh*(((1./(1.+exp(-(V3+63.)/7.8))).^3).*y1.*(V3-70));
x31=Q11./(gT*(VI-V3));
Ca3=0.0085*x31.*(140-V3+dca);
plot (Ca3,x31,'r')
hold on



dca=-55;
dx=-4;
% 
%x3=0.1:0.01:.999;
%V3=dx-50-log(1./x3-1)/0.15;
 V3=-68:1:-26;
x3=1./(1+exp(-0.15*(V3+50-dx)));
% 
gsyn=0.009;
Iapp=+0.00;

VI=30;
VK=-75;
Vs=127*V3/(VI-VK)-(115*VK+VI*12)/(VI-VK); 
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
h3=hinf;                                                                                                                                             
n3=ninf;  
Vhh=-50;
y1=1./(1+exp(10*(V3-Vhh)));   

%Jack 
% r = 127*V3/105 + 8265/105;
% amm=0.1*(50-r)./(exp((50-r)./10)-1);
% bmm=4*exp((25-r)./18);
% minfty = amm./(amm+bmm);

% TTX current
gT=0.01;
% little change from -0.15 as above
xinf=1./(1+exp(-0.15*(V3+50-dx)));  
Ittx=gT*xinf.*(V3-VI);

% leak current 
gL=0.003; VL=-40; 
Ileak=gL*(V3-VL);

%Calcium current 
rho=0.0003; gKCa=0.03; VCa=140;
Ca3=0.0085*x3.*(VCa-V3+dca);
ICa=gKCa*Ca3.*(V3-VK)./(0.5+Ca3);

%very small conntributions 
% Potassium currebt
gK=0.3; 
IK=gK*n3.^4.*(V3-VK);
%Sodium current
gI=4;
INa=-gI*minf.^3.*h3.*(30-V3);

% h-current
gh=0.0000;
Ihh=gh./((1.+exp(-(V3+63.)/7.8))).^3.*y1.*(V3-70);


%all currents  
Q1 =-gsyn*(V3+70) -Iapp -INa  -IK -Ittx -Ileak -Ihh;
 
% solve for Q1 = -0.03*Ca3./(Ca3+0.5).*(V3+75)
Ca2=0.5*Q1./(0.03*(V3-VK)-Q1);

%x22=Ca3./(0.0085*(140-V3+dca));
plot (Ca2,x3,'Color',[0 0.1 0.99],'LineWidth',1)
hold on

% nullcline Ca'=0
Q11=-gsyn*(V3+70)-Iapp+gI*minf.^3.*h3.*(VI-V3) + gK*n3.^4.*(VK-V3) +0.03*Ca3.*(V3+75)./(0.5+Ca3) -gL*(V3-VL) -gh*(((1./(1.+exp(-(V3+63.)/7.8))).^3).*y1.*(V3-70));
x31=Q11./(gT*(VI-V3));
Ca3=0.0085*x31.*(140-V3+dca);
plot (Ca3,x31,'r')
hold on

%print(gcf,'-djpeg','-r600' ,'fredo3a.jpeg');
