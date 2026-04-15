
%clear all 
 table1=xlsread('Dendronotus whole swim.xlsx');


time_1=table1(:,1);
l_si1_1=table1(:,2);
%save('l_si1_1.abf','l_si1_1','-ascii')
l_si2_1=table1(:,3);
r_si3_1=table1(:,4);
% plot (time_1-t1,smoothdata(l_si1_1,'gaussian'),'Color',[1. 0. .0])
% hold on
%plot(time_1-t1,l_si2_1-85, 'Color',[0 0 .5]);
hold on
%plot(time_1-t1,r_si3_1-185, 'Color',[1 0 0]);
%hold on
%axis([0 280 -240 40])





%alpha synape 
%slow m*s
%fast
%beta =  0.01;
%taum=100;
%slow
%taum=1500;


noise=.0;

%Parameters3

Ca_shift1=  -30.8;
Ca_shift2 = -0.6;
x_shift =  -3.;



alpha = 0.002;
beta =  0.0001;
g12=0.0;
g21=0.0;
gelec=.0;
taum=500;

g13=0.01;
alphae=.005;
betae=.002;


% Constant stimuli if amy
Iapp=.008; %1Hz
Iapp=.022; %  2Hz
%Iapp=.115; % 5Hz
%Iapp=.655; % 10Hz
Iapp=0.06;
Iapp=0.001;

t1=10*1000;
t2=230*1000;
      
% H-current      
  % H-current      
    gh    = 0.0000;
    Vhh   = -53;

% Intergration    
t_final=time_1(end)-55   
step=.5;
step=1000*time_1(end)/length(time_1)

% Initial values
V1= 0; V2= 0; Ca1=1.; Ca2=.99; h1 =0; h2 =0; n1 =0; n2 =0; x1 =0.8;
x2 =0.85; y1 =0; y2 =0; s1 =0; s2 =0; m1 =0; m2 =0;
i=1;
tt=1;
se=0;

clear ss1 mm1  clear time vv1 vv2 see ss2 mm2 clear Caa1 Caa2 xx1 xx2

time=zeros(length(time_1),1);vv1=time;vv2=time;ss1=time;ss2=time;mm1=time;mm2=time;Caa1=time;
Caa2=time;xx1=time;xx2=time;

%V3=l_si1_1;
V3=l_si2_1;
tic
for i=1:length(time_1) 
   
    


se=se+step*(alphae*(1-se)/(1+exp(-10*(V3(i)+10)))-betae*se);
see(i)=se;

V1 =V1 +step*(4*((0.1*(50-(127*V1/105+8265/105))/(exp((50 - ...
    (127*V1/105 + 8265/105))/10) - 1))/((0.1*(50 - (127*V1/105 + 8265/105))/(exp((50 - (127*V1/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V1/105 + 8265/105))/18))))^3*h1*(30 - V1) + 0.3*n1^4*(-75 - V1)+0.01*x1*(30-V1) +0.03*Ca1/(.5 + Ca1)*(-75 - V1)...
    +0.003*(-40 - V1)   +gh*((1/(1+exp(-(V1+63)/7.8)))^3)*y1*(+120-V1)...
    -g13*(V1-55)*se-g21*(V1+80)*s2*m2+gelec*(V3(i)-V1)+noise*(rand-0.5));
V2= V2 +step*(4*((0.1*(50-(127*V2/105+8265/105))/(exp((50 - ...
    (127*V2/105 + 8265/105))/10) - 1))/((0.1*(50 - (127*V2/105 + 8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V2/105 + 8265/105))/18))))^3*h2*(30 - V2) + 0.3*n2^4*(-75 - V2)+0.01*x2*(30-V2) +0.03*Ca2/(.5 + Ca2)*(-75 - V2)...
    +0.003*(-40 - V2) +gh*((1/(1+exp(-(V2+63)/7.8)))^3)*y2*(-V2+120)...
    -g13*(V2-55)*se-g12*(V2+80)*s1^2*m1+gelec*(V3(i)-V2)+noise*(rand-0.5));

Ca1=Ca1+step*(0.0001*(0.0085*x1*(140-V1+Ca_shift1)-Ca1));
Ca2=Ca2+step*(0.0001*(0.0085*x2*(140-V2+Ca_shift2)-Ca2));
x1 =x1+step*(((1/(exp(0.15*(-V1-50+x_shift))+1))-x1)/235);
x2 =x2+step*(((1/(exp(0.15*(-V2-50+x_shift))+1))-x2)/235);
h1 =h1+step*(((1-h1)*(0.07*exp((25 - (127*V1/105 + 8265/105))/20))-h1*(1.0/(1 + exp((55 - (127*V1/105 + 8265/105))/10))))/12.5);
h2 =h2+step*(((1-h2)*(0.07*exp((25 - (127*V2/105 + 8265/105))/20))-h2*(1.0/(1 + exp((55 - (127*V2/105 + 8265/105))/10))))/12.5);
n1 =n1+step*(((1-n1)*(0.01*(55 - (127*V1/105 + 8265/105))/(exp((55 - (127*V1/105 + 8265/105))/10) - 1))-n1*(0.125*exp((45 - (127*V1/105 + 8265/105))/80)))/12.5);
n2 =n2+step*(((1-n2)*(0.01*(55 - (127*V2/105 + 8265/105))/(exp((55 - (127*V2/105 + 8265/105))/10) - 1))-n2*(0.125*exp((45 - (127*V2/105 + 8265/105))/80)))/12.5);
y1 =y1+step*(.5*((1/(1+exp(10*(V1-Vhh))))-y1)/(7.1+10.4/(1+exp((V1+68)/2.2))));
y2 =y2+step*(.5*((1/(1+exp(10*(V2-Vhh))))-y2)/(7.1+10.4/(1+exp((V2+68)/2.2))));

s1 =s1+step*(alpha*(1-s1)/(1+exp(-10*(V1+30)))-beta*s1);
s2 =s2+step*(alpha*(1-s2)/(1+exp(-10*(V2+30)))-beta*s2);
m1 =m1+step*((1/(1+exp(-(V1+30)))-m1)/taum);
m2 =m2+step*((1/(1+exp(-(V2+30)))-m2)/taum);

time(i)=(i-1)*step;
vv1(i)=V1;
vv2(i)=V2;
ss1(i)=s1;
ss2(i)=s2;
mm1(i)=m1;
mm2(i)=m2;
Caa1(i)=Ca1;
Caa2(i)=Ca2;
xx1(i)=x1;
xx2(i)=x2;
end  
toc    
    
size(time)

%  
%  [vp1,tp1] = findpeaks(vv1,time,'MinPeakHeight',0);
%  fr1=length(tp1)/(tp1(end)-tp1(1))*1000
%  tp1=tp1/1000;


time=time/1000;
t0=25;
%     
fig = figure(1);
clf 

subplot(3,1,1)
left_color = [0  0 0];
right_color = [0 0 .0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
right_color = [0 0 .0];
 
yyaxis left 
plot(time_1-t0,V3, 'Color',[0.5 0.5 .5]);
hold on
%plot(time_1-t0,vv2, 'Color',[0.5 0.0 .0]);
hold on
xlim([0 t_final]) 
ylim([-55 45]) 
xlabel('Time [sec]','FontSize',14),ylabel(' Si1R V[mV]','FontSize',14)
 
yyaxis right 
plot(time-t0,see,'Color',[0 0 0])
hold on
xlim([0 t_final]) 
ylim([-.01 1])
ylabel('probability s(t)','FontSize',14) 

 left_color = [0  0 0];
right_color = [.0 0.0 .0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

subplot(3,1,2)
left_color = [0  0 0];
right_color = [.0 .0 .0];
 
yyaxis left 
plot(time(1:end)-t0,vv1(1:end),'Color',[0 0  1],'LineWidth',1)
 hold on
 plot(time_1-t0,V3, 'Color',[0.5 0.5 .5]);
hold on
 xlim([0 t_final]) 
 ylim([-80 40])
 %xlabel('Time [sec]','FontSize',14)
 ylabel('Math neuron V_1[mV]','FontSize',14)
 
 yyaxis right 
 plot(time(1:end)-t0,mm1(1:end).*ss1(1:end),'Color',[0 0  .51])
 hold on
 
 xlim([0 t_final]) 
 ylim([-.01 1])
 ylabel('probability s(t)','FontSize',14)

set(fig,'defaultAxesColorOrder',[left_color; right_color]);

  
  
 subplot(3,1,3)
  
 yyaxis left 
 plot(time-t0,vv2,'Color',[0 1 0],'LineWidth',1)
 hold on
 plot(time_1-t0,V3, 'Color',[0.5 0.5 .5]);
hold on

 xlim([0 t_final]) 
 ylim([-80 45])
 %xlabel('Time [sec]','FontSize',14)
 ylabel('Math neuron V_2[mV]','FontSize',14)
 
 yyaxis right 
 plot(time(1:end)-t0,mm2(1:end).*ss2(1:end),'Color',[0 0.5  0])
  hold on
 xlim([0 t_final]) 
 ylim([-.01 1])
 ylabel('probability s(t)','FontSize',14)
  


 