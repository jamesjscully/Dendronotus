% Model 
% V : Membrane Potential(mV)
% t : Time (ms)33
% C_m : Membrane capacitance(pico farad/cm^2)
% g_i : Maximum conductance 
% V_i : Nerst reversal potential(mV) 
% x, n : dimensionless ativation variables
% h : Dimensionless inactivation variable
% m_inf : Instantenous activation function for the fast inward current
% Ca : The dimensionless calcium ion concentration at the innner cell
function Melibe_Fig


global alpha beta alpha1 beta1 g12 g21 gelec Iapp Iapp1 Iapp2 t1 t2  t3 t4  t5 t6 Ca_shift x_shift gh Vhh

beta = 0.001;
alpha= beta;

alpha1 = 0.012;
beta1 =  alpha1/10;

g12=.02;
g21=0;
gelec=0.0002;


Ca_shift = -20;
x_shift =  -4;

% Positive Constant stimuli
Iapp=.013; %1Hz
Iapp2=.04;
Iapp1=.15;

%Iapp=.03; %  2Hz
%Iapp=.115; % 5Hz
%Iapp=.655; % 10Hz


% Negative constant stimuli
%Iapp=-.025; %  

t1=10*1000;
t2=15*1000;

t5=25*1000;
t6=31*1000;

t3=42*1000;
t4=50*1000;


% H-current      
    gh    = 0.0000;
    Vhh   = -53;

% Initial Conditions
%      V Ca  h  n  x y s r m   
s0 = [-43.85 1.1 0 0 1 0   0 0 0 0 0 43.85 1.2 0. 0 1 0  ];
opts = odeset('Stats','off','RelTol',1e-5,'AbsTol',1e-5);
[t,s] = ode15s(@Plant_Model,[0:0.5:70*1000], s0,opts);
t=t/1000;
tmax=70; 
%  
% [vp1,tp1] = findpeaks(s(:,1),t,'MinPeakHeight',-20);
% fr1=length(tp1)/(tp1(end)-tp1(1))*1000
%     

clf 
 figure(4)
 %clf
 subplot(5,1,1)
 plot(t,s(:,1),'b','LineWidth',1.5)
 hold on
 ylim([-55 35]) 

 subplot(5,1,2)
 plot(t,s(:,12),'Color',[0 0.5 0],'LineWidth',1.5)
 hold on
 ylim([-50 -40]) 
 
 subplot(5,1,3)
 plot(t,s(:,7),'Color',[0 0 1],'LineWidth',1.5)
 hold on
 plot(t,s(:,8),'Color',[0 0 .5],'LineWidth',1.5)
 hold on
 ylabel('\alpha-synapse')
title(['alpha and \beta =  ', num2str(alpha), ' and ',  num2str(beta/2), ' and ', num2str(beta/4)] ); 
%title(['alpha and \beta =  ', num2str(alpha), num2str(beta/2)] ); 

 %xlim([25 60])   
 
 subplot(5,1,4)
 plot(t,s(:,9),'Color',[.6 .0 .0],'LineWidth',1.5)
 hold on
 plot(t,s(:,10),'Color',[.8 .5 .8],'LineWidth',1.5)
 hold on
 xlabel('Time'),ylabel('Voltage')
 %xlim([25 30]) 
 title(['alpha and \beta =  ', num2str(alpha), ' and ',  num2str(beta/10) , ' and ',  num2str(beta/20)] ); 

 subplot(5,1,5)
  plot(t,s(:,11),'red','LineWidth',1.5)
 hold on
%   plot(t,movavg(s(:,11)),'red','LineWidth',1.5)
%  hold on
title(['alpha and \beta =  ', num2str(alpha1),'and ', num2str(beta1) ] ); 
 xlabel('Time'),ylabel('Voltage')
 %xlim([25 60])   
 

 %% Second plot: Full 5-cell staggered view
fig=figure(44);
clf
set(gcf, 'Renderer', 'painters') ;           % Use the Painters renderer for smoother lines
set(0, 'DefaultLineLineSmoothing', 'on')    % Turn on line smoothing
set(gcf, 'Position',  [10, 300, 600, 350]); % Larger figure size for 4 cells
set(gcf, 'PaperPositionMode', 'auto');
hold on

% Define offset and scaling factors for staggered plot
offset_step = 110;  % Vertical offset between traces
synscale = 30;     % Base scaling factor for synaptic variables

% Start with highest offset for top trace (Si2L)
offset = 670;       % Starting from the top

% Plot Si2L (top trace)
plot(t,s(:,1) + offset-100,'b','Color', [0, 0, 0.6],'LineWidth',1.5)
 hold on
 
% Decrease offset for Si2R
offset = offset - offset_step;

% Plot Si2R trace
plot(t,s(:,12)*6+ offset+380,'Color',[0 0.0 0.9],'LineWidth',1.5)
 
% Decrease offset for Si3L
offset = offset - offset_step-30;
% Plot Si3L trace

 plot(t,s(:,7)*100 + offset,'Color',[0 0 0.66],'LineWidth',1.5)
 hold on
 plot(t,s(:,8)*120+ offset,'Color',[0 0 .5],'LineWidth',1.)
 hold on

% Decrease offset for Si3R
offset = offset - offset_step-20;

% Plot Si3R trace

% plot(t,s(:,9)*400+offset,'Color',[.0 .0 .6],'LineWidth',1.5)
%  hold on
%  plot(t,s(:,10)*400+offset,'Color',[.0 .0 .5],'LineWidth',1.5)
%  hold on
%  xlabel('Time'),ylabel('Voltage')
%  %xlim([25 30]) 
%  %title(['alpha and \beta =  ', num2str(alpha), ' and ',  num2str(beta/10) , ' and ',  num2str(beta/20)] ); 
% 
% % Decrease offset for Si1
% offset = offset - offset_step;

% Plot Si1R trace
plot(t,s(:,11)*100+ offset+20,'Color',[0.3 0.3 .4],'LineWidth',1.5)
 hold on

% Remove axis ticks
set(gca, 'YTick', [])   % Remove y-axis ticks
set(gca, 'XTick', [])   % Remove x-axis ticks

% % Add L‐shaped scale bar
x_start   = 31;
y_start   = 620;
yh_start   = y_start+5;
time_len  = 5;
volt_h    = 50;
%plot([x_start, x_start], [yh_start, yh_start + volt_h], 'k', 'LineWidth', 2);
plot([x_start, x_start - time_len], [y_start, y_start], 'k', 'LineWidth', 2);
text(x_start - time_len/2, y_start + 15, [num2str(time_len) ' s'], 'HorizontalAlignment', 'center', 'FontSize', 14);
%text(x_start + .3, yh_start + volt_h/2, [num2str(volt_h) ' mV'], 'VerticalAlignment', 'middle', 'FontSize', 14);


% Set axis limits and labels
xlim([5 tmax-12])     % Extended to show labels
ylim([300 720])         % Adjusted to fit all traces
xlabel('Time [s]', 'Fontsize', 20)
ylabel('  Rate      &     Voltage [mV]', 'Fontsize', 20)
box on

% End of stacked‐trace plotting
exportgraphics(fig, 'ZOO_synapses.jpg', 'Resolution', 75,'ContentType', 'image', 'BackgroundColor', 'white')


end

function [qd] = Plant_Model(t,q)
global alpha beta alpha1 beta1 g12 g21 gelec Iapp Iapp1 Iapp2 t1 t2 t3 t4 t5 t6 Ca_shift x_shift gh Vhh
  qd = zeros(17,1);
          
      V1 = q(1);                       % Fast
     Ca1 = q(2);                       % Slow
      h1 = q(3);                       % Fast 
      n1 = q(4);                       % Fast
      x1 = q(5);                       % Slow
      y1 = q(6);         
      s1 = q(7);
      s2 = q(8);
      s3 = q(9);
      s4= q(10); 
      sl = q(11);                       % Fast
     
      V2 = q(12);
      Ca2 =q(13);                       % Slow
      h2 = q(14);                       % Fast 
      n2 = q(15);                       % Fast
      x2 = q(16);                       % Slow
      y2 = q(17);         
%       V3 = q(18);
%       V4 = q(19);

    
qd(1) =Iapp*heaviside(t-t1)*heaviside(t2-t)+Iapp1*heaviside(t-t3)*heaviside(t4-t)+...
      Iapp2*heaviside(t-t5)*heaviside(t6-t)+ 4*((0.1*(50-(127*V1/105+8265/105))/(exp((50 - ...
    (127*V1/105 + 8265/105))/10) - 1))/((0.1*(50 - (127*V1/105 + 8265/105))/(exp((50 -...
    (127*V1/105 + 8265/105))/10) - 1))+(4*exp((25 - (127*V1/105 + 8265/105))/18))))^3*h1*(30- V1)...
    + 0.3*n1^4*(-75 - V1) + 0.01*x1*(30-V1) +0.03*Ca1/(.5 + Ca1)*(-75 - V1) +...
    0.003*(-40 - V1) +gh*((1/(1+exp(-(V1+63)/7.8)))^3)*y1*(+120-V1);  
qd(2)  =0.00003*(0.0085*x1*(140 - V1+Ca_shift)-Ca1);
qd(3)  =((1-h1)*(0.07*exp((25 - (127*V1/105 + 8265/105))/20))-h1*(1.0/(1 + exp((55 - (127*V1/105 + 8265/105))/10))))/12.5;
qd(4)  =((1-n1)*(0.01*(55 - (127*V1/105 + 8265/105))/(exp((55 - (127*V1/105 + 8265/105))/10) - 1))-n1*(0.125*exp((45 - (127*V1/105 + 8265/105))/80)))/12.5;
qd(5)  =((1/(exp(0.15*(-V1-50+x_shift))+1))-x1)/100;
qd(6)  =.5*((1/(1+exp(10*(V1-Vhh))))-y1)/(7.1+10.4/(1+exp((V1+68)/2.2)));
qd(7)  =alpha*(1-s1)/(1+exp(-100*(V1+40)))-beta*s1;
qd(8)  =2*alpha*(1-s2)/(1+exp(-100*(V1+40)))-beta*s2;
qd(9) = .5*alpha*(1-s3)/(1+exp(-100*(V1+40)))-beta*s3;
qd(10) =.1*alpha*(1-s4)/(1+exp(-100*(V1+40)))-beta*s4;
qd(11) =alpha1*sl*(1-sl)/(1+exp(-100*(V1+40)))-beta1*(sl-0.002);

qd(12)= 4*((0.1*(50-(127*V2/105+8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))...
    /((0.1*(50 - (127*V2/105 + 8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V2/105 + 8265/105))/18))))^3*h2*(30 - V2) + 0.3*n2^4*(-75 - V2) +...
    0.01*x2*(30-V2) +0.03*Ca2/(.5 + Ca2)*(-75 - V2) + 0.008*(-40 - V2) +...
    gh*((1/(1+exp(-(V2+63)/7.8)))^3)*y2*(-V2+120)  -g12*(V2+55)*s2;
qd(13) =0.00003*(0.0085*x2*(140 - V2+Ca_shift)-Ca2);
qd(14) =((1-h2)*(0.07*exp((25 - (127*V2/105 + 8265/105))/20))-h2*(1.0/(1 + exp((55 - (127*V2/105 + 8265/105))/10))))/12.5;
qd(15) =((1-n2)*(0.01*(55 - (127*V2/105 + 8265/105))/(exp((55 - (127*V2/105 + 8265/105))/10) - 1))-n2*(0.125*exp((45 - (127*V2/105 + 8265/105))/80)))/12.5;
qd(16) =((1/(exp(0.15*(-V2-50+x_shift))+1))-x2)/100;
qd(17) =.5*((1/(1+exp(10*(V2-Vhh))))-y2)/(7.1+10.4/(1+exp((V2+68)/2.2)));
end


    

