
%close all 
%clear all
    table1=load('Fig5A.txt');
    table8=load('JN2016_Fig5B_130917_01_3351.txt');
%   table2=xlsread('Fig. 5C (130917-01) Sakurai and Katz 2016.xlsx');
%   table3=xlsread('Fig. 5E (150413-01) Sakurai and Katz 2016.xlsx');
%   table4=xlsread('Fig. 5F (150413-01) Sakurai and Katz 2016.xlsx');
%   table7=xlsread('Fig. 5G (130917-01) Sakurai and Katz 2016.xlsx');
%   table5=xlsread('Dendronotus Si3 paper Fig6A (Si2).xlsx');
%   table6=xlsread('Dendronotus Si3 paper Fig6B (Si3).xlsx');
  
%Fig 5A
timea=table1(:,1);
length(timea)
l_si2a=table1(:,2);
r_si2a=table1(:,3);
l_si3a=table1(:,4);
r_si3a=table1(:,5);
figure(1)
clf 
t1=table1(1,1)+12; 
tt1=1000;
plot(timea(1:tt1)-t1,l_si2a(1:tt1)+5,'Linewidth',1.5','Color',[0.0 0 .5]);
hold on
tt2=47079;
plot(timea(tt1+42:tt2)-t1,l_si2a(tt1+42:tt2)+25,'Linewidth',1.5','Color',[0.0 0 .5]);
hold on
tt3=85100;
plot(timea(tt2+42:tt3)-t1,l_si2a(tt2+42:tt3)+100,'Linewidth',1.5','Color',[0.0 0 .5]);
hold on
plot(timea(tt3+42:end)-t1,l_si2a(tt3+42:end)+20,'Linewidth',1.5','Color',[0.0 0 .5]);
hold on
plot(timea-t1,r_si2a-90,'Linewidth',1.5','Color',[0 0 1]);
hold on
plot(timea-t1,l_si3a-190,'Linewidth',1.5','Color',[0.5 0 0]);
hold on
plot(timea-t1,r_si3a-290,'Linewidth',1.5','Color',[1 0 0]);
hold on
axis([0 43 -360 70])
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8. 3.]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8. 3.]);
set(gca,'FontName','Italica');
%title('Fig5A  Si2L / Si2R / Si3L  /  Si3R', 'Fontsize', 11);
%set(gca,'ytick',[])
%axis off 
box on 

% Set axis limits and labels

%xlabel('Time [s]', 'Fontsize', 23)
ylabel('Voltage [mV]', 'Fontsize', 12)
box on
%print(gcf,'-dpng','-r600' ,'den_Fig5A.png');
fig=figure(1)
exportgraphics(fig, 'Fig5A.jpg', 'Resolution', 300,'ContentType', 'image', 'BackgroundColor', 'white')

%Fig 5B
time_5b=table8(:,1);

l_si2_5b=table8(:,2);
size(l_si2_5b)
r_si2_5b=table8(:,3);
size(r_si2_5b)
l_si3_5b=table8(:,4);
r_si3_5b=table8(:,5); 

figure(2)
clf
t1=table8(1,1)+12; 
tt1=36700;
plot(time_5b(1:tt1)-t1,l_si2_5b(1:tt1)+0,'Linewidth',1.5','Color',[0.0 0 .5]);
hold on
tt1=36800;
tt2=72450;
plot(time_5b(tt1:tt2)-t1,l_si2_5b(tt1:tt2)+250,'Linewidth',1.5','Color',[0.0 0 .5]);
hold on
tt2=72590;
tt3=76000;
plot(time_5b(tt2:tt3)-t1,l_si2_5b(tt2:tt3)+0,'Linewidth',1.5','Color',[0.0 0 .5]);
hold on
tt3=75600;
plot(time_5b(tt3:end)-t1,l_si2_5b(tt3:end)+0,'Linewidth',1.5','Color',[0.0 0 .5]);
hold on


plot(time_5b-t1,r_si2_5b-135,'Linewidth',1.5','Color',[0.0 0 1]);
hold on
plot(time_5b-t1,l_si3_5b-230,'Linewidth',1.5','Color',[0.5 0 0]);
hold on
plot(time_5b-t1,r_si3_5b-320,'Linewidth',1.5','Color',[1 0 0]);
hold on
axis([0 35 -380 60])
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [8. 3.]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8. 3.]);
set(gca,'FontName','Italica');
%title('Fig5B Si2L / Si2R / Si3L / Si3R');
set(gca,'ytick',[])
%print(gcf,'-dpdf','-r0' ,'den_swim_curare.pdf');
fig=figure(2)
exportgraphics(fig, 'Fig5B.jpg', 'Resolution', 300,'ContentType', 'image', 'BackgroundColor', 'white')

