clear all 
close all

LW=2;
LW2=4;
MS=6;
MS3=9;
FigFont=16;
C1=3;
Colors
%%

G12PH=3.2;
ThreeDPH=16.4;
PADPH=11.7;
SintPH=10.2;
P09PH=6;
MB19PH=16.1;

CSint=7;
CPAD=8;
CG12=9;
C3D=1;

CSint=9;
CPAD=2;
CG12=1;
C3D=5;
CP09=3;
CMB19=10;

%% Model plots

load G12.mat

figure
plot(t.*1e-3,zeros(size(t)),'k')
hold on
plot(t.*1e-3,d,'Color',Color(3,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('G12','FontWeight','Normal')
set(gca,'fontsize',FigFont);
ylabel('Dipole')
xlabel('Time (Myr)')
%xlim([25 30])
xlim([30 35])
ylim([-2.5 2.5])
box off

load ThreeDCoarse.mat

figure
plot(t.*1e-3,zeros(size(t)),'k')
hold on
plot(t.*1e-3,d,'Color',Color(3,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('3-D','FontWeight','Normal')
set(gca,'fontsize',FigFont);
ylabel('Dipole')
xlabel('Time (Myr)')
%xlim([100 105])
xlim([50 55])
ylim([-2.5 2.5])
box off


load Sint2000.mat

figure
plot(t.*1e-3,zeros(size(t)),'k')
hold on
h2=plot(t.*1e-3,d,'Color',Color(CSint,:),'LineWidth',LW);
load PADM2Mstar
h1=plot(t.*1e-3,d,'Color',Color(CPAD,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('Paleomagnetic reconstructions','FontWeight','Normal')
legend([h1 h2],'PADM2M','Sint-2000')
set(gca,'fontsize',FigFont);
ylabel('Dipole')
xlabel('Time (Myr)')
ylim([-2.5 2.5])
box off

%% G12 SVM

load('Results/G12s_180_Combined.mat')
PH=G12PH;

WarningThresh=.1:.0025:1.3;
[~,WTbreak]=min(abs(WarningThresh-.8));
MCC=zeros(length(WarningThresh),length(PH));

load('G12sb.mat')
[start,~]=FindEvents(t,d,.1,.8);
idx=find(t==start(FullResults(4,1),1));
t=t(1:idx);
d=d(1:idx);
[start,stop]=FindEvents(t,d,.1,.8);
for ii=1:length(PH)
    PN=PosNeg(t,start,stop,PH(ii));
  
    for jj=1:length(WarningThresh)
        WT=WarningThresh(jj);
        [TP,FP,TN,FN]=TestThresh(t,d,WT,PN);
        MCC(jj,ii)=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    end
end
[G12Mtrain,I1]=max(MCC(1:WTbreak,:));

load('G12s.mat')
[start,stop]=FindEvents(t,d,.1,.8);
for ii=1:length(PH)
    PN=PosNeg(t,start,stop,PH(ii));
    
    WT=WarningThresh(I1(ii));
    [TP,FP,TN,FN]=TestThresh(t,d,WT,PN);
    M1(ii)=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    
end
G12Thresh=M1;
load('Results/G12s_180_Combined.mat')

L=unique(tvMCCResults(2,:));
G12L=L;
MCCvL=zeros(length(L),length(PH));

for ii=1:length(PH)
    idx=find(tvMCCResults(1,:)==PH(ii));
    temp=tvMCCResults(:,idx);
    for jj=1:length(L)
        idx=find(temp(2,:)==L(jj));
        MCCvL(jj,ii)=temp(end,idx);
    end
end
G12MCCvL=MCCvL;


%% G12, Short training

load('Results/G12s_5_180_Comb.mat')
PH=G12PH;

WarningThresh=.1:.0025:1.3;
[~,WTbreak]=min(abs(WarningThresh-.8));
MCC=zeros(length(WarningThresh),length(PH));

load('G12sb.mat')
[start,~]=FindEvents(t,d,.1,.8);
idx=find(t==start(FullResults(4,1),1));
t=t(1:idx);
d=d(1:idx);
[start,stop]=FindEvents(t,d,.1,.8);
for ii=1:length(PH)
    PN=PosNeg(t,start,stop,PH(ii));
  
    for jj=1:length(WarningThresh)
        WT=WarningThresh(jj);
        [TP,FP,TN,FN]=TestThresh(t,d,WT,PN);
        MCC(jj,ii)=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    end
end
[M1,I1]=max(MCC(1:WTbreak,:));

load('G12s.mat')
[start,stop]=FindEvents(t,d,.1,.8);
for ii=1:length(PH)
    PN=PosNeg(t,start,stop,PH(ii));
    
    WT=WarningThresh(I1(ii));
    [TP,FP,TN,FN]=TestThresh(t,d,WT,PN);
    M1(ii)=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    
end
G12Thresh5=M1;


load('Results/G12s_5_180_Comb.mat')

L=unique(FullResults(2,:));
G12L5=[.2,L];
MCCvL=zeros(length(L),length(PH));

for ii=1:length(PH)
    idx=find(tvMCCResults(1,:)==PH(ii));
    temp=tvMCCResults(:,idx);
    for jj=1:length(G12L5)
        idx=find(temp(2,:)==G12L5(jj));
        MCCvL(jj,ii)=temp(end,idx);
    end
end
G12MCCvL5=MCCvL;


%% ThreeD CFA

load('Results/ThreeDCoarse_180_CFA_Full.mat')
PH=ThreeDPH;

WarningThresh=.1:.0025:1.3;
[~,WTbreak]=min(abs(WarningThresh-.8));
MCC=zeros(length(WarningThresh),length(PH));

load('ThreeDCoarse.mat')
[start,~]=FindEvents(t,d,.1,.8);
idx=find(t==start(FullResults(4,1),1));
t=t(1:idx);
d=d(1:idx);
[start,stop]=FindEvents(t,d,.1,.8);
for ii=1:length(PH)
    PN=PosNeg(t,start,stop,PH(ii));
  
    for jj=1:length(WarningThresh)
        WT=WarningThresh(jj);
        [TP,FP,TN,FN]=TestThresh(t,d,WT,PN);
        MCC(jj,ii)=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    end
end
[ThreeDMtrain,I1]=max(MCC(1:WTbreak,:));

%load('ThreeDkyr.mat')
load('ThreeDCoarse.mat')
t=t(idx:end);
d=d(idx:end);

[start,stop]=FindEvents(t,d,.1,.8);
for ii=1:length(PH)
    PN=PosNeg(t,start,stop,PH(ii));
    
    WT=WarningThresh(I1(ii));
    [TP,FP,TN,FN]=TestThresh(t,d,WT,PN);
    M1(ii)=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    
end

ThreeDThresh=M1;


%% CFA
L=unique(tvMCCResults(2,:));
ThreeDL=L;

MCCvL=zeros(length(L),length(PH));

for ii=1:length(PH)
    idx=find(tvMCCResults(1,:)==PH(ii));
    temp=tvMCCResults(:,idx);
    for jj=1:length(L)
        idx=find(temp(2,:)==L(jj));
        MCCvL(jj,ii)=temp(end,idx);
    end
end
ThreeDMCCvL=MCCvL;


%% ThreeD CFA short


load('Results/ThreeDCoarse_5_CFA.mat')
PH=ThreeDPH;

WarningThresh=.1:.0025:1.3;
[~,WTbreak]=min(abs(WarningThresh-.8));
MCC=zeros(length(WarningThresh),length(PH));

load('ThreeDCoarse.mat')
[start,~]=FindEvents(t,d,.1,.8);
idx=find(t==start(FullResults(4,1),1));
t=t(1:idx);
d=d(1:idx);
[start,stop]=FindEvents(t,d,.1,.8);
for ii=1:length(PH)
    PN=PosNeg(t,start,stop,PH(ii));
  
    for jj=1:length(WarningThresh)
        WT=WarningThresh(jj);
        [TP,FP,TN,FN]=TestThresh(t,d,WT,PN);
        MCC(jj,ii)=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    end
end
[M1,I1]=max(MCC(1:WTbreak,:));

load('ThreeDCoarse.mat')
t=t(idx:end);
d=d(idx:end);
[start,stop]=FindEvents(t,d,.1,.8);
for ii=1:length(PH)
    PN=PosNeg(t,start,stop,PH(ii));
    
    WT=WarningThresh(I1(ii));
    [TP,FP,TN,FN]=TestThresh(t,d,WT,PN);
    M1(ii)=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    
end

ThreeDThresh5=M1;


%% CFA
L=unique(tvMCCResults(2,:));
L(find(L==0))=[];
ThreeDL5=L;

MCCvL=zeros(length(L),length(PH));

for ii=1:length(PH)
    idx=find(tvMCCResults(1,:)==PH(ii));
    temp=tvMCCResults(:,idx);
    for jj=1:length(L)
        idx=find(temp(2,:)==L(jj));
        MCCvL(jj,ii)=temp(end,idx);
    end
end
ThreeDMCCvL5=MCCvL;

figure
h1=plot(G12L5,G12MCCvL5(:,1),'-o','Color',Color(CG12,:),'MarkerEdgeColor',Color(CG12,:),'MarkerFaceColor',Color(CG12,:),'MarkerSize',MS);
hold on
plot(G12L5,G12Thresh5.*ones(length(G12L5),1),'--','Color',Color(CG12,:),'LineWidth',LW);
h2=plot(ThreeDL5,ThreeDMCCvL5(:,1),'-o','Color',Color(C3D,:),'MarkerEdgeColor',Color(C3D,:),'MarkerFaceColor',Color(C3D,:),'MarkerSize',MS);
plot([0,ThreeDL5],ThreeDThresh.*ones(length(ThreeDL5)+1,1),'--','Color',Color(C3D,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('Short training','FontWeight','Normal')
set(gca,'fontsize',FigFont);
%legend('PH=1.6 kyr}','\rho^{Clim}_{16,1456}','\rho^{Ens}_{16,961}','\rho^{Clim}_{16,961}')%'Location','SouthEast')
legend('G12 SVM','G12 Threshold','3-D SVM','3-D Threshold','Location','southeast')
ylabel('MCC')
%legend([h1 h2],'G12','3-D','Location','east')
xlabel('Search window (kyr)')
xlim([0,max(G12L5)])
ylim([0,1])
box off

%% G12 and ThreeD together

figure
h1=plot(G12L,G12MCCvL(:,1),'-o','Color',Color(CG12,:),'MarkerEdgeColor',Color(CG12,:),'MarkerFaceColor',Color(CG12,:),'MarkerSize',MS);
hold on
plot(G12L,G12Thresh.*ones(length(G12L),1),'--','Color',Color(CG12,:),'LineWidth',LW);
h2=plot(ThreeDL,ThreeDMCCvL(:,1),'-o','Color',Color(C3D,:),'MarkerEdgeColor',Color(C3D,:),'MarkerFaceColor',Color(C3D,:),'MarkerSize',MS);
plot(ThreeDL,ThreeDThresh.*ones(length(ThreeDL),1),'--','Color',Color(C3D,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('Long training','FontWeight','Normal')
set(gca,'fontsize',FigFont);
%legend('PH=1.6 kyr}','\rho^{Clim}_{16,1456}','\rho^{Ens}_{16,961}','\rho^{Clim}_{16,961}')%'Location','SouthEast')
ylabel('MCC')
%legend([h1 h2],'G12','3-D','Location','southeast')
legend('G12 SVM','G12 Threshold','3-D SVM','3-D Threshold','Location','southeast')
xlabel('Search window (kyr)')
xlim([0,max(L)])
ylim([0,1])
box off



%% G12 Training Plot

load('Results/G12TrainingPlot.mat')

figure
yyaxis right
plot(cij(3:end),100.*tACC(3:end),'-o','Color',Color(8,:),'MarkerEdgeColor',Color(8,:),'MarkerFaceColor',Color(8,:),'MarkerSize',MS)
ylabel('Accuracy')
ytickformat('percentage');
ylim(100.*[.9 1.001])
hold on
yyaxis left
plot(cij(3:end),tMCC(3:end),'-o','Color',Color(10,:),'MarkerEdgeColor',Color(10,:),'MarkerFaceColor',Color(10,:),'MarkerSize',MS)
ylim([.7 1.01])
box off
title('SVM training','FontWeight','Normal')
legend('MCC','Accuracy','Location','SouthEast')
set(gca,'fontsize',FigFont);
xlabel('FP/FN penalty ratio')
ylabel('MCC')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%% P09 SVM

load('Results/P09s_180_CFA.mat')
PH=P09PH;

WarningThresh=.1:.0025:1.3;
[~,WTbreak]=min(abs(WarningThresh-.8));
MCC=zeros(length(WarningThresh),length(PH));

load('P09sb.mat')
[start,~]=FindEvents(t,d,.1,.8);
idx=find(t==start(FullResults(4,1),1));
t=t(1:idx);
d=d(1:idx);
[start,stop]=FindEvents(t,d,.1,.8);
for ii=1:length(PH)
    PN=PosNeg(t,start,stop,PH(ii));
  
    for jj=1:length(WarningThresh)
        WT=WarningThresh(jj);
        [TP,FP,TN,FN]=TestThresh(t,d,WT,PN);
        MCC(jj,ii)=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    end
end
[M1,I1]=max(MCC(1:WTbreak,:));

load('P09s.mat')
[start,stop]=FindEvents(t,d,.1,.8);
for ii=1:length(PH)
    PN=PosNeg(t,start,stop,PH(ii));
    
    WT=WarningThresh(I1(ii));
    [TP,FP,TN,FN]=TestThresh(t,d,WT,PN);
    M1(ii)=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    
end
P09Thresh=M1;
load('Results/P09s_180_CFA.mat')

L=unique(tvMCCResults(2,:));
P09L=L;
MCCvL=zeros(length(L),length(PH));

for ii=1:length(PH)
    idx=find(tvMCCResults(1,:)==PH(ii));
    temp=tvMCCResults(:,idx);
    for jj=1:length(L)
        idx=find(temp(2,:)==L(jj));
        MCCvL(jj,ii)=temp(end,idx);
    end
end
P09MCCvL=MCCvL;



%% MB19 SVM

load('Results/MB19s_180_CFA.mat')
PH=MB19PH;

WarningThresh=.1:.0025:1.3;
[~,WTbreak]=min(abs(WarningThresh-.8));
MCC=zeros(length(WarningThresh),length(PH));

load('MB19sb.mat')
[start,~]=FindEvents(t,d,.1,.8);
idx=find(t==start(FullResults(4,1),1));
t=t(1:idx);
d=d(1:idx);
[start,stop]=FindEvents(t,d,.1,.8);
for ii=1:length(PH)
    PN=PosNeg(t,start,stop,PH(ii));
  
    for jj=1:length(WarningThresh)
        WT=WarningThresh(jj);
        [TP,FP,TN,FN]=TestThresh(t,d,WT,PN);
        MCC(jj,ii)=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    end
end
[M1,I1]=max(MCC(1:WTbreak,:));

load('MB19s.mat')
[start,stop]=FindEvents(t,d,.1,.8);
for ii=1:length(PH)
    PN=PosNeg(t,start,stop,PH(ii));
    
    WT=WarningThresh(I1(ii));
    [TP,FP,TN,FN]=TestThresh(t,d,WT,PN);
    M1(ii)=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    
end
MB19Thresh=M1;
load('Results/MB19s_180_CFA.mat')

L=unique(tvMCCResults(2,:));
MB19L=L;
MCCvL=zeros(length(L),length(PH));

for ii=1:length(PH)
    idx=find(tvMCCResults(1,:)==PH(ii));
    temp=tvMCCResults(:,idx);
    for jj=1:length(L)
        idx=find(temp(2,:)==L(jj));
        MCCvL(jj,ii)=temp(end,idx);
    end
end
MB19MCCvL=MCCvL;


%% DW and P09

figure
h1=plot(L,P09MCCvL(:,1),'-o','Color',Color(CP09,:),'MarkerEdgeColor',Color(CP09,:),'MarkerFaceColor',Color(CP09,:),'MarkerSize',MS);
hold on
plot(L,P09Thresh.*ones(length(L),1),'--','Color',Color(CP09,:),'LineWidth',LW);
h3=plot(L,MB19MCCvL(:,1),'-o','Color',Color(CMB19,:),'MarkerEdgeColor',Color(CMB19,:),'MarkerFaceColor',Color(CMB19,:),'MarkerSize',MS);
plot(L,MB19Thresh.*ones(length(L),1),'--','Color',Color(CMB19,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('Stochastic models','FontWeight','Normal')
set(gca,'fontsize',FigFont);
legend([h3 h1],'DW','P09')
%legend('PH=1.6 kyr}','\rho^{Clim}_{16,1456}','\rho^{Ens}_{16,961}','\rho^{Clim}_{16,961}')%'Location','SouthEast')
%legend('P09 SVM','P09 Threshold','DW SVM','DW Threshold','Location','southeast')
ylabel('MCC')
xlabel('Search window (kyr)')
xlim([0,max(L)])
ylim([0,1])
box off


%% G12 and 3-D RBF

load('Results/ThreeDCoarse_180_CFA_rbf.mat')
PH=ThreeDPH;

L=unique(tvMCCResults(2,:));
RBFL=L;
MCCvL=zeros(length(L),length(PH));

ThreeDRBFMCCvL=tvMCCResults(end-1:end,:);

figure
plot(RBFL,ThreeDRBFMCCvL(1,:),'-o','Color',Color(C3D,:),'MarkerEdgeColor',Color(C3D,:),'MarkerFaceColor',Color(C3D,:),'MarkerSize',MS)
hold on
plot(RBFL,ThreeDMtrain.*ones(length(RBFL)),'--','Color',Color(C3D,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('3-D training','FontWeight','Normal')
set(gca,'fontsize',FigFont);
%legend('3-D','Location','southeast')
ylabel('MCC')
xlabel('Search window (kyr)')
xlim([0,max(L)])
ylim([0,1])
box off

figure
plot(RBFL,ThreeDRBFMCCvL(2,:),'-o','Color',Color(C3D,:),'MarkerEdgeColor',Color(C3D,:),'MarkerFaceColor',Color(C3D,:),'MarkerSize',MS)
hold on
plot(RBFL,ThreeDThresh.*ones(length(RBFL)),'--','Color',Color(C3D,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('3-D validation','FontWeight','Normal')
set(gca,'fontsize',FigFont);
%legend('3-D','Location','southeast')
ylabel('MCC')
xlabel('Search window (kyr)')
xlim([0,max(L)])
ylim([0,1])
box off


%% G12 and 3-D LSTM

load('Results/ThreeDCoarse_LSTM.mat')
PH=ThreeDPH;

L=unique(tvMCCResults(2,:));
RBFL=L;
MCCvL=zeros(length(L),length(PH));

ThreeDRBFMCCvL=tvMCCResults(end-1:end,:);


figure
plot(RBFL,ThreeDRBFMCCvL(1,:),'-o','Color',Color(C3D,:),'MarkerEdgeColor',Color(C3D,:),'MarkerFaceColor',Color(C3D,:),'MarkerSize',MS)
hold on
plot(RBFL,ThreeDMtrain.*ones(length(RBFL)),'--','Color',Color(C3D,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('3-D training','FontWeight','Normal')
set(gca,'fontsize',FigFont);
%legend([h1 h2],'G12','3-D','Location','southeast')
ylabel('MCC')
xlabel('Search window (kyr)')
xlim([0,max(L)])
ylim([0,1])
box off

figure
plot(RBFL,ThreeDRBFMCCvL(2,:),'-o','Color',Color(C3D,:),'MarkerEdgeColor',Color(C3D,:),'MarkerFaceColor',Color(C3D,:),'MarkerSize',MS)
hold on
plot(RBFL,ThreeDThresh.*ones(length(RBFL)),'--','Color',Color(C3D,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('3-D validation','FontWeight','Normal')
set(gca,'fontsize',FigFont);
%legend([h1 h2],'G12','3-D','Location','southeast')
legend('3-D LSTM','3-D Threshold','Location','southeast')
ylabel('MCC')
xlabel('Search window (kyr)')
xlim([0,max(L)])
ylim([0,1])
box off

%% Stacked plots

clear all

alpha=.05;
WindowWidth=50;
bins = linspace(0.2,2.35,50);

LW=2;
MS=6;
FigFont=16;
C1=3;
Colors
%%

G12PH=3.2;
ThreeDPH=16.4;
PADPH=11.7;
SintPH=10.2;
P09PH=6;
MB19PH=16.1;

CSint=9;
CPAD=2;
CG12=1;
C3D=5;
CP09=3;
CMB19=4;

%% G12

PH=G12PH;

load G12sb                      % Use long training to make plots

[start,stop]=FindEvents(t,d,.1,.8);
PN=PosNeg(t,start,stop,PH);

tidx=[];
for ii=1:length(start(:,1))
    [~,temp]=min(abs(t-start(ii,1)));
    tidx=[tidx;temp];
end

dt=t(2)-t(1);
G12tsteps=round(WindowWidth/dt);

data=zeros(length(tidx),G12tsteps+1);
for ii=1:length(tidx)
    if tidx(ii)-G12tsteps<=0
        data(ii,:)=-999;
    elseif any(PN(tidx(ii)-G12tsteps:tidx(ii)-1)==2)
        data(ii,:)=-999;
    else
        data(ii,:)=d(tidx(ii)-G12tsteps:tidx(ii));
        if data(ii,1)<0
            data(ii,:) = -1.*data(ii,:);
        end
    end
end

data(find(data(:,1)==-999),:)=[];

G12data=transpose(data);

%PlotDensity(dt.*(-tsteps:1:0),bins,G12data)

%% G12 add more events

PH=G12PH;

load G12s                      % Use long training to make plots

[start,stop]=FindEvents(t,d,.1,.8);
PN=PosNeg(t,start,stop,PH);

tidx=[];
for ii=1:length(start(:,1))
    [~,temp]=min(abs(t-start(ii,1)));
    tidx=[tidx;temp];
end

dt=t(2)-t(1);
G12tsteps=round(WindowWidth/dt);

data=zeros(length(tidx),G12tsteps+1);
for ii=1:length(tidx)
    if tidx(ii)-G12tsteps<=0
        data(ii,:)=-999;
    elseif any(PN(tidx(ii)-G12tsteps:tidx(ii)-1)==2)
        data(ii,:)=-999;
    else
        data(ii,:)=d(tidx(ii)-G12tsteps:tidx(ii));
        if data(ii,1)<0
            data(ii,:) = -1.*data(ii,:);
        end
    end
end

data(find(data(:,1)==-999),:)=[];

G12data=[G12data,transpose(data)];


%% 3D

PH=ThreeDPH;

load ThreeDCoarse
[start,~]=FindEvents(t,d,.1,.8);
idx=find(t==start(360,1));
t=t(1:idx);
d=d(1:idx);

[start,stop]=FindEvents(t,d,.1,.8);
PN=PosNeg(t,start,stop,PH);

tidx=[];
for ii=1:length(start(:,1))
    [~,temp]=min(abs(t-start(ii,1)));
    tidx=[tidx;temp];
end

dt=t(2)-t(1);
ThreeDtsteps=round(WindowWidth/dt);

data=zeros(length(tidx),ThreeDtsteps+1);
for ii=1:length(tidx)
    if tidx(ii)-ThreeDtsteps<=0
        data(ii,:)=-999;
    elseif any(PN(tidx(ii)-ThreeDtsteps:tidx(ii)-1)==2)
        data(ii,:)=-999;
    else
        data(ii,:)=d(tidx(ii)-ThreeDtsteps:tidx(ii));
        if data(ii,1)<0
            data(ii,:) = -1.*data(ii,:);
        end
    end
end

data(find(data(:,1)==-999),:)=[];
        
ThreeDdata=transpose(data);


%% Plots

% [~,~,~,DensityG12] = PlotDensity(dt.*(-G12tsteps:1:0),bins,G12data);
% 
% MaxDensity=max(max(DensityG12));
% caxis([0 MaxDensity])
% colorbar
% set(gca,'fontsize',FigFont);
% title('G12','FontWeight','Normal')
% set(gca,'fontsize',FigFont);
% ylabel('Dipole')
% xlabel('Time (kyr)')
% xlim([-WindowWidth 0])
% 
% [~,~,~,Density3D] = PlotDensity(dt.*(-ThreeDtsteps:1:0),bins,ThreeDdata);
% 
% caxis([0 MaxDensity])
% colorbar
% set(gca,'fontsize',FigFont);
% title('3-D','FontWeight','Normal')
% set(gca,'fontsize',FigFont);
% ylabel('Dipole')
% xlabel('Time (kyr)')
% xlim([-WindowWidth 0])

figure
plot(dt.*(-G12tsteps:1:0),G12data,'Color',[Color(CG12,:) alpha],'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('G12','FontWeight','Normal')
set(gca,'fontsize',FigFont);
ylabel('Dipole')
xlabel('Time (kyr)')
ylim([.1 1.6])
xlim([-WindowWidth 0])
box off 

figure
plot(dt.*(-ThreeDtsteps:1:0),ThreeDdata,'Color',[Color(C3D,:) alpha],'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('3-D','FontWeight','Normal')
set(gca,'fontsize',FigFont);
ylabel('Dipole')
xlabel('Time (kyr)')
xlim([-WindowWidth 0])
ylim([.1 1.6])
box off

%% Paleo stacked

%% PADM2M

PH=PADPH;

load PADM2Mstar                      % Use long training to make plots

[start,stop]=FindEvents(t,d,.1,.8);
PN=PosNeg(t,start,stop,PH);

tidx=[];
for ii=1:length(start(:,1))
    [~,temp]=min(abs(t-start(ii,1)));
    tidx=[tidx;temp];
end

dt=t(2)-t(1);
PADtsteps=round(WindowWidth/dt);

data=zeros(length(tidx),PADtsteps+1);
for ii=1:length(tidx)
    if tidx(ii)-PADtsteps<=0
        data(ii,:)=-999;
    elseif any(PN(tidx(ii)-PADtsteps:tidx(ii)-1)==2)
        data(ii,:)=-999;
    else
        data(ii,:)=d(tidx(ii)-PADtsteps:tidx(ii));
        if data(ii,1)<0
            data(ii,:) = -1.*data(ii,:);
        end
    end
end

data(find(data(:,1)==-999),:)=[];

PADdata=transpose(data);

PH=SintPH;

%% Sint-2000
load Sint2000                     % Use long training to make plots

[start,stop]=FindEvents(t,d,.1,.8);
PN=PosNeg(t,start,stop,PH);

tidx=[];
for ii=1:length(start(:,1))
    [~,temp]=min(abs(t-start(ii,1)));
    tidx=[tidx;temp];
end

dt=t(2)-t(1);
Sinttsteps=round(WindowWidth/dt);

data=zeros(length(tidx),Sinttsteps+1);
for ii=1:length(tidx)
    if tidx(ii)-Sinttsteps<=0
        data(ii,:)=-999;
    elseif any(PN(tidx(ii)-Sinttsteps:tidx(ii)-1)==2)
        data(ii,:)=-999;
    else
        data(ii,:)=d(tidx(ii)-Sinttsteps:tidx(ii));
        if data(ii,1)<0
            data(ii,:) = -1.*data(ii,:);
        end
    end
end

data(find(data(:,1)==-999),:)=[];

Sintdata=transpose(data);

%% paleo plots

figure
plot(dt.*(-PADtsteps:1:0),PADdata,'Color',Color(CPAD,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('PADM2M','FontWeight','Normal')
set(gca,'fontsize',FigFont);
ylabel('Dipole')
xlabel('Time (kyr)')
ylim([.1 1.3])
xlim([-WindowWidth 0])
box off 

figure
plot(dt.*(-Sinttsteps:1:0),Sintdata,'Color',Color(CSint,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('Sint-2000','FontWeight','Normal')
set(gca,'fontsize',FigFont);
ylabel('Dipole')
xlabel('Time (kyr)')
xlim([-WindowWidth 0])
ylim([.1 1.3])
box off

%% Paleo k-fold plots

clear all 

LW=2;
MS=6;
MS2=50;
MS3=9;
at=.5;
k1=.95;
k2=.8;
FigFont=16;
C1=3;
Colors
%%

G12PH=3.2;
ThreeDPH=16.4;
PADPH=11.7;
SintPH=10.2;
P09PH=6;
MB19PH=16.1;

CSint=7;
CPAD=8;
CG12=9;
C3D=1;

CSint=9;
CPAD=2;
CG12=1;
C3D=5;
CP09=3;
CMB19=4;

%% PADM2M

PH=PADPH;

M1=Kfold_Thresh(11.7,.1:.0025:1,'PADM2Mstar',5,'Thresh');


% Cost function adjustment

load('Results/PADM2Mstar_5_fold_CFA.mat')

sigma=std(PartResults(2:end,:));

L=unique(tvMCCResults(2,:));
MCCvL=zeros(length(L),length(PH));

for ii=1:length(PH)
    idx=find(tvMCCResults(1,:)==PH(ii));
    temp=tvMCCResults(:,idx);
    for jj=1:length(L)
        idx=find(temp(2,:)==L(jj));
        MCCvL(jj,ii)=temp(end,idx);
    end
end

Top=max(PartResults(2:end,:));
Bot=min(PartResults(2:end,:));

load('Results/PADM2Mstar_0_Full.mat')


figure
fill([L,fliplr(L)],[Top,fliplr(Bot)],[k1,k1,k1])
hold on
%plot(L,PartResults(2:end,:),'o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerFaceColor',Color(CPAD,:),'MarkerSize',MS)
%scatter(reshape([L;L;L;L;L],25*5,1),reshape(PartResults(2:end,:),25*5,1),MS2,Color(CPAD,:),'filled')
scatter(reshape([L;L;L;L;L],25*5,1),reshape(PartResults(2:end,:),25*5,1),MS2,[k2,k2,k2],'filled')
%alpha(at)
%fill([L,fliplr(L)],[Plus2sig,fliplr(Minus2sig)],[k1,k1,k1])
plot(L,MCCvL(:,1),'-o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerFaceColor',Color(CPAD,:),'MarkerSize',MS)
plot(L,tvMCCResults(end,:),'-d','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerSize',MS3)
%plot(L,MCCvL(:,1),'-o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerFaceColor',Color(CPAD,:),'MarkerSize',MS)
%plot(L,tvMCCResults(end,:),'-o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerSize',MS)
% plot(L,PartResults(2,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(3,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(4,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(5,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(6,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
plot(L,M1(1).*ones(length(L)),'--','Color',Color(CPAD,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('PADM2M','FontWeight','Normal')
set(gca,'fontsize',FigFont);
%legend('PH=1.6 kyr}','\rho^{Clim}_{16,1456}','\rho^{Ens}_{16,961}','\rho^{Clim}_{16,961}')%'Location','SouthEast')
ylabel('MCC')
xlabel('Search window (kyr)')
xlim([min(L),max(L)])
ylim([0,1])
box off



%% Sint-2000

PH=SintPH;

M1=Kfold_Thresh(SintPH,.1:.0025:1,'Sint2000',5,'Thresh');


% Cost function adjustment

load('Results/Sint2000_5_fold_CFA.mat')

sigma=std(PartResults(2:end,:));

L=unique(tvMCCResults(2,:));
MCCvL=zeros(length(L),length(PH));

for ii=1:length(PH)
    idx=find(tvMCCResults(1,:)==PH(ii));
    temp=tvMCCResults(:,idx);
    for jj=1:length(L)
        idx=find(temp(2,:)==L(jj));
        MCCvL(jj,ii)=temp(end,idx);
    end
end

Top=max(PartResults(2:end,:));
Bot=min(PartResults(2:end,:));

load('Results/Sint2000_0_Full.mat')



figure
fill([L,fliplr(L)],[Top,fliplr(Bot)],[k1,k1,k1])
hold on
%plot(L,PartResults(2:end,:),'o','Color',Color(CSint,:),'MarkerEdgeColor',Color(CSint,:),'MarkerFaceColor',Color(CSint,:),'MarkerSize',MS)
%scatter(reshape([L;L;L;L;L],25*5,1),reshape(PartResults(2:end,:),25*5,1),MS2,Color(CSint,:),'filled')
scatter(reshape([L;L;L;L;L],25*5,1),reshape(PartResults(2:end,:),25*5,1),MS2,[k2,k2,k2],'filled')
%alpha(at)
%fill([L,fliplr(L)],[Plus2sig,fliplr(Minus2sig)],[k1,k1,k1])
plot(L,MCCvL(:,1),'-o','Color',Color(CSint,:),'MarkerEdgeColor',Color(CSint,:),'MarkerFaceColor',Color(CSint,:),'MarkerSize',MS)
plot(L,tvMCCResults(end,:),'-d','Color',Color(CSint,:),'MarkerEdgeColor',Color(CSint,:),'MarkerSize',MS3)
%plot(L,MCCvL(:,1),'-o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerFaceColor',Color(CPAD,:),'MarkerSize',MS)
%plot(L,tvMCCResults(end,:),'-o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerSize',MS)
% plot(L,PartResults(2,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(3,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(4,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(5,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(6,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
plot(L,M1(1).*ones(length(L)),'--','Color',Color(CSint,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('Sint-2000','FontWeight','Normal')
set(gca,'fontsize',FigFont);
%legend('PH=1.6 kyr}','\rho^{Clim}_{16,1456}','\rho^{Ens}_{16,961}','\rho^{Clim}_{16,961}')%'Location','SouthEast')
ylabel('MCC')
xlabel('Search window (kyr)')
xlim([min(L),max(L)])
ylim([0,1])
box off




%% PADM2M rbf

PH=PADPH;

M1=Kfold_Thresh(11.7,.1:.0025:1,'PADM2Mstar',5,'Thresh');


% Cost function adjustment

load('Results/PADM2Mstar_5_fold_rbf.mat')

sigma=std(PartResults(2:end,:));

L=unique(tvMCCResults(2,:));
MCCvL=zeros(length(L),length(PH));

for ii=1:length(PH)
    idx=find(tvMCCResults(1,:)==PH(ii));
    temp=tvMCCResults(:,idx);
    for jj=1:length(L)
        idx=find(temp(2,:)==L(jj));
        MCCvL(jj,ii)=temp(end,idx);
    end
end

Top=max(PartResults(2:end,:));
Bot=min(PartResults(2:end,:));

load('Results/PADM2Mstar_0_rbfFull.mat')



figure
fill([L,fliplr(L)],[Top,fliplr(Bot)],[k1,k1,k1])
hold on
%plot(L,PartResults(2:end,:),'o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerFaceColor',Color(CPAD,:),'MarkerSize',MS)
%scatter(reshape([L;L;L;L;L],25*5,1),reshape(PartResults(2:end,:),25*5,1),MS2,Color(CPAD,:),'filled')
scatter(reshape([L;L;L;L;L],25*5,1),reshape(PartResults(2:end,:),25*5,1),MS2,[k2,k2,k2],'filled')
%alpha(at)
%fill([L,fliplr(L)],[Plus2sig,fliplr(Minus2sig)],[k1,k1,k1])
plot(L,MCCvL(:,1),'-o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerFaceColor',Color(CPAD,:),'MarkerSize',MS)
plot(L,tvMCCResults(end,:),'-d','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerSize',MS3)
%plot(L,MCCvL(:,1),'-o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerFaceColor',Color(CPAD,:),'MarkerSize',MS)
%plot(L,tvMCCResults(end,:),'-o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerSize',MS)
% plot(L,PartResults(2,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(3,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(4,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(5,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(6,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
plot(L,M1(1).*ones(length(L)),'--','Color',Color(CPAD,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('PADM2M','FontWeight','Normal')
set(gca,'fontsize',FigFont);
%legend('PH=1.6 kyr}','\rho^{Clim}_{16,1456}','\rho^{Ens}_{16,961}','\rho^{Clim}_{16,961}')%'Location','SouthEast')
ylabel('MCC')
xlabel('Search window (kyr)')
xlim([min(L),max(L)])
ylim([0,1])
box off





%% Sint-2000

PH=SintPH;

M1=Kfold_Thresh(SintPH,.1:.0025:1,'Sint2000',5,'Thresh');


% Cost function adjustment

load('Results/Sint2000_5_fold_rbf.mat')

sigma=std(PartResults(2:end,:));

L=unique(tvMCCResults(2,:));
MCCvL=zeros(length(L),length(PH));

for ii=1:length(PH)
    idx=find(tvMCCResults(1,:)==PH(ii));
    temp=tvMCCResults(:,idx);
    for jj=1:length(L)
        idx=find(temp(2,:)==L(jj));
        MCCvL(jj,ii)=temp(end,idx);
    end
end

Top=max(PartResults(2:end,:));
Bot=min(PartResults(2:end,:));

load('Results/Sint2000_0_rbfFull.mat')



figure
fill([L,fliplr(L)],[Top,fliplr(Bot)],[k1,k1,k1])
hold on
%plot(L,PartResults(2:end,:),'o','Color',Color(CSint,:),'MarkerEdgeColor',Color(CSint,:),'MarkerFaceColor',Color(CSint,:),'MarkerSize',MS)
%scatter(reshape([L;L;L;L;L],25*5,1),reshape(PartResults(2:end,:),25*5,1),MS2,Color(CSint,:),'filled')
scatter(reshape([L;L;L;L;L],25*5,1),reshape(PartResults(2:end,:),25*5,1),MS2,[k2,k2,k2],'filled')
%alpha(at)
%fill([L,fliplr(L)],[Plus2sig,fliplr(Minus2sig)],[k1,k1,k1])
plot(L,MCCvL(:,1),'-o','Color',Color(CSint,:),'MarkerEdgeColor',Color(CSint,:),'MarkerFaceColor',Color(CSint,:),'MarkerSize',MS)
plot(L,tvMCCResults(end,:),'-d','Color',Color(CSint,:),'MarkerEdgeColor',Color(CSint,:),'MarkerSize',MS3)
%plot(L,MCCvL(:,1),'-o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerFaceColor',Color(CPAD,:),'MarkerSize',MS)
%plot(L,tvMCCResults(end,:),'-o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerSize',MS)
% plot(L,PartResults(2,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(3,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(4,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(5,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(6,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
plot(L,M1(1).*ones(length(L)),'--','Color',Color(CSint,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('Sint-2000','FontWeight','Normal')
set(gca,'fontsize',FigFont);
%legend('PH=1.6 kyr}','\rho^{Clim}_{16,1456}','\rho^{Ens}_{16,961}','\rho^{Clim}_{16,961}')%'Location','SouthEast')
ylabel('MCC')
xlabel('Search window (kyr)')
xlim([min(L),max(L)])
ylim([0,1])
box off



%% PADM2M LSTM

PH=PADPH;

M1=Kfold_Thresh(11.7,.1:.0025:1,'PADM2Mstar',5,'Thresh');


% Cost function adjustment

load('Results/PADM2Mstar_5_fold_h50.mat')

sigma=std(PartResults(2:end,:));

L=unique(tvMCCResults(2,:));
MCCvL=zeros(length(L),length(PH));

for ii=1:length(PH)
    idx=find(tvMCCResults(1,:)==PH(ii));
    temp=tvMCCResults(:,idx);
    for jj=1:length(L)
        idx=find(temp(2,:)==L(jj));
        MCCvL(jj,ii)=temp(end,idx);
    end
end

Top=MCCvL'+2*sigma;
idx=find(isnan(Top));
idx=[idx,length(Top)+1];
%Plus2sig(isnan(Plus2sig))=1e1;
Bot=MCCvL'-2*sigma;
%Minus2sig(isnan(Minus2sig))=-1e1;

load('Results/PADM2Mstar_0_Full_LSTM.mat')



Top=max(PartResults(2:end,:));
Bot=min(PartResults(2:end,:));

figure
fill([L,fliplr(L)],[Top,fliplr(Bot)],[k1,k1,k1])
hold on
%plot(L,PartResults(2:end,:),'o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerFaceColor',Color(CPAD,:),'MarkerSize',MS)
%scatter(reshape([L;L;L;L;L],25*5,1),reshape(PartResults(2:end,:),25*5,1),MS2,Color(CPAD,:),'filled')
scatter(reshape([L;L;L;L;L],25*5,1),reshape(PartResults(2:end,:),25*5,1),MS2,[k2,k2,k2],'filled')
%alpha(at)
%fill([L,fliplr(L)],[Plus2sig,fliplr(Minus2sig)],[k1,k1,k1])
plot(L,MCCvL(:,1),'-o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerFaceColor',Color(CPAD,:),'MarkerSize',MS)
plot(L,tvMCCResults(end,:),'-d','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerSize',MS3)
%plot(L,MCCvL(:,1),'-o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerFaceColor',Color(CPAD,:),'MarkerSize',MS)
%plot(L,tvMCCResults(end,:),'-o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerSize',MS)
% plot(L,PartResults(2,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(3,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(4,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(5,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(6,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
plot(L,M1(1).*ones(length(L)),'--','Color',Color(CPAD,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('PADM2M','FontWeight','Normal')
set(gca,'fontsize',FigFont);
%legend('PH=1.6 kyr}','\rho^{Clim}_{16,1456}','\rho^{Ens}_{16,961}','\rho^{Clim}_{16,961}')%'Location','SouthEast')
ylabel('MCC')
xlabel('Search window (kyr)')
xlim([min(L),max(L)])
ylim([0,1])
box off




%% Sint-2000 LSTM

PH=SintPH;

M1=Kfold_Thresh(SintPH,.1:.0025:1,'Sint2000',5,'Thresh');


% Cost function adjustment

load('Results/Sint2000_5_fold_h50.mat')

sigma=std(PartResults(2:end,:));

L=unique(tvMCCResults(2,:));
MCCvL=zeros(length(L),length(PH));

for ii=1:length(PH)
    idx=find(tvMCCResults(1,:)==PH(ii));
    temp=tvMCCResults(:,idx);
    for jj=1:length(L)
        idx=find(temp(2,:)==L(jj));
        MCCvL(jj,ii)=temp(end,idx);
    end
end

Top=MCCvL'+2*sigma;
idx=find(isnan(Top));
idx=[idx,length(Top)+1];
%Plus2sig(isnan(Plus2sig))=1e1;
Bot=MCCvL'-2*sigma;
%Minus2sig(isnan(Minus2sig))=-1e1;

load('Results/Sint2000_0_Full_LSTM.mat')



Top=max(PartResults(2:end,:));
Bot=min(PartResults(2:end,:));

figure
fill([L,fliplr(L)],[Top,fliplr(Bot)],[k1,k1,k1])
hold on
%plot(L,PartResults(2:end,:),'o','Color',Color(CSint,:),'MarkerEdgeColor',Color(CSint,:),'MarkerFaceColor',Color(CSint,:),'MarkerSize',MS)
%scatter(reshape([L;L;L;L;L],25*5,1),reshape(PartResults(2:end,:),25*5,1),MS2,Color(CSint,:),'filled')
scatter(reshape([L;L;L;L;L],25*5,1),reshape(PartResults(2:end,:),25*5,1),MS2,[k2,k2,k2],'filled')
%alpha(at)
%fill([L,fliplr(L)],[Plus2sig,fliplr(Minus2sig)],[k1,k1,k1])
plot(L,MCCvL(:,1),'-o','Color',Color(CSint,:),'MarkerEdgeColor',Color(CSint,:),'MarkerFaceColor',Color(CSint,:),'MarkerSize',MS)
plot(L,tvMCCResults(end,:),'-d','Color',Color(CSint,:),'MarkerEdgeColor',Color(CSint,:),'MarkerSize',MS3)
%plot(L,MCCvL(:,1),'-o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerFaceColor',Color(CPAD,:),'MarkerSize',MS)
%plot(L,tvMCCResults(end,:),'-o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerSize',MS)
% plot(L,PartResults(2,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(3,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(4,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(5,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
% plot(L,PartResults(6,:),'o','MarkerEdgeColor',Color(1,:),'MarkerFaceColor',Color(1,:),'MarkerSize',MS)
plot(L,M1(1).*ones(length(L)),'--','Color',Color(CSint,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('Sint-2000','FontWeight','Normal')
set(gca,'fontsize',FigFont);
%legend('PH=1.6 kyr}','\rho^{Clim}_{16,1456}','\rho^{Ens}_{16,961}','\rho^{Clim}_{16,961}')%'Location','SouthEast')
ylabel('MCC')
xlabel('Search window (kyr)')
xlim([min(L),max(L)])
ylim([0,1])
box off

%% Models k-fold

clear all 

LW=2;
MS=6;
MS2=50;
MS3=9;
at=.5;
k1=.95;
k2=.8;
FigFont=16;
C1=3;
Colors
%%

G12PH=3.2;
ThreeDPH=16.4;
PADPH=11.7;
SintPH=10.2;
P09PH=6;
MB19PH=16.1;

CSint=7;
CPAD=8;
CG12=9;
C3D=1;

CSint=9;
CPAD=2;
CG12=1;
C3D=5;
CP09=3;
CMB19=4;


%%

PH=G12PH;

M1=Kfold_Thresh(3.2,.1:.0025:1,'G12Kfold',5,'Thresh');


% Cost function adjustment

load('Results/G12Kfold_5_fold_lin.mat')

sigma=std(PartResults(2:end,:));

L=unique(tvMCCResults(2,:));
MCCvL=zeros(length(L),length(PH));

for ii=1:length(PH)
    idx=find(tvMCCResults(1,:)==PH(ii));
    temp=tvMCCResults(:,idx);
    for jj=1:length(L)
        idx=find(temp(2,:)==L(jj));
        MCCvL(jj,ii)=temp(end,idx);
    end
end

Top=MCCvL'+2*sigma;
idx=find(isnan(Top));
idx=[idx,length(Top)+1];
%Plus2sig(isnan(Plus2sig))=1e1;
Bot=MCCvL'-2*sigma;
%Minus2sig(isnan(Minus2sig))=-1e1;

load('Results/G12Kfold_0_Full.mat')

Top=max(PartResults(2:end,:));
Bot=min(PartResults(2:end,:));



G12L=L;
G12Top=Top;
G12Bot=Bot;
G12PartResults=PartResults;
G12M1=M1;
G12MCCvL=MCCvL;
G12tvMCCResults=tvMCCResults;

%%

PH=ThreeDPH;

M1=Kfold_Thresh(3.2,.1:.0025:1,'ThreeDKfold',5,'Thresh');


% Cost function adjustment

load('Results/ThreeDKfold_5_fold_lin.mat')

sigma=std(PartResults(2:end,:));

L=unique(tvMCCResults(2,:));
MCCvL=zeros(length(L),length(PH));

for ii=1:length(PH)
    idx=find(tvMCCResults(1,:)==PH(ii));
    temp=tvMCCResults(:,idx);
    for jj=1:length(L)
        idx=find(temp(2,:)==L(jj));
        MCCvL(jj,ii)=temp(end,idx);
    end
end

Top=MCCvL'+2*sigma;
idx=find(isnan(Top));
idx=[idx,length(Top)+1];
%Plus2sig(isnan(Plus2sig))=1e1;
Bot=MCCvL'-2*sigma;
%Minus2sig(isnan(Minus2sig))=-1e1;

load('Results/ThreeDKfold_0_Full.mat')

Top=max(PartResults(2:end,:));
Bot=min(PartResults(2:end,:));



figure
fill([L,fliplr(L)],[Top,fliplr(Bot)],[k1,k1,k1])
hold on
scatter(reshape([L;L;L;L;L],length(L)*5,1),reshape(PartResults(2:end,:),length(L)*5,1),MS2,[k2,k2,k2],'filled')
plot(L,MCCvL(:,1),'-o','Color',Color(C3D,:),'MarkerEdgeColor',Color(C3D,:),'MarkerFaceColor',Color(C3D,:),'MarkerSize',MS)
plot(L,tvMCCResults(end,:),'-d','Color',Color(C3D,:),'MarkerEdgeColor',Color(C3D,:),'MarkerSize',MS3)
plot(L,M1(1).*ones(length(L)),'--','Color',Color(C3D,:),'LineWidth',LW);
fill([G12L,fliplr(G12L)],[G12Top,fliplr(G12Bot)],[k1,k1,k1])
hold on
scatter(reshape([G12L;G12L;G12L;G12L;G12L],length(G12L)*5,1),reshape(G12PartResults(2:end,:),length(G12L)*5,1),MS2,[k2,k2,k2],'filled')
plot(G12L,G12MCCvL(:,1),'-o','Color',Color(CG12,:),'MarkerEdgeColor',Color(CG12,:),'MarkerFaceColor',Color(CG12,:),'MarkerSize',MS)
plot(G12L,G12tvMCCResults(end,:),'-d','Color',Color(CG12,:),'MarkerEdgeColor',Color(CG12,:),'MarkerSize',MS3)
plot(G12L,G12M1(1).*ones(length(G12L)),'--','Color',Color(CG12,:),'LineWidth',LW);
set(gca,'fontsize',FigFont);
title('K-fold cross validation','FontWeight','Normal')
set(gca,'fontsize',FigFont);
ylabel('MCC')
xlabel('Search window (kyr)')
xlim([min(L),max(L)])
ylim([0,1])
box off

%% Autocorrelation plots

clear all 

LW=2;
MS=6;
FigFont=20;
C1=3;
Colors
%%

G12PH=3.2;
ThreeDPH=16.4;
PADPH=11.7;
SintPH=10.2;

CSint=7;
CPAD=8;
CG12=9;
C3D=1;

CSint=9;
CPAD=2;
CG12=1;
C3D=5;
CP09=3;
CMB19=6;

L=100;

load G12.mat
t=t(1:5:length(t));
d=d(1:5:length(d));
dt=t(2)-t(1);
Lags=round(L/dt);
time=dt.*(0:1:Lags);
d=abs(d);
dlag=d(2:end)-d(1:end-1);
aclag=autocorr(dlag,'NumLags',Lags);

figure
stem(time,aclag,'-o','Color',Color(CG12,:),'MarkerEdgeColor',Color(CG12,:),'MarkerFaceColor',Color(CG12,:),'MarkerSize',MS)
hold on
plot(time,zeros(length(time),1),'k');
set(gca,'fontsize',FigFont);
title('G12','FontWeight','Normal')
set(gca,'fontsize',FigFont);
%legend('Dipole','\Delta Dipole')%'Location','SouthEast')
ylabel('Correlation')
xlabel('Lag (kyr)')
xlim([0,max(L)])
ylim([-.2,1])
box off


load Sint2000.mat
%t=t(1:5:length(t));
%d=d(1:5:length(d));
dt=t(2)-t(1);
Lags=round(L/dt);
time=dt.*(0:1:Lags);
d=abs(d);
dlag=d(2:end)-d(1:end-1);
ac=autocorr(d,'NumLags',Lags);
aclag=autocorr(dlag,'NumLags',Lags);

figure
stem(time,aclag,'-o','Color',Color(CSint,:),'MarkerEdgeColor',Color(CSint,:),'MarkerFaceColor',Color(CSint,:),'MarkerSize',MS)
hold on
plot(time,zeros(length(time),1),'k');
set(gca,'fontsize',FigFont);
title('Sint-2000','FontWeight','Normal')
set(gca,'fontsize',FigFont);
%legend('Dipole','\Delta Dipole')%'Location','SouthEast')
ylabel('Correlation')
xlabel('Lag (kyr)')
xlim([0,max(L)])
ylim([-.2,1])
box off

load PADM2Mstar.mat
%t=t(1:5:length(t));
%d=d(1:5:length(d));
dt=t(2)-t(1);
Lags=round(L/dt);
time=dt.*(0:1:Lags);
d=abs(d);
dlag=d(2:end)-d(1:end-1);
ac=autocorr(d,'NumLags',Lags);
aclag=autocorr(dlag,'NumLags',Lags);

figure
stem(time,aclag,'-o','Color',Color(CPAD,:),'MarkerEdgeColor',Color(CPAD,:),'MarkerFaceColor',Color(CPAD,:),'MarkerSize',MS)
hold on
plot(time,zeros(length(time),1),'k');
set(gca,'fontsize',FigFont);
title('PADM2M','FontWeight','Normal')
set(gca,'fontsize',FigFont);
%legend('Dipole','\Delta Dipole')%'Location','SouthEast')
ylabel('Correlation')
xlabel('Lag (kyr)')
xlim([0,max(L)])
ylim([-.2,1])
box off

load ThreeD.mat
t=t(1:23:length(t));
d=d(1:23:length(d));
dt=t(2)-t(1);
Lags=round(L/dt);
time=dt.*(0:1:Lags);
d=abs(d);
dlag=d(2:end)-d(1:end-1);
ac=autocorr(d,'NumLags',Lags);
aclag=autocorr(dlag,'NumLags',Lags);

figure
stem(time,aclag,'-o','Color',Color(C3D,:),'MarkerEdgeColor',Color(C3D,:),'MarkerFaceColor',Color(C3D,:),'MarkerSize',MS)
hold on
plot(time,zeros(length(time),1),'k');
set(gca,'fontsize',FigFont);
title('3-D','FontWeight','Normal')
set(gca,'fontsize',FigFont);
%legend('Dipole','\Delta Dipole')%'Location','SouthEast')
ylabel('Correlation')
xlabel('Lag (kyr)')
xlim([0,max(L)])
ylim([-.2,1])
box off

load MB19.mat
dt=t(2)-t(1);
Lags=round(L/dt);
time=dt.*(0:1:Lags);
d=abs(d);
dlag=d(2:end)-d(1:end-1);
aclag=autocorr(dlag,'NumLags',Lags);

figure
stem(time,aclag,'-o','Color',Color(CMB19,:),'MarkerEdgeColor',Color(CMB19,:),'MarkerFaceColor',Color(CMB19,:),'MarkerSize',MS)
hold on
plot(time,zeros(length(time),1),'k');
set(gca,'fontsize',FigFont);
title('DW','FontWeight','Normal')
set(gca,'fontsize',FigFont);
%legend('Dipole','\Delta Dipole')%'Location','SouthEast')
ylabel('Correlation')
xlabel('Lag (kyr)')
xlim([0,max(L)])
ylim([-.2,1])
box off

load P09.mat
dt=t(2)-t(1);
Lags=round(L/dt);
time=dt.*(0:1:Lags);
d=abs(d);
dlag=d(2:end)-d(1:end-1);
aclag=autocorr(dlag,'NumLags',Lags);

figure
stem(time,aclag,'-o','Color',Color(CP09,:),'MarkerEdgeColor',Color(CP09,:),'MarkerFaceColor',Color(CP09,:),'MarkerSize',MS)
hold on
plot(time,zeros(length(time),1),'k');
set(gca,'fontsize',FigFont);
title('P09','FontWeight','Normal')
set(gca,'fontsize',FigFont);
%legend('Dipole','\Delta Dipole')%'Location','SouthEast')
ylabel('Correlation')
xlabel('Lag (kyr)')
xlim([0,max(L)])
ylim([-.2,1])
box off






