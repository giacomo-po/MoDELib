function plotCells
clc
close all
%clear


%% Load files
fileID=0;
cellSize=1;
P=load(['P/P_' num2str(fileID) '.txt']);
C=load(['C/C_' num2str(fileID) '.txt']);

nCells=size(C,1)

%% Prepare colors
colormap jet
cMap=colormap;
rankMax=max(C(:,6))+1 % one-based
rankColorIDs=round([1:rankMax]/rankMax*size(cMap,1));
cellColors=cMap(rankColorIDs,:);



%% Processor efficiency
rankLoads=zeros(rankMax,1);
for c=1:size(C,1)
    r=C(c,6)+1;
    w=C(c,4)*C(c,5);
    rankLoads(r)=rankLoads(r)+w;
    %color=cellColors(C(c,6)+1,:);
    %disp(color)
    %plotCube(C(c,[1:3]),cellSize,color)
end
%totalLoad=sum(rankLoads);
optimumLoad=mean(rankLoads);
inefficiency=max(rankLoads-optimumLoad)/optimumLoad;

fontSize=18;

figure(2)
clf
bar([1:rankMax]-1,rankLoads)
hold on 
plot(C(:,6),C(:,4).*C(:,5),'ro')
plot([1:rankMax]-1,[1:rankMax]*0+optimumLoad,'g--','Linewidth',2)
xlabel('rank','FontSize',fontSize)
ylabel('weight (number of interactions)','FontSize',fontSize)
grid on
legend('cumulative processor weight','individual cell weight','optimum (average) weight','location','northwest')
set(gca,'FontSize',fontSize)
title(['load imbalance = ' num2str(inefficiency*100) '%'])

%% Plot particles and cells
figure(1)
clf
plot3(P(:,2),P(:,3),P(:,4),'k.')
hold on
for c=1:size(C,1)
    color=cellColors(C(c,6)+1,:);
    plotCube(C(c,[1:3]),cellSize,color)
end
grid on
axis equal


function plotCube(cellID,cellSize,faceColor)

%rank=rank+1;
%rank
%colors(rank,:)
%figure(1)
transparency=1.0;
LW=1;

P=[cellID+[0 0 0];
   cellID+[1 0 0];
   cellID+[1 1 0];
   cellID+[0 1 0];
   cellID+[0 0 0]]*cellSize;
p=plot3(P(:,1)',P(:,2)',P(:,3)','color',faceColor,'linewidth',LW);
%set(p,'facecolor',faceColor,'edgecolor','none','FaceAlpha',transparency);

%p=patch(P(:,1),P(:,2),P(:,3),ones(size(P,1),1));
%set(p,'facecolor',faceColor,'edgecolor','none','FaceAlpha',transparency);

P=[cellID+[0 0 1];
   cellID+[1 0 1];
   cellID+[1 1 1];
   cellID+[0 1 1];
   cellID+[0 0 1]]*cellSize;
p=plot3(P(:,1)',P(:,2)',P(:,3)','color',faceColor,'linewidth',LW);
%p=patch(P(:,1),P(:,2),P(:,3),ones(size(P,1),1));
%set(p,'facecolor',faceColor,'edgecolor','none','FaceAlpha',transparency);

P=[cellID+[0 0 0];
   cellID+[0 0 1]]*cellSize;
p=plot3(P(:,1)',P(:,2)',P(:,3)','color',faceColor,'linewidth',LW);

%p=patch(P(:,1),P(:,2),P(:,3),ones(size(P,1),1));
%set(p,'facecolor',faceColor,'edgecolor','none','FaceAlpha',transparency);

P=[cellID+[1 0 0];
   cellID+[1 0 1]]*cellSize;
p=plot3(P(:,1)',P(:,2)',P(:,3)','color',faceColor,'linewidth',LW);

%p=patch(P(:,1),P(:,2),P(:,3),ones(size(P,1),1));
%set(p,'facecolor',faceColor,'edgecolor','none','FaceAlpha',transparency);

P=[cellID+[1 1 0];
   cellID+[1 1 1]]*cellSize;
p=plot3(P(:,1)',P(:,2)',P(:,3)','color',faceColor,'linewidth',LW);

%p=patch(P(:,1),P(:,2),P(:,3),ones(size(P,1),1));
%set(p,'facecolor',faceColor,'edgecolor','none','FaceAlpha',transparency);

P=[cellID+[0 1 0];
   cellID+[0 1 1]]*cellSize;
p=plot3(P(:,1)',P(:,2)',P(:,3)','color',faceColor,'linewidth',LW);

%p=patch(P(:,1),P(:,2),P(:,3),ones(size(P,1),1));
%set(p,'facecolor',faceColor,'edgecolor','none','FaceAlpha',transparency);

%set(gca,'Alpha',0.5)
%alpha(0.8);
