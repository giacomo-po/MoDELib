clc
clear all
close all

U=load('D/D_0.txt');

%% Plot displacement
comp={'x_1','x_2','x_3','u_1','u_2','u_3'};
figure(1)
hold on
axis equal

clrCol=4;
uMax=max(U(:,clrCol))
uMin=min(U(:,clrCol))
clrID=round((U(:,clrCol)-uMin)/(uMax-uMin)*size(colormap,1));

%def=1;

for f=1:3:size(U,1)
%    plot3(U(f:f+2,1),U(f:f+2,2),U(f:f+2,3),'Color',[0.5 0.5 0.5])
    fill3(U(f:f+2,1),U(f:f+2,2),U(f:f+2,3),clrID(f:f+2))
end
grid on
xlabel('X_1','FontSize',14)
ylabel('X_2','FontSize',14)
zlabel('X_3','FontSize',14)
title(comp{clrCol},'FontSize',14)

