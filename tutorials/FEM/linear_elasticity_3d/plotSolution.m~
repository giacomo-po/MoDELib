clc
clear all
close all

U=load('U/U_0.txt');
S=load('S/S_0.txt');

%% Plot displacement 
figure(1)
hold on
axis equal

clrCol=4;
uMax=max(U(:,clrCol))
uMin=min(U(:,clrCol))
clrID=round((U(:,clrCol)-uMin)/(uMax-uMin)*size(colormap,1));

def=1;

for f=1:3:size(S,1)
fill3(U(f:f+2,1)+def*U(f:f+2,4),U(f:f+2,2)+def*U(f:f+2,5),U(f:f+2,3)+def*U(f:f+2,6),clrID(f:f+2))
plot3(U(f:f+2,1),U(f:f+2,2),U(f:f+2,3),'k')
end
grid on
xlabel('X_1','FontSize',14)
ylabel('X_2','FontSize',14)
zlabel('X_3','FontSize',14)

%% Plot stress 
figure(2)
hold on
axis equal

clrCol=9;
sMax=max(S(:,clrCol))
sMin=min(S(:,clrCol))
clrID=round((S(:,clrCol)-sMin)/(sMax-sMin)*size(colormap,1));

for f=1:3:size(S,1)
fill3(S(f:f+2,1),S(f:f+2,2),S(f:f+2,3),clrID(f:f+2))
end
grid on
xlabel('X_1','FontSize',14)
ylabel('X_2','FontSize',14)
zlabel('X_3','FontSize',14)


