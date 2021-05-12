clc
clear all
close all

U=load('U/U_0.txt');
S=load('S/S_0.txt');

%% Plot displacement
comp={'x_1','x_2','x_3','u_1','u_2','u_3'};
figure(1)
hold on
axis equal

clrCol=1;
uMax=max(U(:,clrCol))
uMin=min(U(:,clrCol))
caxis([uMin uMax])
%clrID=round((U(:,clrCol)-uMin)/(uMax-uMin)*size(colormap,1));

def=1;



for f=1:3:size(S,1)
%    plot3(U(f:f+2,1),U(f:f+2,2),U(f:f+2,3),'Color',[0.5 0.5 0.5])
%    fill3(U(f:f+2,1)+def*U(f:f+2,4),U(f:f+2,2)+def*U(f:f+2,5),U(f:f+2,3)+def*U(f:f+2,6),clrID(f:f+2))
     fill3(U(f:f+2,1)+def*U(f:f+2,4),U(f:f+2,2)+def*U(f:f+2,5),U(f:f+2,3)+def*U(f:f+2,6),U(f:f+2,clrCol))

end
grid on
xlabel('X_1','FontSize',14)
ylabel('X_2','FontSize',14)
zlabel('X_3','FontSize',14)
colorbar
title(comp{clrCol},'FontSize',14)

%return
%% Plot stress
comp={'x_1','x_2','x_3','\sigma_{11}','\sigma_{22}','\sigma_{33}','\sigma_{12}','\sigma_{23}','\sigma_{13}'};
figure(2)
clf
hold on
axis equal

clrCol=7;
s=S(:,clrCol);

slip=[0.000e+00 5.774e-01 8.165e-01];
n=[-8.165e-01  4.714e-01 -3.333e-01];

%slip=[1 0 0];
%n=[0 1 0];

s=s*0;
for i=1:3
    for j=1:3
        col=tensor2voigt(i,j)+3;
        s=s+slip(i)*S(:,col)*n(j);
    end
end

sMax=max(s)
sMin=min(s)
caxis([sMin sMax])
%clrID=round((s-sMin)/(sMax-sMin)*size(colormap,1));

for f=1:3:size(S,1)
%    patch(S(f:f+2,1),S(f:f+2,2),S(f:f+2,3),clrID(f:f+2))
    patch(S(f:f+2,1),S(f:f+2,2),S(f:f+2,3),s(f:f+2))
end
grid on
xlabel('X_1','FontSize',14)
ylabel('X_2','FontSize',14)
zlabel('X_3','FontSize',14)
colorbar
title(comp{clrCol},'FontSize',14)


function col=tensor2voigt(i,j)
ids=[1 4 5;
     4 2 6;
     5 6 3];
 col=ids(i,j);
end

