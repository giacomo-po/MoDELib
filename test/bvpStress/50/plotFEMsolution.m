clc
clear all
close all

fileID=1070;
U=load(['U/U_' num2str(fileID) '.txt']);
S=load(['S/S_' num2str(fileID) '.txt']);
Y=load(['Y/Y_' num2str(fileID) '.txt']);
Z=load(['Z/Z_' num2str(fileID) '.txt']);

%% Plot displacement
comp={'x_1','x_2','x_3','u_1','u_2','u_3'};
figure(1)
hold on
axis equal

%plot3(Y(:,1),Y(:,2),Y(:,3),'g.')


clrCol=6;
uMax=max(U(:,clrCol))
uMin=min(U(:,clrCol))
clrID=round((U(:,clrCol)-uMin)/(uMax-uMin)*size(colormap,1));

def=200;



for f=1:3:size(S,1)
    %    plot3(U(f:f+2,1),U(f:f+2,2),U(f:f+2,3),'Color',[0.5 0.5 0.5])
    fill3(U(f:f+2,1)+def*U(f:f+2,4),U(f:f+2,2)+def*U(f:f+2,5),U(f:f+2,3)+def*U(f:f+2,6),clrID(f:f+2))
end
grid on
xlabel('X_1','FontSize',14)
ylabel('X_2','FontSize',14)
zlabel('X_3','FontSize',14)
title(comp{clrCol},'FontSize',14)

%return
%% Plot stress
comp={'x_1','x_2','x_3','\sigma_{11}','\sigma_{22}','\sigma_{33}','\sigma_{12}','\sigma_{23}','\sigma_{13}'};
figure(2)
hold on
axis equal

clrCol=6;
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
title(comp{clrCol},'FontSize',14)

%%
%%%%%%%%%%%%%%%%%%%%%
figure(3)
clf
hold on
R=1000;
for f=1:3:size(S,1)
    theta=angle(S(f:f+2,1)+i*S(f:f+2,2));
    dTheta=theta-theta(1);
    id=find(abs(dTheta)>pi);
    if length(id)>0
        for k=id
            theta(k)=theta(k)+sign(theta(1))*2*pi;
        end
    end
    fill(R*theta,S(f:f+2,3),clrID(f:f+2))
    %drawnow
end

%%
figure(4)
id=find(S(:,3)==4000);
plot3(S(id,1),S(id,2),S(id,clrCol),'.')
xlabel('X1')
ylabel('X2')

H=4000;

figure(5)
clf
hold on

curCol=6;
id=find(abs(Y(:,3)-H)<0.001);
plot3(Y(id,1),Y(id,2),Y(id,curCol),'rx')
% quiver3(Y(id,1),Y(id,2),Y(id,3)*0,Y(id,7),Y(id,8),Y(id,9),0.01,'k')
%axis equal


%hold on
%plot3(Y(id,7),Y(id,8),Y(id,curCol),'gs')

id=find(abs(U(:,3)-H)<0.001);
plot3(U(id,1),U(id,2),U(id,curCol),'.')
%xlabel('X1')
%ylabel('X2')
%plot3(9.018e+02, 3.205e+02, 0,'ms')


plotCol=4;


figure(6)
%plot3(Z(:,1),Z(:,2),Z(:,3),'.')
%xlabel('X1')
%ylabel('X2')
R=1000;
theta=angle(Z(:,1)+i*Z(:,2));
plot3(R*theta,Z(:,3),Z(:,plotCol),'.')


figure(7)
%plot3(Z(:,1),Z(:,2),Z(:,3),'.')
%xlabel('X1')
%ylabel('X2')
id=find(abs(Z(:,3)-H)<0.001);
plot3(Z(id,1),Z(id,2),Z(id,plotCol),'.')

 %figure(8)
 id=find(abs(Y(:,3)-H)<0.001);

 %quiver3(Y(id,1),Y(id,2),Y(id,3),Y(id,7),Y(id,8),Y(id,9),0,'k')
for n=1:size(Y,1)
aa=norm(Y(n,[7:9]));
if abs(aa-1)>0.001
    format long
    aa
error('norm os s verctor different than 1')
end
end

% plot3(Y(:,1),Y(:,2),Y(:,3),'rx')
% hold on
% plot3(Y(:,7),Y(:,8),Y(:,9),'b.')
% norm(Y(:,1)-Y(:,7))
% norm(Y(:,2)-Y(:,8))
% norm(Y(:,3)-Y(:,9))
