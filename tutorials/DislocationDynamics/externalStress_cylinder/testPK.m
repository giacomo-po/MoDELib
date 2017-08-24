clear all
clc
close all

sf=1000;
Ph=load('P_hr0/P_0.txt');
Pc=load('P_cr0/P_0.txt');
P=load('P/P_0.txt');

V=load('V/V_0.txt');

quiver3(Ph(:,4),Ph(:,5),Ph(:,6),sf*Ph(:,7),sf*Ph(:,8),sf*Ph(:,9),0)
hold on
quiver3(Pc(:,4),Pc(:,5),Pc(:,6),sf*Pc(:,7),sf*Pc(:,8),sf*Pc(:,9),0,'r')
quiver3(P(:,4),P(:,5),P(:,6),sf*P(:,7),sf*P(:,8),sf*P(:,9),0,'k')

plot3(V(:,2),V(:,3),V(:,4),'mo')


u=[1:1000]/1000;
figure(2)
plot(u,-6*u+6*u.^2)