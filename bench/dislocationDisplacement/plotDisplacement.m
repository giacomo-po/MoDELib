close all
clear all
clc

fontSize=14;

Dn=load('D/D_0.txt');
Da=load('D/D_1.txt');

origin=find(~isfinite(Da(:,4)));
Da(origin,[4:6])=0;
%Da(origin,[4:6])=0;

df=1;
useOffset=0;

offX=mean(Da(:,4)-Dn(:,4))
offY=mean(Da(:,5)-Dn(:,5))

figure(1)
clf
plot(Da(:,1)+df*Da(:,4),Da(:,2)+df*Da(:,5),'o')
hold on
plot(Dn(:,1)+df*(Dn(:,4)+offX*useOffset),Dn(:,2)+df*(Dn(:,5)+offY*useOffset),'rx')
axis equal 
grid on
xlabel('x_1 / b','Fontsize',fontSize)
ylabel('x_2 / b','Fontsize',fontSize)
legend('analytical','numerical')
set(gca,'FontSize',fontSize)

figure(2)
plot(Da(:,1),Da(:,5)-Dn(:,5),'.')