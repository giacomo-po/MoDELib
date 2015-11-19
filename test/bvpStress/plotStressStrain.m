clear all
close all
clc

epsilonDot=1.0e-11;
Burgers=0.2556e-9; % Burgers vector for Cu [m]
L1=2000 % radius of cylinder (units of Burgers vector)
L2=2000
H=4000
R=1000
%A=L1*L2;
A=pi*R^2;
V=H*A;

files=[50];
clrs=rand(length(files),3);


fontSize=16;
for k=1:length(files)
    file=files(k)
F=load(['./' num2str(file) '/F/F_0.txt']);

runID=F(:,1);
time=F(:,2);
dt=F(:,3);

dc=18;
u3=F(:,3+dc+1);
theta3=F(:,3+dc+2);
f3=F(:,3+dc+3);
t3=F(:,3+dc+4);
u3_DD=F(:,3+dc+5);
f3_DD=F(:,3+dc+6);
pdr=F(:,13:21)/V;

figure(1)
hold on
plot(runID,dt,'Color',clrs(k,:))
legend(num2str(files'))

figure(2)
hold on
plot(u3/H,f3/A,'Color',clrs(k,:),'Linewidth',2)
legend(num2str(files'),'Location','NorthWest')
%plot(u3/H,f3_DD/A,'--','Color',clrs(k,:))
v=axis;
axis([0 v(2) 0 v(4)])

figure(3)
hold on
plot(runID,f3/A,'Color',clrs(k,:))
xlabel('runID','FontSize',fontSize)
ylabel('\sigma','FontSize',fontSize)
set(gca,'FontSize',fontSize)
legend(num2str(files'))

figure(4)
hold on
plot(runID,u3_DD,'Color',clrs(k,:))
plot(runID,epsilonDot*time*H,'--','Color',clrs(k,:))
plot(runID,epsilonDot*time*H-u3_DD,'-.','Color',clrs(k,:))
%plot(u3/H,u3_DD,'Color',clrs(k,:))
%plot(u3/H,epsilonDot*time*H,'--','Color',clrs(k,:))
%plot(u3/H,epsilonDot*time*H-u3_DD,'-.','Color',clrs(k,:))
%plot(u3/H,[0;diff(epsilonDot*time*H-u3_DD)],'-.','Color',clrs(k,:))

xlabel('runID','FontSize',fontSize)
ylabel('u_{DD}','FontSize',fontSize)
set(gca,'FontSize',fontSize)
legend(num2str(files'))

figure(5)
hold on
plot(runID,pdr(:,9),'Color',clrs(k,:))

figure(6)
hold on
plot(runID,(f3-f3_DD)/A,'Color',clrs(k,:))
xlabel('runID','FontSize',fontSize)
ylabel('\sigma_{FEM}','FontSize',fontSize)
set(gca,'FontSize',fontSize)
legend(num2str(files'))
end

