clear all
close all
clc

fontSize=16;
F=load('F_0.txt');

runID=F(:,1);
dt=F(:,2);
ddLength=F(:,3);
 
if size(F,2)>=24 % using remote stress field
	externalstrain33=F(:,end-12);
	externalstress33=F(:,end-3);
	figure(10)
	plot(externalstrain33,externalstress33);
end

u3=F(:,end-2);
stress33=F(:,end-1);
stress33DD=F(:,end);

 

L=1000; % cube side length
V=L^3;
A=L^2;

e33=u3/L;
s33=stress33;

figure(1)
subplot(2,1,1)
plot(-e33*100,-s33)
v=axis;
axis([0 max(-e33*100) 0 v(4)])
xlabel('\epsilon_{33} [%]','Fontsize',fontSize)
ylabel('\sigma_{33}/\mu','Fontsize',fontSize)
grid on
set(gca,'Fontsize',fontSize)

subplot(2,1,2)
plot(-e33*100,ddLength)
xlabel('\epsilon_{33} [%]','Fontsize',fontSize)
ylabel('line length / b','Fontsize',fontSize)
v=axis;
axis([0 max(-e33*100) 0 v(4)])
grid on
set(gca,'Fontsize',fontSize)

figure(2)
subplot(2,1,1)
plot(runID,F(:,6:end)/V)
v=axis;
axis([0 max(runID) v(3) v(4)])
grid on
xlabel('DD step','Fontsize',fontSize)
ylabel('d\epsilon^P/dt [mu/B]','Fontsize',fontSize)
legend('\epsilon^P_{11}','\epsilon^P_{12}','\epsilon^P_{13}','\epsilon^P_{22}','\epsilon^P_{23}','\epsilon^P_{33}','location','northwest')
set(gca,'Fontsize',fontSize)

subplot(2,1,2)
plot(runID,cumsum(F(:,6:end).*repmat(dt,1,size(F(:,6:end),2))/V))
v=axis;
axis([0 max(runID) v(3) v(4)])
grid on
xlabel('DD step','Fontsize',fontSize)
ylabel('\epsilon^P','Fontsize',fontSize)
legend('d\epsilon^P_{11}/dt','d\epsilon^P_{12}/dt','d\epsilon^P_{13}/dt','d\epsilon^P_{22}/dt','d\epsilon^P_{23}/dt','d\epsilon^P_{33}/dt','location','northwest')
set(gca,'Fontsize',fontSize)

%return 
figure(3)
subplot(2,1,1)
plot(runID,-s33)
grid on

subplot(2,1,2)
plot(runID,dt)
axis([0 max(runID) 0 1.5*max(dt)])
grid on
xlabel('DD step')
ylabel('dt')

t=cumsum(dt);
t=t-t(1);
figure(4)
plot(t,u3)
hold on
plot(t,-1.0e-9*t*1000,'r')
grid on
