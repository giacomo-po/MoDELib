clear all
close all
clc

Burgers=0.2851e-9; % Burgers vector for Al [m]
R=0.5e-6/2/Burgers % radius of cylinder (units of Burgers vector)
H=4*R;    % height of cylinder (units of Burgers vector)
A=pi*R^2;
V=H*A;
J=pi*R^4/2;

fontSize=16;
F=load('F_0.txt');

runID=F(:,1);
time=F(:,2);
dt=F(:,3);
ddLength=F(:,4);
ddImmobLength=F(:,5);
u3=F(:,6);
theta3=F(:,7); % rotation angle [rad]
f3=F(:,8); % force
t3=F(:,9); % torque
pdr=F(:,10:end)/V; % plastic distortion rate
pd=cumsum(pdr.*repmat(dt,1,size(pdr,2))); % plastic distortion



e33=u3/H;
s33=f3/A;

gamma=theta3*R/H;
tau=t3*R/J;


%return

%%
system('mkdir tensionMovie')
figure(1)
for k=1%:size(F,1)
    k
    clf
    plot(e33*100,s33,'Linewidth',2)
    hold on
    v=axis;
    axis([0 max(e33*100) 0 v(4)])
    plot(e33(k)*100,s33(k),'ro','Linewidth',2)
    xlabel('\epsilon [%]','Fontsize',fontSize)
    ylabel('\sigma / \mu [-]','Fontsize',fontSize)
    grid on
    set(gca,'Fontsize',fontSize)
    print(gcf,'-djpeg100',['./tensionMovie/image_' num2str(runID(k))])
end

system('mkdir torsionMovie')
figure(2)
for k=1%:size(F,1)
    k
    clf
    plot(gamma*100,tau,'Linewidth',2)
    %set(gca,'Fontsize',fontSize)
    %plot(e33*100,s33,'Linewidth',2)
    hold on
    v=axis;
    axis([0 max(gamma*100) 0 v(4)])
    plot(gamma(k)*100,tau(k),'ro','Linewidth',2)
    xlabel('\gamma [%]','Fontsize',fontSize)
    ylabel('\tau / \mu [-]','Fontsize',fontSize)
    grid on
    set(gca,'Fontsize',fontSize)
    print(gcf,'-djpeg100',['./torsionMovie/image_' num2str(runID(k))])
end

%%
figure(3)
hold on
plot(e33*100,ddLength/V/Burgers^2,'r','Linewidth',2)
plot(e33*100,ddImmobLength/V/Burgers^2,'b','Linewidth',2)
plot(e33*100,(ddLength-ddImmobLength)/V/Burgers^2,'k','Linewidth',2)
xlabel('\epsilon (or \gamma) [%]','Fontsize',fontSize)
ylabel('dislocaiton density [m^{-2}]','Fontsize',fontSize)
v=axis;
axis([0 max(e33*100) 0 v(4)])
grid on
legend('total','immobile+boundary','mobile','Location','NorthWest')
set(gca,'Fontsize',fontSize)

figure(4)
subplot(2,1,1)
plot(runID,pdr)
v=axis;
axis([0 max(runID) v(3) v(4)])
grid on
xlabel('DD step','Fontsize',fontSize)
ylabel('d\beta^P/dt [mu/B]','Fontsize',fontSize)
legend('d\beta^P_{11}/dt','d\beta^P_{12}/dt','d\beta^P_{13}/dt','d\beta^P_{21}/dt','d\beta^P_{22}/dt','d\beta^P_{23}/dt','d\beta^P_{31}/dt','d\beta^P_{32}/dt','d\beta^P_{33}/dt','location','northwest')
set(gca,'Fontsize',fontSize)

subplot(2,1,2)
plot(runID,pd)
v=axis;
axis([0 max(runID) v(3) v(4)])
grid on
xlabel('DD step','Fontsize',fontSize)
ylabel('\beta^P','Fontsize',fontSize)
legend('\beta^P_{11}','\beta^P_{12}','\beta^P_{13}','\beta^P_{21}','\beta^P_{22}','\beta^P_{23}','\beta^P_{31}','\beta^P_{32}','\beta^P_{33}','location','northwest')
set(gca,'Fontsize',fontSize)

%return
figure(5)
subplot(2,1,1)
plot(runID,s33)
grid on

subplot(2,1,2)
plot(runID,dt)
axis([0 max(runID) 0 1.5*max(dt)])
grid on
xlabel('DD step')
ylabel('dt')

%t=cumsum(dt);
%t=t-t(1);
figure(6)
plot(runID,u3,'Linewidth',2)
hold on
plot(runID,1.0e-9*time*H,'r','Linewidth',1)
grid on
xlabel('DD step','Fontsize',fontSize)
ylabel('u_3 / b','Fontsize',fontSize)
legend('measured','prescribed','Location','NorthWest')
set(gca,'Fontsize',fontSize)

figure(7)
%plot(t,theta3)
plot(runID,theta3,'b','Linewidth',2)
hold on
plot(runID,1.0e-9*time*H/R,'r','Linewidth',1)
grid on
xlabel('DD step ','Fontsize',fontSize)
ylabel('theta_3 [rad]','Fontsize',fontSize)
legend('measured','prescribed','Location','NorthWest')
set(gca,'Fontsize',fontSize)

%return

%%

