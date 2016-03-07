clear all
close all
clc

fontSize=16;

Burgers=0.2556e-9; % Burgers vector for Cu [m]
L1=2000 % radius of cylinder (units of Burgers vector)
L2=2000
H=4000
%R=1000
A=L1*L2;
%A=pi*R^2;
V=H*A;

%eDot=1.0e-11;
%gDot=1.0e-11;
dotTheta3=4e-11;

fontSize=16;
F=load('F_0.txt');

runID=F(:,1);
time=F(:,2);
dt=F(:,3);

dc=2;
%u3=F(:,3+dc+1);
theta3=F(:,3+dc+1);
%f3=F(:,3+dc+3);
t3=F(:,3+dc+2);

figure(1)
plot(time,theta3,'r')
hold on
plot(time,dotTheta3*time,'b')


figure(2)
plot(theta3,t3)
grid on
xlabel('twist angle [rad]','FontSize',fontSize)
ylabel('torque [\mu b^3]','FontSize',fontSize)

