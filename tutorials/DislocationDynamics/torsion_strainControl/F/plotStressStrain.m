clear all
close all
clc

fontSize=16;

Burgers=0.2556e-9; % Burgers vector for Cu [m]
L1=2000 % side length of pillar
L2=2000 % side length of pillar
H=4000  % heigth of pillar
A=L1*L2;
V=H*A;

dotTheta3=1e-9;

fontSize=16;
F=load('F_0.txt');

runID=F(:,1);
time=F(:,2);
dt=F(:,3);


theta3=F(:,end-1);
t3=F(:,end);

figure(1)
plot(time,theta3,'r')
hold on
plot(time,dotTheta3*time,'b')


figure(2)
plot(theta3,t3)
grid on
xlabel('twist angle [rad]','FontSize',fontSize)
ylabel('torque [\mu b^3]','FontSize',fontSize)

