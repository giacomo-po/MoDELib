clear all
close all
clc

% Burgers=0.2556e-9; % Burgers vector for Cu [m]
% L1=2000 % radius of cylinder (units of Burgers vector)
% L2=2000
% H=4000
% R=1000
% %A=L1*L2;
% A=pi*R^2;
% V=H*A;
% 
% eDot=1.0e-11;
% gDot=1.0e-11;


fontSize=16;
F=load('F_0.txt');

runID=F(:,1);
time=F(:,2);
dt=F(:,3);

dc=4;
u3=F(:,dc);
%theta3=F(:,3+dc+2);
f3=F(:,dc+1);
%t3=F(:,3+dc+4);

figure(1)
plot(runID,dt,'r')

figure(2)
plot(u3,f3)
grid on

% figure(3)
% plot(runID,f3/A)
