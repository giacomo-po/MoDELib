close all
clear all
clc

X=[-1.0,-1.0,0.0;
 0.0,-1.0,0.0;
 0.0, 0.0,0.0;
-1.0, 0.0,0.0;
-1.0, 0.0,0.0;
 0.0, 0.0,0.0;
 0.0, 1.0,0.0;
-1.0, 1.0,0.0;
 0.0, 0.0,0.0;
 1.0, 0.0,0.0;
 1.0, 1.0,0.0;
 0.0, 1.0,0.0];


P0=load('P/P_0.txt');
P1=load('P/P_1.txt');
P2=load('P/P_2.txt');
P3=load('P/P_3.txt');


figure(1)
clf
hold on
plot(P0(:,1),P0(:,2),'o')
plot(P1(:,1),P1(:,2),'r.')
plot(P2(:,1),P2(:,2),'g.')
plot(P3(:,1),P3(:,2),'kx')
plot(X(:,1),X(:,2),'ms')

text(X(:,1),X(:,2),num2str([0:size(X,1)-1]'),'FontSize',24)


axis equal