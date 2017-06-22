close all
clear all
clc
P0=load('P/P_0.txt');
P1=load('P/P_1.txt');
P2=load('P/P_2.txt');

figure(1)
clf
hold on
plot(P0(:,1),P0(:,2),'o')
plot(P1(:,1),P1(:,2),'r.')
plot(P2(:,1),P2(:,2),'g.')

axis equal