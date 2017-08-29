close all
clear all
clc

F0=load('F_0.txt');
F1=load('F_1.txt');
F2=load('F_2.txt'); % no openMP
F3=load('F_3.txt'); % no openMP
F4=load('F_4.txt'); % no openMP
F5=load('F_5.txt'); % no openMP
F6=load('F_10.txt'); % no openMP
F7=load('F_11.txt'); % no openMP

col=3;
figure(1)
plot(F0(:,col));
hold on
plot(F1(:,col),'r');
plot(F2(:,col),'m');
plot(F2(:,col),'k');

figure(2)
%plot(F0(:,col)-F1(:,col))
hold on
%plot(F2(:,col)-F3(:,col),'m')
%plot(F4(:,col)-F5(:,col),'k')
plot(F6(:,col)-F7(:,col),'k')
id=find(abs(F6(:,col)-F7(:,col))>1e-12);
plot(id,id*0,'or')