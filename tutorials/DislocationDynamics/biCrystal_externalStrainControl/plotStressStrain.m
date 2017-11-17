close all
clear all
clc

F=load('F/F_0.txt');

e33=F(:,end-9);
s33=F(:,end);

figure(1)
plot(e33,s33,'x-')