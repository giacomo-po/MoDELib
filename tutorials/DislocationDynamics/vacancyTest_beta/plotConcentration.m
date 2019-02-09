clear all
close all
clc

data=load('concentration.txt');

figure(1)
plot(data(:,1),data(:,4))
grid on