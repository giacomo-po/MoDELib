clear all
clc
close all

m=10;
s=0.1;

system(['./logNormalTest ' num2str(log(m)) ' ' num2str(s)])

data=load('probability.txt');

x=data(:,1);
y=[1:size(data,1)]'/size(data,1);
clear data

figure(1)
plot(x,y)

figure(2)
plot(x(1:end-1),diff(y)./diff(x))

figure(3)
plot(log(x(1:end-1)),diff(y)./diff(x))