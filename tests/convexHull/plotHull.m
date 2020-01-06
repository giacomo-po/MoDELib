clear all
close all
clc
fontSize=16;

%N=300;
%inPoints=rand(N,2);
%dlmwrite('inPoints.txt', inPoints,' ');

system('./hull ')

inPoints=load('inPoints.txt');
outPoints=load('outPoints.txt');
hullSize=size(outPoints,1)

figure(1)
plot(inPoints(:,1),inPoints(:,2),'xb')
hold on
plot(outPoints(:,1),outPoints(:,2),'r')

