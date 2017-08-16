close all
clear all
clc

fileID=0;
P=load(['P/P_' num2str(fileID) '.txt']);

figure(1)
plot3(P(:,4),P(:,5),P(:,6),'r.')
hold on
quiver3(P(:,4),P(:,5),P(:,6),P(:,7),P(:,8),P(:,9),10)