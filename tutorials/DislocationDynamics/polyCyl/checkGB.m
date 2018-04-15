clear all
close all
clc

figure(1)
clf
hold on
N=42;
for f=1:N
P=load(['file' num2str(f) '.txt']);
patch(P(:,1),P(:,2),P(:,3),'g')
end
axis equal