clear all
close all
clc

for n=0:400
D=load(['Z/Z_' num2str(n) '.txt']);
LN(n+1)=sum(D(:,3));
end

plot(LN)