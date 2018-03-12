clear all
close all
clc

F=load('F/F_0.txt');
s33=F(:,21);
e33=F(:,12);
figure(1)
plot(e33,s33)
%plot(s33)

for n=0:F(end,1)
D=load(['Z/Z_' num2str(n) '.txt']);
notTouchID=find(D(:,3)==0);
LNa(n+1)=sum(abs(D(notTouchID,4)));
LN(n+1)=sum(D(notTouchID,4));
end

figure(2)
clf
hold on
plot(e33,LN)
plot(e33,LNa)

