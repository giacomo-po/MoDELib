close all
clear all
clc

dim=3;


%% Load
P=load('P/P_0.txt');
D=load('D/D_0.txt');

%% Plot
p=round(rand(1)*(size(P,1)/dim));

figure(1)
clf
hold on
grid on
axis equal
plot3(P(dim*p+1,1:dim+1),P(dim*p+2,1:dim+1),P(dim*p+3,1:dim+1),'g')
plot3(P(dim*p+1,[1 3]),P(dim*p+2,[1 3]),P(dim*p+3,[1 3]),'g')
plot3(P(dim*p+1,[1 4]),P(dim*p+2,[1 4]),P(dim*p+3,[1 4]),'g')
plot3(P(dim*p+1,[2 4]),P(dim*p+2,[2 4]),P(dim*p+3,[2 4]),'g')
Pc=mean(P(dim*p+[1:dim],1:dim+1)')';
plot3(Pc(1),Pc(2),Pc(3),'ro')
set(gca,'View',[45 25])

for i=1:dim+1
    F1=[mean(P(dim*p+1,setdiff(1:dim+1,i))) mean(P(dim*p+2,setdiff(1:dim+1,i))) mean(P(dim*p+3,setdiff(1:dim+1,i)))];
    quiver3(F1(1),F1(2),F1(3),D(dim*p+1,i),D(dim*p+2,i),D(dim*p+3,i),0.001,'b')
end
