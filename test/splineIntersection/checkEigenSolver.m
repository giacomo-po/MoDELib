clear all
close all
clc

AB=load('AB.txt');
eigSize=size(AB,2);
A=AB(1:eigSize,:);
B=AB(eigSize+1:2*eigSize,:);
clear AB

[V,D]=eig(A,B);

for k=1:eigSize
s(k,1)=D(k,k);
t(k,1)=V(2,k)/V(1,k);
end

[t s]

%return

C=inv(B)*A;

[V1,D1]=eig(A,B);
for k=1:eigSize
s1(k,1)=D1(k,k);
t1(k,1)=V1(2,k)/V(1,k);
end

[t1 s1]

COEFFs=load('coeffFile.txt');
coeff1=COEFFs(1:2,:);
coeff2=COEFFs(3:4,:);

u=[0:0.001:1]';

u1=0.6647;
u2=1;

P1=zeros(length(u),2);
P2=zeros(length(u),2);
P1c=[0 0];
P2c=[0 0];
for k=1:size(COEFFs,2)
P1=P1+u.^(k-1)*coeff1(:,k)';
P2=P2+u.^(k-1)*coeff2(:,k)';
P1c=P1c+u1^(k-1)*coeff1(:,k)';
P2c=P2c+u2^(k-1)*coeff2(:,k)';
end
figure(1)
clf
plot(P1(:,1),P1(:,2))
hold on
plot(P2(:,1),P2(:,2),'r')
plot(P1c(:,1),P1c(:,2),'x')
plot(P2c(:,1),P2c(:,2),'ro')
grid on
axis equal