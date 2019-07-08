close all
clear all
clc

a=rand(2,1)
b=rand(2,1);
t=(b-a)/norm(b-a);

N=50;
Ax=ones(N,N)*a(1);
Ay=ones(N,N)*a(2);
Bx=ones(N,N)*b(1);
By=ones(N,N)*b(2);
Tx=ones(N,N)*t(1);
Ty=ones(N,N)*t(2);

x=[0:N-1]/(N-1);
y=[0:N-1]/(N-1);
[X,Y]=meshgrid(x,y);

num=(Bx-X).*Tx+(By-Y).*Ty+sqrt((Bx-X).^2+(By-Y).^2);
den=(Ax-X).*Tx+(Ay-Y).*Ty+sqrt((Ax-X).^2+(Ay-Y).^2);
f=log(num./den);

figure(1)
clf
hold on
pcolor(X,Y,f)
plot([a(1) b(1)],[a(2) b(2)],'k')
axis equal

num1=-(Ax-X).*Tx-(Ay-Y).*Ty+sqrt((Ax-X).^2+(Ay-Y).^2);
den1=-(Bx-X).*Tx-(By-Y).*Ty+sqrt((Bx-X).^2+(By-Y).^2);
f1=log(num1./den1);

figure(2)
clf
hold on
pcolor(X,Y,f1)
plot([a(1) b(1)],[a(2) b(2)],'k')
axis equal

% figure(3)
% clf
% hold on
% surf(X,Y,f-f1)
% plot([a(1) b(1)],[a(2) b(2)],'k')
% axis equal