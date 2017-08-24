%clear all
%close all
%clc
function main
clc
%[t,x]=ode45(@overDamped,[0 100],1.6);

Np=200
x=rand(1)*10
f=pkForce(x,0,-1);
B=1;
dt=100;
t=0;
for k=1:Np

    v=1/B*pkForce(x(end),0,-1);
    dt=5/abs(v);
    dx=v*dt
    x=[x;x(end)+dx];
t=[t;t(end)+dt];
end
figure(1)
xx=[-10:0.1:10];
fx=pkForce(xx,xx*0,-1);
plotyy(x,t,xx,fx)

function xDot=overDamped(t,x) 
B=1;
b=-1;
%xi=;
fx=pkForce(x,0,b);
xDot=1/B*fx;

function fx=pkForce(x,y,b)
b1=1;
mu=1;
nu=0.3;
%x=xV(1);
%y=xV(2);
D=mu*b1/(2*pi*(1-nu));
R2=(x.^2+y.^2+b1^2).^2;
%sxx=-D*y.*(3*x.^2+y.^2)./R2;
%syy=D*y.*(x.^2-y.^2)./R2;
sxy=D*x.*(x.^2-y.^2)./R2;
%S=[sxx sxy 0;
%   sxy syy 0;
%   0    0  0];
%b1=[b 0 0]';
%xi=[0 0 1]';
%f=cross(S*b,[0 0 1]');
fx=sxy*b;