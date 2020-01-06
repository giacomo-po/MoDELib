clear all
close all
clc



L1=10;
L2=20;

nx=10;
ny=10;

X=(rand(3,1)-[1;1;1]*0.5)*30
X(3)=abs(X(3));
d=-10;
X=[1*L1/2;L2/2;d];

%X(3)=0;

P=[0  0  0;
   L1 0  0;
   L1 L2 0;
   0  L2 0];

figure(1)
grid on
axis equal
hold on
plot3([P(:,1);P(1,1)],[P(:,2);P(1,2)],[P(:,3);P(1,3)])
plot3(X(1),X(2),X(3),'ro')

dx=L1/nx;
dy=L2/ny;
dA=dx*dy;
S=0;
Acheck=0;
for i=1:nx
for j=1:ny
    X1=[-dx/2+i*dx;-dy/2+j*dy;0];
    plot3(X1(1),X1(2),X1(3),'.')
    r=X-X1;
    S=S+r(3)/norm(r)^3*dA;
    Acheck=Acheck+dA;
end
end

Acheck=Acheck/L1/L2

% syms x1 y1 real
% X1=[x1 y1 0];
% R=X-X1;
% 
% ix=int(R(3)/norm(R)^3,x1,0,L1);
% ixy=int(ix,y1,0,L2)

I=0;
P=P';
q=[0 1 1];
q=q/norm(q);
for i=1:4
    i1=i+1;
if i1==5
    i1=1;
end
    %[i i1]
    Rn=P(:,i)-X;
    Rm=P(:,i1)-X;
    Rn=Rn/norm(Rn);
    Rm=Rm/norm(Rm);
    
    RdR=dot(Rn,Rm);
    RcR=cross(Rn,Rm);
    RcR=RcR/norm(RcR);
    A=sqrt((1-RdR)/(1+RdR));
    
    num=dot(q,RcR)*A;
    den=1-dot(q,Rn)-dot(q,cross(RcR,Rn))*A;
    
    I=I+2*atan(num/den);
%        I=I+2*atan(num/den);

%    I=I+2*atan(A);

end



Irect=4*asin(L1*L2/sqrt((L1^2+4*d^2)*(L2^2+4*d^2)));
[S I Irect]

quiver3(X(1),X(2),X(3),q(1),q(2),q(3))