clear all
clc

Ap=eye(3);
Af=[0 1 1;
   1 0 1;
   1 1 0]/2;

A=Ap;



%Q=[12 -3 4;4 12 -3;-3 4 12];
%p=13;
Q=[5 0 0;0 4 -3; 0 3 4]
p=5
R=Q/p;
A1=R*A;
s=[0 1 0]';
n=[0 0 1]';
U=eye(3)+s*n';
%U=[1 0 1;0 1 0;0 1 1];
det(U)

T=eye(3)-U*inv(A)*R'*A
sigma=2/det(T)
X01=inv(T)
X0=A*inv(T)
n
det(T)
return
%inv(A)*Q*A
DSCL1=R/p;

%return
%A=eye(3);

%R=[12 -3 4;4 12 -3;-3 4 12]/13;
det(R)
%R'*R
%A1=R*A;
CSL=[-1  3  1;
     -1  0 -3;
     -1 -1  0];

DSCL=inv(CSL)';

figure(1)
clf
hold on
np=2
for i=[-np:np]
    for j=[-np:np]
        for k=[-np:np]
            P=i*A(:,1)+j*A(:,2)+k*A(:,3);
            P1=i*A1(:,1)+j*A1(:,2)+k*A1(:,3);
            %PC=i*C(:,1)+j*C(:,2)+k*C(:,3);

            plot3(P(1),P(2),P(3),'bo','Linewidth',2)
            plot3(P1(1),P1(2),P1(3),'rx','Linewidth',2)
            %plot3(PC(1),PC(2),PC(3),'g.','Linewidth',2)
        end
    end
end

np=4
for i=[-np:np]
    for j=[-np:np]
        for k=[-np:np]
            %P=i*A(:,1)+j*A(:,2)+k*A(:,3);
            %P1=i*A1(:,1)+j*A1(:,2)+k*A1(:,3);
            PD=i*DSCL(:,1)+j*DSCL(:,2)+k*DSCL(:,3);
            PD1=i*DSCL1(:,1)+j*DSCL1(:,2)+k*DSCL1(:,3);

            %plot3(P(1),P(2),P(3),'bo')
            %plot3(P1(1),P1(2),P1(3),'rx')
            plot3(PD(1),PD(2),PD(3),'k.','Linewidth',1)
                        plot3(PD1(1),PD1(2),PD1(3),'m.','Linewidth',1)

        end
    end
end

%return

np=2
for i=[-np:np]
    for j=[-np:np]
        for k=[-np:np]
            %P=i*A(:,1)+j*A(:,2)+k*A(:,3);
            %P1=i*A1(:,1)+j*A1(:,2)+k*A1(:,3);
            PC=i*CSL(:,1)+j*CSL(:,2)+k*CSL(:,3);

            %plot3(P(1),P(2),P(3),'bo','Linewidth',1)
            %plot3(P1(1),P1(2),P1(3),'rx','Linewidth',1)
            plot3(PC(1),PC(2),PC(3),'gs')
        end
    end
end

% A*n+R*A*m=C*k
% let R=Q/p, with Q integer matrix and p integer
% p*A*n+Q*A*m=p*C*k
% let A*B=Q*A
% B=inv(A)*Q*A
% p*A*n+A*B*m=p*C*k
% if B is integer, then B*m=M is integer
% p*A*n+A*M=p*C*k
% 
axis equal