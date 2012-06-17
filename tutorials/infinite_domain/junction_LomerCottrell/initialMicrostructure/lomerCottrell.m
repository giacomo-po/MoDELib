%function FCC_loop_generator
clc
clear
close all


addpath('../../../../matlab/')
L=1000;
[N,V,R]=thompson_tetrahedron(L/2/sqrt(6));

delta=4;
beta=2;
A=V(:,4);
B=V(:,3);
P0=0.5*(A+B)
plot3(A(1),A(2),A(3),'sm')
plot3(B(1),B(2),B(3),'sg')
plot3(P0(1),P0(2),P0(3),'sk')

AB=B-A;

theta=30; % deg


Pd(:,1)=P0-angleAxis(N(:,delta),theta*pi/180)*AB/2;
Pd(:,2)=P0+angleAxis(N(:,delta),theta*pi/180)*AB/2;






for k=2:size(Pd,2)
quiver3(Pd(1,k-1),Pd(2,k-1),Pd(3,k-1),[Pd(1,k)-Pd(1,k-1)],[Pd(2,k)-Pd(2,k-1)],[Pd(3,k)-Pd(3,k-1)],0,'m','Linewidth',2)
end


Pb(:,1)=P0-angleAxis(N(:,beta),-theta*pi/180)*AB/2;
Pb(:,2)=P0+angleAxis(N(:,beta),-theta*pi/180)*AB/2;


angle=acosd(dot((Pd(:,2)-Pd(:,1))/L,(Pb(:,2)-Pb(:,1))/L))

for k=2:size(Pb,2)
quiver3(Pb(1,k-1),Pb(2,k-1),Pb(3,k-1),[Pb(1,k)-Pb(1,k-1)],[Pb(2,k)-Pb(2,k-1)],[Pb(3,k)-Pb(3,k-1)],0,'m','Linewidth',2)
end

file_V = fopen('../V/V_0.txt','w');
nodeID=0;
nodeformat='%i %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %i \n';
for n=1:size(Pd,2)
fprintf(file_V,nodeformat, [nodeID Pd(:,n)' zeros(1,3) 0]);
nodeID=nodeID+1;
end
for n=1:size(Pb,2)
fprintf(file_V,nodeformat, [nodeID Pb(:,n)' zeros(1,3) 0]);
nodeID=nodeID+1;
end
fclose(file_V)

file_E = fopen('../E/E_0.txt','w');
linkformat='%i %i %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %i \n';
Bd=[1 -1 0]/sqrt(2);
nodeID=0;
for n=1:size(Pd,2)-1
fprintf(file_E,linkformat, [nodeID nodeID+1 Bd zeros(1,3) 0]);
nodeID=nodeID+1;
end
nodeID=size(Pd,2);
Bb=[-1 0 -1]/sqrt(2);
for n=1:size(Pb,2)-1
fprintf(file_E,linkformat, [nodeID nodeID+1 Bb zeros(1,3) 0]);
nodeID=nodeID+1;
end
fclose(file_E);

