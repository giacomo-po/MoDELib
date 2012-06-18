%function FCC_loop_generator
clc
clear
close all

Ra=2000;
H=6*Ra;
r=0.9*Ra;

addpath('../../../../../matlab/')

[N,V,R]=thompson_tetrahedron(1);

delta=4;
beta=2;

np=11;
theta=[0:np-1]/np*2*pi;
P0=[Ra Ra H/3]';
x=P0(1)+r*cos(theta);
y=P0(2)+r*sin(theta);
z=P0(3)-((x-P0(1))*N(1,delta)+(y-P0(2))*N(2,delta))/N(3,delta);
Pd=[x' y' z'];

figure(1)
hold on
plot3(x,y,z)


P0=[Ra Ra H*2/3]';
x=P0(1)+r*cos(theta);
y=P0(2)+r*sin(theta);
z=P0(3)-((x-P0(1))*N(1,beta)+(y-P0(2))*N(2,beta))/N(3,beta);
Pb=[x' y' z'];
plot3(x,y,z)


file_V = fopen('../V/V_0.txt','w');
nodeID=0;
nodeformat='%i %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %i \n';
for n=1:size(Pd,1)
fprintf(file_V,nodeformat, [nodeID Pd(n,:) zeros(1,3) 0]);
nodeID=nodeID+1;
end
for n=1:size(Pb,1)
fprintf(file_V,nodeformat, [nodeID Pb(n,:) zeros(1,3) 0]);
nodeID=nodeID+1;
end
fclose(file_V)

file_E = fopen('../E/E_0.txt','w');
linkformat='%i %i %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %i \n';
Bd=[1 -1 0]/sqrt(2);
nodeID=0;
for n=0:np-1
fprintf(file_E,linkformat, [nodeID rem(nodeID+1,np) Bd zeros(1,3) 0]);
nodeID=nodeID+1;
end

nodeID=0;
Bb=[-1 0 -1]/sqrt(2);
for n=0:np-1
fprintf(file_E,linkformat, [nodeID+np rem(nodeID+1,np)+np  Bb zeros(1,3) 0]);
nodeID=nodeID+1;
end
fclose(file_E);

