clc
close all
clear all

A=[0 1 1;
   1 0 1;
   1 1 0]/sqrt(2);

v1=[1 -1 0];
v1=v1/norm(v1);
v3=[1 1 1];
v3=v3/norm(v3);
C2G=[v1;cross(v3,v1); v3]

A=C2G*A;  
%invA=inv(A);

g1=A(:,1)-A(:,3); % primitive lattice vector on the 111 plane
g2=A(:,2)-A(:,3); % primitive lattice vector on the 111 plane
G=[g1 g2];

%% Initial configuration
b=[1 0 0];
%r=200;
np=11;
r=1900;

%return
theta=[0:np-1]/np*2*pi;

XL=[r*cos(theta);r*sin(theta);zeros(1,np)];
XL=G*round(inv(G'*G)*G'*XL);

%XL(:,[1])=XL(:,[1])*0.75;

figure(1)
hold on
plot3([XL(1,:) XL(1,1)],[XL(2,:) XL(2,1)],[XL(3,:) XL(3,1)])
axis equal

% Write inout files E and V
file_V = fopen('V/V_0.txt','w');
file_E = fopen('E/E_0.txt','w');

nodeformat='%i %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %i \n';
for n=1:np
fprintf(file_V,nodeformat, [n-1 XL(:,n)' zeros(1,3) 0]);
end

linkformat='%i %i %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %i \n';
for n=1:np-1
fprintf(file_E,linkformat, [n-1 n b zeros(1,3) 0]);
end
fprintf(file_E,linkformat, [n 0 b zeros(1,3) 0]);

fclose(file_V);
fclose(file_E);