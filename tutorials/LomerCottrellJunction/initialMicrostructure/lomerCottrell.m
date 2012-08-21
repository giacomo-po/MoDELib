
L=500;    % initila length on the straigth segments [b]
theta=50; % angle between dislocations and common line[deg]

figure(1)
addpath('../../../../matlab/')
[N,V,R]=thompson_tetrahedron(L/2/sqrt(6));

Bd=[1 -1 0]/sqrt(2);    % Burgers vector on delta plane
Bb=[-1 0 -1]/sqrt(2);   % Burgers vector on beta plane

delta=4;   % the column into N corresponding to the delta plane
beta=2;    % the column into N corresponding to the  beta plane
A=V(:,4);  % the A-vertex of the thompson tetrahedron
B=V(:,3);  % the B-vertex of the thompson tetrahedron
P0=0.5*(A+B)
plot3(A(1),A(2),A(3),'sm')
plot3(B(1),B(2),B(3),'sg')
plot3(P0(1),P0(2),P0(3),'sk')


AB=B-A;

Pd(:,1)=P0+angleAxis(N(:,delta),-theta*pi/180)*AB/2;
Pd(:,2)=P0-angleAxis(N(:,delta),-theta*pi/180)*AB/2;



for k=2:size(Pd,2)
quiver3(Pd(1,k-1),Pd(2,k-1),Pd(3,k-1),[Pd(1,k)-Pd(1,k-1)],[Pd(2,k)-Pd(2,k-1)],[Pd(3,k)-Pd(3,k-1)],0,'m','Linewidth',2)
end


Pb(:,1)=P0+angleAxis(N(:,beta),theta*pi/180)*AB/2;
Pb(:,2)=P0-angleAxis(N(:,beta),theta*pi/180)*AB/2;


for k=2:size(Pb,2)
quiver3(Pb(1,k-1),Pb(2,k-1),Pb(3,k-1),[Pb(1,k)-Pb(1,k-1)],[Pb(2,k)-Pb(2,k-1)],[Pb(3,k)-Pb(3,k-1)],0,'m','Linewidth',2)
end

file_V = fopen('../V/V_0.txt','w');
nodeID=0;
nodeformat='%i %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %i \n';
fprintf(file_V,nodeformat, [0 Pd(:,1)' zeros(1,3) 0]);
fprintf(file_V,nodeformat, [1 Pd(:,2)' zeros(1,3) 0]);
fprintf(file_V,nodeformat, [2 Pb(:,1)' zeros(1,3) 0]);
fprintf(file_V,nodeformat, [3 Pb(:,2)' zeros(1,3) 0]);
fclose(file_V)

file_E = fopen('../E/E_0.txt','w');
linkformat='%i %i %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %i \n';

fprintf(file_E,linkformat, [0 1 Bd zeros(1,3) 0]);
fprintf(file_E,linkformat, [2 3 Bb zeros(1,3) 0]);

fclose(file_E);

