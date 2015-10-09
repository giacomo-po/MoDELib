v1=[1 -1 0];
v1=v1/norm(v1);
v3=[1 1 1];
v3=v3/norm(v3);
C2G=[v1;cross(v3,v1); v3]

A=[0 1 1;
   1 0 1;
   1 1 0]/sqrt(2);

A=C2G*A;
invA=inv(A);

meshID=0;
N=load(['N/N_' num2str(meshID)]);

N1=N(:,2:end)';

% P=A*n
% n=invA*P
N1=A*round(invA*N1);
N=[N(:,1) N1'];

