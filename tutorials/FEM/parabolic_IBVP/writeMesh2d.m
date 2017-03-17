clear all
close all
clc

Lx=2;
Ly=2;
nx=51;
ny=51;
dx=Lx/(nx-1);
dy=Ly/(ny-1);
Px=[-Lx/2:dx:Lx/2];
Py=[-Ly/2:dy:Ly/2];

%% Write node file N/N_2.txt
file_N = fopen('N/N_2.txt','w');
nodeformat='%i %1.15e %1.15e\n';
n=0;
P=[];
for j=1:ny
for i=1:nx
    P=[P;[Px(i) Py(j)]];
fprintf(file_N,nodeformat, [n Px(i) Py(j)]);
n=n+1;
end
end
fclose(file_N)

%% Write element file T/T_2.txt
figure(1)
hold on
axis equal

file_T = fopen('T/T_2.txt','w');
nodeformat='%i %i %i %i %i \n';
region=0;
n=0;
E=[];
for j=1:ny-1
for i=1:nx-1
    E=[E;[i+(j-1)*nx i+1+(j-1)*nx i+j*nx]];
fprintf(file_T,nodeformat, [n E(end,:)-1  region]);
n=n+1;
patch(P(E(end,:),1),P(E(end,:),2),'g')
    E=[E;[i+j*nx i+1+(j-1)*nx i+1+j*nx]];
fprintf(file_T,nodeformat, [n E(end,:)-1 region]);
n=n+1;
patch(P(E(end,:),1),P(E(end,:),2),'m')

end
end
% for n=1:np-1
% fprintf(file_T,nodeformat, [n-1 n-1 n region]);
% end
fclose(file_T)
