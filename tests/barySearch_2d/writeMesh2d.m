clear all
close all
clc

meshID=0;
Lx=1;
Ly=1;
nx=11;
ny=11;
dx=Lx/(nx-1);
dy=Ly/(ny-1);
Px=[-Lx/2:dx:Lx/2];
Py=[-Ly/2:dy:Ly/2];

nRegions=4;
%regionClrs={'g','m','b'}

%% Write node file N/N_x.txt
file_N = fopen(['N/N_' num2str(meshID) '.txt'],'w');
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

%% Write element file T/T_x.txt
figure(1)
hold on
axis equal

file_T = fopen(['T/T_' num2str(meshID) '.txt'],'w');
nodeformat='%i %i %i %i %i \n';

n=0;
E=[];
for j=1:ny-1
for i=1:nx-1
    if i<nx/2
    E=[E;[i+(j-1)*nx i+1+(j-1)*nx i+j*nx]];
    else
    E=[E;[i+j*nx i+(j-1)*nx i+1+j*nx]];
    end
    Pm=mean(P(E(end,:),:));
    if((Pm(2)<0 || abs(Pm(2)/Pm(1))<1) && Pm(1)<0)
    region=1;
    elseif ((Pm(2)<0 || abs(Pm(2)/Pm(1))<1) && Pm(1)>0)
            region=2;
    else
        region=3;
    end
fprintf(file_T,nodeformat, [n E(end,:)-1  region]);
n=n+1;
patch(P(E(end,:),1),P(E(end,:),2),scalar2color(region,0,nRegions))

if i<nx/2
    E=[E;[i+j*nx i+1+(j-1)*nx i+1+j*nx]];
else
        E=[E;[i+(j-1)*nx i+1+(j-1)*nx i+1+j*nx]];

end

    Pm=mean(P(E(end,:),:));
    if((Pm(2)<0 || abs(Pm(2)/Pm(1))<1) && Pm(1)<0)
    region=1;
    elseif ((Pm(2)<0 || abs(Pm(2)/Pm(1))<1) && Pm(1)>0)
            region=2;
    else
        region=3;
    end
fprintf(file_T,nodeformat, [n E(end,:)-1 region]);
n=n+1;
patch(P(E(end,:),1),P(E(end,:),2),scalar2color(region,0,nRegions))

end
end
% for n=1:np-1
% fprintf(file_T,nodeformat, [n-1 n-1 n region]);
% end
fclose(file_T)
