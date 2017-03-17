L=2;
np=101;
dx=L/(np-1);
P=[-L/2:dx:L/2];

%% Write node file N/N_0.txt
file_N = fopen('N/N_0.txt','w');
nodeformat='%i %1.15e\n';
for n=1:np
fprintf(file_N,nodeformat, [n-1 P(n)]);
end
fclose(file_N)

%% Write element file T/T_0.txt
file_T = fopen('T/T_0.txt','w');
nodeformat='%i %i %i %i \n';
region=0;
for n=1:np-1
fprintf(file_T,nodeformat, [n-1 n-1 n region]);
end
fclose(file_T)