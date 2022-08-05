clc
clear
close all

format long

alignToSlipSystem0=0;
lattice='fcc';

if(alignToSlipSystem0)
switch lattice
    case 'bcc'
        s=[5.773502691896258e-01  5.773502691896258e-01 -5.773502691896258e-01];
        n=[8.660254037844385e-01 0.000000000000000e+00 8.660254037844385e-01];
        n=n/norm(n);
    case 'fcc'
        s=[ 0.000000000000000e+00 7.071067811865475e-01 7.071067811865475e-01];
        n=[-7.071067811865476e-01  7.071067811865476e-01 -7.071067811865476e-01];
        n=n/norm(n);
    otherwise
end
    C2G1=[s;cross(n,s);n]
    else
    C2G1=eye(3)
end

switch lattice
    case 'bcc'
        A=1/sqrt(3)*[-1 1 1;1 -1 1;1 1 -1]; % bcc
        material='W';
    case 'fcc'
        A=1/sqrt(2)*[0 1 1;1 0 1;1 1 0]; % fcc
        material='Cu';
    otherwise
end
A=C2G1*A;


g=0.0;

%Fs=1000*eye(3);
Fs=diag([1000 1000 4000]);

F12=[1 g 0;
    0 1 0;
    0 0 1];

F31=[1 0 0;
    0 1 0;
    g 0 1];

F23=[1 0 0;
    0 1 g;
    0 0 1];

F=F12*F23*F31

Ba=F*Fs;
Ma=inv(A)*Ba
M=round(Ma)
B=A*M

x0=[0 0 0];


fID=fopen('polycrystal.txt','w');
fprintf(fID,['materialFile=../../../MaterialsLibrary/' material '.txt; \n']);
fprintf(fID,'enablePartials=0; \n');
fprintf(fID,'absoluteTemperature = 300; # [K] simulation temperature \n');
fprintf(fID,'meshFile=../../../MeshLibrary/unitCube.msh; \n');
fprintf(fID,'C2G1=');
fprintf(fID,'%1.15f ',C2G1(1,:));
fprintf(fID,'\n ');
fprintf(fID,'%1.15f ',C2G1(2,:));
fprintf(fID,'\n ');
fprintf(fID,'%1.15f ',C2G1(3,:));
fprintf(fID,';\n ');
fprintf(fID,'\n ');
fprintf(fID,'\n ');

fprintf(fID,'A=');
fprintf(fID,'%1.15f ',B(1,:));
fprintf(fID,'\n ');
fprintf(fID,'%1.15f ',B(2,:));
fprintf(fID,'\n ');
fprintf(fID,'%1.15f ',B(3,:));
fprintf(fID,';\n ');
fprintf(fID,'\n ');
fprintf(fID,'\n ');

fprintf(fID,'x0=');
fprintf(fID,'%1.15f ',x0(1,:));
fprintf(fID,';\n ');

fprintf(fID,'periodicFaceIDs=');
fprintf(fID,'%d ',[0 1 2 3 4 5]);
fprintf(fID,';\n ');


 


