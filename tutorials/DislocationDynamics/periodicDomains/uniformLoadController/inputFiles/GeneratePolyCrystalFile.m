clc
clear
close all

format long

alignToSlipSystem0=1;
lattice='fcc';

if(alignToSlipSystem0)
switch lattice
    case 'bcc'
        s=[1 1 -1];
        n=[1 0 1];
    case 'fcc'
        s=[0 1 1];
        n=[-1 1 -1];
        n=n/norm(n);
    otherwise
end
        s=s/norm(s);
        n=n/norm(n);
        C2G1=[s;cross(n,s);n]
        %C2G1=[s;n;cross(s,n)]

    else
    C2G1=eye(3)
end

switch lattice
    case 'bcc'
        A=1/sqrt(3)*[-1 1 1;1 -1 1;1 1 -1]; % bcc
        material='W';
    case 'fcc'
        A=1/sqrt(2)*[0 1 1;1 0 1;1 1 0]; % fcc
%        material='AlMg15';
%        material='Al';
        material='Ni';
    otherwise
end
A=C2G1*A;


g=0.0;

%Fs=1000*eye(3);
Fs=diag([1000 1000 1000]);

F12=[1 g 0;
    0 1 0;
    0 0 1];

F31=[1 0 g;
    0 1 0;
    0 0 1];

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

fprintf(fID,'\n#################\n');
fprintf(fID,'# GlidePlaneNoise\n');
fprintf(fID,'gridSize=');
fprintf(fID,'%d ',[256 256]);
fprintf(fID,'; # size of grid on the glide plane\n ');
fprintf(fID,'gridSpacing_SI=');
fprintf(fID,'%d ',[1 1]*1e-10);
fprintf(fID,'; # [m] spacing of grid on the glide plane\n ');

fprintf(fID,'solidSolutionNoiseMode=0; # 0=no noise, 1= read noise, 2=compute noise\n');
fprintf(fID,'solidSolutionNoiseFile_xz=../../../NoiseLibrary/noise_xz.vtk;\n');
fprintf(fID,'solidSolutionNoiseFile_yz=../../../NoiseLibrary/noise_yz.vtk;\n');
fprintf(fID,'stackingFaultNoiseMode=0; # 0=no noise\n');



