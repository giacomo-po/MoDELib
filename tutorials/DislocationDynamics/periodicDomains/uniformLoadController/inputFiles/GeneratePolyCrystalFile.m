clc
clear
close all

format long

A=1/sqrt(3)*[-1 1 1;1 -1 1;1 1 -1];


g=0.1;

Fs=1000*eye(3);

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

C2G1=eye(3);
x0=[0 0 0];


fID=fopen('polycrystal.txt','w');
fprintf(fID,'materialFile=../../MaterialsLibrary/W.txt; \n');
fprintf(fID,'meshFile=../../MeshLibrary/unitMesh.msh; \n');
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
fprintf(fID,'%d ',[0 1 2]);
fprintf(fID,';\n ');


