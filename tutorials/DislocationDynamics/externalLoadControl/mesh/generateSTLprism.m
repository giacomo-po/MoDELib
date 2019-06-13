% This file is part of MODEL, the Mechanics Of Defect Evolution Library.
%
% Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
% model is distributed without any warranty under the
% GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.


clc
close all
clear all

MODEL_DIR='../../../../';
addpath([MODEL_DIR '/matlab/']);

%% Define output file name
meshID=0;
filename='cube'; % this creates file cube.stl
nElements=1e4;

%% Size and position of the cube
L1=4000; % the side length of the cube, in units of Burgers vector
L2=4000; % the side length of the cube, in units of Burgers vector
L3=4000; % the side length of the cube, in units of Burgers vector

%% Compute element size
V=L1*L2*L3;
averageElementVolume=V/nElements;

% coordinates of the 8 vertices of the cube.
% The base of the cube is at z=0. The cube is centered in x and y
v0=[0 0 0]-[L1 L2 L3]/2*0;
v1=[L1 0 0]-[L1 L2 L3]/2*0;
v2=[L1 L2 0]-[L1 L2 L3]/2*0;
v3=[0 L2 0]-[L1 L2 L3]/2*0;
v4=[0 0 L3]-[L1 L2 L3]/2*0;
v5=[L1 0 L3]-[L1 L2 L3]/2*0;
v6=[L1 L2 L3]-[L1 L2 L3]/2*0;
v7=[0 L2 L3]-[L1 L2 L3]/2*0;
v8=[L1/2 L2/2 0]-[L1 L2 L3]/2*0;
v9=[L1/2 L2/2 L3]-[L1 L2 L3]/2*0;
v10=[L1/2 0 L3/2]-[L1 L2 L3]/2*0;
v11=[L1/2 L2 L3/2]-[L1 L2 L3]/2*0;
v12=[0 L2/2 L3/2]-[L1 L2 L3]/2*0;
v13=[L1 L2/2 L3/2]-[L1 L2 L3]/2*0;


% Create 12 facets by splitting each face of the cube in two triangles
Facets={};
% Bottom facets
Facets=addFacet(Facets,v8,v0,v3);
Facets=addFacet(Facets,v8,v1,v0);
Facets=addFacet(Facets,v8,v2,v1);
Facets=addFacet(Facets,v8,v3,v2);
% Top facets
Facets=addFacet(Facets,v9,v4,v5);
Facets=addFacet(Facets,v9,v5,v6);
Facets=addFacet(Facets,v9,v6,v7);
Facets=addFacet(Facets,v9,v7,v4);
% Left facets
Facets=addFacet(Facets,v10,v4,v0);
Facets=addFacet(Facets,v10,v0,v1);
Facets=addFacet(Facets,v10,v1,v5);
Facets=addFacet(Facets,v10,v5,v4);
% Right facets
Facets=addFacet(Facets,v11,v2,v3);
Facets=addFacet(Facets,v11,v3,v7);
Facets=addFacet(Facets,v11,v7,v6);
Facets=addFacet(Facets,v11,v6,v2);
% Back facets
Facets=addFacet(Facets,v12,v0,v4);
Facets=addFacet(Facets,v12,v4,v7);
Facets=addFacet(Facets,v12,v7,v3);
Facets=addFacet(Facets,v12,v3,v0);
% Front facets
Facets=addFacet(Facets,v13,v1,v2);
Facets=addFacet(Facets,v13,v2,v6);
Facets=addFacet(Facets,v13,v6,v5);
Facets=addFacet(Facets,v13,v5,v1);


%% plot facets
figure(1)
clf
plotFacets(Facets)
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')

%% write stl file
writeSTL(Facets,filename)

%% Run Tetgen
system([MODEL_DIR '/scripts/tetgenSTL.sh ' filename ' ' num2str(averageElementVolume)]);

%% Create T and N files and clean tetgent output
system([MODEL_DIR '/scripts/tetgen2TN.sh ' filename ' ' num2str(meshID)]);

%%
B=[2.474873734152911e+02 0.000000000000000e+00 0.000000000000000e+00,0.000000000000000e+00 2.474873734152916e+02 0.000000000000000e+00
0.000000000000000e+00 2.474873734152911e+02 0.000000000000000e+00,0.000000000000000e+00 4.000000000000000e+03 3.752512626584709e+03
0.000000000000000e+00 4.000000000000000e+03 3.752512626584709e+03,2.474873734152910e+02 4.000000000000000e+03 4.000000000000000e+03
2.474873734152903e+02 4.000000000000000e+03 4.000000000000000e+03,4.000000000000000e+03 2.474873734152910e+02 4.000000000000000e+03
4.000000000000000e+03 2.474873734152903e+02 4.000000000000000e+03,4.000000000000000e+03 0.000000000000000e+00 3.752512626584709e+03
4.000000000000000e+03 0.000000000000000e+00 3.752512626584709e+03,2.474873734152916e+02 0.000000000000000e+00 0.000000000000000e+00]

for k=1:size(B,1)
plot3([B(k,1) B(k,4)],[B(k,2) B(k,5)],[B(k,3) B(k,6)],'r','Linewidth',2)
end

B=[2.964795672342894e+03 0.000000000000000e+00 4.000000000000000e+03,0.000000000000000e+00 0.000000000000000e+00 1.035204327657105e+03
0.000000000000000e+00 0.000000000000000e+00 1.035204327657106e+03,0.000000000000000e+00 1.035204327657105e+03 0.000000000000000e+00
0.000000000000000e+00 1.035204327657106e+03 0.000000000000000e+00,2.964795672342895e+03 4.000000000000000e+03 0.000000000000000e+00
2.964795672342895e+03 4.000000000000000e+03 0.000000000000000e+00,4.000000000000000e+03 4.000000000000000e+03 1.035204327657105e+03
4.000000000000000e+03 4.000000000000000e+03 1.035204327657105e+03,4.000000000000000e+03 1.035204327657105e+03 4.000000000000000e+03
4.000000000000000e+03 1.035204327657106e+03 4.000000000000000e+03,2.964795672342895e+03 0.000000000000000e+00 4.000000000000000e+03]

for k=1:size(B,1)
plot3([B(k,1) B(k,4)],[B(k,2) B(k,5)],[B(k,3) B(k,6)],'b','Linewidth',2)
end

plot3([P0(1) P2(1)],[P0(2) P2(2)],[P0(3) P2(3)],'m','Linewidth',4)