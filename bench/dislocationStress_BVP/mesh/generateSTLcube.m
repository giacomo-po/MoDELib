% This file is part of MODEL, the Mechanics Of Defect Evolution Library.
%
% Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
% model is distributed without any warranty under the
% GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.


clc
close all
clear all

MODEL_DIR='../../..';
addpath([MODEL_DIR '/matlab/']);

filename='cube'; % this creates file cube.stl
meshID=0; % creates ../N/N_0.txt and ../T/T_0.txt
L=400; % the side length of the cube, in units of Burgers vector
Lz=40*6*sqrt(3);
V=L^3;  % volume of the cube
nElements=5e4; % target number of mesh elements

% coordinates of the 8 vertices of the cube, centered at origin
v0=[0 0 0]-[L L Lz]/2;
v1=[L 0 0]-[L L Lz]/2;
v2=[L L 0]-[L L Lz]/2;
v3=[0 L 0]-[L L Lz]/2;
v4=[0 0 Lz]-[L L Lz]/2;
v5=[L 0 Lz]-[L L Lz]/2;
v6=[L L Lz]-[L L Lz]/2;
v7=[0 L Lz]-[L L Lz]/2;

% Create 12 facets by splitting each face of the cube in two triangles
Facets={};
Facets=addFacet(Facets,v0,v3,v1);
Facets=addFacet(Facets,v1,v3,v2);

Facets=addFacet(Facets,v4,v5,v7);
Facets=addFacet(Facets,v5,v6,v7);

Facets=addFacet(Facets,v0,v1,v5);
Facets=addFacet(Facets,v0,v5,v4);

Facets=addFacet(Facets,v0,v7,v3);
Facets=addFacet(Facets,v0,v4,v7);

Facets=addFacet(Facets,v3,v7,v2);
Facets=addFacet(Facets,v7,v6,v2);

Facets=addFacet(Facets,v1,v2,v6);
Facets=addFacet(Facets,v1,v6,v5);

%% plot facets
figure(1)
clf
plotFacets(Facets)
grid on

%% write stl file
writeSTL(Facets,filename)

%% Run Tetgen 
averageElementVolume=V/nElements;
system([MODEL_DIR '/scripts/tetgenSTL.sh ' filename ' ' num2str(averageElementVolume)]);

%% Create T and N files and clean tetgent output
system([MODEL_DIR '/scripts/tetgen2TN.sh ' filename ' ' num2str(meshID)]);

