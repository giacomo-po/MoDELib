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
nElements=2e4;

%% Size and position of the cube
L1=2000; % the side length of the cube, in units of Burgers vector
L2=2000; % the side length of the cube, in units of Burgers vector
L3=4000; % the side length of the cube, in units of Burgers vector

%% Compute element size
V=L1*L2*L3;
averageElementVolume=V/nElements;

% coordinates of the 8 vertices of the cube.
% The base of the cube is at z=0. The cube is centered in x and y
v0=[0 0 0]-[L1 L2 0]/2;
v1=[L1 0 0]-[L1 L2 0]/2;
v2=[L1 L2 0]-[L1 L2 0]/2;
v3=[0 L2 0]-[L1 L2 0]/2;
v4=[0 0 L3]-[L1 L2 0]/2;
v5=[L1 0 L3]-[L1 L2 0]/2;
v6=[L1 L2 L3]-[L1 L2 0]/2;
v7=[0 L2 L3]-[L1 L2 0]/2;

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
system([MODEL_DIR '/scripts/tetgenSTL.sh ' filename ' ' num2str(averageElementVolume)]);

%% Create T and N files and clean tetgent output
system([MODEL_DIR '/scripts/tetgen2TN.sh ' filename ' ' num2str(meshID)]);

