% This file is part of MODEL, the Mechanics Of Defect Evolution Library.
%
% Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
% model is distributed without any warranty under the
% GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.

clc
close all
clear all
MODEL_DIR='../../../..';

addpath([MODEL_DIR '/matlab/']);
filename='cylinder'; % this creates file cylinder.stl
meshID=1;
targetElements=1e3;



R=4000; % radius of cylinder (units of Burgers vector)
H=3*R;    % height of cylinder (units of Burgers vector)
V=H*pi*R^2;
x0=0;     % offset of cylinder axis 
y0=0;     % offset of cylinder axis

averageElementVolume=V/targetElements;


np=30;    % number of circumferential points
theta=[0:np-1]/np*2*pi;
x=R*cos(theta)+x0;
y=R*sin(theta)+y0;

Facets={};

for k=1:np
    k1=rem(k,np)+1;
    
    % triangles on the bottom surface
    v1=[x(k)  y(k) 0];
    v2=[x0 y0 0];
    v3=[x(k1)  y(k1) 0];
    Facets=addFacet(Facets,v1,v2,v3);
    
    % triangles on the top surface
    v1=[x(k1)  y(k1) H];
    v2=[x0 y0 H];
    v3=[x(k)  y(k) H];
    Facets=addFacet(Facets,v1,v2,v3); 
    
    % triangles on the lateral surface
    v1=[x(k)  y(k) 0];
    v2=[x(k1) y(k1) 0];
    v3=[x(k)  y(k) H];
    Facets=addFacet(Facets,v1,v2,v3);

    v1=[x(k1)  y(k1) H];
    v2=[x(k)  y(k) H];
    v3=[x(k1)  y(k1) 0];
    Facets=addFacet(Facets,v1,v2,v3);
end

plotFacets(Facets)
grid on

%% write stl file
writeSTL(Facets,filename) % creates file mesh.stl

%% Run Tetgen 
system([MODEL_DIR '/scripts/tetgenSTL.sh ' filename ' ' num2str(averageElementVolume)]);

%% Create T and N files and clean tetgent output
system([MODEL_DIR '/scripts/tetgen2TN.sh ' filename ' ' num2str(meshID)]);


