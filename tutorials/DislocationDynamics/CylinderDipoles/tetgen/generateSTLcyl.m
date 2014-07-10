% This file is part of MODEL, the Mechanics Of Defect Evolution Library.
%
% Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
% model is distributed without any warranty under the
% GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.

clc
close all
clear all
addpath('../../../../matlab/');


R=2127/2; % radius of cylinder (units of Burgers vector)
H=6*R;    % height of cylinder (units of Burgers vector)

np=30;    % number of circumferential points

theta=[0:np-1]/np*2*pi;
x=R*cos(theta)+R;
y=R*sin(theta)+R;

Facets={};

for k=1:np
    k1=rem(k,np)+1;
    
    % triangles on the bottom surface
    v1=[x(k)  y(k) 0];
    v2=[R R 0];
    v3=[x(k1)  y(k1) 0];
    Facets=addFacet(Facets,v1,v2,v3);
    
    % triangles on the top surface
    v1=[x(k1)  y(k1) H];
    v2=[R R H];
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
writeSTL(Facets,'mesh') % creates file mesh.stl
