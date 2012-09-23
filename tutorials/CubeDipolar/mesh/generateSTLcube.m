clc

addpath('../../../matlab/');

L=4000; % the side length of the cube, in units of b

% coordinates of the 8 vertices of the cube
v0=[0 0 0];
v1=[L 0 0];
v2=[L L 0];
v3=[0 L 0];
v4=[0 0 L];
v5=[L 0 L];
v6=[L L L];
v7=[0 L L];

% Create 12 facets structure by splitting each face of the cube in 2
% triangles
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

% plot facets and make sure that normals point outward
figure(1)
clf
plotFacets(Facets)
grid on

% write stl file
writeSTL(Facets,['cube_' num2str(L)])
