%function generateSTLcyl
clc
close all
clear all
addpath('../../../matlab/');


L=2000;


v0=[0 0 0];
v1=[L 0 0];
v2=[L L 0];
v3=[0 L 0];
v4=[0 0 L];
v5=[L 0 L];
v6=[L L L];
v7=[0 L L];


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


plotFacets(Facets)
grid on
writeSTL(Facets,['cube_' num2str(L)])



return

H=6*R;

np=30;

theta=[0:np-1]/np*2*pi;
x=R*cos(theta)+R;
y=R*sin(theta)+R;

Facets={};
    Facets=addFacet(Facets,v1,v2,v3);

for k=1:np
    k1=rem(k,np)+1;
    
    v1=[x(k)  y(k) 0];
    v2=[x(k1) y(k1) 0];
    v3=[x(k)  y(k) H];
    Facets=addFacet(Facets,v1,v2,v3);
        
    v1=[x(k)  y(k) 0];
    v2=[R R 0];
    v3=[x(k1)  y(k1) 0];
    Facets=addFacet(Facets,v1,v2,v3);
    
    v1=[x(k1)  y(k1) H];
    v2=[x(k)  y(k) H];
    v3=[x(k1)  y(k1) 0];
    Facets=addFacet(Facets,v1,v2,v3);
    
    v1=[x(k1)  y(k1) H];
    v2=[R R H];
    v3=[x(k)  y(k) H];
    Facets=addFacet(Facets,v1,v2,v3);
    
    
end



