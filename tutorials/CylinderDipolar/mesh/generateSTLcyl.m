%function generateSTLcyl
clc
close all
clear all
addpath('../../../../matlab/');


R=2000;
H=6*R;

np=30;

theta=[0:np-1]/np*2*pi;
x=R*cos(theta)+R;
y=R*sin(theta)+R;

Facets={};

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



plotFacets(Facets)
grid on
writeSTL(Facets,'cylinder_2k_12k')
