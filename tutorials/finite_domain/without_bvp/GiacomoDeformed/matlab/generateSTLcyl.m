function generateSTLcyl
addpath('../../../../../matlab');


R=2000;
H=6*R;

np=40;

theta=[0:np-1]/np*2*pi;
x=R*cos(theta)+R;
y=R*sin(theta)+R;

%facetID=1;
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

writeSTL(Facets,'cylinder_2k_12k')

function Facets=addFacet(Facets,v1,v2,v3)
%facetID=;
facetID=length(Facets)+1;

Facets{facetID}.v1=v1;
Facets{facetID}.v2=v2;
Facets{facetID}.v3=v3;
n=cross(v1-v3,v2-v1);
n=n/norm(n)
Facets{facetID}.normal=n;

function plotFacets(Facets)
figure(1)
clf
hold on
axis equal
for k=1:length(Facets)
v1=Facets{k}.v1;
v2=Facets{k}.v2;
v3=Facets{k}.v3;
n=Facets{k}.normal;
plot3([v1(1) v2(1) v3(1)],[v1(2) v2(2) v3(2)],[v1(3) v2(3) v3(3)],'-o')
quiver3((v1(1)+v2(1)+v3(1))/3,(v1(2)+v2(2)+v3(2))/3,(v1(3)+v2(3)+v3(3))/3,n(1),n(2),n(3),1000)
end
