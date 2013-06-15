function Facets=addFacet(Facets,v1,v2,v3)

facetID=length(Facets)+1;

Facets{facetID}.v1=v1;
Facets{facetID}.v2=v2;
Facets{facetID}.v3=v3;
n=cross(v1-v3,v2-v1);
n=n/norm(n);
Facets{facetID}.normal=n;
facetID=length(Facets)+1;
