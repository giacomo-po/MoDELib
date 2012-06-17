function writeSTL(Facets,filename)

facetformat=' facet normal  %1.15e %1.15e %1.15e \n';
vertexformat='   vertex  %1.15e %1.15e %1.15e \n';


STLfile = fopen([filename '.stl'],'w');

fprintf(STLfile,['solid ' filename '\n']);




for k=1:length(Facets)
    fprintf(STLfile,facetformat,Facets{k}.normal);
     fprintf(STLfile,'  outer loop \n');
     fprintf(STLfile,vertexformat,Facets{k}.v1);
     fprintf(STLfile,vertexformat,Facets{k}.v2);
     fprintf(STLfile,vertexformat,Facets{k}.v3);
     fprintf(STLfile,'  endloop \n');
 fprintf(STLfile,' endfacet \n');

end

fprintf(STLfile,'endsolid');


fclose(STLfile);