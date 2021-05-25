clear all
close all
clc


points=load("points.txt");

patches=load("edges.txt");

poly=load("poly.txt");
poly3D=load("poly3D.txt");
center=poly(1,:);
poly=poly(2:end,:);





figure(1)
clf
hold on
axis equal
plot(points(:,2),points(:,3),'o')
plot3(poly3D(:,1),poly3D(:,2),poly3D(:,3),'ko')
plot([poly(:,1);poly(1,1)],[poly(:,2);poly(1,2)],'x-r','Linewidth',2)

%text(points(:,2),points(:,3),num2str(points(:,1)))
%for n=1:size(poly,1)
%plot([center(1) poly(n,1)],[center(2) poly(n,2)],'r--')
%end
grid on

patchIDs=unique(patches(:,3));
patchSize=length(patchIDs)
for p=1:length(patchIDs)
    patchID=patchIDs(p);
    rows=find(patches(:,3)==patchID);
    edges=patches(rows,[1 2]);
    sourceRows=[];
    sinkRows=[];
    for r=1:size(edges,1)
        sourceRows=[sourceRows;find(points(:,1)==edges(r,1))];
        sinkRows=[sinkRows;find(points(:,1)==edges(r,2))];
    end
    patch(points(sourceRows,2),points(sourceRows,3),'b','FaceAlpha',0.5);
    drawnow
end
