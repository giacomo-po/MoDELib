clear all
close all
clc


system(['rm edges_*.txt'])
system(['rm lastPatch_*.txt'])
system(['rm points_*.txt'])
system(['rm poly.txt'])
system(['./microstructureGenerator'])

poly=load("poly.txt");

%%
for stepID=0:4


points=load(['points_' num2str(stepID) '.txt']);

patches=load(['edges_' num2str(stepID) '.txt']);

lastPatch=load(['lastPatch_' num2str(stepID) '.txt']);


%poly3D=load("poly3D.txt");
%center=poly(1,:);
%poly=poly(2:end,:);





figure(stepID+1)
clf
hold on
axis equal
text(poly(:,1),poly(:,2),num2str([0:size(poly,1)-1]'))
if ~isempty(points)
plot(points(:,1),points(:,2),'ms','Linewidth',2)
end
%plot3(poly3D(:,1),poly3D(:,2),poly3D(:,3),'ko')
plot([poly(:,1);poly(1,1)],[poly(:,2);poly(1,2)],'x-r','Linewidth',2)

%text(points(:,2),points(:,3),num2str(points(:,1)))
%for n=1:size(poly,1)
%plot([center(1) poly(n,1)],[center(2) poly(n,2)],'r--')
%end
if(~isempty(patches))
patchIDs=unique(patches(:,1));
%patchSize=length(patchIDs)
for p=1:length(patchIDs)
    patchID=patchIDs(p);
    rows=find(patches(:,1)==patchID);
    edges=patches(rows,[2:5]);
    %sourceRows=[];
    %sinkRows=[];
    %for r=1:size(edges,1)
    %    sourceRows=[sourceRows;find(points(:,1)==edges(r,1))];
    %    sinkRows=[sinkRows;find(points(:,1)==edges(r,2))];
    %end
    patch(edges(:,1),edges(:,2),'b','FaceAlpha',0.5);
    drawnow
end
end

patch(lastPatch(:,2),lastPatch(:,3),'g','FaceAlpha',0.5);


grid on
end
% 

