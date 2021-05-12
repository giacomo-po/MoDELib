clear all
close all
clc

boundaryPoints=[0 0;
     1 0
     0.75 0.5
     1 1;
     0 1]*100;
%boundaryPoints=ginput(10);
%boundaryPoints=rand(50,2);
save('boundaryPoints.txt', 'boundaryPoints', '-ascii', '-double', '-tabs');


internalPoints=[];
internalPoints=[
    50 50;
 0.2 0.2;
 50 50.1
];
internalPoints=rand(150,2)*50;
save('internalPoints.txt', 'internalPoints', '-ascii', '-double', '-tabs');

system('./tri')
plotMesh('vertices_mid.txt','triangles_mid.txt',boundaryPoints,internalPoints)
%plotMesh('vertices_out.txt','triangles_out.txt',boundaryPoints,internalPoints)

function plotMesh(vertexFile,trianglesFile,boundaryPoints,internalPoints)
vertices=load(vertexFile);
triangles=load(trianglesFile);

figure
clf
hold on
IDs=[1:size(boundaryPoints,1)]';
IDs=[IDs IDs+1];
IDs(end,2)=1;
quiver(boundaryPoints(IDs(:,1),1),boundaryPoints(IDs(:,1),2),boundaryPoints(IDs(:,2),1)-boundaryPoints(IDs(:,1),1),boundaryPoints(IDs(:,2),2)-boundaryPoints(IDs(:,1),2),0,'k','Linewidth',2)

if(~isempty(internalPoints))
plot(internalPoints(:,1),internalPoints(:,2),'mo','Linewidth',2)
end
plot(vertices(:,1),vertices(:,2),'xb')
for t=1:size(triangles,1)
    vIDs=triangles(t,:)+1;
    patch(vertices(vIDs,1),vertices(vIDs,2),'g','FaceALpha',0.2)
end
axis equal
grid on
end