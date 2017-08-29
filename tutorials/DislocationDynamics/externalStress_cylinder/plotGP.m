clear all
close all
clc

meshID=0;
N=load(['N/N_' num2str(meshID) '.txt']);
T=load(['T/T_' num2str(meshID) '.txt']);

stepID=0;
G=load(['G/G_' num2str(stepID) '.txt']);

figure(1)
clf
hold on
plot3(N(:,2),N(:,3),N(:,4),'k.')

for k=1:size(G,1)
id1=G(k,2)+1;
id2=G(k,3)+1;
    %plot3([N(id1,2) N(id2,2)],[N(id1,3) N(id2,3)],[[N(id1,4) N(id2,4)]],'b','Linewidth',1)
    X=G(k,[4:6]);
    plot3(X(1),X(2),X(3),'ro')
end

return
%%
figure(1)
id1=422+1;
id2=69+1;
plot3([N(id1,2) N(id2,2)],[N(id1,3) N(id2,3)],[[N(id1,4) N(id2,4)]],'m','Linewidth',2)