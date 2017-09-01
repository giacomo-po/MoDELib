clear all
close all
clc

meshID=2;
N=load(['N/N_' num2str(meshID) '.txt']);
T=load(['T/T_' num2str(meshID) '.txt']);

stepID=0;
G=load(['G/G_' num2str(stepID) '.txt']);

figure(1)
clf
hold on
plot3(N(:,2),N(:,3),N(:,4),'k.')

%for k=1:size(G,1)
%id1=G(k,2)+1;
%id2=G(k,3)+1;
%    %plot3([N(id1,2) N(id2,2)],[N(id1,3) N(id2,3)],[[N(id1,4) N(id2,4)]],'b','Linewidth',1)
%    X=G(k,[4:6]);
%    plot3(X(1),X(2),X(3),'ro')
%end

return
%%
figure(1)
id=[1383 597]+1;
id1=id(1);
id2=id(2);
plot3([N(id1,2) N(id2,2)],[N(id1,3) N(id2,3)],[[N(id1,4) N(id2,4)]],'r','Linewidth',5)
%plot3(P0(:,5),P0(:,6),P0(:,7),'k-x','Linewidth',2)
%plot3(P1(:,5),P1(:,6),P1(:,7),'b-x','Linewidth',2)
%plot3(P2(:,5),P2(:,6),P2(:,7),'r-x','Linewidth',2)
plot3(P3(:,5),P3(:,6),P3(:,7),'b-x','Linewidth',2)
xlabel('X')
ylabel('Y')
zlabel('Z')