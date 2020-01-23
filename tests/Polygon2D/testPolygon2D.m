clear all
close all
clc

system('rm polyPoints.txt');
system('rm testPoints.txt');
system('rm results.txt');
system('rm nodes.txt');

%N0=4;
%poly0=rand(N0,2);
h=1.5;
polyPoints=[0 0;
       2 0;
       2 h;
       1 h;
       1 1;
       4 1;
       4 h;
       3 h;
       3 0;
       5 0;
       5 3;
       0 3];
   
%    polyPoints=[0 0;
%        3 0;
%        3 2;
%        2 2;
%        2 -1;
%        1 -1;
%        1 2;
%        0 2];
   
%    polyPoints=[0 0;
%        2 2;
%        2 0;
%        0 2];

patchPoints=[0.5 -1.2;
             5 0;
             5 3
             0.5 3];
         
%         patchPoints=[1.5 1.00 ; %//This is the original issue
%             3.5 1.00 ;
%             2.5 2.5];
         
 %  patchPoints=flipud(patchPoints);

   
% polyPoints=[0 0;
%        2 0;
%        2 2;
%        0 2;
%        0 0;
%        2 0;
%        2 2;
%        0 1];
  % polyPoints=flipud(polyPoints);
   
  
testPoints=[];
dlmwrite('polyPoints.txt',polyPoints,' ');
dlmwrite('testPoints.txt',testPoints,' ');
dlmwrite('patchPoints.txt',patchPoints,' ');

system('./test')
results=load('results.txt');
edges=load('edges.txt');
nodes=load('nodes.txt');

%check=max(abs(results(:,[1 2])-testPoints))

wnCol=3;

clrs={'b','g','r','m'};
clrsLab=[1 -1 2 -2];

%figure(1)
%plot(results(:,wnCol))

figure(2)
clf
hold on
axis equal
plot([polyPoints(:,1);polyPoints(1,1)],[polyPoints(:,2);polyPoints(1,2)],'b--','Linewidth',2);
patch([patchPoints(:,1)],[patchPoints(:,2)],'m','FaceAlpha',0.2,'Linewidth',2);
text(nodes(:,2),nodes(:,3),num2str(nodes(:,1)),'FontSize',16)

% zeroIDs=find(results(:,wnCol)==0);
% plot(testPoints(zeroIDs,1),testPoints(zeroIDs,2),'k.','Linewidth',2);
% for k=1:2
% plusIDs=find(results(:,wnCol)==k);
% minusIDs=find(results(:,wnCol)==-k);
% plot(testPoints(plusIDs,1),testPoints(plusIDs,2),[clrs{2*k-1} '.'],'Linewidth',2);
% plot(testPoints(minusIDs,1),testPoints(minusIDs,2),[clrs{2*k} '.'],'Linewidth',2);
% end

for k=1:size(edges,1)
    wn=edges(k,5);
    if wn==0
        clr='k';
    elseif wn==1
                clr='b';
                    elseif wn==-1
                clr='g';    
    elseif wn==2
                clr='r';
                    elseif wn==-2
                clr='m';
    end
quiver(edges(k,1),edges(k,2),edges(k,3)-edges(k,1),edges(k,4)-edges(k,2),0,'Linewidth',2,'Color',clr)
text(0.5*(edges(k,1)+edges(k,3)),0.5*(edges(k,2)+edges(k,4)),num2str(wn),'FontSize',16)
end
