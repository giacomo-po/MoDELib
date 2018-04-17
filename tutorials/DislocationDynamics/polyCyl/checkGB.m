clear all
close all
clc

X=[7.767e+02 -6.955e+02  7.693e+02];


figure(1)
clf
hold on
N=221;
for f=1:N
%R=load(['rb' num2str(f) '.txt']);
%for r=1:size(R,1)
%plot3([R(r,1) R(r,4)],[R(r,2) R(r,5)],[R(r,3) R(r,6)],'bx-','Linewidth',2)
%end


%try
G=load(['gb' num2str(f) '.txt']);
patch(G(:,1),G(:,2),G(:,3),'g')
%end
for g=1:size(G,1)
gx=norm(G(g,:)-X);
    if gx<1
    gx
        G(g,:)
    
end
   % plot3([R(r,1) R(r,4)],[R(r,2) R(r,5)],[R(r,3) R(r,6)],'bx-','Linewidth',2)
end

end

axis equal


%%
plot3(X(:,1),X(:,2),X(:,3),'r.','Linewidth',2)
%plot3(X(1,1),X(1,2),X(1,3),'ks','Linewidth',2)

%for r=1:size(X,1)
%plot3([X(r,1) X(r,4)],[X(r,2) X(r,5)],[X(r,3) X(r,6)],'r','Linewidth',2)
%plot3([X(r,1)],[X(r,2) X(r,5)],[X(r,3) X(r,6)],'r','Linewidth',2)
%drawnow
%pause(0.1)

%end


