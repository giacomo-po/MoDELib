close all
clear all
clc

SIGMA=[3 5]

for s=1:length(SIGMA)
    sigma=SIGMA(s);
    
    fileName=['SimpleCubic_sigma_' num2str(sigma)];
M=load([fileName '.txt']);
A=M(1:3,:);
B=M(4:6,:);
CSL=M(7:9,:);

figure(2*s)
clf
hold on
np=2
for i=[-np:np]
    for j=[-np:np]
        for k=[-np:np]
            P0=i*A(:,1)+j*A(:,2)+k*A(:,3);
            P1=i*B(:,1)+j*B(:,2)+k*B(:,3);
            P2=i*CSL(:,1)+j*CSL(:,2)+k*CSL(:,3);
            
            plot3(P0(1),P0(2),P0(3),'bo','Linewidth',1)
            plot3(P1(1),P1(2),P1(3),'rx','Linewidth',2)
            plot3(P2(1),P2(2),P2(3),'g.','Linewidth',2)
        end
    end
end
title(fileName)
axis equal

end

