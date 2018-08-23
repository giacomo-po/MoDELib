clear all
close all
clc

        P=[  0 0
            0 1
            %1 0.5
            0.5 -0.5
%                        1.8 0.8
            1 1
            1 0];
        
        P=ginput(5);
        
        figure(1)
        plot([P(:,1);P(1,1)],[P(:,2);P(1,2)],'k','Linewidth',2)
        fill(P(:,1),P(:,2),'r')
        axis equal

pgon2 = polyshape(P,'Simplify',false)
A2 = area(pgon2)