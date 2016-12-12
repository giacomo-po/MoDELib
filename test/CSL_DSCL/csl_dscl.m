A=eye(3);

R=[12 -3 4;4 12 -3;-3 4 12]/13;
det(R)
R'*R
A1=R*A;
C=[-1  3  1;
   -1  0 -3;
   -1 -1  0];
figure(1)
clf
hold on
np=2
for i=[-np:np]
    for j=[-np:np]
        for k=[-np:np]
            %P=i*A(:,1)+j*A(:,2)+k*A(:,3);
            %P1=i*A1(:,1)+j*A1(:,2)+k*A1(:,3);
            PC=i*C(:,1)+j*C(:,2)+k*C(:,3);

            %plot3(P(1),P(2),P(3),'bo')
            %plot3(P1(1),P1(2),P1(3),'rx')
            plot3(PC(1),PC(2),PC(3),'g.','Linewidth',2)
        end
    end
end

np=2
for i=[-np:np]
    for j=[-np:np]
        for k=[-np:np]
            P=i*A(:,1)+j*A(:,2)+k*A(:,3);
            P1=i*A1(:,1)+j*A1(:,2)+k*A1(:,3);
            %PC=i*C(:,1)+j*C(:,2)+k*C(:,3);

            plot3(P(1),P(2),P(3),'bo','Linewidth',1)
            plot3(P1(1),P1(2),P1(3),'rx','Linewidth',1)
            %plot3(PC(1),PC(2),PC(3),'gs')
        end
    end
end


axis equal