filename='S/S_1.txt';
S=load(filename);

[Y1,Z1]=meshgrid([0:1:101],[0:1:101]);

c=4;
s=griddata(S(:,2),S(:,3),S(:,c),Y1,Z1);

figure(200)
surf(Y1,Z1,s)


figure(201)
contour(Y1,Z1,s,20)
axis equal


%for c=1:6

%end