filename='S/S_1.txt';
S=load(filename);

%[Y1,Z1]=meshgrid([0:1:101],[0:1:101]);

%c=4;
%s=griddata(S(:,2),S(:,3),S(:,c),Y1,Z1);

figure(200)
quiver3(S(:,1),S(:,2),S(:,3),S(:,4),S(:,4),S(:,4))



figure(201)
plot(S(:,1),S(:,4:9),'.-')
%for c=1:6

%end