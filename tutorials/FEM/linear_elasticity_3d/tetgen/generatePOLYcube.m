clc
close all
clear all

filename='cube';
modelDir='~/Documents/model/';

L=1001;     % side length of the cube
L1=900;  % side length of flat punch

%% Define vertices
% Coordinates of the outer cube
% The base of the cube is at z=0. The cube is centered in x and y
P(1,:)=[0 0 0]*L-[1 1 0]*L/2;
P(2,:)=[1 0 0]*L-[1 1 0]*L/2;
P(3,:)=[1 1 0]*L-[1 1 0]*L/2;
P(4,:)=[0 1 0]*L-[1 1 0]*L/2;
P(5,:)=[0 0 1]*L-[1 1 0]*L/2;
P(6,:)=[1 0 1]*L-[1 1 0]*L/2;
P(7,:)=[1 1 1]*L-[1 1 0]*L/2;
P(8,:)=[0 1 1]*L-[1 1 0]*L/2;

% Coordinates of the inner cube
P( 9,:)=[0   0 0]-[L1 L1 0]/2;
P(10,:)=[L1  0 0]-[L1 L1 0]/2;
P(11,:)=[L1 L1 0]-[L1 L1 0]/2;
P(12,:)=[0  L1 0]-[L1 L1 0]/2;
P(13,:)=[0   0 L]-[L1 L1 0]/2;
P(14,:)=[L1  0 L]-[L1 L1 0]/2;
P(15,:)=[L1 L1 L]-[L1 L1 0]/2;
P(16,:)=[0  L1 L]-[L1 L1 0]/2;


figure(1)
clf
plot3(P(:,1),P(:,2),P(:,3),'ro','Linewidth',2)
hold on
text(P(:,1),P(:,2),P(:,3),num2str([1:size(P,1)]'),'FontSize',16)
axis equal
xlabel('x')
ylabel('y')
grid on


%% Define facets
% Facets normal to x
Fx=[ 1  5  8  4;
     2  6  7  3;
     9 13 16 12;
    10 14 15 11];
% Facets normal to y
Fy=[ 1  2  6  5;
     4  3  7  8;
     9 10 14 13;
    12 11 15 16];
% Facets normal to z
Fz=[1 2 10 9;
    2 3 11 10;
    3 4 12 11;
    4 1 9 12;
    9 10 11 12];
Fz=[Fz;Fz+4];

 % Inclined facets
Fxy=[1 9 13 5;
    2 10 14 6;
    3 11 15 7;
    4 12 16 8];

F=[Fx;Fy;Fz;Fxy];
for f=1:size(F,1)
    vID=F(f,:);
    fill3(P(vID,1),P(vID,2),P(vID,3),'g','Facealpha',0.2)
    drawnow
    pause(0.1)
end

%% Define Regions
% Each row in the following matrix correponds to a region.
% Each row contains the vertices that bound that region.
RP=[9 10 11 12 13 14 15 16;
    1 2 9 10 5 6 13 14;
    2 10 11 3 6 14 15 7;
    4 12 11 3 8 16 15 7;
    1 4 12 9 8 5 16 13];

% The rows of the matrix Xm are the barycenters of each region
for r=1:size(RP,1)
Xm(r,:)=mean(P(RP(r,:),:));
end
plot3(Xm(:,1),Xm(:,2),Xm(:,3),'bx','Linewidth',2)
text(Xm(:,1),Xm(:,2),Xm(:,3),num2str([1:size(Xm,1)]'),'FontSize',16)

%% Write .poly file
polyFile = fopen('cube.poly','w');

% Part 1- the node list.
pointFormat='%i %d %d %d \n';
fprintf(polyFile,'# Part 1 - the node list.\n');
fprintf(polyFile,'%i 3 0 0 \n',size(P,1));  % number of nodes
for k=1:size(P,1)
fprintf(polyFile,pointFormat,[k P(k,:)]);
end
fprintf(polyFile,'\n');
fprintf(polyFile,'\n');

% Part 2- the facet list.
fprintf(polyFile,'# Part 2 - the facet list.\n');
fprintf(polyFile,'%i 0 \n',size(F,1)); % number of facets
for r=1:size(F,1)
fprintf(polyFile,'1 \n'); 
fprintf(polyFile,'4   %d %d %d %d \n',F(r,:));
end

% Part 3- the hole list.
fprintf(polyFile,'# Part 3 - the hole list.\n');
fprintf(polyFile,'0 \n\n');

% Part 4- the region list.
fprintf(polyFile,'# Part 4 - the region list.\n');
fprintf(polyFile,'%i \n',size(Xm,1));

vL=((L-L1)/2)^3/6/sqrt(2);
meshSize=ones(size(Xm,1),1)*vL;
meshSize(1)=vL*10;
for r=1:size(Xm,1)
fprintf(polyFile,'%i %d %d %d %d \n',[r Xm(r,:) meshSize(r)]);
end

fclose(polyFile);

%% Run Tetgen 
system([modelDir 'scripts/tetgenPOLY.sh ' filename]);

%% Create T and N files and clean tetgent output
meshID=1;
system([modelDir 'scripts/tetgen2TN.sh ' filename ' ' num2str(meshID)]);
