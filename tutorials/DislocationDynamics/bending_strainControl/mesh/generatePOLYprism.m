clc
close all
clear all

MODEL_DIR='../../../..';

%% Define output file name
meshID=0; % creates ../N/N_1.txt and ../T/T_1.txt
filename='prism'; % this creates file cube.poly
nElements=2e4;

L1=2000; % the side length of the cube, in units of Burgers vector
L2=6000; % the side length of the cube, in units of Burgers vector
L3=2000; % the side length of the cube, in units of Burgers vector

V=L1*L2*L3;
averageElementVolume=V/nElements;


%% Define vertices
% Coordinates of the outer cube
% The base of the cube is at z=0. The cube is centered in x and y
P(1,:)=[0 0 0]-[L1 L2 L3]/2;
P(2,:)=[L1 0 0]-[L1 L2 L3]/2;
P(3,:)=[L1 L2 0]-[L1 L2 L3]/2;
P(4,:)=[0 L2 0]-[L1 L2 L3]/2;
P(5,:)=[0 0 L3]-[L1 L2 L3]/2;
P(6,:)=[L1 0 L3]-[L1 L2 L3]/2;
P(7,:)=[L1 L2 L3]-[L1 L2 L3]/2;
P(8,:)=[0 L2 L3]-[L1 L2 L3]/2;
P(9,:)=[0 0 L3/2]-[L1 L2 L3]/2;
P(10,:)=[L1 0 L3/2]-[L1 L2 L3]/2;
P(11,:)=[0 L2 L3/2]-[L1 L2 L3]/2;
P(12,:)=[L1 L2 L3/2]-[L1 L2 L3]/2;

figure(1)
clf
plot3(P(:,1),P(:,2),P(:,3),'ro','Linewidth',2)
hold on
text(P(:,1),P(:,2),P(:,3),num2str([1:size(P,1)]'),'FontSize',16)
axis equal
xlabel('x')
ylabel('y')
grid on

%return

%% Define facets
% Facets normal to x
Fx=[ 1  2 10 9;
     9 10 6 5;
     4 3 12 11;
    11 12 7 8];
% Facets normal to y
Fy=[ 1  9 11 4;
     11 9 5 8;
    2 10 12 3;
     10 6 7 12];
% Facets normal to z
Fz=[ 1  2  3  4;
     9 10 12 11;
    5 6 7 8];

F=[Fx;Fy;Fz];
for f=1:size(F,1)
    vID=F(f,:);
    fill3(P(vID,1),P(vID,2),P(vID,3),'g','Facealpha',0.2)
    drawnow
    pause(0.1)
end

%return

%% Define Regions
% Each row in the following matrix correponds to a region.
% Each row contains the vertices that bound that region.
RP=[1 2 3 4 9 10 12 11;
    1 10 12 11 5 6 7 8];

% The rows of the matrix Xm are the barycenters of each region
for r=1:size(RP,1)
Xm(r,:)=mean(P(RP(r,:),:));
end
plot3(Xm(:,1),Xm(:,2),Xm(:,3),'bx','Linewidth',2)
text(Xm(:,1),Xm(:,2),Xm(:,3),num2str([1:size(Xm,1)]'),'FontSize',16)

%% Write .poly file
polyFile = fopen([filename '.poly'],'w');

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

%vL=t^3/6/sqrt(2);
meshSize=ones(size(Xm,1),1)*averageElementVolume;
for r=1:size(Xm,1)
%fprintf(polyFile,'%i %d %d %d %d \n',[r Xm(r,:) meshSize(r)]);
fprintf(polyFile,'%i %d %d %d %i %d \n',[r Xm(r,:) r meshSize(r)]);
end

fclose(polyFile);

%% Run Tetgen 
system([MODEL_DIR '/scripts/tetgenPOLY.sh ' filename]);

%% Create T and N files and clean tetgent output
system([MODEL_DIR '/scripts/tetgen2TN.sh ' filename ' ' num2str(meshID)]);
