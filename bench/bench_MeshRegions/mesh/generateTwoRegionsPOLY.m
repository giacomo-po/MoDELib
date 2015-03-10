clc
close all
clear all

modelDir='../../../';


%% Define output file name
meshID=2;
filename='bicrystal'; % this creates file cube.poly

L=1000;     % side length of the cube
Lp=300;     % side length of flat punch
t=50;       % tickness of the boundary layer
%L1=L-2*t;   % side length of the inner cube

%% Define vertices
% Coordinates of the outer cube
% The base of the cube is at z=0. The cube is centered in x and y
P(1,:)=[0 0 0]*L;%-[1 1 0]*L/2;
P(2,:)=[1 0 0]*L;%-[1 1 0]*L/2;
P(3,:)=[1 1 0]*L;%-[1 1 0]*L/2;
P(4,:)=[0 1 0]*L;%-[1 1 0]*L/2;
P(5,:)=[0 0 1]*L;%-[1 1 0]*L/2;
P(6,:)=[1 0 1]*L;%-[1 1 0]*L/2;
P(7,:)=[1 1 1]*L;%-[1 1 0]*L/2;
P(8,:)=[0 1 1]*L;%-[1 1 0]*L/2;

% Coordinates of mid-plane (GB interface)
P(9,:)=[0 0 0.5]*L;%-[1 1 0]*L/2;
P(10,:)=[1 0 0.5]*L;%-[1 1 0]*L/2;
P(11,:)=[1 1 0.5]*L;%-[1 1 0]*L/2;
P(12,:)=[0 1 0.5]*L;%-[1 1 0]*L/2;

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
Fx=[ 9  5  8  12;
     1  9  12 4;
     10 6  7  11;
     2  10 11 3];
     
% Facets normal to y
Fy=[ 1  2  10 9;
     4  3  11 12;
     9  10 6  5;
     12 11 7  8];
% Facets normal to z
Fz=[ 1  2  3  4;
     5  6  7  8;
     9  10 11 12];


F=[Fx;Fy;Fz];
for f=1:size(F,1)
    vID=F(f,:);
    if P(vID,3)<=L/2
       fill3(P(vID,1),P(vID,2),P(vID,3),'r','Facealpha',0.2);
       drawnow
       pause(0.1)
    else
        fill3(P(vID,1),P(vID,2),P(vID,3),'b','Facealpha',0.2);
        drawnow
        pause(0.1)
    end
end

%% Define Regions
% Each row in the following matrix correponds to a region.
% Each row contains the vertices that bound that region.
RP=[9  10  6   5  12  11  7   8;
    1  2   10  9  4   3   11  12];

%---------- DO NOT CHANGE BELOW:  -----------------------
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

vL=1e+5
%vL=t^3/6/sqrt(2);
meshSize=ones(size(Xm,1),1)*vL;
meshSize(1)=vL; % increase volume element in central cube
for r=1:size(Xm,1)
fprintf(polyFile,'%i %d %d %d %i %d \n',[r Xm(r,:) r meshSize(r)]);
end

fclose(polyFile);

%% Run Tetgen 
system([modelDir 'scripts/tetgenPOLY.sh ' filename]);

%% Create T and N files and clean tetgent output
system([modelDir 'scripts/tetgen2TN.sh ' filename ' ' num2str(meshID)]);
