clc
close all
clear all

MODEL_DIR='../../../..';

%% Define output file name
meshID=3; % creates ../N/N_1.txt and ../T/T_1.txt
targetElements=1e4;

filename='bicrystalCyl'; % this creates file cube.poly
t=400;      	% Film thickness (of columnar grains)
H=200;      	% Height of one of the crystals
X=H/cos(degtorad(30));


V=X*1.5*H*t;		%Volume of the triple crystal



t=t/2;

P(1,:)=[0,-t,0];
P(2,:)=[X/2,-t,-H];
P(3,:)=[X/2,-t,0];
P(4,:)=[X/2,-t,H];
P(5,:)=[-X,-t,H];
P(6,:)=[-X,-t,-H];

P(7,:)=[0,t,0];
P(8,:)=[X/2,t,-H];
P(9,:)=[X/2,t,0];
P(10,:)=[X/2,t,H];
P(11,:)=[-X,t,H];
P(12,:)=[-X,t,-H];
P(13,:)=[-X,-t,0];
P(14,:)=[-X,t,0];


figure(1)
clf
plot3(P(:,1),P(:,2),P(:,3),'ro','Linewidth',2)
hold on
text(P(:,1),P(:,2),P(:,3),num2str([1:size(P,1)]'),'FontSize',16)
axis equal
xlabel('x')
ylabel('y')
grid on

FT(1,:)=[1 2 3];
FT(2,:)=[3 4 1];
FT(3,:)=[1 4 5];
FT(4,:)=[1 5 13];
FT(5,:)=[1 6 2];


FT(1+5,:)=[1+6 2+6 3+6];
FT(2+5,:)=[3+6 4+6 1+6];
FT(3+5,:)=[1+6 4+6 5+6];
FT(4+5,:)=[1+6 5+6 14];
FT(5+5,:)=[1+6 6+6 2+6];


%%% 			BOTTOM LEFT 					internal facet (grain buondary)
FT(11,:)=[12 6 1];
FT(12,:)=[12 1 7];
%%% 				TOP LEFT 						internal facet (grain buondary)
FT(13,:)=[7 1 5];
FT(14,:)=[7 5 11];
%%% 			MIDDLE RIGHT 					internal facet (grain buondary)
FT(15,:)=[7 1 3];
FT(16,:)=[3 7 9];


%%% 			BOTTOM SQUARE
FT(17,:)=[12 6 2];
FT(18,:)=[12 2 8];
%%% 			RIGHT SQUARE
FT(19,:)=[2 8 3];
FT(20,:)=[8 3 9];
FT(21,:)=[9 3 4];
FT(22,:)=[9 4 10];
%%% 			TOP SQUARE
FT(23,:)=[10 4 5];
FT(24,:)=[5 10 11];
%%% 			LEFT SQUARE
FT(25,:)=[11 14 5];
FT(26,:)=[14 5 13];
FT(27,:)=[13 14 12];
FT(28,:)=[12 13 6];
FT(29,:)=[1 13 6];
FT(30,:)=[7 14 12];


for f=1:size(FT,1)
    vID=FT(f,:);
    fill3(P(vID,1),P(vID,2),P(vID,3),'g','Facealpha',0.2)
    drawnow
    %pause(0.1)
end

%% Define Regions
% Each row in the following matrix correponds to a region.
% Each row contains the vertices that bound that region.
RP=[	[1 5 13 6 7   11  14  12]
			[6 2 3  1 6+6 2+6 3+6 1+6]
			[1 3 4  5 1+6 3+6 4+6 5+6];];

size(RP);
% The rows of the matrix Xm are the barycenters of each region
for r=1:size(RP,1)
Xm(r,:)=mean(P(RP(r,:),:));
end
plot3(Xm(:,1),Xm(:,2),Xm(:,3),'bx','Linewidth',2)
text(Xm(:,1),Xm(:,2),Xm(:,3),num2str([1:size(Xm,1)]'),'FontSize',16)

%% Write .poly file
polyFile = fopen([filename '.poly'],'w');

% Part 1- the node list.
pointFormat='%i %.15e %.15e %.15e \n';
fprintf(polyFile,'# Part 1 - the node list.\n');
fprintf(polyFile,'%i 3 0 0 \n',size(P,1));  % number of nodes
for k=1:size(P,1)
fprintf(polyFile,pointFormat,[k P(k,:)]);
end
fprintf(polyFile,'\n');
fprintf(polyFile,'\n');

% Part 2- the facet list.
fprintf(polyFile,'# Part 2 - the facet list.\n');
fprintf(polyFile,'%i 0 \n',size(FT,1)); % number of facets
for r=1:size(FT,1)
fprintf(polyFile,'1 \n'); 
fprintf(polyFile,'3   %d %d %d \n',FT(r,:));
end

% Part 3- the hole list.
fprintf(polyFile,'# Part 3 - the hole list.\n');
fprintf(polyFile,'0 \n\n');

% Part 4- the region list.
fprintf(polyFile,'# Part 4 - the region list.\n');
fprintf(polyFile,'%i \n',size(Xm,1));

meshSize=ones(size(Xm,1),1)*V/targetElements;
for r=1:size(Xm,1)
fprintf(polyFile,'%i %.15e %.15e %.15e %i %.15e \n',[r Xm(r,:) r meshSize(2)]);
end

fclose(polyFile);

%% Run Tetgen 
system([MODEL_DIR '/scripts/tetgenPOLY.sh ' filename]);

%% Create T and N files and clean tetgent output
system([MODEL_DIR '/scripts/tetgen2TN.sh ' filename ' ' num2str(meshID)]);


%FT
