clc
close all
clear all

MODEL_DIR='../../../..';

%% Define output file name
meshID=1; % creates ../N/N_meshID.txt and ../T/T_meshID.txt
targetElements=1e4;

filename='verticalBicrystalCyl'; % this creates file bicrystalCyl.poly

R=1000;      % radius of cyl
H=4*R;       % height of cyl
A=pi*R^2;    % cross section area
V=A*H;       % volume
np=40;       % number of points along circumference
N=[0 3 1]';  % normal to GB plane
N=N/norm(N);

%% Define mesh vertices 
% bottom plane
for k=1:np
    theta=2*pi/np*(k-1);
    P(k,:)=[R*cos(theta) R*sin(theta) -H/2];
end

% top plane
for k=1:np
    theta=2*pi/np*(k-1);
    P(np+k,:)=[R*cos(theta) R*sin(theta) +H/2];
end


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

% Grain Boundary triangular facets
FB(1,:)=[np*3/4+1 np*7/4+1 np/4+1+np];
FB(2,:)=[np*5/4+1 np/4+1 np*3/4+1];


% inclined GB facets
for k=1:np
FI(k,:)=[k k+1 np+k];
end
%FI
% top facets
for k=1:np
FT(k,:)=[np+k k+1 k+np+1];
if k==np
FT(k,:)=[1 k+1 k];
end
end
%
%FT



FB(3,:)=[np/4 np/4+1 np*3/4+1];
FB(4,:)=[np/4+1 np/4+2 np*3/4+1];
FB(5,:)=[np*5/4 np*5/4+1 np*7/4+1];
FB(6,:)=[np*5/4+1 np*5/4+2 np*7/4+1];

F3=[FB;FI;FT];

for f=1:size(F3,1)
    vID=F3(f,:);
    fill3(P(vID,1),P(vID,2),P(vID,3),'g','Facealpha',0.2)
    drawnow
    pause(0.1)
end

% Top and bottom facets
%for k=1:np
%FL(k,:)=[k k+1 np-k np-k+1];
%    if k==np
%        FL(k,:)=[k 1 np np-1];
%    end
%end




% lower lateral facets
%for k=1:np/4
%	j=k*2;
%	FL(k,:)=[j-1 j np-j+1 np-j];
%end
% upper lateral facets
%for k=1:np/4
%	j=np+k*2;
%	FL(k,:)=[np+j-1 np+j-2 2*np-j-1 2*np-j];
%end

for(k=1:np/4-1)
n=k-1;
FL(k,:)=[np-n k k+1 np-k];
FL(n+np/4,:)=[np/4+2+n np/4+2+k np*3/4-n np*3/4+1-n];
end



%FL(5,:)=[7 8 15 16];
%FL(6,:)=[8 9 14 15];
%FL(7,:)=[9 10 13 14];
%FL(8,:)=[10 11 12 13];

for(i=1:np/4-1)
FL(np/2-2+i,:)=[2*np+1-i i+np 1+i+np 2*np-i];
FL(np*3/4-3+i,:)=[np+np/4+1+i np+np/4+2+i np+np*3/4+1-i np+np*3/4+2-i];
end

%FL(9,:)=[20+np 1+np 2+np 19+np];
%FL(10,:)=[19+np 2+np 3+np 18+np];
%FL(11,:)=[18+np 3+np 4+np 17+np];
%FL(12,:)=[17+np 4+np 5+np 16+np];
%FL(13,:)=[7+np 8+np 15+np 16+np];
%FL(14,:)=[8+np 9+np 14+np 15+np];
%FL(15,:)=[9+np 10+np 13+np 14+np];
%FL(16,:)=[10+np 11+np 12+np 13+np];



FL
% upper lateral facets

F4=[FL];

for f=1:size(F4,1)
    vID=F4(f,:);
    fill3(P(vID,1),P(vID,2),P(vID,3),'g','Facealpha',0.2)
    drawnow
    pause(0.1)
end

%% Define Regions
% Each row in the following matrix correponds to a region.
% Each row contains the vertices that bound that region.
%RP=[	[1:6 16:26 36:40]
%	[6:16  26:36];];
RP=[	[1:np/4+1 np*3/4+1:np*5/4+1 np*7/4+1:2*np]
	[np/4+1:np*3/4+1  np*5/4+1:np*7/4+1];];
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
fprintf(polyFile,'%i 0 \n',size(F3,1)+size(F4,1)); % number of facets
for r=1:size(F3,1)
fprintf(polyFile,'1 \n'); 
fprintf(polyFile,'3   %d %d %d \n',F3(r,:));
end
for r=1:size(F4,1)
fprintf(polyFile,'1 \n'); 
fprintf(polyFile,'4   %d %d %d %d \n',F4(r,:));
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

%% Print C2G1 and C2G2 (paste in DDinput.txt)
format long

v1=cross([1 0 0]',N) % vector on GB plane in grain 1
v2=[v1(1) v1(2) -v1(3)]';

c=dot(v1,v2);   % cos(Theta)
s=sqrt(1-c^2);  % sin(Theta)

C2G1=eye(3)

C2G2=[c 0  s;
      0 1  0;
      -s 0  c]

meshID
filename
