clc
close all
clear all

MODEL_DIR='../../../..';
addpath([MODEL_DIR '/matlab'])

%% Define output file name
meshID=0; % creates ../N/N_1.txt and ../T/T_1.txt
targetElements=1e4;

filename='bicrystalCyl'; % this creates file cube.poly

R=1000;        % radius of cyl
H=4*R;         % height of cyl
V=pi*R^2*H;    % volume
np=20;         % number of points along circumference
A=[0 1 1;1 0 1;1 1 0]/sqrt(2); % matrix of primitive lattice vectors for FCC

STGBtype=4;

switch STGBtype
    case 1
        a=[0 1 0];     % tilt axis
        N=[  1 0 5]';    % normal to GB plane. STGB is sigma=13, theta=22.6199 deg
        N1=[-1 0 5]';
    case 2
        a=[0 1 0];     % tilt axis
        N=[  1 0 3]';    % normal to GB plane. STGB is sigma=5, theta=36.8699 deg
        N1=[-1 0 3]';
    case 3
        a=[0 1 0];     % tilt axis
        N=[  2 0 3]';     % normal to GB plane. STGB is sigma=13, theta=67.3801 deg
        N1=[-2 0 3]';     % normal to GB plane. STGB is sigma=13, theta=67.3801 deg
        
    case 4
        a=[1 -1 0]
        N=[1 1 1]';     % normal to GB plane. STGB is sigma=13, theta=67.3801 deg
        N1=[1 1 -1]';     % normal to GB plane. STGB is sigma=13, theta=67.3801 deg
        
    otherwise
        error('GB not implemented')
end
a=a/norm(a);
N=N/norm(N);
N1=N1/norm(N1);

if(dot(a,N)~=0)
    error('plane normal and tilt axis must be orthogonal')
end

r=reciprocalLatticeDirection(N,A);
d=1/norm(r);    % spacing of planes
h=0*d;          % offset of the GB plane, in units of N
P0=h*N;        % a point on the GB plane
%P0=[0 0 0]'; % a point on the GB plane
%P0=A*round(inv(A)*P0); % snap P0 to lattice


%return

%% Define mesh vertices
% bottom plane
for k=1:np
    theta=2*pi/np*(k-1);
    P(k,:)=[R*cos(theta) R*sin(theta) -H/2];
end
% GB plane
for k=1:np
    theta=2*pi/np*(k-1);
    P(np+k,:)=[R*cos(theta) R*sin(theta) (dot(P0,N)-R*cos(theta)*N(1)-R*sin(theta)*N(2))/N(3)];
    P1(k,:)=[R*cos(theta) R*sin(theta) (dot(P0,N1)-R*cos(theta)*N1(1)-R*sin(theta)*N1(2))/N1(3)];
end
% top plane
for k=1:np
    theta=2*pi/np*(k-1);
    P(2*np+k,:)=[R*cos(theta) R*sin(theta) +H/2];
end
P(3*np+1,:)=[0 0 -H/2];
P(3*np+2,:)=[0 0 dot(P0,N)/N(3)];
P(3*np+3,:)=[0 0 +H/2];

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

% bottom triangular facets
for k=1:np
    FB(k,:)=[k k+1 3*np+1];
    if k==np
        FB(k,:)=[k   1 3*np+1];
    end
end
% inclined GB facets
for k=1:np
    FI(k,:)=[k+np k+1+np 3*np+2];
    if k==np
        FI(k,:)=[k+np   1+np 3*np+2];
    end
end
% top facets
for k=1:np
    FT(k,:)=[k+2*np k+1+2*np 3*np+3];
    if k==np
        FT(k,:)=[k+2*np   1+2*np 3*np+3];
    end
end
F3=[FB;FI;FT];

for f=1:size(F3,1)
    vID=F3(f,:);
    fill3(P(vID,1),P(vID,2),P(vID,3),'g','Facealpha',0.2)
    drawnow
    pause(0.01)
end

% lower lateral facets
for k=1:np
    FL(k,:)=[k k+1 k+1+np k+np];
    if k==np
        FL(k,:)=[k 1 1+np 2*np];
    end
end
% upper lateral facets
for k=1:np
    FU(k,:)=[k+np k+1+np k+1+2*np k+2*np];
    if k==np
        FU(k,:)=[k+np np+1 1+2*np 3*np];
    end
end
F4=[FL;FU];

for f=1:size(F4,1)
    vID=F4(f,:);
    fill3(P(vID,1),P(vID,2),P(vID,3),'g','Facealpha',0.2)
    drawnow
    pause(0.01)
end

plot3(P1([1:end 1],1),P1([1:end 1],2),P1([1:end 1],3),'m','Linewidth',2)


%% Define Regions
% Each row in the following matrix correponds to a region.
% Each row contains the vertices that bound that region.
RP=[[1:2*np];
    [np+1:3*np]];

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
    fprintf(polyFile,'%i %.15e %.15e %.15e %i %.15e \n',[r Xm(r,:) r meshSize(r)]);
end

fclose(polyFile);

%% Run Tetgen
system([MODEL_DIR '/scripts/tetgenPOLY.sh ' filename]);

%% Create T and N files and clean tetgent output
system([MODEL_DIR '/scripts/tetgen2TN.sh ' filename ' ' num2str(meshID)]);

%% Print C2G1 and C2G2 (paste in DDinput.txt)
format long

%v1=cross([0 1 0]',N) % vector on GB plane in grain 1
%v2=[v1(1) v1(2) -v1(3)]';

f=R*0.75;
quiver3(P0(1),P0(2),P0(3),N(1)*f ,N(2)*f,N(3)*f,'Linewidth',4)
quiver3(P0(1),P0(2),P0(3),N1(1)*f,N1(2)*f,N1(3)*f,'Linewidth',4)
quiver3(P0(1),P0(2),P0(3),a(1)*f,a(2)*f,a(3)*f,'Linewidth',4)

c=dot(N,N1);   % cos(Theta)
%c=dot(v1,v2);   % cos(Theta)
s=sqrt(1-c^2);  % sin(Theta)
theta=atan2(s,c);
theta_deg=atan2(s,c)*180/pi
C2G1=eye(3);
C2G1=integerC2G(C2G1)
C2G2=angleAxis(a',theta)
%C2G2=[c 0  s;
%    0 1  0;
%    -s 0  c];
C2G2=integerC2G(C2G2)
