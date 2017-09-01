function writeBicrystalPoly(meshID,useRegions,targetElements,filename,R,H,np,N,P0,A)

MODEL_DIR='../../../..';

P0=A*round(inv(A)*P0); % snap P0 to lattice
N=N/norm(N);
s=2*R*sin(pi/np);
V=np/4*H*s^2*cot(pi/np)

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
    pause(0.1)
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
    pause(0.1)
end

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
if ~useRegions
T=load(['../T/T_' num2str(meshID) '.txt']);
T(:,end)=T(:,end)*0;
dlmwrite(['../T/T_' num2str(meshID) '.txt'],T,'delimiter','\t')
end

%% Print C2G1 and C2G2 (paste in DDinput.txt)
format long

v1=cross([0 1 0]',N) % vector on GB plane in grain 1
v2=[v1(1) v1(2) -v1(3)]';

c=dot(v1,v2);   % cos(Theta)
s=sqrt(1-c^2);  % sin(Theta)

C2G1=eye(3)

C2G2=[c 0  s;
      0 1  0;
      -s 0  c]