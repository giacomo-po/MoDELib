function generateDipolarMicrostructure_cylinder(targetDensity)

%clc
%close all
%clear all

Burgers=0.2851e-9; % Burgers vector for Al [m]

R=0.5e-6/2/Burgers;   % cylinder radius [in units of b]
H=4*R;      % cylinder height [in units of b]
V=pi*R^2*H; % volume of cylinder
x0=0;       % offset of cylinder axis
y0=0;       % offset of cylinder axis


lMin=R/2;           % minimum size of loop [b]
lMax=R;             % max size of loop [b]
%targetDensity=1e15; % target dislocaiton density [m^-2]
fractionSessile=0.5;

% Matrix of lattice vectors
A=[0 1 1;
   1 0 1;
   1 1 0]/sqrt(2);
invA=inv(A);


% Plane normals for FCC (unsigned)
N =[-1    1   -1    1;
     1   -1   -1    1;
    -1   -1    1    1]/sqrt(3);

% Burgers vectors for FCC (unsigned)
B=[0  1  1;
    0  1 -1;
    1  0  1;
    1  0 -1;
    1  1  0;
    1 -1  0];

%% initialize variables
dislocationLength=0;
dislocationDensity=dislocationLength/V/Burgers^2; % current density [m^-2]

figure

k=1; % counter
rowB0=floor(rand(1)*size(B,1));  % a random row into B

while dislocationDensity<targetDensity
    s=sign(rand(1)-0.5);        % random sign of Burgers vector
    rowB=rem(rowB0,size(B,1))+1;
    b=s*B(rowB,:)/norm(B(rowB,:));    % the Burgers vector of the current loop
    
    r=rand(1)*R;        % random radius
    phi=rand(1)*2*pi;   % random angle
    z=(rand(1))*H;      % random height
    P1=[r*cos(phi) r*sin(phi) z]'; % the first point of the loop
    
    Nc=[];
    for n=1:size(N,2)       % loop over plane normals
        if abs(dot(b,N(:,n)))<0.001
            Nc=[Nc N(:,n)];
        end
    end
    if size(Nc,2)~=2
        error('something went wrong')
    end
    
    d1=cross(b,Nc(:,1))'; % "edge direction" on plane 1
    u=rand(1);
    a1=lMin*(1-u)+lMax*u;   % random size of loop on plane 1
    a1=round(a1/sqrt(3))*sqrt(3);
    
    d2=cross(b,Nc(:,2))'; % "edge direction" on plane 2
    u=rand(1);
    a2=lMin*(1-u)+lMax*u;   % random size of loop on plane 2
    a2=round(a2/sqrt(3))*sqrt(3);

    % Generate a sessile loop by overwriting d2
    if dislocationDensity<fractionSessile*targetDensity
        for kk=1:size(B,1)
            if abs(dot(B(kk,:),b))<1e-7
                d2=B(kk,:)'/norm(B(kk,:));
            end
        end
    a2=lMin*(1-u)+lMax*u;   % random size of loop on plane 2
    a2=round(a2);
    end
    
    if dot(P1,d1)<0     % correct d1 to point inside the cylinder
        d1=-d1;
    end
    if dot(P1,d2)<0     % correct d2 to point inside the cylinder
        d2=-d2;
    end
    P(:,1)=P1+[x0 y0 0]'; % shift the axis to match mesh file
    P=A*round(invA*P);

    
    P(:,2)=P(:,1)+a1*d1;    % second point of the loop
    P(:,3)=P(:,2)+a2*d2;    % third point of the loop
    P(:,4)=P(:,3)-a1*d1;    % fourth point of the loop
    P=A*round(invA*P);
    
    fR=0.95;
    if (   norm(P([1 2],1)-[x0 y0]')<fR*R && P(3,1)>0 && P(3,1)<H ...
            && norm(P([1 2],2)-[x0 y0]')<fR*R && P(3,2)>0 && P(3,2)<H ...
            && norm(P([1 2],3)-[x0 y0]')<fR*R && P(3,3)>0 && P(3,3)<H ...
            && norm(P([1 2],4)-[x0 y0]')<fR*R && P(3,4)>0 && P(3,4)<H) % all points are inside
        
        %rowB0
        rowB
        % plot loop
        hold on
        plot3([P(1,:) P(1,1)],[P(2,:) P(2,1)],[P(3,:) P(3,1)],'color',[0.5 0.5 0.5]+b/2,'Linewidth',3)
        axis equal
        set(gca,'View',[-34 26])
        
        % store loop for later
        dipolarLoops{k}.b=b;
        dipolarLoops{k}.P=P;
        dipolarLoops{k}.L=2*(a1+a2);
        
        % update dislocationLength and dislocationDensity
        dislocationLength=dislocationLength+dipolarLoops{k}.L;
        dislocationDensity=dislocationLength/V/Burgers^2;
        rowB0=rowB0+1;
        k=k+1;
    end
    
end


dislocationDensity

%return
%% write input files
file_V = fopen('V/V_0.txt','w');
file_E = fopen('E/E_0.txt','w');

nodeID=0;
nodeformat='%i %1.15e %1.15e %1.15e %i \n';
linkformat='%i %i %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %i \n';

for k=1:length(dipolarLoops)
    fprintf(file_V,nodeformat, [nodeID   dipolarLoops{k}.P(:,1)' 0]);
    fprintf(file_V,nodeformat, [nodeID+1 dipolarLoops{k}.P(:,2)' 0]);
    fprintf(file_V,nodeformat, [nodeID+2 dipolarLoops{k}.P(:,3)' 0]);
    fprintf(file_V,nodeformat, [nodeID+3 dipolarLoops{k}.P(:,4)' 0]);
    
    fprintf(file_E,linkformat, [nodeID   nodeID+1 dipolarLoops{k}.b zeros(1,3) 0]);
    fprintf(file_E,linkformat, [nodeID+1 nodeID+2 dipolarLoops{k}.b zeros(1,3) 0]);
    fprintf(file_E,linkformat, [nodeID+2 nodeID+3 dipolarLoops{k}.b zeros(1,3) 0]);
    fprintf(file_E,linkformat, [nodeID+3 nodeID   dipolarLoops{k}.b zeros(1,3) 0]);
    
    nodeID=nodeID+4;
end
fclose(file_V);
fclose(file_E);
