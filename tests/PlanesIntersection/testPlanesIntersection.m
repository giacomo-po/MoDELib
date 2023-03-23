clc
close all
clear all
addpath('../../matlab')

dim=3;
numPlanes=3;
tol=1e-7;
forceDegenerate=0;

%N=rand(dim,numPlanes);
%P=rand(dim,numPlanes);

N=[ 5.773502691896257e-01 -5.773502691896257e-01  5.773502691896257e-01
 5.773502691896257e-01  5.773502691896257e-01 -5.773502691896257e-01
-5.773502691896257e-01  5.773502691896257e-01  5.773502691896257e-01];

P=[ 2.045895620233077e+02 -3.323401871576773e+02  3.379970414071697e+02
 2.045895620233077e+02  3.323401871576773e+02 -3.379970414071697e+02
-2.045895620233077e+02  3.323401871576773e+02  3.379970414071697e+02];

%N(1:3,1:3)=[1 cos(120*pi/180)  cos(240*pi/180)
%   0 sin(120*pi/180) sin(240*pi/180)
%   0 0 0];

%P(1:3,1:3)=0;
%return

%P=[zeros(size(N)) rand(3,1)];

if forceDegenerate
N(:,end)=N(:,1);
P(:,end)=P(:,1);
end

%x=rand(dim,1);
x=[8.138799051457162e+02 8.053946237714779e+02 1.005505842847271e+03]';

inputFileName=[pwd '/build/inputFile.txt'];
outputFileName=[pwd '/build/outputFile.txt'];

fID=fopen(inputFileName,'w');
fprintf(fID,['dim=' num2str(dim) '; \n']);
fprintf(fID,['numPlanes=' num2str(numPlanes) '; \n']);
printMatrix(fID,'N',N);
printMatrix(fID,'P',P);
printMatrix(fID,'x',x);
fclose(fID)

%% 
system(['build/test ' inputFileName ' ' outputFileName])
outData=load(outputFileName)
success=outData(1);
snappedPoint=outData(2:end);


%return



%% Plot planes as circles
theta=[0:100]/100*2*pi;
figure(1)
hold on
for k=1:numPlanes
    plotCircle3(P(:,k),N(:,k),2000,theta,tol);
end


%% Sep up solution in MATLAB FOR COMPARISON
b=zeros(numPlanes,1);
for k=1:numPlanes
b(k)=dot(N(:,k),(x-P(:,k)));
end

[U,S,V] = reducedSVD(N,tol);
lam=V*diag(1./diag(S).^2)*V'*b;
x1=x-N*lam;

plot3(x(1),x(2),x(3),'sm') % point to snap
plot3(x1(1),x1(2),x1(3),'og') % point snapped (MATLAB)
plot3(snappedPoint(1),snappedPoint(2),snappedPoint(3),'xk') % point snapped (c++)


for k=1:numPlanes 
    plot3([x(1) x1(1)],[x(2) x1(2)],[x(3) x1(3)])
end

axis equal

function printMatrix(fID,name,m)
fprintf(fID,[name '=']);
for k=1:size(m,1)
fprintf(fID,'%1.15f',m(k,:));
if k<size(m,1)
fprintf(fID,'\n');
else
    fprintf(fID,';\n');
end
end
end

function [Ur,Sr,Vr]=reducedSVD(A,tol)
[U,S,V] = svd(A);
r=sum(abs(diag(S))>tol)
Ur=U(:,[1:r]);
Sr=S([1:r],[1:r]);
Vr=V(:,[1:r]);
end

function plotCircle3(P,n,r,theta,tol)
n=n/norm(n);
b=cross(n,rand(size(n)));

while norm(b) < tol
b=cross(n,rand(size(n)));
end
b=b/norm(b);
x=zeros(size(P,1),length(theta));
for(t=1:length(theta))
    x(:,t)=P+angleAxis(n,theta(t))*b*r;
end
    plot3(x(1,:),x(2,:),x(3,:))
end