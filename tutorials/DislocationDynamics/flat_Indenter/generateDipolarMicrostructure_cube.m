clc
close all
clear all

Burgers=0.2556e-9; % for Cu [m]
targetDensity=2e14; % target initial dislocaiton density [m^-2]

L=1000;
V=L^3;

lMin=L/6;
lMax=L/3;

%alpha=0;
%c=cos(alpha);
%s=sin(alpha);
%Rot=[c -s 0;s c 0;0 0 1];

N =[-1    1   -1    1;
     1   -1   -1    1;
    -1   -1    1    1];
%N=Rot*N;
    

B=[0  1  1;
   0  1 -1;
   1  0  1;
   1  0 -1;
   1  1  0;
   1 -1  0];

%B=(Rot*B')';


totalLength=0;
density=totalLength/V/Burgers^2;
k=1;

figure(1)
hold on


while density<targetDensity
r=ceil(rand(1)*size(B,1));  % a random row into B
b=B(r,:)/norm(B(r,:)); % normalize b
bs=round(rand(1))*2-1; % random sign of b
b=b*bs;
Nc=[];
for n=1:size(N,2)
if abs(dot(b,N(:,n)))<0.001
Nc=[Nc N(:,n)/norm(N(:,n))];
end
end
if size(Nc,2)~=2
error('something went wrong')
end
d1=cross(b,Nc(:,1))';
d2=cross(b,Nc(:,2))';
P1=rand(3,1)*L;
if dot(P1,d1)<0
    d1=-d1;
end
if dot(P1,d2)<0
    d2=-d2;
end
P(:,1)=P1;
u=rand(1);
a1=lMin*(1-u)+lMax*u;
u=rand(1);
a2=lMin*(1-u)+lMax*u;
P(:,2)=P(:,1)+a1*d1;
P(:,3)=P(:,2)+a2*d2;
P(:,4)=P(:,3)-a1*d1;

if ( length(find(P<0))==0 && length(find(P>L))==0)
figure(1)
hold on
plot3(P(1,:),P(2,:),P(3,:),'color',[0.5 0.5 0.5]+b/2)

dipolarLoops{k}.b=b;
dipolarLoops{k}.P=P;
dipolarLoops{k}.L=2*(a1+a2);
totalLength=totalLength+dipolarLoops{k}.L;
density=totalLength/V/Burgers^2;
k=k+1;
end

end
density



file_V = fopen('V/V_0.txt','w');
file_E = fopen('E/E_0.txt','w');

nodeID=0;
nodeformat='%i %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %i \n';
linkformat='%i %i %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %i \n';

for k=1:length(dipolarLoops)
fprintf(file_V,nodeformat, [nodeID   dipolarLoops{k}.P(:,1)' zeros(1,3) 0]);
fprintf(file_V,nodeformat, [nodeID+1 dipolarLoops{k}.P(:,2)' zeros(1,3) 0]);
fprintf(file_V,nodeformat, [nodeID+2 dipolarLoops{k}.P(:,3)' zeros(1,3) 0]);
fprintf(file_V,nodeformat, [nodeID+3 dipolarLoops{k}.P(:,4)' zeros(1,3) 0]);

fprintf(file_E,linkformat, [nodeID   nodeID+1 dipolarLoops{k}.b zeros(1,3) 0]);
fprintf(file_E,linkformat, [nodeID+1 nodeID+2 dipolarLoops{k}.b zeros(1,3) 0]);
fprintf(file_E,linkformat, [nodeID+2 nodeID+3 dipolarLoops{k}.b zeros(1,3) 0]);
fprintf(file_E,linkformat, [nodeID+3 nodeID   dipolarLoops{k}.b zeros(1,3) 0]);

nodeID=nodeID+4;
end
fclose(file_V);
fclose(file_E);
grid on

