% This file is part of MODEL, the Mechanics Of Defect Evolution Library.
%
% Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
% model is distributed without any warranty under the
% GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
clc
close all
clear all

L=4000;             % units of Burgers vector
V=L^3;              % units of Burgers vector
rhoMax=5e13;        % [m^-2]
Burgers=0.2556e-9;  % [m]

%% Generate dipoles in cube with density rhoMax
lMin=L/20;
lMax=L;

N =[-0.5774    0.5774   -0.5774    0.5774;
     0.5774   -0.5774   -0.5774    0.5774;
    -0.5774   -0.5774    0.5774    0.5774];

B=[0  1  1;
   0  1 -1;
   1  0  1;
   1  0 -1;
   1  1  0;
   1 -1  0];

totalLength=0;
density=totalLength/V/Burgers^2;
k=1;

figure(1)
hold on
while density<rhoMax
r=ceil(rand(1)*size(B,1));  % a random row into B

b=B(r,:)/norm(B(r,:));
Nc=[];
for n=1:size(N,2)
if abs(dot(b,N(:,n)))<0.001
Nc=[Nc N(:,n)];
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


%% Write E and V files
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

%% Run DD code
N=10;   % number of subdivisions
system(['./DDomp ' num2str(L/2) ' ' num2str(10)]);

%% Plot PK for different cell sizes
for n=1:N
P=load(['P/P_' num2str(n-1) '.txt']);
size(P)
pk_noMult{n}=P(:,[5:7]);
end

for n=1:N
P=load(['P/P_' num2str(n-1+N) '.txt']);
size(P)
pk_Mult{n}=P(:,[5:7]);
end

% compute norms
for n=1:N
pk_noMult_norm(n)=norm(pk_noMult{n}-pk_noMult{N})/norm(pk_noMult{N});
pk_Mult_norm(n)=norm(pk_Mult{n}-pk_Mult{N})/norm(pk_Mult{N});
end
% make sure that for cell size==systemSize, then multipole expansion does
% nothing
if(norm(pk_noMult{N}-pk_Mult{N})/norm(pk_noMult{N})>1.0e-7)
error('at cellSize==systemSize, PK with and without multipole must be the same.')
end


%% Plot
fontsize=14;
figure(1)
clf
hold on
grid on
cellSizes=[1:N]/N*L/2;
plot(cellSizes,pk_noMult_norm,'b')
plot(cellSizes,pk_Mult_norm,'r')
ylabel('error in PK force','FontSize',fontsize)
xlabel('cell size / b','FontSize',fontsize)
legend('nearest neighbor','nearest neighbor + multipole')
set(gca,'FontSize',fontsize)
print(gcf,'-depsc','error_PK')