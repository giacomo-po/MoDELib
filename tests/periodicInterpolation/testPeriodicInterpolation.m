clear all
close all
clc

system('rm input.txt');

%% Data points for gamma surface
A=[1 0.5;  
  0 sqrt(3)/2];

N=[2 2];
D=[1 1];

f=[0 0 0
  0.5 sqrt(3)/6 42;    % ISF
  0.25 sqrt(3)/12 182; % USF
  0.75 sqrt(3)/12 182; % USF
  0.5  sqrt(3)/3 182; % USF
  ]; 

% df=[0.5 sqrt(3)/6  1 0 0; % ISF
%     0.5 sqrt(3)/6 0 1 0; % ISF
% %     0.25 sqrt(3)/12 0.25  sqrt(3)/12 0;
% %     0.75 sqrt(3)/12 0.25 -sqrt(3)/12 0;
%     0.25 sqrt(3)/12 1 0 0;
%     0.75 sqrt(3)/12 1 0 0;
%         0.25 sqrt(3)/12 0 1 0;
%     0.75 sqrt(3)/12 0 1 0;
%     ];

df=[0.25 sqrt(3)/12 -0.5 sqrt(3)/2 0; % USF
    0.75 sqrt(3)/12 0.5 sqrt(3)/2 0; % USF
    0.5  sqrt(3)/3 1 0 0; % USF
    ];


%% Write input file
fid=fopen('input.txt','w')
printMatrixToFile(fid,A,'A');
printMatrixToFile(fid,N,'N');
printMatrixToFile(fid,D,'D');
printMatrixToFile(fid,f,'f');
printMatrixToFile(fid,df,'df');
fclose(fid)

%% run code
system('./test')

%% Load data and plot
data=load('output.txt')

np=400;
[X,Y] = meshgrid([0:1/(np-1):1],[0:sqrt(3)/(np-1):sqrt(3)]);

f=zeros(size(X));
for i=1:size(data,1)
k=data(i,[1 2]);
S=data(i,3);
C=data(i,4);
    f=f+S*sin(k(1)*X+k(2)*Y)+C*cos(k(1)*X+k(2)*Y);
end

fMax=max(max(f));

B=inv(A');

figure
clf
hold on
surf(X,Y,f,'edgecolor','none')
plot3([0 A(1,1)],[0 A(2,1)],[fMax fMax]+1,'m','Linewidth',2)
plot3([0 A(1,2)],[0 A(2,2)],[fMax fMax]+1,'m','Linewidth',2)
plot3([0 B(1,1)],[0 B(2,1)],[fMax fMax]+1,'g','Linewidth',2)
plot3([0 B(1,2)],[0 B(2,2)],[fMax fMax]+1,'g','Linewidth',2)
grid on
xlabel('x')
ylabel('y')
h = get(gca,'DataAspectRatio') 
if h(3)==1
      set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))])
else
      set(gca,'DataAspectRatio',[1 1 h(3)])
end
colormap jet
colorbar

f1=functionCut(data,[1 0],1,np);
f2=functionCut(data,[0.5 sqrt(3)/6],sqrt(3),np);
figure(2)
hold on
plot([0:(np-1)]/np,f1,'Linewidth',2)
plot([0:(np-1)]/np,f2,'Linewidth',2)
axis([0 1 0 max(max(f1),max(f2))])
grid on
xlabel('reaction coordinate')
ylabel('\gamma-surface [mJ/m^2]')
set(gca,'Fontsize',16)