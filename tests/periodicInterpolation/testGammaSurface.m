clear all
close all
clc

system('rm input.txt');

%% Data points for gamma surface
A=[1 0.5;
    0 sqrt(3)/2];
B=2*pi*inv(A');

gamma=0;
gammaPrime=1;
structure=gamma

fid=fopen('input.txt','w')
printMatrixToFile(fid,A,'A');

switch structure
    
    case gamma
        N=[2 2];
        D=[1 1];
        
        f=[0 0 0
            0.5 sqrt(3)/6 42;    % ISF
            0.25 sqrt(3)/12 182; % USF
            0.75 sqrt(3)/12 182; % USF
            0.5  sqrt(3)/3 182; % USF
            ];
        
        df=[0.25 sqrt(3)/12 -0.5 sqrt(3)/2 0; % USF
            0.75 sqrt(3)/12 0.5 sqrt(3)/2 0; % USF
            0.5  sqrt(3)/3 1 0 0; % USF
            ];
        
        printMatrixToFile(fid,N,'N');
        printMatrixToFile(fid,D,'D');
        readWaveVectors=0;
        
        
    case gammaPrime
        
        N=[3 3];
        D=[2 2];
        
        
        waveVec_ctr=0;
        
        for i=0:N(1)-1
            for j=0:N(2)-1
                % if( ~(i==(N(1)-1) & j==(N(2)-1) ) )
                waveVec_ctr= waveVec_ctr + 1;
                waveVec(waveVec_ctr,:)=[i/D(1) j/D(2)];
                % end
            end
        end
        
        %extra_points of waveVec
        
        % % N1=N2=3 case
        waveVec_ctr=waveVec_ctr+1;
        waveVec(waveVec_ctr,:)= [1/D(1) -1/D(2)];
        
        APB=175;
        SISF=10;
        CESF=270;
        CISF=230;
        SESF=75;
        
        
        
        
        f=[0 0 0;
            1 sqrt(3)/3 SISF;
            1 0 APB;
            0.5 sqrt(3)/2 APB;
            1.5 sqrt(3)/2 APB;
            0.5 sqrt(3)/6 CESF;
            1 sqrt(3)*2/3 CESF;
            1.5 sqrt(3)/6 CESF; % ---- till here cell 1
            2 2*sqrt(3)/3 SESF;
            2 sqrt(3)/3 CISF;
            1.5 5*sqrt(3)/6 CISF;
            2.5 5*sqrt(3)/6 CISF;
            ]; %  12 conditions
        
        
        df=[0 0 1 0 0; %PC
            0 0 0 1 0; %PC
            1 sqrt(3)/3 1 0 0; %SISF
            1 sqrt(3)/3 0 1 0; %SISF
            1 0 1 0 0; %APB
            0.5 sqrt(3)/2 0.5 sqrt(3)/2 0; %APB
            1.5 sqrt(3)/2 -0.5 sqrt(3)/2 0 %APB
            ]; % 7 conditions
        
        printMatrixToFile(fid,waveVec,'waveVecInt');
        readWaveVectors=1;
        
end
printMatrixToFile(fid,f,'f');
printMatrixToFile(fid,df,'df');
fclose(fid)

%% run code
system(['./test ' num2str(readWaveVectors)])

%% Load data and plot
data=load('output.txt')

%% Plot Wave vectors
figure(1)
hold on;
grid on
axis equal
plot(data(:,1),data(:,2),'ro')

%% Construct GammaSurface and plot it
np=400;
[X,Y] = meshgrid([0:np]/np*2*D(1),[0:np]/np*sqrt(3)/2*2*D(2));
f=zeros(size(X));
for i=1:size(data,1)
    k=data(i,[1 2]);
    S=data(i,3);
    C=data(i,4);
    f=f+S*sin(k(1)*X+k(2)*Y)+C*cos(k(1)*X+k(2)*Y);
end
fMax=max(max(f));

figure(2)
clf
hold on
surf(X,Y,f,'edgecolor','none')
plot3([0 D(1)*A(1,1)],[0 D(1)*A(2,1)],[fMax fMax]+1,'m','Linewidth',2)
plot3([0 D(2)*A(1,2)],[0 D(2)*A(2,2)],[fMax fMax]+1,'m','Linewidth',2)
plot3([D(2)*A(1,2) D(1)*A(1,1)],[D(2)*A(2,2) D(1)*A(2,1)],[fMax fMax]+1,'m','Linewidth',2)
plot3([D(2)*A(1,2) D(2)*A(1,2)+D(1)*A(1,1)],[D(2)*A(2,2)  D(2)*A(2,2)],[fMax fMax]+1,'m','Linewidth',2)
plot3([D(1)*A(1,1) D(2)*A(1,2)+D(1)*A(1,1)],[D(1)*A(2,1)  D(2)*A(2,2)],[fMax fMax]+1,'m','Linewidth',2)
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
h=colorbar
h.Label.String = 'mJ/m^2';
print('gammaSurface_gammaPrime', '-dpng', '-r300');

%% Plot cuts of GammaSurface
f1=functionCut(data,A(:,1),D(1),np);
f2=functionCut(data,(A(:,1)+A(:,2))/2,D(2)*sqrt(3),np);
figure(3)
hold on
plot([0:(np-1)]/np,f1,'Linewidth',2)
plot([0:(np-1)]/np,f2,'Linewidth',2)
axis([0 1 0 max(max(f1),max(f2))])
grid on
xlabel('reaction coordinate')
ylabel('\gamma-surface [mJ/m^2]')
set(gca,'Fontsize',16)


