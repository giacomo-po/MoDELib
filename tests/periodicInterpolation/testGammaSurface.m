clear all
close all
clc

system('rm input.txt');
system('rm output.txt');

%% Data points for gamma surface


%structure='gamma';
%structure='gammaPrime';
%structure='hcpBasal';
%structure='hcpPrismatic';
structure='hcpPrismaticEAM';

fid=fopen('input.txt','w');

switch structure
    
    case 'gamma'
        
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
        
        df=[0.25 sqrt(3)/12 -0.5 sqrt(3)/2 0; % USF
            0.75 sqrt(3)/12 0.5 sqrt(3)/2 0; % USF
            %            0.5  sqrt(3)/3 1 0 0; % USF
            ];
        
        printMatrixToFile(fid,N,'N');
        printMatrixToFile(fid,D,'D');
        readWaveVectors=0;
        
        cutDirs=[A(:,1) A(:,2) A(:,1)+A(:,2)];
        cutDirsLabels={'a_1', 'a_2', 'a_1+a_2'};
        cutDirsStyles={'k-', 'm--', 'b-'};

    case 'gammaPrime'
        
        A=[1 0.5;
            0 sqrt(3)/2];
        
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
        
        cutDirs=2*[A(:,1) A(:,2) A(:,1)+A(:,2)];
        cutDirsLabels={'a_1', 'a_2', 'a_1+a_2'};
        cutDirsStyles={'k-', 'm--', 'b-'};

        
    case 'hcpBasal'
        
        A=[1 0.5;
            0 sqrt(3)/2];
        
        D=[1 1];
        
        ISF=213;
        USF=261;
        
        waveVec=[0  0
            1  0
            0  1
            1  1
            1  2
            2  1
            1  -1
            ];
        
        
        f=[0 0 0
            0.5 sqrt(3)/6 ISF;    % ISF
            0.25 sqrt(3)/12 USF; % USF
            0.75 sqrt(3)/12 USF; % USF
            0.5  sqrt(3)/3 USF; % USF
            1 sqrt(3)/3 653
            ];
        
        df=[
            0 0 1 0 0
            0 0 0 1 0
            1 sqrt(3)/3 1 0 0
            1 sqrt(3)/3 0 1 0
            0.5 0 1 0 0
            0.25 sqrt(3)/4 0.5 sqrt(3)/2 0
            0.75 sqrt(3)/4 -0.5 sqrt(3)/2 0
            ];
        
        printMatrixToFile(fid,waveVec,'waveVecInt');
        readWaveVectors=1;
        
        figure(3)
        clf
        xx = [0.05:0.05:0.95];
        yy = [30.79 93.21 174.8 242.99 256.51 225.92 216.43 272.143 351.82 448.77 ...
            540.92 616.76 653.3 648.6 585.4 461.8 303.7 155.1 43];
        plot(xx,yy,'o')
        
        cutDirs=[A(:,1) A(:,2) A(:,1)+A(:,2)];
        cutDirsLabels={'a_1', 'a_2', 'a_1+a_2'};
        cutDirsStyles={'k-', 'm--', 'b-'};

    case 'hcpPrismatic'
        c=sqrt(8/3);
        A=[1 0;
            0 c];
        
        D=[1 1];
        
        waveVec=[
            0  0
            1  0
            0  1
            1 1
            1 -1
            2 0];
        
        
        h=c/4; % cannot be 0 or c/2
        hv=645;
        
        f=[0 0 0;
            1/2 0 211
            0 c/2 1350
            1/2 c/2 700
            1/4 h hv
            3/4 h hv
            ];
        
        df=[0 0 1 0 0
            0 0 0 1 0
            1/2 0 1 0 0
            1/2 c/2 1 0 0
            1/2 c/2 0 1 0
            ];
        
        printMatrixToFile(fid,waveVec,'waveVecInt');
        readWaveVectors=1;
        
        
        
        figure(3)
        clf
        xx = [0.05:0.05:0.95];
        
        yy = [12.23 48.38 93.21 139.05 173.68 205.26 220.54 220.07 215.45 210.866 ...
            216.47 222.07 221.56 203.22 174.19 138.54 92.7 48.9 12.22];
        
        plot(xx,yy,'ko')
        
        cutDirs=[A(:,1) A(:,2) A(:,1)+A(:,2)];
        cutDirsLabels={'a', 'c', 'a+c'}
        cutDirsStyles={'k-', 'm-', 'b-'};
        
         case 'hcpPrismaticEAM'
        c=sqrt(8/3);
        A=[1 0;
            0 c];
        
        D=[1 1];
        
        waveVec=[
            0  0
            1  0
            0  1
            1 1
            1 -1
            2 0
            0 2
            1 2
            1 -2
%            2 1
%            2 -1
%            3 0
            ];
        
        hx=1/4;
        hy=c/4; % cannot be 0 or c/2
        hv=700;
        
        f=[0 0 0;
            1/2 0.14*c 135
                        1/2 0.86*c 135
            0 c/3 960
            0 2/3*c 960
            0 1/2*c 840
            1/2 c/2 680
            hx hy hv
            1-hx hy hv
%            3/4 c-h hv
            hx c-hy hv
%            2/3 0 250
%            1/3 0 250
%            1/2 0 274
            ];
        
        df=[
            0 0 1 0 0
            0 0 0 1 0
            1/2 0.14*c 1 0 0
            1/2 0.14*c 0 1 0
           1/2 0.86*c 1 0 0
           1/2 0.86*c 0 1 0
            0 c/2 1 0 0
%             1/2 c/2 1 0 0
%             1/2 c/2 0 1 0
%            1/2 0 1 0 0
            ];
        
        printMatrixToFile(fid,waveVec,'waveVecInt');
        readWaveVectors=1;
        
        
        
        figure(3)
        clf
        xa = [0.201878788527690,0.321291819234938,0.423731452215650,0.548842582376464,0.656964545162808,0.816080292819118,0.992206234047144,1.23642025850583,1.48622725143186,1.63953451209906,1.87227874596026,2.12768922045225,2.36595283109042,2.53041974594575,2.66074011647416,2.79104471735465,2.93838579455341,3.07443054728195,3.18789842071799]/3.2;
    ya = [22.9071986111293,59.6945018983517,107.448938255163,148.217420577113,196.964687685115,226.756180907518,252.553353721915,271.338289754895,272.182082932366,272.073959135295,272.950911775655,264.845794286642,231.835978748330,199.864988767684,151.975853334769,101.098369618318,51.2071459664202,12.2797700407907,1.73156214876628];
        
        plot(xa,ya,'ko')
        
        xc = [0.122347197446092,0.241775997801275,0.355653882083610,0.458230185346419,0.549457598645899,0.686180446235612,0.788793545343601,0.902787073710788,1.02809795274544,1.13633555961663,1.24452585754402,1.34704959531372,1.46089594030018,1.55769003932045,1.69980084995435,1.83046289618801,1.93832728805810,2.01776426125209, ...
        2.13122687813882,2.25599107604509,2.38644811685560,2.47718141451981,2.59631584811355,2.73251829732144,2.83470035938589,2.94826285070954,3.06757075043056,3.18118055069801,3.27212411033468,3.38579173264456,3.50508911926696,3.60167295631476,3.73229820670324,3.84010477653090,3.95920767082878,4.07255464363066, ...
        4.16307242277310,4.27631426458875,4.39522792311142,4.51410478578890,4.61023655959593,4.72901354783649,4.83077508595603,4.92681749842476,5.04561551286258,5.17006957435946,5.27187842142281,5.37946421617939,5.49285849792507,5.56091504185987]/5.56091504185987;
    yc = [12.9920324724865,52.7676840432443,116.453407225604,190.106862039720,264.763003634988,354.334603901374,434.960871377073,520.561148638691,599.182042552087,669.843863739349,731.540640076004,795.232933945003,852.941960560292,906.676378827493,941.458308552414,958.316719262762,958.254297739689,950.239374177120, ...
        935.231925893051,910.256745977226,888.266629001616,869.287857712778,853.281007990878,844.237114900601,843.181862626336,847.100620137990,863.965601534978,876.849403897238,897.715276388871,921.556355790760,936.429104998724,950.318879485455,960.204477534221,949.184778971519,927.201232682548,890.279230319221, ...
        830.459699155402,773.615374901839,715.771649210445,650.955110857469,579.178901216189,495.436157067490,414.691617232285,325.981433984305,246.223154213653,162.477124721634,90.6976297370342,37.8410552048384,9.88409769211680,1.10483565637153];
    hold on
    plot(xc, yc, 'mo')
        
        cutDirs=[A(:,1) A(:,2) A(:,1)+A(:,2)];
        cutDirsLabels={'a', 'c', 'a+c'}
        cutDirsStyles={'k-', 'm-', 'b-'};
end

B=2*pi*inv(A');

printMatrixToFile(fid,A,'A');
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
plot(data(:,1)/2/pi,data(:,2)/2/pi,'bo','Linewidth',2)
plot(-data(:,1)/2/pi,-data(:,2)/2/pi,'mx','Linewidth',2)
plot([-max(data(:,1)) max(data(:,1))]/2/pi,[0 0],'k')
plot([0 0],[-max(data(:,2)) max(data(:,2))]/2/pi,'k')
title('wave vectors / 2\pi')
set(gca,'Fontsize',16)
print([structure '_waveVectors'], '-dpng', '-r300');


%% Construct GammaSurface and plot it
np=400;
%[X,Y] = meshgrid([0:np]/np*2*D(1),[0:np]/np*sqrt(3)/2*2*D(2));
[X,Y] = meshgrid([0:np]/np*sum(cutDirs(1,[1 2])),[0:np]/np*sum(cutDirs(2,[1 2])));
fs=zeros(size(X));
for i=1:size(data,1)
    k=data(i,[1 2]);
    S=data(i,3);
    C=data(i,4);
    fs=fs+S*sin(k(1)*X+k(2)*Y)+C*cos(k(1)*X+k(2)*Y);
end
fMax=max(max(fs));

figure(2)
hold on
surf(X,Y,fs,'edgecolor','none')
plot3([0 D(1)*A(1,1)],[0 D(1)*A(2,1)],[fMax fMax]+1,'m','Linewidth',2)
plot3([0 D(2)*A(1,2)],[0 D(2)*A(2,2)],[fMax fMax]+1,'m','Linewidth',2)
plot3([D(2)*A(1,2) D(1)*A(1,1)],[D(2)*A(2,2) D(1)*A(2,1)],[fMax fMax]+1,'m','Linewidth',2)
plot3([D(2)*A(1,2) D(2)*A(1,2)+D(1)*A(1,1)],[D(2)*A(2,2)  D(2)*A(2,2)],[fMax fMax]+1,'m','Linewidth',2)
plot3([D(1)*A(1,1) D(2)*A(1,2)+D(1)*A(1,1)],[D(1)*A(2,1)  D(2)*A(2,2)],[fMax fMax]+1,'m','Linewidth',2)
plot3(f(:,1),f(:,2),f(:,1)*0+fMax+1,'ok','Linewidth',2)
%plot3(df(:,1),df(:,2),df(:,1)*0+fMax+1,'xg','Linewidth',2)
quiver3(df(:,1),df(:,2),df(:,1)*0+fMax+1,df(:,3),df(:,4),df(:,1)*0,0.2,'k','Linewidth',2)
grid on
xlabel('x')
ylabel('y')
h = get(gca,'DataAspectRatio');
if h(3)==1
    set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))])
else
    set(gca,'DataAspectRatio',[1 1 h(3)])
end
colormap jet
h=colorbar
h.Label.String = 'mJ/m^2';
axis image
print([structure '_gammaSurface'], '-dpng', '-r300');

%% Plot cuts of GammaSurface
figure(3)
hold on
cutPlots=[];
for c=1:size(cutDirs,2)
    fcut=functionCut(data,cutDirs(:,c),np);
    cutPlots=[cutPlots plot([0:(np-1)]/np,fcut,cutDirsStyles{c},'Linewidth',2)];
end
grid on
xlabel('reaction coordinate')
ylabel('\gamma-surface [mJ/m^2]')
legend(cutPlots,cutDirsLabels)
set(gca,'Fontsize',16)
print([structure '_gammaSurfaceCuts'], '-dpng', '-r300');


