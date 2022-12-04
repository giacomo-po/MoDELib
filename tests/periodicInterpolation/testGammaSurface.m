clear all
close all
clc
fontSize=16;

%% Select structure
%structure='gamma';
structure='gammaPrime';
%structure='hcpBasal';
%structure='hcpBasalEAM';
%structure='hcpPrismatic';
%structure='hcpPrismaticEAM';


%% Input data
system('rm build/input.txt');
system('rm build/output.txt');
system('rm build/points.txt');
system('rm build/data.txt');

fid=fopen('build/input.txt','w');

switch structure
    
    case 'gamma'
        
        cellScale=1;
        A=[1 -0.5;
            0 sqrt(3)/2];
        rotSymm=3;
        mirSymm=[];
        
        ISF=142;
        USF=863;
        MSF=2800;
        
        waveVec=[0  0
            0  1
            1  1
            ]/cellScale;
        
        
        f=[0 0 0
            0.5 sqrt(3)/6 ISF;    % ISF
            0.25 sqrt(3)/12 USF; % USF
            1 sqrt(3)/3 MSF
            ];
        
        printMatrixToFile(fid,waveVec,'waveVecInt');
        
        cutDirs={A(:,1),A(:,1)+A(:,2),2*A(:,1)+A(:,2) };
        cutDirsLabels={'a_1', 'a_2', 'a_1+a_2'};
        cutDirsStyles={'m-', 'k--', 'k-.'};
        x1Lab='a/2\langle110\rangle';
        y1Lab='a/2\langle112\rangle';
        y1TicksFactor=sqrt(3); % length of 1/2a<112> in units of b
        
    case 'gammaPrime'
        cellScale=2;

        A=[1 -0.5;
            0 sqrt(3)/2];
        rotSymm=3;
        mirSymm=[];
        
%          APB=175;
%          SISF=10;
%          CESF=270;
%         CISF=230;
%         SESF=75;
        
         CESF=1350;
         SISF=74;
         APB=162;
        SESF=1800;
        CISF=191;



        
        f=[0 0 0;
            0.5 sqrt(3)/6 CISF;
            1 sqrt(3)/3 SESF;
            1.5 sqrt(3)/2 APB;
            2 2*sqrt(3)/3 SISF;
            2.5 5*sqrt(3)/6 CESF;
            ]; %  6 conditions
        
figure(1)
hold on
textH=1.5*max([CESF,SISF,APB,SESF,CISF]);
text(0.5, sqrt(3)/6,textH, 'CISF','Color','b','FontSize',16)
text(1, sqrt(3)/3,textH, 'SESF','Color','b','FontSize',16)
text(1.5, sqrt(3)/2,textH, 'APB','Color','b','FontSize',16)
text(2, 2*sqrt(3)/3,textH, 'SISF','Color','b','FontSize',16)
text(2.5, 5*sqrt(3)/6,textH,'CESF','Color','b','FontSize',16)

        waveVec=[0  0
            0  1
            1  1
            2 0
            ]/cellScale;
        
        
        
        printMatrixToFile(fid,waveVec,'waveVecInt');
        
%        cutDirs={A(:,1),A(:,2) A(:,1)+A(:,2)};
                cutDirs={A(:,1),A(:,1)+A(:,2),2*A(:,1)+A(:,2) };

        cutDirsLabels={'a_1', 'a_2', 'a_1+a_2'};
        cutDirsStyles={'k-', 'm--', 'k.-'};
%        cutDirs={A(:,1),A(:,2) A(:,1)+A(:,2)};
        cutDirsLabels={'a_1', 'a_2', 'a_1+a_2'};
        cutDirsStyles={'m-', 'k--', 'k-.'};
        x1Lab='a/2\langle110\rangle';
        y1Lab='a/2\langle112\rangle';
        y1TicksFactor=sqrt(3); % length of 1/2a<112> in units of b
        
        
    case 'hcpBasal'
        
        cellScale=1;
        A=[1 0.5;
            0 sqrt(3)/2];
        rotSymm=3;
        mirSymm=[];
        
        
        ISF=213;
        USF=261;
        
        waveVec=[0  0
            0  1
            1  -1
            ]/cellScale;
        
        
        f=[0 0 0
            0.5 sqrt(3)/6 ISF;    % ISF
            0.25 sqrt(3)/12 USF; % USF
            1 sqrt(3)/3 653
            ];
        
        printMatrixToFile(fid,waveVec,'waveVecInt');
        
        figure(2)
        clf
        xx = [0.05:0.05:0.95];
        yy = [30.79 93.21 174.8 242.99 256.51 225.92 216.43 272.143 351.82 448.77 ...
            540.92 616.76 653.3 648.6 585.4 461.8 303.7 155.1 43];
        plot(xx,yy,'o')
        
        cutDirs={A(:,1),A(:,2) A(:,1)+A(:,2)};
        cutDirsLabels={'a_1', 'a_2', 'a_1+a_2'};
        cutDirsStyles={'m-', 'k--', 'k-.'};
        x1Lab='a/2\langle110\rangle';
        y1Lab='a/2\langle112\rangle';
        y1TicksFactor=sqrt(3); % length of 1/2a<112> in units of b
        
    case 'hcpBasalEAM'
        
        cellScale=1;
        A=[1 0.5;
            0 sqrt(3)/2];
        rotSymm=3;
        mirSymm=[];
        
        
        ISF=198;
        USF=323;
        
        waveVec=[0  0
            0  1
            1  -1
            ]/cellScale;
        
        
        f=[0 0 0
            0.5 sqrt(3)/6 ISF;    % ISF
            0.25 sqrt(3)/12 USF; % USF
            1 sqrt(3)/3 371
            ];
        
        printMatrixToFile(fid,waveVec,'waveVecInt');
        
        figure(2)
        clf
        hold on
        xx = [0.05:0.05:0.95];
        yy = [59.9 192.1 292.0 325.3 287.2 214.0 202.6 267.3 406.1 472.7 ...
            453.7 407.1 371.9 378.5 411.8 428.0 382.3 214.9 54.2];
        plot(xx, yy ,'bo')
        
        %        cutDirs=[A(:,1) A(:,2) A(:,1)+A(:,2)];
        cutDirs={A(:,1),A(:,2) A(:,1)+A(:,2)};
        cutDirsLabels={'a_1', 'a_2', 'a_1+a_2'};
        cutDirsStyles={'k-', 'm--', 'b-'};
        
        cutDirsLabels={'a_1', 'a_2', 'a_1+a_2'};
        cutDirsStyles={'m-', 'k--', 'k-.'};
        x1Lab='a/2\langle110\rangle';
        y1Lab='a/2\langle112\rangle';
        y1TicksFactor=sqrt(3); % length of 1/2a<112> in units of b
        
        
        
    case 'hcpPrismatic'
        cellScale=1;
        c=sqrt(8/3);
        A=[1 0;
            0 c];
        rotSymm=1;
        mirSymm=[1 0
            0 1]; % each row is the normal vector to a plane of symm
        
        
        waveVec=[
            0  0
            1  0
            0  1
            1 1
            2 0
            ]/cellScale;
        
               
        f=[0 0 0;
            1/2 0 211
            0 c/2 1350
            1/4 c/4 645
            1/2 c/2 700
            ];
        
        printMatrixToFile(fid,waveVec,'waveVecInt');
        
        
        
        figure(2)
        clf
        xx = [0.05:0.05:0.95];
        
        yy = [12.23 48.38 93.21 139.05 173.68 205.26 220.54 220.07 215.45 210.866 ...
            216.47 222.07 221.56 203.22 174.19 138.54 92.7 48.9 12.22];
        
        plot(xx,yy,'ko')
        
        
        
        cutDirs={A(:,1),A(:,2) A(:,1)+A(:,2)};
        cutDirsLabels={'a', 'c', 'a+c'}
        cutDirsStyles={'m-', 'k--', 'k-.'};
        
        x1Lab='a';
        y1Lab='c';
        y1TicksFactor=c; % length of 1/2a<112> in units of b
        
        
    case 'hcpPrismaticEAM'

        cellScale=1;
        c=sqrt(8/3);
        A=[1 0;
            0 c];
        rotSymm=1;
        mirSymm=[1 0
            0 1]; % each row is the normal vector to a plane of symm
        
        
        waveVec=[
            0  0
            1  0
            0  1
            1 1
            2 0
            0 2
            1 2
            %            2 1
            ]/cellScale;
        
        
        hx=1/8;
        hy=c/4; % cannot be 0 or c/2
        hv=900;
        
        
        L=4.3292*0.1273;
        
        
        f=[ 0 0 0;
            0.5 0 272.5;
            0 c/2 840;
            0.5 c/2 704;
            0.5*L 0.14*L*c 289.4;
            0 c/3 960;
            0.5 0.14*c 135;
            ];
        
        
        printMatrixToFile(fid,waveVec,'waveVecInt');
        
        
        figure(2)
        clf
        xa = [0.201878788527690,0.321291819234938,0.423731452215650,0.548842582376464,0.656964545162808,0.816080292819118,0.992206234047144,1.23642025850583,1.48622725143186,1.63953451209906,1.87227874596026,2.12768922045225,2.36595283109042,2.53041974594575,2.66074011647416,2.79104471735465,2.93838579455341,3.07443054728195,3.18789842071799]/3.2;
        ya = [22.9071986111293,59.6945018983517,107.448938255163,148.217420577113,196.964687685115,226.756180907518,252.553353721915,271.338289754895,272.182082932366,272.073959135295,272.950911775655,264.845794286642,231.835978748330,199.864988767684,151.975853334769,101.098369618318,51.2071459664202,12.2797700407907,1.73156214876628];
        hold on
        plot(xa, ya ,'ko')
        
        xc = [0.122347197446092,0.241775997801275,0.355653882083610,0.458230185346419,0.549457598645899,0.686180446235612,0.788793545343601,0.902787073710788,1.02809795274544,1.13633555961663,1.24452585754402,1.34704959531372,1.46089594030018,1.55769003932045,1.69980084995435,1.83046289618801,1.93832728805810,2.01776426125209, ...
            2.13122687813882,2.25599107604509,2.38644811685560,2.47718141451981,2.59631584811355,2.73251829732144,2.83470035938589,2.94826285070954,3.06757075043056,3.18118055069801,3.27212411033468,3.38579173264456,3.50508911926696,3.60167295631476,3.73229820670324,3.84010477653090,3.95920767082878,4.07255464363066, ...
            4.16307242277310,4.27631426458875,4.39522792311142,4.51410478578890,4.61023655959593,4.72901354783649,4.83077508595603,4.92681749842476,5.04561551286258,5.17006957435946,5.27187842142281,5.37946421617939,5.49285849792507,5.56091504185987]/5.56091504185987;
        yc = [12.9920324724865,52.7676840432443,116.453407225604,190.106862039720,264.763003634988,354.334603901374,434.960871377073,520.561148638691,599.182042552087,669.843863739349,731.540640076004,795.232933945003,852.941960560292,906.676378827493,941.458308552414,958.316719262762,958.254297739689,950.239374177120, ...
            935.231925893051,910.256745977226,888.266629001616,869.287857712778,853.281007990878,844.237114900601,843.181862626336,847.100620137990,863.965601534978,876.849403897238,897.715276388871,921.556355790760,936.429104998724,950.318879485455,960.204477534221,949.184778971519,927.201232682548,890.279230319221, ...
            830.459699155402,773.615374901839,715.771649210445,650.955110857469,579.178901216189,495.436157067490,414.691617232285,325.981433984305,246.223154213653,162.477124721634,90.6976297370342,37.8410552048384,9.88409769211680,1.10483565637153];
        hold on
        plot(xc, yc, 'ks')
        
        xac = [0.0257791300462219,0.184847568758729,0.383913091181738,0.543365257993977,0.714055927234007,0.958338286833744,1.11156144271197,1.27035128431097,1.49148905729258,1.65025787269433,1.80910027978645,2.00228898007931,2.25237457011884,2.38898177362370,2.63358478273143,2.82139077652941,3.02054566029071, ...
            3.23663714593337,3.40154561093085,3.60603589224323,3.78771800609341,3.98067016166725,4.21330926454222,4.38332183892108,4.58733903079543,4.77420410226666,4.98381426260836,5.17070036027683,5.35195143708348,5.55042296943430,5.76606239183618,5.94208320207797,6.10132510691776,6.30595731506155,6.48208851283888,...
            6.66943193029745,6.86226318523713,7.03783193223815,7.20200973688134,7.37769412796721,7.53086997490163]/7.53086997490163;
        yac = [2.09060626929080,22.9170546410882,92.5301942407240,186.073117511879,257.694916017136,289.429361278768,278.383379969509,246.415675332183,195.485626122618,159.533457107245,137.526913415037,169.290926766546,224.928872953140,292.585919140268,385.083446167117,472.633246154603,559.180359360939,627.787526836104,...
            679.490288794446,701.286570514198,704.169787811505,691.108576909985,668.063207660110,611.186029973351,543.331863187045,452.576882256928,368.781572615114,282.011056063044,203.212753610387,160.264774530195,143.205957877349,149.080808801511,202.779088292196,251.470504563765,278.263793472674,278.155377143126, ...
            242.183496067836,162.392362863987,75.6349876851955,17.7584085606053,1.25261759925980];
        plot(xac, yac ,'k^')
        
        cutDirs={A(:,1),A(:,2),[[0.5;0.14*c] [0.5;0.86*c] A(:,2)]};

        cutDirsStyles={'k-', 'm-', 'b-'};
        
        cutDirsLabels={'a', 'c', 'path'}
        cutDirsStyles={'m-', 'k--', 'k-.'};
        
        x1Lab='a';
        y1Lab='c';
        y1TicksFactor=c; % length of 1/2a<112> in units of b
        
end

B=2*pi*inv(A');

printMatrixToFile(fid,A,'A');
printMatrixToFile(fid,f,'f');
printMatrixToFile(fid,mirSymm,'m');
fclose(fid)

%% Plot wave vectors
figure(3)
clf
hold on

quiver(0,0,B(1,1),B(2,1),0,'k','Linewidth',2,'MaxHeadSize',0.5)
quiver(0,0,B(1,2),B(2,2),0,'k','Linewidth',2,'MaxHeadSize',0.5)
text(0.5*B(1,1),0.5*B(2,1)+0.7,'b_1','FontSize',fontSize)
text(0.5*B(2,1)-0.8,0.5*B(2,2),'b_2','FontSize',fontSize)
nMax=max(waveVec)+1;
for i=-nMax:nMax
    for j=-nMax:nMax
        k1=B*[i;j];
        plot(k1(1),k1(2),'ko','Linewidth',1)
    end
end

for i=1:size(waveVec,1)
    n=waveVec(i,:)';
    k=B*n;
    plot(k(1),k(2),'mo','Linewidth',2)
    kNorm=norm(k);
    np=1000;
    phi=2*pi*[0:np]/np;
    plot(kNorm*cos(phi),kNorm*sin(phi),'k--')
end
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'Visible','off')
axis equal
axis image
set(gca,'FontSize',fontSize)
print(['build/' structure '_waveVectors'], '-depsc');

%return
%% Run and plot gamma-surface
figure(1)
%clf
hold on
fMax=0;
plot3(f(:,1),f(:,2),f(:,1)*0+fMax+1,'om','Linewidth',2)

fid=fopen('build/points.txt','w');
np=400;
printMatrixToFile(fid,cellScale*[0:np-1]'/(np-1)*(abs(A(1,1))+abs(A(1,2))),'x');
printMatrixToFile(fid,cellScale*[0:np-1]'/(np-1)*(abs(A(2,1))+abs(A(2,2))),'y');
fclose(fid)

system(['build/test build ' num2str(rotSymm) ' 1'])
%return
wsc=load([pwd '/build/output.txt']); % wave vector, S and C coeffs in each row
data=load([pwd '/build/data.txt']);
X=reshape(data(:,1),np,np);
Y=reshape(data(:,2),np,np);
G=reshape(data(:,3),np,np);
fMax=max(data(:,3));

%return
%surf(X,Y,G,'edgecolor','none')
contour(X,Y,G,30)

xlabel(x1Lab)
ylabel(y1Lab)
yticks([0:0.25:1]*y1TicksFactor*0.9999) % *0.999 makes sure upper limit is included
yticklabels(num2str([0:0.25:1]'))
h = get(gca,'DataAspectRatio');
if h(3)==1
    set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))])
else
    set(gca,'DataAspectRatio',[1 1 h(3)])
end
axis image
colormap jet
h=colorbar
h.Label.String = 'mJ/m^2';
h.Label.FontSize = fontSize;
view([0,90])
%plot3([0 A(1,1)],[0 A(2,1)],[fMax fMax]+1,'k','Linewidth',2)
%plot3([A(1,1) A(1,1)+A(1,2)],[A(2,1) A(2,1)+A(2,2)],[fMax fMax]+1,'k','Linewidth',2)
%plot3([A(1,1)+A(1,2) A(1,2)],[A(2,1)+A(2,2) A(2,2)],[fMax fMax]+1,'k','Linewidth',2)
%plot3([0 A(1,2)],[0 A(2,2)],[fMax fMax]+1,'k','Linewidth',2)
plot3(f(:,1),f(:,2),f(:,1)*0+fMax+2,'ok','Linewidth',2)
%set(gca,'FontSize',fontSize)
%print(['build/' structure '_gammaSurface'], '-depsc');

%return
%% Plot cuts of GammaSurface
np=100;
figure(2)
hold on
cutPlots=[];
for c=1:size(cutDirs,2)
    fid=fopen('build/points.txt','w');
    cpts=[[0;0] cutDirs{c}];
    pts=[[0;0]];
    rc=[0];
    for k=1:size(cutDirs{c},2)
        dx=(cpts(:,k+1)-cpts(:,k))*cellScale;
        if length(rc)==0
            rc=norm(dx)*[1:np-1]/(np-1);
            pts=[dx*[1:np-1]/(np-1)];
        else
            rc=[rc rc(end)+norm(dx)*[1:np-1]/(np-1)];
            pts=[pts pts(:,end)+dx*[1:np-1]/(np-1)];
        end
    end
    rc=rc/rc(end);
    
    printMatrixToFile(fid,pts(1,:)','x');
    printMatrixToFile(fid,pts(2,:)','y');
    fclose(fid)
    
    figure(1)
    plot3(pts(1,:),pts(2,:),pts(1,:)*0+fMax+1,cutDirsStyles{c},'Linewidth',2)
    figure(2)
    system(['build/test build ' num2str(rotSymm) ' 0'])
    wsc=load('build/output.txt'); % wave vector, S and C coeffs in each row
    data=load('build/data.txt');
    cutPlots=[cutPlots plot(rc,data(:,3),cutDirsStyles{c},'Linewidth',2)];
end
grid on
xlabel('reaction coordinate')
ylabel('misfit energy [mJ/m^2]')
legend(cutPlots,cutDirsLabels)
set(gca,'Fontsize',fontSize)
print(['build/' structure '_gammaSurfaceCuts'], '-depsc');



