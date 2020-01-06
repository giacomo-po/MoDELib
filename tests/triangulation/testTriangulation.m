clear all
close all
clc

pointGenerator=6
switch pointGenerator
    case 0
N=30;
[x, y, dt] = simple_polygon(N);
P=[x y];
%P=P(1:end-1,:);
    case 1
        N=10;
        P=rand(N,2);
    case 2
P=[0.267224995715998   0.233928342494905
   0.261454542297731   0.603237361263948
   0.180668194442003   0.707105522792741
   0.068144352785810   0.623433948227880
   0.076800032913210   0.490713519607755
   0.867352151215692   0.510910106571687
   0.930827138816622   0.326255597187166
   0.754828309559499   0.262780609586236
   %0.492272679028383   0.268551063004503
   0.125848886968473   0.271436289713636
  -0.110739703180445   0.456090799098157
  -0.102084023053045   0.885989578758996
   0.596140840557176   0.946579339650792
   1.049121433891081   0.883104352049863
   1.098170287946344   0.332026050605432
   1.028924846927149   0.069470420074315
   0.821188523869562   0.049273833110383
   0.694238548667703   0.031962472855584
   0.679812415122038   0.626319174937013
   0.469190865355318   0.620548721518747
   0.457649958518785   0.031962472855584];

    case 3
        P=[  0 0
            0 4
            4 4
            4 1
            1 1
            1 2
            3 2
            3 0
];
            x=P(:,1);
    y=P(:,2);
    N=length(x);
            case 4
        P=[  0 0
            0 1
            1 0
            1 1];
                    case 5
        P=[  0 0
            0 1
            1 1
            -0.5 0.5
            1 0.5
            1 0];
                            case 6
        P=[  4.5 0 
            4.5 4
            0 4
            0 1
            4 1
            4 2
            1 2
            1 3
            2 3
            2 0
            3 0
            3 3
            3.5 3
            3.5 0
            ];
        P=flipud(P);

end
    x=P(:,1);
    y=P(:,2);
    N=length(x);

%% write input points
fid=fopen('I/I_0.txt','w');
for n=1:N
    fprintf(fid,'%i %e %e\n',[n x(n) y(n)]);
end
fclose(fid);

    plot([x;x(1)],[y;y(1)],'k','Linewidth',2)

%return
%% Call code
system('./test');

%% Read output



%% plot results
X=load("P/P_0.txt");
x=X(:,1);
y=X(:,2);
figure(1)
clf
hold on
%for(k=1:N)
%    k1=k+1;
%    if k==N
%        k1=1;
%    end
    plot([x;x(1)],[y;y(1)],'k','Linewidth',2)
%end

%return
T=load('T/T_0.txt');
X=load('X/X_0.txt');

%return
text(P(:,1),P(:,2),num2str([0:size(P,1)-1]'),'Fontsize',16,'color','r')
clr=colormap;
loopIDs=unique(T(:,1));
maxID=max(loopIDs)+1;
for k=1:length(loopIDs)
    loopID=loopIDs(k);
    u=ceil(loopID/maxID*size(colormap,1)+1);
    Y=X(find(X(:,1)==loopID),[2:3]);

    ids=T(find(T(:,1)==loopID),[2:4])+1; % one based
    for r=1:size(ids,1)
    patch(Y(ids(r,:),1),Y(ids(r,:),2),clr(u,:),'FaceAlpha',0.2)
    pause(0.5)
    end
%

end


    return
%return
for(t=1:size(O,1))
    ids=O(t,:)+1; % one-based
    C=[mean([x(ids) y(ids)]) 0];
    n=cross([x(2)-x(1) y(2)-y(1) 0],[x(3)-x(2) y(3)-y(2) 0]);
    n=n/norm(n);
    if(dot(n,[0 0 1])<0)
    patch(x(ids),y(ids),'g','FaceAlpha',0.5)
    else
            patch(x(ids),y(ids),clr(u,:),'FaceAlpha',0.5)

    end
    %quiver3(C(1),C(2))
end
%text(x,y,num2str([1:N]'),'FontSize',16,'Color','magenta')
axis equal


% 
% figure(2)
% C=load('C/C_0.txt');
% for f=[1:3:size(C,1)]
%         patch(C(f:f+2,1),C(f:f+2,2),'b','FaceAlpha',0.5)
% end
% axis equal
