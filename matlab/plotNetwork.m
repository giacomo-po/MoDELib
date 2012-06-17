clc
close all
FS=16


step=191;

V=load(['../GiacomoCylinderDipolar/V/V_' num2str(step) '.txt']);
E=load(['../GiacomoCylinderDipolar/E/E_' num2str(step) '.txt']);
alpha=0.5;
%for e=[1 2 5 7 9 11 13 15]
figure(1)
hold on
axis equal
% sourceID sinkID bx by bz nx ny nz sourceTf sinkTf snID
for e=1:size(E,1)
source=E(e,1);
sink=E(e,2);
source_row=find(V(:,1)==source);
sink_row=find(V(:,1)==sink);
sourceTf=E(e,9);
sinkTf=E(e,10);
H=[ V(source_row,[2:4]);
    V(source_row,[5:7])*sourceTf;
    V(sink_row,[2:4]);
   -V(sink_row,[5:7])*sinkTf];

tubePlotter(H,alpha,'b-',1);

%text(H(1,1),H(1,2),H(1,3),num2str(source),'Fontsize',FS)
%text(H(3,1),H(3,2),H(3,3),num2str(sink),'Fontsize',FS)

boxSize=1000;
plot3([0 boxSize],[0 0],[0 0],'k')
plot3([0 0],[0 boxSize],[0 0],'k')
plot3([0 0],[0 0],[0 boxSize],'k')
end


return

%%
H2=[1.181e+04 0.000e+00 9.344e+03 -1.483e+02;
1.037e+04 0.000e+00 8.237e+03 -2.124e+02;
2.000e+04 0.000e+00 1.980e+04 -6.636e+02];
tubePlotter(H2',0,'m',2);
xlabel('x')
ylabel('y')



%return
%% 
%sourceV=[46  899 747 803 879 926 945 955];
%sinkV =[841 790 762 644 495 47  511 945];
sourceV=[2069];
sinkV = [3990];

for n=1:length(sourceV)
source=sourceV(n);
sink=sinkV(n);
source_row=find(V(:,1)==source)
sink_row=find(V(:,1)==sink)
sourceTf=E(e,9);
sinkTf=E(e,10);
H1=[ V(source_row,[2:4]);
    V(source_row,[5:7])*sourceTf;
    V(sink_row,[2:4]);
   -V(sink_row,[5:7])*sinkTf];

tubePlotter(H1,alpha,'m',2);
text(H1(1,1),H1(1,2),H1(1,3),num2str(source),'Fontsize',FS)
text(H1(3,1),H1(3,2),H1(3,3),num2str(sink),'Fontsize',FS)

end



% source=278;
% sink=14;
% source_row=find(V(:,1)==source)
% sink_row=find(V(:,1)==sink)
% sourceTf=E(e,9);
% sinkTf=E(e,10);
% H2=[ V(source_row,[2:4]);
%     V(source_row,[5:7])*sourceTf;
%     V(sink_row,[2:4]);
%    -V(sink_row,[5:7])*sinkTf];
% 
% tubePlotter(H2,alpha,'m',2);
%%

%H1L=H1'
%H2L=H2'

%B1=H2BZplanar(H1,0.5)
%B2=H2BZplanar(H2,0.5)
%CubicSplineIntersection(B1,B2)

