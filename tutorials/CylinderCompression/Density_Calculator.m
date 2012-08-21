clc
close all

b = 0.25e-9;
D = 1e-6;
area = pi/4*D^2;
height = 4*D;
volume = area*height;
% sourceID sinkID bx by bz nx ny nz sourceTf sinkTf snID

for i=1:303;
V=load(['./V/V_' num2str(i-1) '.txt']);
E=load(['./E/E_' num2str(i-1) '.txt']);

for e=1:size(E,1)
source=E(e,1);
sink=E(e,2);
source_row=find(V(:,1)==source);
sink_row=find(V(:,1)==sink);
p1 = [V(source_row,2),V(source_row,3),V(source_row,4)];
p2 = [V(sink_row,2),V(sink_row,3),V(sink_row,4)];
L(e) = (norm(p2-p1));
% sourceTf=E(e,9);
% sinkTf=E(e,10);
end
density=sum(L)*b/volume;
figure(1)
plot(i,density,'+r'),hold on
figure(2)
plot(i,size(E,1),'+r'),hold on
end
 