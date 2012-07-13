%%
close all,clear all,clc
% n,disp,load,ep11,ep12,ep13,ep22,ep23,ep33,dt
SS = load('./S/S_0.txt');
Disp = SS(21:end,2)./SS(21:end,1);
F = SS(21:end,3);

dt = SS(1:end,10);
% ep11 = SS(21:end,4).*dt;
% ep12 = SS(21:end,5).*dt;
% ep13 = SS(21:end,6).*dt;
% ep22 = SS(21:end,7).*dt;
% ep23 = SS(21:end,8).*dt;
ep33 = SS(1:end,9).*dt;
ep33_rate = SS(1:end,9);

mu  = 48e9; b = 0.25e-9;
f1 = mu*b^2;
D = 1e-6;
area = pi/4*D^2;
height = 4*D;
volume = area*height;
%%
figure(1)
%  subplot(2,2,1),plot(1e2*b*Disp/height/dt,'.b')%displ
%  subplot(2,2,2),
plot(f1*F/area/1e6,'.b')%load
%%
figure(2)
plot(1e2*b*Disp/height,f1*F/area/1e6,'.r')
xlabel('strain in %'),ylabel('stress [MPa]')
title('stress-strain curve')
%%
figure(3)
clf
% plot(cumsum(ep11)*b^3*1e2/volume,'r'),hold on
%  plot(cumsum(ep12)*b^3*1e2/volume,'b'),hold on
%  plot(cumsum(ep13)*b^3*1e2/volume,'g'),hold on
%  plot(cumsum(ep22)*b^3*1e2/volume,'k'),hold on
%  plot(cumsum(ep23)*b^3*1e2/volume,'m'),hold on

ep33(find(abs(ep33*b^3*1e2/volume)>1))=NaN

plot((ep33(1:end))*b^3*1e2/volume,'c'),hold on
ylabel('plastic strain %')
title('plastic strain')

return
%%
for i=1:size(dir('./E'))-4;
    i
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
density(i)=sum(L)*b/volume;
% figure(20)
% plot(i,size(E,1),'+r'),hold on
end
figure(10)
plot(density,'+r'),hold on

%%
E = (F*f1/area/1e9)./(Disp*b/height);