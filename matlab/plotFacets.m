function plotFacets(Facets)
figure(1)
clf
hold on
axis equal
for k=1:length(Facets)
v1=Facets{k}.v1;
v2=Facets{k}.v2;
v3=Facets{k}.v3;
n=Facets{k}.normal;
%plot3([v1(1) v2(1) v3(1) v1(1)],[v1(2) v2(2) v3(2) v1(2)],[v1(3) v2(3) v3(3) v1(3)],'-o')
fill3([v1(1) v2(1) v3(1) v1(1)],[v1(2) v2(2) v3(2) v1(2)],[v1(3) v2(3) v3(3) v1(3)],'g','Facealpha',0.2)
quiver3((v1(1)+v2(1)+v3(1))/3,(v1(2)+v2(2)+v3(2))/3,(v1(3)+v2(3)+v3(3))/3,n(1),n(2),n(3),1000,'r')
end
%axis equal

