clear all
close all
clc

system('rm poly0.txt');
system('rm poly1.txt');
system('rm result.txt');

%N0=4;
%poly0=rand(N0,2);
h=1.5;
poly0=[0 0;
       2 0;
       2 h;
       1 h;
       1 1;
       4 1;
       4 h;
       3 h;
       3 0;
       5 0;
       5 3;
       0 3];
   poly0=flipud(poly0);
   
   A0=area(poly0)

dlmwrite('poly0.txt',poly0,' ');
%writeGPCfile('poly0.txt',poly0)

%N1=4;
%poly1=rand(N1,2);
poly1=[-1 -1;
        6 -1;
        6 4;
       -1 4];
dlmwrite('poly1.txt',poly1,' ');
%writeGPCfile('poly1.txt',poly1)


figure(1)
clf
hold on
axis equal
patch(poly0(:,1),poly0(:,2),'g','Facealpha',0.3);
patch(poly1(:,1),poly1(:,2),'b','Facealpha',0.3);

system('./clip')
readGPCfile('result.txt')

%fid = fopen('result.txt') ;
%data = textscan(fid,'%f %f','HeaderLines',2) ;
%data = cell2mat(data)
%fclose(fid);
%data=load('result.txt');

%if(size(data,1))
%patch(data(:,1),data(:,2),'r','Facealpha',0.3);
%end

function writeGPCfile(filename,P)
fid=fopen(filename,'w');
fprintf(fid,'%i\n',1) % number of polygons
fprintf(fid,'%i\n',size(P,1))  % number of points
for i=1:size(P,1)
    fprintf(fid,'%i %i\n',P(i,:))
end
end

function readGPCfile(filename,P)
fid=fopen(filename,'r');
nPoly = cell2mat(textscan(fid,'%d',1))
nPoly=nPoly(1);
skipLines=1;
for p=1:nPoly
    nV=cell2mat(textscan(fid,'%d',1,'HeaderLines',0));
nV=nV(1);
%skipLines=skipLines+1;
P=cell2mat(textscan(fid,'%f %f',nV,'HeaderLines',0))
%skipLines=skipLines+nV;
patch(P(:,1),P(:,2),'r','Facealpha',0.3);

end

end

function A=area(P)
P0=P(1,:);
A=0;
for k=1:size(P,1)
if k==size(P,1)
    k1=1;
else
    k1=k+1
end
s=P(k,:)-P0;
c=P(k1,:)-P(k,:);
A=A+0.5*(s(1)*c(2)-s(2)*c(1));
end
end