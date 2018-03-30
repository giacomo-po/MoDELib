clear all
close all
clc

N=200;
%P=rand(N,2);
[x, y, dt] = simple_polygon(N);
%return

%% write input points
fid=fopen('I/I_0.txt','w')
for n=1:N
    fprintf(fid,'%i %e %e\n',[n x(n) y(n)])
end
fclose(fid)

%% Call code
system('./test')

%% Read output
O=load('O/O_0.txt');

%% plot results
figure(1)
clf
hold on
for(k=1:N)
    k1=k+1;
    if k==N
        k1=1;
    end
    plot(x([k k1]),y([k k1]),'k','Linewidth',2)
end

for(t=1:size(O,1))
    ids=O(t,:)+1; % one-based
    patch(x(ids),y(ids),'g','FaceAlpha',0.5)
end
axis equal