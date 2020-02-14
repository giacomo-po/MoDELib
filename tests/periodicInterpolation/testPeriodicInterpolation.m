clear all
close all
clc


system('./test')

data=load('output.txt')

[X,Y] = meshgrid([0:0.01:1],[0:0.01:sqrt(3)]);

f=zeros(size(X));
for i=1:size(data,1)
k=data(i,[1 2]);
S=data(i,3);
C=data(i,4);
    f=f+S*sin(k(1)*X+k(2)*Y)+C*cos(k(1)*X+k(2)*Y);
end

figure
clf
hold on
surf(X,Y,f)
grid on
xlabel('x')
ylabel('y')
