close all
clear all
clf

n=7;

figure(1)
hold on
grid on
plot([0 1 0 0], [0 0 1 0],'k','Linewidth',2)

for i=1:n
plot([i/n i/n],[0 1-i/n])
plot([0 1-i/n],[i/n i/n])
plot([i/n 0],[0 i/n])
end

for i=0:n-1
    x=1/n*(i+1/3);

    for j=0:n-1-i
y=1/n*(j+1/3);
plot(x,y,'or')

    end
end

for i=0:n-2
    x1=1/n*(i+2/3);

    for j=0:n-2-i
y1=1/n*(j+2/3);
plot(x1,y1,'om')

    end
end
axis equal 
grid on

data=load('output.txt');
plot(data(:,1),data(:,2),'kx')