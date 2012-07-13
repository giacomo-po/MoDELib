clear all;

data = load('triSegInter.txt');

x1 = data(:,2);
y1 = data(:,3);
z1 = data(:,4);

x2 = data(:,5);
y2 = data(:,6);
z2 = data(:,7);

nn = size(x1);

n = nn(1);

figure(1);
hold on

for i = 1 : n
    X = [x1(i) x2(i)];
    Y = [y1(i) y2(i)];
    Z = [z1(i) z2(i)];
    
    plot3(X,Y,Z);
end

