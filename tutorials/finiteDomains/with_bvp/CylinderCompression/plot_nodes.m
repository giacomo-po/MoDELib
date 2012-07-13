clear all;

data = load('dis_nodes.txt');

nn = size(data);
n = nn(1);

figure(1);
hold on;

for i = 1 : n
    plot3(data(i,1),data(i,2),data(i,3),'o');
end