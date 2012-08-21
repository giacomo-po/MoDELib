clear all;

data = load ('triangles_gp.txt');
gp = load ('gp_gp.txt');

nn = size(data);

n = nn(1);

figure(1);
hold on;

for i = 1 : n
    %X(1) = data(i,1);     X(2) = data(i,2);
    %Y(1) = data(i,4);     Y(2) = data(i,5);
    %Z(1) = data(i,7);     Z(2) = data(i,8);
    
    X(1) = data(i,1);     X(2) = data(i,4);
    Y(1) = data(i,2);     Y(2) = data(i,5);
    Z(1) = data(i,3);     Z(2) = data(i,6);
    
    plot3(X,Y,Z);
    
    
    %X(1) = data(i,1);     X(2) = data(i,3);
    %Y(1) = data(i,4);     Y(2) = data(i,6);
    %Z(1) = data(i,7);     Z(2) = data(i,9);
    
    X(1) = data(i,1);     X(2) = data(i,7);
    Y(1) = data(i,2);     Y(2) = data(i,8);
    Z(1) = data(i,3);     Z(2) = data(i,9);
    
    plot3(X,Y,Z);
    
    %X(1) = data(i,3);     X(2) = data(i,2);
    %Y(1) = data(i,6);     Y(2) = data(i,5);
    %Z(1) = data(i,9);     Z(2) = data(i,8);
    
    X(1) = data(i,7);     X(2) = data(i,4);
    Y(1) = data(i,8);     Y(2) = data(i,5);
    Z(1) = data(i,9);     Z(2) = data(i,6);
    
    plot3(X,Y,Z);
end

nn = size(gp);
n= nn(1);

for i = 1 : n
    plot3(gp(i,1),gp(i,2),gp(i,3),'rx');
end
