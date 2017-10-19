syms x1 x2 x3 A1 A2 A3 B1 B2 B3 u real

x=[x1 x2 x3]';
y=[A1 A2 A3]'+u*[B1-A1 B2-A2 B3-A3]';

R=x-y;
r=norm(R);

int(R/r^3,u,0,1)
