syms h f0 f1 f2 real

M3=[1  0     0;
    1 -h     h^2;
    1 -2*h 4*h^2];

F=[f0 f1 f2]';

A=M3^-1*F;

I=simple([h^1/1 h^2/2 h^3/3]*A)