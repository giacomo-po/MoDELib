
c2=1; % the shear wave speed
fc2=0.7*c2; % the fraction of c2 that dislocation velocity can reach

np=100;
x=[0:3*np]/np;  % x= |Fpk|/B

v=fc2*(1-exp(-x/fc2));

figure(1)
plot(x,v)

figure(2)
plot(x(1:end-1),diff(v)./diff(x))

