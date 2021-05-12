a=3.230	# A
c=1.601*a	# A
nx=1
nz=1
surface=nx*nz*a*c
E06=-2379.6662

# conversions Ry => eV
Ry=13.605805
# Conversion eV.A^2 -> J.m^2
Jm2=1.6021895e1

factor=Ry/surface
show var

# Stacking fault energy at saddle point
g=(-2379.64998849-E06)*Ry/surface
gJm2=g*Jm2

factor=Ry*Jm2/surface
show var


set xlabel "x/a"
set ylabel "z/a"
set zlabel "E  (J/m^2)"

set hidden3d

set contour base
set cntrparam cubicspline
set cntrparam levels 10
unset clabel

splot "energy.new" u 3:4:(($5-E06)*Ry*Jm2/surface) t "" w lp
