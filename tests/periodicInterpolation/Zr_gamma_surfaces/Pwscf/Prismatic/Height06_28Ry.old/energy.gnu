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

splot "energy.new" u 3:4:(($5-E06)*Ry*Jm2/surface) t "" w p,\
      "interpole.res" u 1:2:3 w l
pause -1

set title '<1120> direction  -  z = 0'
set ylabel "E  (J/m^2)"
plot "energy.new" u 3:($2==0?($5-E06)*Ry*Jm2/surface:1/0) t '6' w lp


set term wxt enh 2
set ylabel "stress  (kbar)"
plot "energy.new" u 3:($2==0?($6):1/0) t 's11' w lp,\
     "energy.new" u 3:($2==0?($7):1/0) t 's22' w lp,\
     "energy.new" u 3:($2==0?($8):1/0) t 's33' w lp,\
     "energy.new" u 3:($2==0?($9):1/0) t 's23' w lp,\
     "energy.new" u 3:($2==0?($10):1/0) t 's13' w lp,\
     "energy.new" u 3:($2==0?($11):1/0) t 's12' w lp

set term wxt enh 3
set title '<0001> direction  -  x = 0.5'
set ylabel "E  (J/m^2)"
plot "energy.new" u 4:($1==0.5?($5-E06)*Ry*Jm2/surface:1/0) t '6' w lp
