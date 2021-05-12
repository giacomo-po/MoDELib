a=3.234	# A
nx=6
ny=7
surface=nx*ny*a**2*sqrt(3.)/2.
E0=-0.267511463945e+05

factor=1./surface

show var

set xlabel "x"
set ylabel "y"
set zlabel "E  (eV/A^2)"
splot "energy.res" u 3:4:(($5-E0)/surface) w l
pause -1
