a=3.230		# A
nx=1
ny=1
surface=nx*ny*a**2*sqrt(3.)/2.
E10=-1983.0550

# conversions Ry => eV
Ry=13.605805
# conversion ev.A^2 => J.m^-2
Jm2=1.6021895e1

# Intrinsic stable stacking fault
gI=(-1983.04616703-E10)*Ry/surface	# eV/A^2
gI_Jm2=gI*Jm2				# J/m^2

show var

set style line 1 lt 1 lc -1
set style line 2 lt 1 lc -1
set style line 3 lt 1 lc -1

set hidden3d 
set contour base 
set cntrparam cubicspline
set cntrparam levels 10
unset clabel

set xlabel "x"
set ylabel "y"
set zlabel "E  (eV/A^2)"


splot "energy.res" u 3:4:(($5-E10)*Ry/surface) t "" w l
pause -1

