a=3.234	# A
c=5.168	# A
coa = c/a
nx=6
nz=4
surface=nx*nz*a*c
E20=-0.127386411403e+05

# Conversion eV.A^2 => J.m^-2
Jm2=1.6021895e1

factor=Jm2/surface
show var

set xlabel "[1~2{.6-}10]" 
set ylabel "[0001]"
set zlabel "E  (J/m^2)" offset 2,8

set hidden3d offset -1

set contour base
set cntrparam cubicspline
set cntrparam levels 10
unset clabel

set xrange [0:1] 
unset xtics
set yrange [0:c/a] 
unset ytics
set size ratio c/a
set ticslevel 0.
zmin=-0.4 ; zmax=1.4
set zrange [zmin:zmax] ; set ztics 0.5 ; set mztics 2

xpartial=0.5
ypartial=0.146*coa

set arrow from 0,0,zmin to xpartial,ypartial,zmin head filled lt 4 lw 8
set arrow from xpartial,ypartial,zmin to 1., 0.,zmin head filled lt 4 lw 8

# view is 60 rot_x, 30 rot_z, 1 scale, 1 scale_z
set view 60, 30, 1.0, 1.4
splot "energy.res" u 3:4:(($5-E20)*Jm2/surface) t "" w l lt 2 lw 2
pause -1
show view


set output 'energy.eps'
set term post eps enh color solid "Times-Roman" 20
replot
set output
set term GNUTERM


# Gamma line in <1120> direction
set title '<1120> direction  -  z = 0'
set ylabel "E  (eV/A^2)"
plot "energy.res" u 3:($2==0?($5-E20)/surface:1/0) t '' w lp

# Gamma line in <0001> direction
set term wxt 3
E0=0.017
a2=-1.
a4=1.
a6=1.
E(x) = E0 + a2*x**2 + a4*x**4 + a6*x**6
fit [0.:0.45] E(x) "energy.res" u 4:($1==0.5?($5-E20)/surface:1/0) via E0, a2, a4, a6
set title '<0001> direction  -  x = 0.5'
set yrange [0:0.045]
plot "energy.res" u 4:($1==0.5?($5-E20)/surface:1/0) t '6' w lp,\
     E(x)

x0 = sqrt( ( -2.*a4 + sqrt(4.*a4**2-12.*a2*a6) )/( 6.*a6 ) )
x0_coa = x0/coa
Emin = E(x0)
show var
