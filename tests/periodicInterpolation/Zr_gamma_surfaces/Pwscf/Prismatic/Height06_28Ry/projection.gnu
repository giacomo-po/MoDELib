a=3.2300	# A
c=1.601*a	# A
nx=1
nz=1
surface=nx*nz*a*c
E06=-2379.6662	# Ry

# conversions Ry => eV
Ry=13.605805
# Conversion eV.A^2 -> J.m^2
Jm2=1.6021895e1

set samples 1000, 1000
set isosamples 1000, 1000
set xlabel "1/3 [1~2{.7-}10]"
set ylabel "[0001]"
set clabel "E  (J/m^2)"

set contour base
set cntrparam cubicspline 
set cntrparam levels 10
unset clabel

set pm3d explicit
set pm3d scansbackward
set pm3d at b interpolate 5,5

set view map
set xrange [0:1]
set yrange [0:c/a]
set size ratio c/a

set zrange [0:]

splot "energy.new" u 3:4:(($5-E06)*Ry*Jm2/surface) t '' w pm3d


set output 'projection.eps'
set term post eps enh "Times-Roman" 20
replot
set output
set term x11


