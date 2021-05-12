a=3.234	# A
c=5.168	# A
nx=6
nz=4
surface=nx*nz*a*c
E20=-0.127386411403e+05


set samples 1000, 1000
set isosamples 1000, 1000
set xlabel "x"
set ylabel "z"
set zlabel "E  (eV/A^2)"

set contour base
set cntrparam cubicspline 
set cntrparam levels 10
unset clabel

set pm3d explicit
set pm3d scansbackward
set pm3d at b interpolate 5,5

set view map

splot "energy.res" u 3:4:(($5-E20)/surface) t '' w pm3d


set output 'projection.eps'
set term post eps enh "Times-Roman" 20
replot
set output
set term x11


