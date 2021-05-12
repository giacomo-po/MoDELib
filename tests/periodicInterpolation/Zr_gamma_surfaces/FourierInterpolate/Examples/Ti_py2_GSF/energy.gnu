set term wxt enh

set label "{/=28 {/Times-bold (b)} (~2{.6-}112) plane}" at screen 0.31,0.90 left
#set label "{/=28 (~2{.6-}112) plane}" at screen 0.31,0.85 left

# Paramètres du titane
a = 2.936	# Å
c = 4.648
coa = c/a

# conversions Ry => eV
Ry=13.605805
# Conversion eV.A^2 -> J.m^2
Jm2=1.6021895e1


# Énergie de référence
E0 = -3266.29654634	# Ry

# Longueur des vecteurs de périodicité et surface 
lx = sqrt( 1. + coa**2 )
ly = sqrt(3)
S = a**2*lx*ly


set style line 1 lt 1 lc -1
set style line 2 lt 1 lc -1
set style line 3 lt 1 lc -1

set contour base 
set cntrparam cubicspline
unset colorbox
set cbrange [0.:2.5]
set cntrparam levels incremental 0.,0.200, 2.5
unset clabel

set grid

set xlabel "{/=28 1/3 [2~1{.6-}~1{.6-}3]}" offset 3, 4
set xrange [0:lx]
set xtics ("0" 0., "1/2" lx/2., "1" lx)
set ylabel "{/=28 1/2 [01~1{.6-}0]}" offset 6, -1
set ytics ("0" ly, "1/2" ly/2., "1" 0.) offset 0.,-0.5
set yrange [0:ly]
set zlabel "{/=28 {/Symbol g}}  (J m^{-2})" offset 4,7
zmin=-1 ; zmax=2.4
set zrange [zmin:zmax] ; set ztics 0., 1., zmax ; set mztics 2
set xyplane at zmin

set size ratio -1

set hidden3d offset 0 front
set style line 100 lt 5 lw 0.5
set palette rgbformula 30,31,32

set nokey

# Dissociation of the <c+a> dislocation
xPartial=0.45*lx ; yPartial=0.
set arrow from 0,0,zmin to xPartial,yPartial,zmin head filled lw 2 lc 7 front
set arrow from xPartial,yPartial,zmin to lx,0,zmin head filled lw 2 lc 7 front

#         view is 65 rot_x, 301 rot_z, 1 scale, 1.4 scale_z
set view 65, 301, 1, 1.4

splot "energy.interpole" u 1:2:(($3-E0)*Ry/S*Jm2) every 5:5 w l ls 1 nocontour,\
      "energy.interpole" u 1:2:(($3-E0)*Ry/S*Jm2) t '' w pm3d at bs nohidden3d ls 1

#pause -1
#show view



set output 'energy.eps'
set term post eps enh color solid "Times-Roman" 24
replot
set output
set term wxt


