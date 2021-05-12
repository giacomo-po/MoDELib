set term wxt enh

set label "{/=28 {/Times-bold (c)} EAM - basal plane}" at screen 0.01,0.87 left

# conversion ev.A^2 => J.m^-2
Jm2=1.6021895e1

xPartial=.5000000001
yPartial=-.288675134565945


set style line 1 lt 1 lc -1
set style line 2 lt 1 lc -1
set style line 3 lt 1 lc -1

set contour base 
set cntrparam cubicspline
#set cntrparam levels 20
unset clabel

set grid

set xlabel "{/=28 a/3 [1~2{.6-}10]}"
set ylabel "{/=28 a/2 [10~1{.6-}0]}" offset -6.5,3.7 
set zlabel "{/=28 {/Symbol g}}  (J m^{-2})" offset 8,7.5

set xrange [-0.5:1]
set xtics ("0" 1, "1/2" 0.5, "1" 0., "3/2" -0.5)
set yrange [-sqrt(3.)/2.:0.]
set ytics ("0" 0., "1/3" -sqrt(3.)/6., "2/3" -sqrt(3.)/3.,  "1" -sqrt(3.)/2.)
set size ratio sqrt(2.)/2./1.5
set ticslevel 0.
zmin=-0.4 ; zmax=0.8
set zrange [zmin:zmax] ; set ztics ("0" 0, "0.25" 0.25, "0.5" 0.5, "0.75" 0.75)



set hidden3d offset 0 front
set style line 100 lt 5 lw 0.5
set palette gray
#set palette rgbformula 30,31,32	# black-blue-violet-yellow-white
#set palette rgbformula 21,22,23	# black-red-yellow-white
set palette rgbformula 34,35,36	# black-red-yellow-white
#set palette rgbformula 3,11,6	# green-red-violet
#set palette rgbformula 23,28,3	# green-blue-white
#set palette rgbformula 33,13,10	# blue-green-yellow-red


set arrow from 1,0,zmin to 0.5,-sqrt(3.)/2.,zmin nohead front
set arrow from 0,0,zmin to -0.5,-sqrt(3.)/2.,zmin nohead front

set arrow from 1,0,zmin to xPartial,yPartial,zmin head filled lw 2 lc 3 front
set arrow from xPartial,yPartial,zmin to 0,0,zmin head filled lw 2 lc 3 front
set arrow from 1,0,zmin to -0.5,-sqrt(3.)/2.,zmin nohead lw 2 lc 3 lt 2 front

#view is 71 rot_x, 184 rot_z, 1 scale, 1 scale_z
set view 71, 184, 1.4, 1

unset colorbox
set cbrange [0:0.75]
set cntrparam levels incremental 0., 0.050, 0.75

splot "interpole.res" u 1:2:($3*Jm2) t '' w pm3d at bs ls 1,\
      "interpole.res" u 1:2:($3*Jm2) every 5:5 t '' w l ls 1 nocontour
pause -1

set output 'interpole.eps'
set term post eps enh color "Times-Roman" 24
replot
set output
set term x11

