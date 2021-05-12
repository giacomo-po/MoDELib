set term wxt enh

set label "{/=28 {/Times-Bold (c)} EAM - prism plane}" at screen 0.51,0.95 left
coa=1.59802102659245516

# conversion ev.A^2 => J.m^-2
#Jm2=1.6021895e1

xPartial=.5000000000
yPartial=0.146*coa


set style line 1 lt 1 lc -1
set style line 2 lt 1 lc -1
set style line 3 lt 1 lc -1

set contour base 
set cntrparam cubicspline
#set cntrparam levels 20
unset clabel

set grid

set zlabel "{/=28 {/Symbol g}}  (J m^{-2})" offset 4,6

set xlabel "{/=28 a/3 [1~2{.6-}10]}" offset 5, -1.
set xrange [0.:1.]
set xtics 0.5 offset 0, -0.5
set mxtics 2

set ylabel "{/=28 c [0001]}" offset 5, 4
set yrange [0.:coa]
set ytics ("0" 0., "" 0.25*coa 1, "0.5" 0.5*coa, "" 0.75*coa 1, "1" coa) offset 1,-0.5
set size ratio 1.
set ticslevel 0.
zmin=-1.4 ; zmax=1.4
set zrange [zmin:zmax] ; set ztics ("0" 0, "" 0.25 1, "0.5" 0.5, "" 0.75 1, "1" 1., "" 1.25 1)



set hidden3d offset 0 front
set style line 100 lt 5 lw 0.5
set palette gray
set palette rgbformula 30,31,32
set palette rgbformula 21,22,23
set palette rgbformula 34,35,36

set arrow from xPartial,yPartial,zmin to 1,0,zmin head filled lw 2 lc 3 front
set arrow from 0,0,zmin to xPartial,yPartial,zmin head filled lw 2 lc 3 front

#view is 71 rot_x, 184 rot_z, 1 scale, 1 scale_z
set view 59, 31, 1, 1.4

unset colorbox
set cbrange [0.:1.4]
set cntrparam levels incremental 0.,0.075, 1.4


splot "interpole.res" u 1:2:($3) t '' w pm3d at bs ls 1,\
      "interpole.res" u 1:2:($3) every 5:5 t '' w l ls 1 nocontour
pause -1

set output 'interpole.eps'
set term post eps enh color "Times-Roman" 24
replot
set output
set term wxt enh

