
set xrange [0.:sqrt(2.)]
set xtics ("S" 0., "H" sqrt(2.)/3., "E" 2.*sqrt(2.)/3., "S" sqrt(2.))
set grid x

set para
plot "w_all.interpole" u 2:3 t '' w l,\
     "../w_hard-split.dat" u ((1.-$1)*sqrt(2.)/3.):4 t '' w p pt 1 lc 1,\
     2.*sqrt(2.)/3.,0. t '' w p pt 1 lc 1,\
     t*sqrt(2.),0. t '' w l lt 1 lc -1
