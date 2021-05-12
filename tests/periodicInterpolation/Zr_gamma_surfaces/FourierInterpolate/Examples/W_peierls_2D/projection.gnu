# Vecteurs de périodicité
u1x = sqrt(6.)/3. ; u1y = 0.
u2x = sqrt(6.)/6. ; u2y = sqrt(2.)/2.

# Origine
Ox = 0. ; Oy = 0.

# Positions easy et hard
xE = sqrt(6.)/6. ; yE =  sqrt(2.)/6.
xH = sqrt(6.)/3. ; yH =  sqrt(2.)/3.

n=2 ; m=2
set xrange [Ox:Ox+n*u1x+m*u2x] ; set noxtics
set yrange [Oy:Oy+n*u1y+m*u2y] ; set noytics
set size ratio -1

zmax=250
# Bords de boîte
set arrow from Ox,Oy,zmax to Ox+u1x,Oy+u1y,zmax                 nohead lt -1
set arrow from Ox,Oy,zmax to Ox+u2x,Oy+u2y,zmax                 nohead lt -1
set arrow from Ox+u1x,Oy+u1y,zmax to Ox+u1x+u2x,Oy+u1y+u2y,zmax nohead lt -1
set arrow from Ox+u2x,Oy+u2y,zmax to Ox+u1x+u2x,Oy+u1y+u2y,zmax nohead lt -1
# Bord de triangle
set arrow from Ox+u1x,Oy+u1y,zmax to Ox+u2x,Oy+u2y,zmax         nohead lt -1

set label "E" at xE,yE,zmax point pt 1 ps 3 offset 2
set label "H" at xH,yH,zmax point pt 1 ps 3 offset 2

set key bottom

set contour base
set cntrparam cubicspline 
set cntrparam levels 10
unset clabel
set style line 1 lt 1 lc -1

set view map

set para
splot "w_all.interpole" u 1:2:3 t "" w pm3d at s ls 1,\
     "w_peierls.symmetrized" u 1:2:3 t 'peierls' w p nocontour,\
     "w_hard-split.symmetrized" u 1:2:3 t 'hard-split' w p nocontour,\
      xE,yE,0. t '' w p pt 1 ps 3,\
      xH,yH,0. t '' w p pt 1 ps 3

show var
