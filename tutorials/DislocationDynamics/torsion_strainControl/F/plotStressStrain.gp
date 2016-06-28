# =============================================================================
# This is a gnuplot input file for plotting the torque vs. angle
#
# Call this by: gnuplot plotStressStrain.gp
# At any time, press 'e' to refresh the plot (e.g. if data is changed).
#
# Stefan Sandfeld (stefan.sandfeld@fau.de), 2016
# =============================================================================

Burgers=0.2556e-9; # Burgers vector for Cu [m]
L1=2000 # side length of pillar
L2=2000 # side length of pillar
H=4000  # heigth of pillar
A=L1*L2;
V=H*A;

dotTheta3=1e-9;

# open window 1 and plot torsion angle vs. time
set term x11 0
set xlabel 'time'
set ylabel 'torsion angle'
plot 'F_0.txt' u 2:(dotTheta3*$2) t 'dottheta3*t' w l,\
     ''        u 2:6              t 'theta3' w l

# open window 2 and plot torque vs angle
set term x11 1
set xlabel 'twist angle [rad]'
set ylabel 'torque [\mu b^3]'
plot 'F_0.txt' u 6:7 t'' w l

pause -1

