# sudo /usr/local/bin/python3 -m pip install PyQt5
# sudo /usr/local/bin/python3 -m pip install matplotlib
# sudo /usr/local/bin/python3 -m pip install numpy
import sys, string, os
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, '../../../../python')
from modlibUtils import *
plt.rcParams['text.usetex'] = True

# main code
materialFile='inputFiles/'+getStringInFile('inputFiles/polycrystal.txt','materialFile')
print('materialFile='+materialFile)
mu_SI=getValueInFile(materialFile,'mu0_SI')/1e6
print('mu='+ str(mu_SI) + ' MPa')
nu=getValueInFile(materialFile,'nu')
print('nu='+ str(nu))
E_SI=2*mu_SI*(1.0+nu)
print('E='+ str(E_SI) + ' MPa')

rho_SI=getValueInFile(materialFile,'rho_SI');    #[kg/m^3]
b_SI=getValueInFile(materialFile,'b_SI');  #[m]
v_dd2SI=np.sqrt(mu_SI/rho_SI);
t_dd2SI=b_SI/v_dd2SI;
T=getValueInFile('inputFiles/polycrystal.txt','absoluteTemperature')

F,Flabels=readFfile('./F')
runID=getFarray(F,Flabels,'runID')
e_33=getFarray(F,Flabels,'e_33')
s_33=getFarray(F,Flabels,'s_33 [mu]')*mu_SI
rhoG=getFarray(F,Flabels,'glissile density [m^-2]')
rhoS=getFarray(F,Flabels,'sessile density [m^-2]')

#ax11.plot(e_33, s_33)

minIdx = np.where(s_33 == np.amin(s_33))[0][0]
#idx=np.arange(minIdx,len(e_33),1,dtype=int)
e_33=e_33-e_33[minIdx]
s_33=s_33-s_33[minIdx]

hardIDs=e_33>0.0025
e_33_hard=e_33[hardIDs]
s_33_hard=s_33[hardIDs]
slope, intercept = np.polyfit(e_33, s_33, deg=1)


fig1 = plt.figure()
ax11=plt.subplot(2,1,1)
ax11.plot(e_33, s_33)
e_33_E=np.linspace(0,0.0005,2);
ax11.plot(e_33_E,E_SI*e_33_E)
ax11.plot(e_33_E+0.002,E_SI*e_33_E)
ax11.plot(e_33_hard,slope*e_33_hard+intercept)

ax11.grid()
ax11.set(xlabel='$\epsilon_{33} [-]$', ylabel='$\sigma_{33} [MPa]$')

ax12=plt.subplot(2,1,2)
ax12.plot(e_33, rhoG,label='glissile')
ax12.plot(e_33, rhoS,label='sessile')
ax12.grid()
ax12.legend()
ax12.set(xlabel='$\epsilon_{33} [-]$', ylabel='density $[m^{-2}]$')

#plt.show()
fig1.savefig("fig1.pdf", bbox_inches='tight')




