# sudo /usr/local/bin/python3 -m pip install PyQt5
# sudo /usr/local/bin/python3 -m pip install matplotlib
# sudo /usr/local/bin/python3 -m pip install numpy
import sys, string, os
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, '../../../../python')
from readEVL import *

def getStringInFile(fileName,variable):
    with open(fileName) as f:
        datafile = f.readlines()
    found = False  # This isn't really necessary
    for line in datafile:
#        print(line)
        if variable in line:
            # found = True # Not necessary
            foundEqual=line.find('=');
            foundSemiCol=line.find(';');
            return line[foundEqual+1:foundSemiCol]
    return 'Not Found'  # Because you finished the search without finding

def getValueInFile(fileName,variable):
    return float(getStringInFile(fileName,variable))

def readFfile(folder):
    F=np.loadtxt(folder +'/F_0.txt');
    with open('./F/F_labels.txt') as f:
        lines = f.readlines()
        for idx in range(len(lines)):
            lines[idx] = lines[idx].rstrip()
    return F,lines
    
def getFarray(F,Flabels,label):
    k=0;
    for line in Flabels:
        if line==label:
            return F[:,k]
        k=k+1
    return np.zeros(shape=(0,0))


# main code
materialFile='inputFiles/'+getStringInFile('inputFiles/polycrystal.txt','materialFile')
mu_SI=getValueInFile(materialFile,'mu0_SI')
print(mu_SI)
rho_SI=getValueInFile(materialFile,'rho_SI');    #[kg/m^3]
b_SI=getValueInFile(materialFile,'b_SI');  #[m]
v_dd2SI=np.sqrt(mu_SI/rho_SI);
t_dd2SI=b_SI/v_dd2SI;
#D0v=getValueInFile(materialFile,'D0v_SI')
#Ufv=getValueInFile(materialFile,'Ufv_eV')
#Umv=getValueInFile(materialFile,'Umv_eV')
#d=1500*b_SI;
#sigma=getValueInFile('./inputFiles/uniformExternalLoadController.txt','initialStress')*mu_SI;
#omega=1.0/np.sqrt(2.0)*b_SI**3;
#kb_eV=8.617e-5; # eV/K
#kb_SI=1.38e-23; # J/K
T=getValueInFile('inputFiles/polycrystal.txt','absoluteTemperature')

#eDotNH=D0v/d**2*np.exp(-(Ufv+Umv)/kb_eV/T)*(np.exp(sigma*omega/kb_SI/T)-1)
#eDotNH=D0v/d**2*np.exp(-(Ufv+Umv)/kb_eV/T)*(np.sinh(sigma*omega/kb_SI/T))

F,Flabels=readFfile('./F')
runID=getFarray(F,Flabels,'runID')
time=getFarray(F,Flabels,'time [b/cs]')*t_dd2SI
trbetaP=getFarray(F,Flabels,'tr(betaP)')
#s33=getFarray(F,Flabels,'s_33 [mu]')

fig1 = plt.figure()
ax11=plt.subplot(2,1,1)
ax11.plot(time, trbetaP,label='tr(betaP)')
#ax11.plot(time, eDotNH*time*100,label='NH')
ax11.grid()

#ax12=plt.subplot(2,1,2)
#ax12.plot(time, normBp*100,label='Great White')

#ax12.grid()

#ax11.grid()
#ax11.legend()
#ax11.set_xlim((0, 180))
#ax11.set_ylim((0, 400))
#plt.xlabel('time [sec]')
#plt.ylabel('loop radius [$\AA$]')
plt.xlabel('time [sec]')
plt.ylabel('tr(betaP) [cs/b]')
plt.show()
#fig1.savefig("fig1.pdf", bbox_inches='tight')




