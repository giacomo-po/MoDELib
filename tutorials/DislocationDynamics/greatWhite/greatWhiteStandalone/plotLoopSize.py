# sudo /usr/local/bin/python3 -m pip install PyQt5
# sudo /usr/local/bin/python3 -m pip install matplotlib
# sudo /usr/local/bin/python3 -m pip install numpy
import sys, string, os
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, '../../../../python')
from readEVL import *

def getLoopRadius(folderName):
    F=np.loadtxt(folderName+'/F/F_0.txt');
    print(folderName+'/F/F_0.txt has size ' + str(np.shape(F)));
    R=np.empty(np.size(F[:,0]))
    n=0;
    for runID in F[:,0]:
        #print(int(runID))
        evl=readEVLtxt(folderName+'/evl/evl_'+str(int(runID)))
        c=np.mean(evl.nodes, axis=0)
        C=np.tile(c,[evl.nodes.shape[0],1])
        nodesC=evl.nodes-C;
        nodesR=np.sqrt(np.sum(np.square(nodesC),axis=1))
        R[n]=np.mean(nodesR);
        n=n+1;
    return F[:,1], R;

# main code
mu_SI=27e9;     #[Pa]
rho_SI=2700;    #[kg/m^3]
b_SI=0.286e-9;  #[m]
v_dd2SI=np.sqrt(mu_SI/rho_SI);
t_dd2SI=b_SI/v_dd2SI;

# experimental data
expData=np.array([[0.434, 329.651],
                 [15.029, 320.349],
                 [29.913, 290.698],
                 [44.942, 268.023],
                 [60.116, 249.419],
                 [90.029, 198.837],
                 [105.202, 168.023],
                 [120.087, 99.419]]);
ypos=np.array([329.651, 340.116, 310.465, 287.791, 268.023, 218.605,188.372,119.186]);
yneg=np.array([329.651, 299.419, 272.093, 248.837, 229.651, 179.651,148.837,77.907]);
print('experimental data has size '+ str(np.shape(expData)));

# simulation data data
t,R=getLoopRadius('.')

fig1 = plt.figure()
ax11=plt.subplot(1,1,1)
ax11.errorbar(expData[:,0], expData[:,1], yerr=[expData[:,1]-yneg, ypos-expData[:,1]], fmt='o', capthick=2,label='Silcox & Whelan, T=470$\pm$10K')
#ax11.plot(F[:,1]*t_dd2SI, R*b_SI*1e10,label='Great White')
ax11.plot(t*t_dd2SI, R*b_SI*1e10,label='Great White')
ax11.grid()
ax11.legend()
ax11.set_xlim((0, 180))
ax11.set_ylim((0, 400))
plt.xlabel('time [sec]')
plt.ylabel('loop radius [$\AA$]')
plt.show()
fig1.savefig("fig1.pdf", bbox_inches='tight')




