import sys, string, os
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, '../../python')
from modlibUtils import *


curDir=os.getcwd()


os.system('rm -r F/*.*')
os.system('rm -r evl/*.*')
microstrucureFile=getStringInFile('inputFiles/initialMicrostructure.txt','microstructureFile')
microstrucureFileName=curDir+'/inputFiles/'+microstrucureFile

if microstrucureFile=='periodicDipole.txt':
    figFile='periodicDipole.pdf'
    setInputVariable(microstrucureFileName,'periodicDipoleExitFaceIDs',str(1)) # 1=edge, 0=screw
    X=np.linspace(0, 20000, num=390)
    for x in X:
        setInputVariable(microstrucureFileName,'periodicDipoleGlideSteps',str(x))
        os.system('../../tools/microstructureGenerator/build/microstructureGenerator '+ curDir)
        os.system('../../tools/DDomp/build/DDomp '+ curDir)
elif microstrucureFile=='periodicLoops.txt':
    figFile='periodicLoops.pdf'
    X=np.linspace(0, 25, num=49)
    for x in X:
        setInputVariable(microstrucureFileName,'periodicLoopRadii_SI',str(x*1e-8))
        os.system('../../tools/microstructureGenerator/build/microstructureGenerator '+ curDir)
        os.system('../../tools/DDomp/build/DDomp '+ curDir)
else:
    print('unsupported option')
        
F,Flabels=readFfile('./F')
E=getFarray(F,Flabels,'elastic energy [mu b^3]')

fig1 = plt.figure()
ax1=plt.subplot(1,1,1)
ax1.plot(X, E,label='Ee',marker ='o')
plt.show()
fig1.savefig(figFile, bbox_inches='tight')
