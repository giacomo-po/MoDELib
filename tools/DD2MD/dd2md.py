# conda deactivate
# cd build
# cmake ..
# make -j6
# cd ..
# /usr/local/bin/python3.10 dd2md.py
import sys
import matplotlib.pyplot as plt
import numpy as np
sys.path.append("./build")
import DD2MD
#from DD2MD import *

ad=DD2MD.AtomDisplacer("../../tutorials/DislocationDynamics/periodicDomains/uniformLoadController/")
#ad.readConfiguration(567)
ad.readMicrostructure()
#ad.writeConfiguration(0) # optional

#ad.solidAngle(100,100,100)
print(ad.dislocationPlasticDisplacement(100,100,100))

#n1=100;
#n2=100;
#y=500
#x=np.linspace(0, 1000, num=n1)
#z=np.linspace(0, 1000, num=n2)
#sa=np.empty([n1, n2])
#for i1 in range(0,x.size):
#    for i2 in range(0,z.size):
#        sa[i2,i1]=ad.solidAngle(x[i1],y,z[i2])
##plt.imshow(sa, cmap='hot', interpolation='nearest')
#plt.imshow(sa, origin='lower',cmap='jet')
#plt.colorbar()
#plt.show()

