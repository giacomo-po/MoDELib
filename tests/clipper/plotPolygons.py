import sys, string, os
import matplotlib.pyplot as plt
import numpy as np

def getPolyFromFile(fileName):
    P=np.loadtxt(fileName);
    print(fileName + ' has size ' + str(np.shape(P)));
    return P;

# main code
subject=getPolyFromFile('subject.txt')
clipper=getPolyFromFile('clipper.txt')
result=getPolyFromFile('result.txt')

fig1 = plt.figure()
ax11=plt.subplot(1,1,1)
#ax11.plot(F[:,1]*t_dd2SI, R*b_SI*1e10,label='Great White')
ax11.plot(subject[:,0], subject[:,1],'b',label='subject')
ax11.plot(clipper[:,0], clipper[:,1],'r',label='clipper')
ax11.plot(result[:,0], result[:,1],'g',label='result')

ax11.grid()
ax11.legend()
plt.show()




