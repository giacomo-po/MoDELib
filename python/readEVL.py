import sys, string, os
import numpy as np
# evl file

class EVL:
    nodes=np.empty([0,0])

def readEVLtxt(filename):
    evlFile = open(filename+'.txt', "r")
    numNodes=int(evlFile.readline().rstrip())
    numLoops=int(evlFile.readline().rstrip())
    numLinks=int(evlFile.readline().rstrip())
    evl=EVL();
    evl.nodes=np.empty([numNodes, 3])
    for k in np.arange(numNodes):
        data=np.fromstring(evlFile.readline().rstrip(), sep=' ')
        evl.nodes[k,:]=data[1:4]
    return evl



