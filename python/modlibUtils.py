import sys, string, os, fileinput
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
            if line[0:foundEqual].strip()==variable:
                return line[foundEqual+1:foundSemiCol].strip()
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

def setInputVariable(fileName,variable,newVal):
    with fileinput.FileInput(fileName, inplace=True) as file:
        for line in file:
            if variable in line:
                foundEqual=line.find('=');
                foundSemiCol=line.find(';');
#                if line[0:foundEqual-1].strip()==variable:
                oldVal=line[foundEqual+1:foundSemiCol]
                line = line.replace(oldVal,newVal)
            sys.stdout.write(line)

