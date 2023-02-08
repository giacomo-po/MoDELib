#!/bin/python3

# Contact : hyunsol@g.clemson.edu

import numpy as np
import re

def GeneratePolyCrystalFile(simPath,boxSizes,flag):
    #Math variables
    I = np.identity(3)

    #Preset variables
    #modelibDirPath     = '/home/anon/Documents/code_local/MoDELib2/MoDELib/'
    polycrystalDataName = 'polycrystal.txt'
    alignToSlipSystem0 = 1
    # temp = 300
    lattice            = 'fcc'
    # material           = 'AlMg15'
    gridSpacing_SI     = np.mat([1, 1])*1e-10 # in meter

    # if flag = 1 is passed, pull the temperature value from the previous polycrystal.txt
    if flag == 1:
        with open(simPath+'inputFiles/polycrystal.txt', 'r') as o:
            for line in o:
                tempvaldummy = re.match("absoluteTemperature.*",line)
                #print(tempvaldummy) # for debugging
                matdummy = re.match("materialFile.*",line)
                if not tempvaldummy == None:
                    tempval = tempvaldummy.group(0)
                if not matdummy == None:
                    matType = matdummy.group(0)
        tempval = tempval.split(';')
        matType = matType.split('/')
        tempvalinit = []
        matTypeInit = []
        for i in tempval:
            tempvalinit.append(i.split())
        # string filtering, find the material type string in polycrystal.txt
        for j in matType:
            if re.match("AlMg.*",j):
                matTypeInit = j.split()
        for k in matTypeInit:
            matTypeInit = k.split('.')
            if re.match("AlMg.*",k):
                matTypeInit = matTypeInit[0]
                break

        tempval = [ int(new_string) for new_string in tempvalinit[0] if new_string.isdigit() ]
        temp = tempval[0]
        material = matTypeInit

    # if flag = 0, you declare temperature and material manually.
    elif flag == 0:
        temp      = 300
        material  = 'AlMg5'

    if alignToSlipSystem0:
        match lattice:
            case 'bcc':
                s = [1, 1, -1]
                n = [1, 0, 1]
            case 'fcc':
                s = [0, 1, 1]
                n = [-1, 1, -1]
            case _:
                exit("Crystal structure is not valid")
        s = s/np.linalg.norm(s)
        n = n/np.linalg.norm(n)
        C2G1 = np.mat( [s,np.cross(n,s),n] )
    else:
        C2G1 = I

    match lattice:
        case 'bcc':
            A = np.mat([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]) * 1/np.sqrt(3) 
        case 'fcc':
            A = np.mat([[0, 1, 1], [1, 0, 1], [1, 1, 0]]) * 1/np.sqrt(2)
        case _:
            exit("Crystal structure is not valid")

    A = np.matmul(C2G1, A)

    g = 0

    Fs = np.diag([1, 1, 10])*boxSizes
    #Fs = np.diag([100, 100, 1000])

    F12 = [[1, g, 0],[0, 1, 0],[0, 0, 1]]
    F31 = [[1, 0, g],[0, 1, 0],[0, 0, 1]]
    F23 = [[1, 0, 0],[0, 1, g],[0, 0, 1]]

    F = np.matmul(F12, F23)
    F = np.matmul(F, F31)

    Ba = np.matmul(F, Fs)

    Ma = np.matmul(np.linalg.inv(A), Ba)
    M = np.round(Ma)
    B = np.matmul(A, M)

    x0 = np.mat([0, 0, 0])

    #f = open(simPath+'inputFiles/polycrystal3.txt','w')
    f = open(simPath+'inputFiles/'+polycrystalDataName,'w')
    f.write('materialFile=../../../MaterialsLibrary/'+ material +'.txt; \n')
    f.write('enablePartials=0; \n')
    f.write('absoluteTemperature = ' + str(temp) + '; # [K] simulation temperature \n')
    f.write('meshFile=../../../MeshLibrary/unitCube.msh; \n')
    f.write('C2G1=')
    C2G1string = np.array2string(C2G1, formatter={'float_kind':lambda C2G1: '%.15f' % C2G1})
    f.write( re.sub('\[|\]','',C2G1string) )
    f.write(';\n');
    f.write('\n');
    f.write('\n');

    f.write('A=')
    Bstring = np.array2string(B, formatter={'float_kind':lambda B: '%.15f' % B})
    f.write( re.sub('\[|\]','',Bstring) )
    f.write(';\n');
    f.write('\n');
    f.write('\n');

    f.write('x0=');
    x0string = np.array2string(x0, formatter={'float_kind':lambda x0: '%.15f' % x0})
    f.write( re.sub('\[|\]','',x0string) )
    f.write(';\n');

    f.write('periodicFaceIDs=');
    f.write('0 1 2 3 4 5');
    f.write(';\n');

    f.write('\n#################\n');
    f.write('# GlidePlaneNoise\n');
    f.write('gridSize=');
    f.write('256 256 ');
    f.write('; # size of grid on the glide plane\n');
    f.write('gridSpacing_SI=');
    gridSpacing_SIstring = np.array2string(gridSpacing_SI, formatter={'float_kind':lambda gridSpacing_SI: '%.15f' % gridSpacing_SI})
    f.write( re.sub('\[|\]','',gridSpacing_SIstring) )
    f.write('; # [m] spacing of grid on the glide plane\n');
    f.write('solidSolutionNoiseMode=1; # 0=no noise, 1= read noise, 2=compute noise\n');
    f.write('solidSolutionNoiseFile_xz=../../../NoiseLibrary/noise_xz.vtk;\n');
    f.write('solidSolutionNoiseFile_yz=../../../NoiseLibrary/noise_yz.vtk;\n');
    f.write('stackingFaultNoiseMode=0; # 0=no noise\n');
    f.write('spreadLstress_A=1.0;# [AA] spreading length for stresses\n');
    f.write('a_cai_A=3.0; # [AA] spreading length for non-singular theory\n');
    f.write('seed=1234;\n');
    f.close()

if __name__ == '__main__':

    print(" ")
    print(" ")
    print("*** WARNING ***")
    print("*** If you are running this script directly, you need to change the directory path ***")
    print(" ")
    print(" ")

    # evl and F directory path
    simPath = "/mnt/696e72d4-caf1-4ed3-a5f1-7122fe930356/Document/code_local/MoDELib2/MoDELib/tutorials/DislocationDynamics/periodicDomains/uniformLoadController/";
    boxSizes = 100;
    
    # flag = 1, read values from the polycrystal.txt
    # flag = 0, use manually assigned values to generate polycrystal.txt
    flag = 1;
    GeneratePolyCrystalFile(simPath,boxSizes,flag)
