#!/bin/python3

import numpy as np
import re, sys, os
from dataclasses import dataclass

@dataclass
class PolyCrystalFile:
    args : list = None   # user argument
    polyOutPath : str = None
    modelibPath : str = None    
    dataDict : dict = None
    polycrystalDataName : str = 'polycrystal.txt' # polycrystal file default name

    def generate(self, lattice: str=None, partialMode: int=0, temp: int=0, material: str='Al', boxSizeX: int=100, boxSizeY: int=100, boxSizeZ: int=1000, seed: int=1234, SFNoiseMode: int=0, SSNoiseMode: int=1, gridSize_x: int=256, gridSize_y: int=256):
        # If you are running this script independently, take user arguments for the system condition
        if __name__ == '__main__':
                polyOutPath = self.polyOutPath
                partialMode = int(sys.argv[1])
                boxSizeX = int(sys.argv[2])
                boxSizeY = int(sys.argv[3])
                boxSizeZ = int(sys.argv[4])
                gridSize_x = int(sys.argv[5])
                gridSize_y = int(sys.argv[6])
                material = sys.argv[7]
                temp = sys.argv[8]
                lattice = sys.argv[9]
                seed = sys.argv[10]
                SSNoiseMode = sys.argv[11]
                SFNoiseMode = sys.argv[12]
                modelibPath = sys.argv[13] + '/'
                matLibPath = modelibPath + 'tutorials/DislocationDynamics/MaterialsLibrary' + '/'
                inputFilePath = modelibPath + 'tutorials/DislocationDynamics/periodicDomains/uniformLoadController/inputFiles' + '/'
                print(f"[ System Info ]")
                print(f"PartialEnabled = {partialMode} \nMaterial = {material} \nTemperature = {temp} \nCrystal Type = {lattice} \nseed = {seed} \nSolidSolutionNoise={SSNoiseMode} \nStackingFaultNoise={SFNoiseMode} \nMoDELibdir={modelibPath}")
        # If you are using the code as a library, it jumps to this line to set the parameters
        else:
            matLibPath = self.modelibPath + 'tutorials/DislocationDynamics/MaterialsLibrary' + '/'
            inputFilePath = self.modelibPath + 'tutorials/DislocationDynamics/periodicDomains/uniformLoadController/inputFiles' + '/'
            polyOutPath = self.polyOutPath+'inputFiles/'
            if not os.path.exists(polyOutPath):
                os.makedirs(polyOutPath)
            for key, value in self.dataDict.items():
                match key:
                    case 'T':
                        temp = str(value)
                    case 'C':
                        material = 'AlMg'+str(value)
                    case 'DL':
                        boxSizeX = value[0]
                        boxSizeY = value[1]
                        boxSizeZ = value[2]
                    case 'S':
                        seed = str(value)
                    case 'SF':
                        SFNoiseMode = str(value)
                    case 'SS':
                        SSNoiseMode = str(value)
                    case 'GS':
                        gridSize_x = value[0]
                        gridSize_y = value[1]
            if not material and not temp and not lattice and not seed:  
                raise ValueError("The system is ill-defiend")

        #Math variables
        I = np.identity(3)

        #Preset variables
        alignToSlipSystem0 = 1

        gridSizes = np.mat([gridSize_x, gridSize_y])
        burgersVector = float(self.extNumber(matLibPath, material, 'b_SI')[0])
        boxSize2D = boxSizeX * boxSizeY
        gridSpacing_SI = np.mat([boxSizeX, boxSizeY]) *burgersVector / gridSizes

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

        #Fs = np.diag([1, 1, 10])*boxSizes
        Fs = np.diag([boxSizeX, boxSizeY, boxSizeZ])
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
        

        f = open(polyOutPath + self.polycrystalDataName, 'w')
        f.write('materialFile=../../../MaterialsLibrary/'+ material +'.txt; \n')
        f.write(f'enablePartials={partialMode}; \n')
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
        # Chnage this part later
        #gridSizes = str(gridSizes).strip('[]')
        f.write(f'{str(gridSizes).strip("[]")}');
        f.write('; # size of grid on the glide plane\n');
        f.write('gridSpacing_SI=');
        gridSpacing_SIstring = np.array2string(gridSpacing_SI, formatter={'float_kind':lambda gridSpacing_SI: '%.2e' % gridSpacing_SI})
        f.write( re.sub('\[|\]','',gridSpacing_SIstring) )
        f.write('; # [m] spacing of grid on the glide plane\n');
        f.write(f'solidSolutionNoiseMode={SSNoiseMode}; # 0=no noise, 1= read noise, 2=compute noise\n');
        f.write('solidSolutionNoiseFile_xz=../../../NoiseLibrary/noise_xz.vtk;\n');
        f.write('solidSolutionNoiseFile_yz=../../../NoiseLibrary/noise_yz.vtk;\n');
        f.write(f'stackingFaultNoiseMode={SFNoiseMode}; # 0=no noise\n');
        f.write('spreadLstress_A=1.0;# [AA] spreading length for stresses\n');
        f.write('a_cai_A=3.0; # [AA] spreading length for non-singular theory\n');
        f.write('seed='+str(seed)+';\n');
        f.close()

    # extract numbers from a pattern matched line in a text file
    def extNumber(self, dirPath, txtFileName, pattern):
        with open(dirPath+txtFileName+'.txt', 'r') as o:
            for line in o:
                txtMatchdummy = re.match(pattern+".*", line)
                if not txtMatchdummy == None:
                    txtMatch = txtMatchdummy.group(0)
        txtMatch = re.findall(r'[-+]?\d*\.\d+[eE][-+]?\d+|\d+', txtMatch)  # Match all scientific notation or regular number
        # returns a list
        return txtMatch;


if __name__ == '__main__':
    if len(sys.argv) != 14 :
        print("***************** Error: expected 10 arguments *********************")
        print("Usage: ./GeneratePolyCrystalFile.py partialMode Boxsize_x Boxsize_y Boxsize_z gridSize_x gridSize_y Material Temperature Lattice Seed SolidSolutionNoiseMode StackingFaultNoiseMode ${MoDELibPath}")
        print("")
        print("Example: ./GeneratePolyCrystalFile.py 1 100 100 1000 256 256 AlMg5 0 fcc 1234 2 0 ~/{MoDelibPath}")
        print("Example: ./GeneratePolyCrystalFile.py 1 300 300 1000 512 512 W 1000 bcc 1234 2 0 ~/{MoDelibPath}")
        print("*******************************************************************")
    else:
        # Generate polycrystal.txt file in the working directory
        polyOutPath='./'

        # create an object and generate txt file
        Crystal = PolyCrystalFile(sys.argv, polyOutPath=polyOutPath)
        Crystal.generate()

