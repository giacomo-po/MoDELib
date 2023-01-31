#!/usr/bin/python3

import ddpy

folderName = "./"
#print("instantiating MicrostructureGeneratorInterface")
#MG = ddpy.MicrostructureGeneratorInterface( folderName)
#print("instantiated MicrostructureGeneratorInterface")

print("instantiating DefectiveCrystalInterface")
DC = ddpy.DefectiveCrystalInterface(folderName) # crashes?
print(f"getPlasticStrains:")
plasticStrains = DC.getPlasticStrains()
print(f"plasticStrains: {plasticStrains}")
print(f"getResolvedShearStrains:")
rsStrains = DC.getResolvedShearStrains()
print(f"rsStrains : {rsStrains}")
print(f"getResolvedShearStresses:")
rsStresses = DC.getResolvedShearStresses()
print(f"rsStresses : {rsStresses}")
