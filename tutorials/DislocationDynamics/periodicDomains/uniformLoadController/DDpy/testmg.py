#!/usr/bin/python3

import ddpy

folderName = "./"
print("instantiating MicrostructureGeneratorInterface")
MG = ddpy.MicrostructureGeneratorInterface( folderName)
print("instantiated MicrostructureGeneratorInterface")

#print("instantiating DefectiveCrystalInterface")
#DC = ddpy.DefectiveCrystalInterface(folderName) # crashes?
#print("instantiated DefectiveCrystalInterface")
#print(f"currentStep: {DC.getCurrentStep()}")
#print(f"setCurrentStep(0)")
#DC.setCurrentStep(0)
#print(f"currentStep: {DC.getCurrentStep()}")
#print(f"getEndingStep: {DC.getEndingStep()}")
#print(f"setEndingStep(100)")
#DC.setEndingStep(100)
#print(f"getEndingStep()")
#print(f"running DC.runGlideSteps()")
#DC.runGlideSteps()
#print(f"finished DC.runGlideSteps()")
#print(f"printSlipSystemNormals:")
#DC.printSlipSystemNormals()
#print(f"printSlipSystemBurgersVectors:")
#DC.printSlipSystemBurgersVectors()
