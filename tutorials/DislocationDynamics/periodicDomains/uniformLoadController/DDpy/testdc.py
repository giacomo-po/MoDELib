#!/usr/bin/python3

import ddpy

folderName = "./"
#print("instantiating MicrostructureGeneratorInterface")
#MG = ddpy.MicrostructureGeneratorInterface( folderName)
#print("instantiated MicrostructureGeneratorInterface")

print("instantiating DefectiveCrystalInterface")
DC = ddpy.DefectiveCrystalInterface(folderName) # crashes?
print("instantiated DefectiveCrystalInterface")
dcCurrentStep = DC.getCurrentStep()
print(f"currentStep: {dcCurrentStep}")
if dcCurrentStep > 100:
    print(f"setCurrentStep(100)")
    DC.setCurrentStep(100)
    dcCurrentStep = DC.getCurrentStep()
    print(f"currentStep: {dcCurrentStep}")
dcEndingStep = DC.getEndingStep()
print(f"getEndingStep: {dcEndingStep}")
if dcEndingStep <= 200:
    print(f"setEndingStep(200)")
    DC.setEndingStep( 200)
dcEndingStep = DC.getEndingStep()
print(f"getEndingStep(): {dcEndingStep}")
print(f"running DC.runGlideSteps()")
DC.runGlideSteps()
print(f"finished DC.runGlideSteps()")
dcCurrentStep = DC.getCurrentStep()
print(f"currentStep: {dcCurrentStep}")
print(f"printSlipSystemNormals:")
DC.printSlipSystemNormals()
print(f"printSlipSystemBurgersVectors:")
DC.printSlipSystemBurgersVectors()
