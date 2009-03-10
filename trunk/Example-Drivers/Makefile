# $Revision: 1.6 $, $Date: 2009/01/20 15:24:01 $
# $Name: Cornish-Gilliflower $
include ../Makefile
Modules += $(ModuleFlag)$(NetcdfHome)/include $(ModuleFlag)$(CodeDir) $(ModuleFlag)$(IntegDir)

rtObjs =  $(CodeDir)/characterUtils.o $(CodeDir)/numericUtilities.o \
          $(CodeDir)/ErrorMessages.o $(CodeDir)/userInterface_Unix.o \
          $(CodeDir)/$(multipleProcCode) $(CodeDir)/scatteringPhaseFunctions.o \
          $(CodeDir)/inversePhaseFunctions.o $(CodeDir)/opticalProperties.o \
          $(CodeDir)/RandomNumbersForMC.o $(CodeDir)/monteCarloIllumination.o \
          $(CodeDir)/surfaceProperties.o $(IntegDir)/monteCarloRadiativeTransfer.o

Targets = planeParallel monteCarloDriver

monteCarloDriver: monteCarloDriver.o $(rtObjs)
monteCarloDriver.o: monteCarloDriver.f95 $(rtObjs)

planeParallel: planeParallel.o $(rtObjs)
planeParallel.o: planeParallel.f95 $(rtObjs)

clean::
	-/bin/rm *.mod *.o $(Targets)
