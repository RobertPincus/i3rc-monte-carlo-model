# $Revision$, $Date$
# $URL$

include ../Makefile
Modules += $(ModuleFlag)$(CodeDir)

uses =  $(CodeDir)/characterUtils.o $(CodeDir)/numericUtilities.o \
        $(CodeDir)/ErrorMessages.o $(CodeDir)/userInterface_Unix.o \
        $(CodeDir)/RandomNumbersForMC.o \
        $(CodeDir)/scatteringPhaseFunctions.o $(CodeDir)/inversePhaseFunctions.o \
        $(CodeDir)/opticalProperties.o $(CodeDir)/monteCarloIllumination.o \
        $(CodeDir)/surfaceProperties.o $(CodeDir)/$(multipleProcCode)

monteCarloRadiativeTransfer.o: monteCarloRadiativeTransfer.f95 $(uses)

clean::
	-/bin/rm *.mod *.o
