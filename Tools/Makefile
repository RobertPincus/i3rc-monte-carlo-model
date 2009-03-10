# $Revision: 1.6 $, $Date: 2009/01/20 15:24:01 $
# $Name: Cornish-Gilliflower $
include ../Makefile
Modules += $(ModuleFlag)$(CodeDir)

mieObjs =  $(CodeDir)/characterUtils.o $(CodeDir)/numericUtilities.o \
           $(CodeDir)/ErrorMessages.o $(CodeDir)/userInterface_Unix.o \
           $(CodeDir)/scatteringPhaseFunctions.o $(CodeDir)/$(multipleProcCode)
           
optObjs = $(mieObjs) $(CodeDir)/inversePhaseFunctions.o $(CodeDir)/opticalProperties.o

all: MakeMieTable PhysicalPropertiesToDomain OpticalPropertiesToDomain

# Object files
#
MakeMieTable.o: MakeMieTable.f95  $(mieObjs)

PhysicalPropertiesToDomain.o: PhysicalPropertiesToDomain.f95 $(optObjs)

OpticalPropertiesToDomain.o: OpticalPropertiesToDomain.f95 $(optObjs)

# Executables
#
MakeMieTable: MakeMieTable.o mieindsub.o RefractiveIndex-IceAndWater.o $(mieObjs)

PhysicalPropertiesToDomain: PhysicalPropertiesToDomain.o $(optObjs)

OpticalPropertiesToDomain: OpticalPropertiesToDomain.o $(optObjs)

clean::
	-/bin/rm *.mod *.o MakeMieTable PhysicalPropertiesToDomain OpticalPropertiesToDomain
