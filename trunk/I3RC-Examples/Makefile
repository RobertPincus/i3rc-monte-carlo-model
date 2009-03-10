# $Revision: 1.12 $, $Date: 2006/07/18 22:21:22 $
# $Name: Cornish-Gilliflower $
# ------------------------------------------------------------------
# Programs for creating domain files for the I3RC test cases
# Run the radiative transfer using the driver in Example-Drivers

include ../make.common 
F95Flags = $(F95Flags_Default) 
Libs     = $(Libs_Default)

%.o: %.f95
	$(Compile.F95) $<

rtObjs = $(CodeDir)/characterUtils.o $(CodeDir)/numericUtilities.o \
         $(CodeDir)/ErrorMessages.o $(CodeDir)/userInterface_Unix.o \
         $(CodeDir)/inversePhaseFunctions.o $(CodeDir)/scatteringPhaseFunctions.o \
         $(CodeDir)/opticalProperties.o 

Targets = i3rcStepCloud i3rcRadarCloud i3rcLandsatCloud

all: $(Targets)

i3rcStepCloud: i3rcStepCloud.f95  $(rtObjs)

i3rcRadarCloud: i3rcRadarCloud.f95  $(rtObjs)

i3rcLandsatCloud: i3rcLandsatCloud.f95  $(rtObjs)

test: all
	./i3rcStepCloud
	./i3rcRadarCloud
	./i3rcLandsatCloud

clean::
	-/bin/rm *.mod *.o $(Targets)
