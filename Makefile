# $Revision$, $Date$
# $URL$

### Instalation directory

I3RC_Home = /Users/robert/Codes/I3RC
CodeDir    = $(I3RC_Home)/Code
IntegDir   = $(I3RC_Home)/Integrators
 
### Netcdf-specific entries

NetcdfHome = /Users/robert/Codes/netcdf
Netcdf_IncludeDir = $(NetcdfHome)/include 
NetcdfLibs = -L$(NetcdfHome)/lib -lnetcdf

### General compiler macros 

ModuleFlag = -I
Modules     = $(ModuleFlag)$(Netcdf_IncludeDir)
Libs        = $(NetcdfLibs) 
Compile.F95 = $(F95) $(F95Flags) -c
Compile.F77 = $(F77) $(FFlags)   -c
Link.F95    = $(F95) $(F95Flags) 

### MPI-specific entries, if any, for use with a specific compiler
### If your environment supports it, you can specify the mpif90 wraper as the 
###    compiler and specify the compiler-specific flags
#
# Set UseMPI to "yes" and specify the appropriate values for the 
#   MPI_Dir, MPI_IncludeDir, and MPI_LibDir directories to enable 
#   the use of multiple processors
#
UseMPI=no
ifeq ($(UseMPI),yes)  
  MPI_Dir = NULL
  MPI_IncludeDir=$(MPI_Dir)/include
  MPI_LibDir=$(MPI_Dir)/lib
  Libs += $(MPI_LibDir) -lmpi 
  # MPI and non-MPI versions live in $(CodeDir)
  multipleProcCode = multipleProcesses_mpi.o
  F95Flags += $(ModuleFlag)$(MPI_IncludeDir)  
else 
  multipleProcCode = multipleProcesses_nompi.o
endif

#### Compiler-specific entries (may override marcos above) 

# Macros are available for ifort, g95, xlf, absoft, mpif90
compiler=mpif90
debug=no

ifeq ($(compiler),ifort)
  #
  # Intel fortran 
  #
  F95         = ifort
  F77         = ifort
  Compile.F95 += -free -Tf
  ifeq ($(debug),no) 
    # Optimization flags. 
    F95Flags = -O3 -heap-arrays -diag-disable vec -ipo $(Modules) 
    FFlags   = -fast 
  else
    # Debugging flags
    F95Flags = -g -O2 -heap-arrays -C -traceback -diag-disable vec $(Modules) 
    FFLAGS   = -g -O2 -heap-arrays -C -traceback -diag-disable vec 
  endif
endif

ifeq ($(compiler),g95)
  #
  # GNU fortran 
  #
  F95         = g95
  F77         = g77
  ifeq ($(debug),no) 
    # Optimization flags. 
    F95Flags = -O2 -std=f95 $(Modules) 
    FFLAGS   = -O2
  else
    # Debugging flags
    F95Flags = -g -std=f95 -fmodule-private -fimplicit-none -Wall $(Modules)
    FFLAGS   = -g                           -fimplicit-none -Wall 
  endif
endif

ifeq ($(compiler),xlf)
  #
  # IBM xlf 8 for Mac 
  #
  F95         = f95
  F77         = f77
  ifeq ($(debug),no) 
    # Optimization flags. 
    F95Flags = -O3 -qlanglvl=95std $(Modules) 
    FFLAGS   = -O3
  else
    # Debugging flags
    F95Flags = -g -O2 -C -qlanglvl=95std -qmaxmem=-1  $(Modules)
    FFLAGS   = -g                      
  endif
endif

ifeq ($(compiler),absoft)
  #
  # absoft fortran 
  #
  F95        = f95
  F77        = f77
  ModuleFlag = -p
  Libs      += -lU77
  ifeq ($(debug),no) 
    # Optimization flags. 
    F95Flags = -O3  -stack_size 0x10000000 -cpu:host $(Modules)  
    FFlags   = -O3  -stack_size 0x10000000 -cpu:host 
  else 
    ##Debugging flags
    F95Flags = -g $(Modules) 
    FFLAGS   = -g 
  endif
endif

ifeq ($(compiler),mpif90)
  #
  # mpif90 wrapper - compiler-specific flags need to be specified by hand  
  #
  F95         = mpif90
  F77         = mpif77
  Compile.F95 += -free -Tf
  multipleProcCode = multipleProcesses_mpi.o
  ifeq ($(debug),no) 
    # Optimization flags. 
    F95Flags = -O3 -heap-arrays -diag-disable vec -ipo $(Modules) 
    FFlags   = -fast 
  else
    # Debugging flags
    F95Flags = -g -O2 -heap-arrays -C -traceback -diag-disable vec $(Modules) 
    FFLAGS   = -g -O2 -heap-arrays -C -traceback -diag-disable vec 
  endif
endif

### General rules - should not require editing

# Rules for bullding object files
%.o: %.f
	$(Compile.F77) $<

%.o: %.f95
	$(Compile.F95) $<

# Rule to build executables	- some compilers are happier if everything 
#   is an object file 
%: %.o
	$(Link.F95) -o $@ $^ $(Libs)

