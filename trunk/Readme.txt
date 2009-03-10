Quick start guide to the I3RC Community Monte Carlo Radiative Transfer Model
$Revision: 1.11 $, $Date: 2009/03/09 19:17:21 $
$Name: Cornish-Gilliflower $

This directory contains the I3RC community Monte Carlo model for solving 
problems in solar radiative transfer in three-dimensionally variable 
atmospheres. 

*** Changes since the last release (Bramley, July 2006) 
*) Released under the terms of the GNU public license. NASA has a 
separate license with similar rights. 
*) Support for surface bidirectional reflectance functions via module 
SurfaceProperties and changes to monteCarloRadiativeTransfer. An example 
module that computes Lambertian reflectance is included. 
*) Optional use of message passing via MPI in a parallel-capable  
monteCarloDriver program.
*) New variance reduction methods (see namelists for monteCarloDriver). 
*) Revised plane-parallel solver - now shares most namelists with generic 
solver, and can save problem domain to a file. 
*) Templates/ subdirectory had been removed. 

*** Credit where credit is due

If you publish results obtained with this code, please refer to it as the I3RC 
Community Monte Carlo model and reference the I3RC overview paper by Bob 
Cahalan (http://dx.doi.org://10.1175/BAMS-86-9-1275) and the paper by 
Robert Pincus and Frank Evans (submitted to J. Atmos. Sci in March 2009; 
current status available from http://www.psd.esrl.noaa.gov/people/robert.pincus/Papers/3DRT-Models/). 

*** Quick start: An example

Before you can solve a radiative transfer problem you have to define the 
(three-dimensionally varying) properties of the atmosphere. These are often 
specified in terms of the concentration and sizes of cloud and aerosol problems, 
for which the optical properties must first be computed. That is, one must 
*) compute the single scattering properties of cloud and/or aerosol particles, 
probably as a function of particle size, at appropriate wavelength, then
*) describe the three-dimensional distribution of particles within the domain, then
*) compute the radiative transfer. 

We've provided programs that correspond to each step: 
*) Tools/MakeMieTable creates a table of single scattering properties at a given 
wavelength for a size distribution of spheres as a function of size
*) Tools/PhysicalPropertieToDomain reads several kinds of formatted ASCII files
describing concentration, drop sizes or numbers, etc., combines them with the
phaseFunctionTables, and produces a file describing the domain
  (Tools/OpticalPropertiesToDomain can be used if the optical properties, rather 
than the physical properties, are available). 
*) Example-Drivers/monteCarloDriver reads the domain, computes the radiative
transfer, and writes out the results. monteCarloDriver-mpi is almost identical 
but divides the computation across mulitple CPUs for speed. 

If what you need are fluxes or intensities at the domain boundaries or heating 
rates within the domain, you can almost certainly use our programs to solve 
your problem and won't need to program anything yourself. 

In the language of the Programmer's Guide, the three steps correspond to the 
creation of three objects: 
*) a phaseFunctionTable (from module scatteringPhaseFunctions), 
*) a domain (from module opticalProperties), and 
*) an integrator (from module monteCarloRadiativeTransfer) which processes a set
of photons (from module monteCarloIllumination).
Problems can also be solved by creating phaseFunctionTables and domains using
calls from Fortran code, saving them using the write functions included in the
modules, then using our driver programs to compute the radiative transfer.

*** What's in the package

The subdirectories contain the code framework for the I3RC community Monte Carlo
radiative transfer model (Code/); a forward Monte Carlo radiative transfer solver
(Integrators/); tools to build tables of phase functions using Mie theory and to
convert ASCII files to the binary versions used by the code framework (Tools/);
and two general purpose drivers that use the code to solve radiative transfer
problems (Example-Drivers). Programs to set up the three I3RC "Phase 1" cases
(I3RC-Examples/) are included as programming examples.

The example radiative transfer solver computes fluxes at the top and bottom of
the domain and absorption through the domain (in units of flux per volume). It
will also compute intensity if supplied with a set of directions, though this is
of course much slower. See the example programs and their example namelists for
information on how to specify the problem you want to solve, and the User's
Guide for more detailed information.

*** Building the code 

This code relies on the netCDF libraries, which are available at
http://www.unidata.ucar.edu/software/netcdf/. These must be built before
compiling the I3RC code. In my experience, it's easiest to build the Fortran
interface to netCDF using the same compiler with which you'll build the I3RC
code. You can build the I3RC code without netCDF if you remove or replace all
the read_ and write_ subroutines from the modules in the Code/ directory. If you
simply remove the subroutines you won't be be able to store information for
later use, so I don't recommend this if you can possibly avoid it.

This code conforms to the ANSI Fortran 95 standard but it really stresses
compilers, and many fail while building the code. Some will simply dump core,
often while compiling the file scatteringPhaseFunctions.f95 in the Code/
subdirectory. If this happens be sure you have the most up-to-date version of
the compiler available. Some compilers (those from Portland Group, and at least
some versions of the Sun Forte compilers) simply don't work on this code. Most
of the development work has been done on Power PC Macs running xlf 8.1 on system
10.3 and 10.4, and on Intel Linux boxes running the Intel Fortran compiler
(versions 8, 9, and 10). 
 
The code builds under the g95 compiler, which is available for a wide range of
systems (Windows, Solaris, HP-UX, Mac) from http://g95.org/. Unfortunately, as
of this writing the g95 compiler produces code that runs very slowly (4-6 times
more slowly than other compilers) because it spends a lot of time managing
temporary memory. This will hopefully improve in future.
 
Compilation options, including the compiler name and compilation flags, are set in 
the Makefile in the top level directory.This is also where the location of the 
(required) netcdf files is set, as well as the (optional) location of the 
MPI message-passing libraries and include files.  This information is
used by the Makefiles in each subdirectory. We've provided production
and debugging settings for the Intel ifort and g95 compilers. If you're
using a new platform add the compiler and flag definitions following
these examples. I'd appreciate copies of working configurations for other platforms.

In the examples we have provided most of the parameters are specified using
namelists. The name of the namelist file must be supplied at run time. Many Unix
systems support reading arguments from the command line; on platforms that don't
the file name is read from standard in. You can choose which behavior you want
by commenting out the approriate subroutines in Code/userInterface_Unix.f95
before compiling.

Once make.common has been edited, type "./Build" in the top level directory. This
will build everything in all the subdirectories in the proper order. Note that
the directories must be built in order (Code/, Integrators/, Example-Drivers/;
Tools/ must be built after Code/) because the dependencies are set in the directories 
themselves.

The subdirectory Tools/ contains general purpose programs to build "domain"
files that describe the 3D distribution of optical properties within a domain.
Typically one would first build one or more phase function tables describing the
single scattering properties of clouds, aerols, etc. using MakeMieTable. We've
provided example namelists for clouds and aerosol, though you may want to tweak
these. Program PhysicalPropertiesToDomain reads an ASCII description of the 3D
liquid water content field, combines it with a phase function table, and creates
the final domain. Alternatively, you can use the program
OpticalPropertiesToDomain to convert files similar to those used by SHDOM to the
internal format. ASCII file formats are described in the example programs. The
domain files can then used as input to the program
Example-Drivers/monteCarloDriver. 

The programs in I3RC-Examples/ will build files describing the I3RC Phase 1 test
cases using the files provided by the I3RC (available in I3RC-Examples/Data).

This code makes extensive use of dynamic memory. If program seems to crash 
without explanation be sure the shell is not imposing limits on the amount 
of memory each process can request (i.e. "unlimit stacksize" when using tcsh). 

Please see the I3RC Community Model Programmer's Guide for more information. 

Best - 

Robert Pincus, University of Colorado/NOAA Earth System Research Lab
Robert.Pincus@colorado.edu
