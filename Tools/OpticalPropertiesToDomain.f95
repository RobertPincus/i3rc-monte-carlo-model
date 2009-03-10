! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

program opticalPropertiesToDomain
 ! $Revision$, $Date$
 ! $URL$
 
 ! This program converts optical properties from an SHDOM-like ascii file
 ! to an I3RC Monte Carlo optical property domain object, and writes this
 ! object to a file.  The ascii input file of optical properties is almost
 ! the same as an SHDOM tabulated phase function file, but with nz+1 height
 ! levels, where nz is the number of grid cells in the vertical.  This
 ! difference is due to SHDOM representing optical properties at grid points,
 ! while the I3RC Monte Carlo code represents properties in uniform grid cells.
 ! 
 ! The input file format has a 4+NUMPHASE header lines, where NUMPHASE is the
 ! number of phase functions: the first line must begin with 'T'; the second
 ! has the number of grid cells in each dimension; the third has the grid
 ! spacing and Z levels; the fourth has the number of phase functions; and the
 ! following "lines" have the phase functions containing the degree of the
 ! Legendre expansion and then the Legendre coefficients.   There is then one
 ! line for every grid cell containing the grid cell X,Y,Z indices, the
 ! temperature, the extinction, single scattering albedo, and phase function
 ! index (1 to NUMPHASE) to the table in the header.  The format is as follows:
 !     T
 !     Nx  Ny  Nz
 !     delX  delY  Z(1) ... Z(Nz+1)
 !     numphase
 !     NumL  Chi1 ...  ChiL       (for each phase function)
 !       . . .
 !     IX IY IZ  Temp Extinct Albedo  Iphase
 !       . . .
 ! where IX, IY, and IZ are the grid cell spatial indices (1,...,Nx, etc.),
 ! Temp is the grid cell temperature in Kelvin, Extinct is the extinction
 ! in units inverse of the grid spacing units, Albedo is the single
 ! scattering albedo, and Iphase is the phase function table index.
 ! NumL is the Legendre series degree, and Chi(l) are the coefficients of the
 ! Legendre expansion of the phase function.  The Legendre coefficients are
 ! normalized so that Chi0 (which is not listed) is always 1. The Legendre
 ! coefficients are defined so that Chi1 is 3*g, where g is the asymmetry
 ! parameter [e.g. a Henyey-Greenstein phase function is Chi(l)=(2*l+1)g^l].
 ! NOTE that this is not the convention used by the I3RC code; this program 
 ! converts the moments that are read in by 2*l + 1.
 ! The Legendre coefficients do not actually need to be all on one line, but
 ! can continue across multiple lines.  The first Z level is the bottom of the
 ! lowest grid cell and the last is the top of the highest grid cell.
 ! 
 !    Frank Evans    University of Colorado     December 2005


   ! Modules from the I3RC community Monte Carlo model framework
  use ErrorMessages
  use userInterface
  use scatteringPhaseFunctions
  use opticalProperties
  
  implicit none
  
   ! Input parameters
  character(len=256)   :: propFileName     = "", &
                          outputFileName   = ""
  namelist /fileNames/ PropFileName, outputFileName
  
  
   ! Local variables
  character(len = 256) ::  namelistFileName = ""
  integer              :: nX, nY, nZ
  integer              :: nPhaseEntries, iostat, i, j, k, l, n
  real                 :: deltaX, deltaY, tempTemp, x, y, z, chi
  real,    allocatable :: zPosition(:)
  integer              :: maxLegCoefs
  integer, allocatable :: nLegCoef(:)
  real,    allocatable :: legendreCoefficients(:,:)
  real,    allocatable :: extinct(:,:,:), ssa(:,:,:)
  integer, allocatable :: phaseFuncIndex(:,:,:)
  character(len=1)     :: junk

   ! I3RC Monte Carlo code derived type variables
  type(phaseFunction), allocatable :: ppVector(:)
  type(phaseFunctionTable)         :: phaseFuncTable
  type(domain)                     :: thisDomain
  type(ErrorMessage)               :: status


  ! -----------------------------------------
  ! Get the input variables from the namelist file
  !
  ! NOTE: The namelist contains only two arguments - 
  !   the names of the input and output files. On systems 
  !   that support reading arguments from the command line 
  !   via getarg it would make sense to get both args that 
  !   way instead of clumsily, though the namelist. 
  !
  namelistFileName = GetOneArgument() 
  open (unit=1, file=trim(namelistFileName), status='old')
  read (1,nml = fileNames)
  close (1)
  if(len_trim(propFileName) == 0 .or. len_trim(outputFileName) == 0) &
    stop "Must provide both input and output file names"

  ! -----------------------------------------
  !  Read the SHDOM style optical properties file
  !
   ! Open the property file and read the header
  open (3, file=trim(PropFileName), status='old')
  read (3, '(A1)') junk
  if (trim(junk) /= "T") then
    print *, "This doesn't look like an tabulated phase function property file."
    stop
  endif
  read (3,*) nX, nY, nZ
  allocate (zPosition(nZ+1))
  read (3,*) deltaX, deltaY, zPosition(1:nZ+1)
  read (3,*) nPhaseEntries

   ! Read in the phase functions
    ! First read the number of Legendre terms for all phase functions
  allocate (nLegCoef(nPhaseEntries))
  do i = 1, nPhaseEntries
    read (3,*) nLegCoef(i), (chi, l=1,nLegCoef(i))
  end do
  close(3)
  maxLegCoefs = maxval(nLegCoef(:))
  
    ! Then allocate the array and read the Legendre coefficients
    !   (Reopen the file and read it again from the beginning - 
    !    the values we have already read won't have changed.) 
  allocate (legendreCoefficients(maxLegCoefs,nPhaseEntries))
  open (3, file=trim(PropFileName), status='old')
  read (3, '(A1)') junk
  read (3,*) nX, nY, nZ
  read (3,*) deltaX, deltaY, zPosition(1:nZ+1)
  read (3,*) nPhaseEntries

   ! Read in the phase functions
    ! First read the number of Legendre terms for all phase functions
  do i = 1, nPhaseEntries
    read (3,*) nLegCoef(i), legendreCoefficients(1:nLegCoef(i),i)
  end do
  print *, "Read ", nPhaseEntries, " phase functions"

   ! Read in the properties at each grid point
  allocate (extinct(nX,nY,nZ), ssa(nX,nY,nZ), phaseFuncIndex(nX,nY,nZ))
  extinct(:,:,:) = 0.0  ;  ssa(:,:,:) = 0.0  ; phaseFuncIndex(:,:,:) = 1
  readProps: do
    read (3, *, iostat=iostat) i, j, k, tempTemp, &
          extinct(i,j,k), ssa(i,j,k), phaseFuncIndex(i,j,k)
    if (phaseFuncIndex(i,j,k) < 1 .or. phaseFuncIndex(i,j,k) > nPhaseEntries) then
      print *, 'Phase function index out of range: ',i,j,k,phaseFuncIndex(i,j,k)
      stop
    endif
    if (iostat /= 0) exit readProps
  end do readProps
  close (3)

  ! -----------------------------------------
  !  Do the calls to the I3RC Monte Carlo routines to 
  !    create and write out the domain. 
  !
   ! Create the domain
  thisDomain = new_Domain (deltaX * (/ (i, i = 0, nX) /), &
                           deltaY * (/ (i, i = 0, nY) /), &
                           zPosition, status)
  call printStatus(status)

   ! Make all the phase function objects (converting definition of coeffs)
  allocate (ppVector(nPhaseEntries))
  do i = 1, nPhaseEntries
    n = nLegCoef(i)
    legendreCoefficients(1:n,i) = legendreCoefficients(1:n,i) / (/ (2*l+1, l = 1, n) /)
    ppVector(i) = new_PhaseFunction( legendreCoefficients(1:n,i), status=status)
    call printStatus(status)
  end do
  deallocate (nLegCoef, legendreCoefficients)

   ! Create the phase function table
  phaseFuncTable = new_PhaseFunctionTable (ppVector(1:nPhaseEntries), &
                               key = real( (/ (i, i = 1, nPhaseEntries) /) ), &
                               status = status)
  call printStatus(status)

   ! Add the optical properties to the domain
  call addOpticalComponent (thisDomain, "mixture", &
                            extinct(:,:,:), ssa(:,:,:), phaseFuncIndex(:,:,:), &
                            phaseFuncTable, status = status)
  call printStatus(status)
  deallocate (extinct, ssa, phaseFuncIndex)
  call finalize_PhaseFunctionTable (phaseFuncTable)

  call write_Domain(thisDomain, trim(outputFileName), status)
  call printStatus(status)
  
  call finalize_Domain(thisDomain)
end program opticalPropertiesToDomain
