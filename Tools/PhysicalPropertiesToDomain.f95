! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

program ParticleFileToDomain
 ! $Revision: 1.10 $, $Date: 2009/03/09 19:17:22 $
 ! $Name: Cornish-Gilliflower $
 !
 ! Reads an ASCII file describing the three-dimensional distribution
 !   of clouds and/or aerosols and writes the description to a
 !   "domain" object from the I3RC community model. Profiles of 
 !   molecular absoprtion can be added and Rayleigh scattering can 
 !   be included. 
 ! 
 ! Input parameters are specified with a namelist.
 !
 ! The input file of particle properties (mass content and effective
 ! radius) may be specified in three different ways:
 !   1) One parameter LWC file.  Header:
 !        1                [format type]
 !        nX nY nZ         [number of X, Y, Z grid cells]
 !        deltaX deltaY    [X and Y grid cell size in km]
 !        Zlevels(1:nZ+1)  [increasing heights of cell boundaries in km]
 !        Temps(1:nZ+1)    [temperatures of boundaries (K)]
 !      One line per grid cell with:  iX iY iZ  LWC 
 !         iX,iY,iZ are indices from 1 to nX,nY,nZ, and
 !         LWC is cloud liquid water content (g/m^3).
 !         The effective radius is obtained from LWC using
 !           Reff = 100* (LWC*0.75*1.3889/(3.14159*DropNumConc))**(1/3)
 !   2) Two parameter LWC file.  Same header as 1 parameter LWC file,
 !      but with "2" starting the first line. 
 !      One line per grid cell with:  iX iY iZ  LWC Reff
 !         Reff is the effective radius (micron).
 !   3) Multicomponent particle properties file. Header:
 !        3                [format type]  
 !        nX nY nZ         [number of X, Y, Z grid cells]
 !        deltaX deltaY    [X and Y grid spacing in km]   
 !        Zlevels(1:nZ+1)  [increasing heights of cell boundaries in km]
 !        Temps(1:nZ+1)    [temperatures of boundaries (K)]
 !      One line per grid cell with:  
 !        iX iY iZ  Numcomp Type1 Mass1 Reff1 ... TypeN MassN ReffN
 !      Numcomp is the number of particle components in this cell, 
 !      Type? are the type numbers of the components, 
 !      Mass? are the mass contents [g/m^3] of the components, and 
 !      Reff? are the effective radii [microns] of the components.  
 !      The type number refers to the scattering table number 
 !      (in the range 1 to number of scattering tables input).
 !
 ! The scattering properties for each component are specified in
 ! files containing I3RC Monte Carlo phase function table objects
 ! (e.g. written by MakeMieTable).
 ! 
 ! An input file of molecular absorption extinction profile may be input.
 ! The Zlevels must be the same as the profile made by combining the
 ! levels in the particle file with the other levels specified in the
 ! namelist file.  Format (three lines): 
 !        nZ               [number of Z grid cells]
 !        Zlevels(1:nZ+1)  [increasing heights of cell boundaries in km]
 !        GasExt(1:nZ)     [molecular extinction in km^-1 for each cell]
 !
 ! Molecular Rayleigh scattering may be included (if RayleighWavelength > 0).
 ! The specified temperature profile is used to calculate the pressure
 ! profile with the hypsometric equation.  The Rayleigh extinction
 ! profile, which is proportional to air density, is calculated from
 ! the temperature and pressure profiles and the wavelength.  The average 
 ! extinction in each layer is calculated, assuming an exponential decay
 ! in air density.
 !

 !    Frank Evans    University of Colorado     December 2005
 !      Modified (for organization) by Robert Pincus, Climate Diagnostics Center


   ! Modules from the I3RC community Monte Carlo model framework
  use ErrorMessages
  use CharacterUtils
  use numericUtilities
  use scatteringPhaseFunctions
  use opticalProperties
  use UserInterface
  
  implicit none

   ! Input parameters
  character(len=256)   :: NamelistFileName = ""
  character(len=256)   :: ParticleFileName = ""
  integer              :: numScatTables
  integer, parameter   :: maxNumComps=5
  character(len=256)   :: ScatTableFiles(maxNumComps) = ""
  character(len=256)   :: MolecAbsFileName = ""
  
  integer, parameter   :: maxOtherLevels = 20
  integer              :: numOtherLevels = 0
  real                 :: OtherHeights(maxOtherLevels) = 0, &
                          OtherTemps(maxOtherLevels) = 0
                          
  real                 :: RayleighWavelength = 0
  real                 :: DropNumConc = 0.
  character(len=256)   :: outputFileName = ""
  
  character(len=512)   :: tableDescription
  
  namelist /fileNames/ ParticleFileName, ScatTableFiles, MolecAbsFileName, &
                       outputFileName
  namelist /profile/  OtherHeights, OtherTemps
  namelist /physicalProperties/  DropNumConc, RayleighWavelength
  
   ! Local variables
  integer              :: nX, nY, nZp, nZt
  integer              :: i, j, k, n, iscat, ix, iy, iz, il
  integer              :: maxNretab, izLevelBase
  real                 :: deltaX, deltaY, x, y, z, f
  real, allocatable    :: Zpar(:), TempPar(:)
  integer, allocatable :: nComp(:,:,:), ptype(:,:,:,:)
  real, allocatable    :: MassCont(:,:,:,:), Reff(:,:,:,:)
  integer, allocatable :: Nretab(:)
  real, allocatable    :: ReffTable(:,:), ExtinctTable(:,:), ssaTable(:,:)
  real, allocatable    :: Zlevels(:), Temp(:)
  real, allocatable    :: GasExt(:), RaylExt(:)
  real, allocatable    :: ssaProf(:)
  integer, allocatable :: phaseIndex(:)
  real                 :: LegendreCoefs(2)
  real,    allocatable :: extinct(:,:,:,:), ssa(:,:,:,:)
  integer, allocatable :: phaseFuncIndex(:,:,:,:)

   ! I3RC Monte Carlo code derived type variables
  type(phaseFunction)                   :: phaseFuncObject
  type(phaseFunction)                   :: PhaseFuncs(2)
  type(phaseFunctionTable)              :: OnePhaseFuncTable
  type(phaseFunctionTable), allocatable :: phaseFuncTables(:)
  type(domain)                          :: thisDomain
  type(ErrorMessage)                    :: status


  ! -----------------------------------------
  ! Get the input variables from the namelist file
  !
  
  namelistFileName = getOneArgument()
  open (unit=1, file=trim(namelistFileName), status='OLD')
  read (1, nml = fileNames)
  read (1, nml = profile)
  read (1, nml = physicalProperties)
  close (1)
  
  !
  ! Check input arguments as much as possible 
  !   Some checks can't be done until you know what kind of file you're reading
  !
  if(len_trim(ParticleFileName) == 0)   stop "Must specify particle file name." 
  if(len_trim(outputFileName) == 0) stop "Must specify output file name." 
  if(all(len_trim(ScatTableFiles) == 0)) &
    stop "Must specify as many scattering tables as there are components"

  if(DropNumConc < 0) stop "DropNumConc must be positive" 
  if(RayleighWavelength < 0) stop "RayleighWavelength must be non-negative." 
  
  
  numOtherLevels = count(otherTemps(:) > 0) 
  if(numOtherLevels == maxOtherLevels) &
    print *, "Read only the first ", maxOtherLevels, " other levels." 
  numScatTables = count(len_trim(ScatTableFiles) > 0)
  if(numScatTables == maxNumComps) &
    print *, "Read only the first ", maxNumComps, " scattering tables." 
  if(len_trim(outputFileName) == 0) stop "Must specify output file." 

  ! Read in the particle properties file
  call read_particle_file_size (ParticleFileName, nX, nY, nZp)
  allocate (Zpar(nZp+1), TempPar(nZp+1))
  allocate (nComp(nX,nY,nZp), ptype(numScatTables,nX,nY,nZp))
  allocate (MassCont(numScatTables,nX,nY,nZp), Reff(numScatTables,nX,nY,nZp))
  call read_particle_file (ParticleFileName, nX, nY, nZp, numScatTables, &
                           DropNumConc,  deltaX, deltaY, Zpar, TempPar, & 
                           nComp, ptype, MassCont, Reff)

   ! Combine the particle levels and extra levels
  nZt = nZp + numOtherLevels
  allocate (Zlevels(nZt+1), Temp(nZt+1))
  call organize_levels (nZp, Zpar, TempPar, &
                        numOtherLevels,     &
                        OtherHeights(:numOtherLevels), OtherTemps(:numOtherLevels), &
                        nZt, Zlevels, Temp, izLevelBase)

   ! Read the molecular absorption extinction file if there is one
  allocate (GasExt(nZt))
  GasExt(:) = 0
  if(len_trim(MolecAbsFileName) > 0) &
    call read_molec_abs_file (MolecAbsFileName, nZt, Zlevels, GasExt)

   ! Calculate the molecular Rayleigh scattering extinction profile
  allocate (RaylExt(nZt))
  RaylExt(:) = 0
  if(RayleighWavelength > 0.) &
    call rayleigh_extinct (nzt, Zlevels, Temp, RayleighWavelength, RaylExt)


  ! -----------------------------------------
  !  Read in the scattering tables
  !
  allocate (phaseFuncTables(numScatTables), Nretab(numScatTables))
   ! Read the scattering tables and get the number of entries in each table
  do i = 1, numScatTables
    call read_PhaseFunctionTable(fileName = ScatTableFiles(i), &
                                 table = phaseFuncTables(i),   &
                                 status = status) 
    call printStatus(status)
    call getInfo_PhaseFunctionTable (phaseFuncTables(i), nEntries=Nretab(i), &
                                     tabledescription = tableDescription,    &
                                     status=status)
    call printStatus(status)
  enddo

   ! Get the effective radius, extinction, and single scattering albedo for
   !   all the entries in each scattering table
  maxNretab = maxval(Nretab(:))
  allocate (ReffTable(maxNretab,numScatTables))
  allocate (ExtinctTable(maxNretab,numScatTables), ssaTable(maxNretab,numScatTables))
  do i = 1, numScatTables
    call getInfo_PhaseFunctionTable (phaseFuncTables(i), &
                                     key                    = ReffTable(1:Nretab(i),i), &
                                     extinction             = ExtinctTable(1:Nretab(i),i), &
                                     singleScatteringAlbedo = ssaTable(1:Nretab(i),i), &
                                     status=status)
    call printStatus(status)
  enddo



  ! -----------------------------------------
  !  Use the effective radius for each particle type in each grid cell 
  !  to index into the scattering tables to calculate the extinction, 
  !  single scattering albedo, and phase function index fields.
  !  
  allocate (extinct(nX,nY,nZp,numScatTables))
  allocate (ssa(nX,nY,nZp,numScatTables))
  allocate (phaseFuncIndex(nX,nY,nZp,numScatTables))

  extinct(:,:,:,:) = 0.0
  ssa(:,:,:,:) = 0.0
  phaseFuncIndex(:,:,:,:) = 1
  do iz = 1, nZp
   do iy = 1, nY
    do ix = 1, nX
      do k = 1, nComp(ix,iy,iz)
        ! Get the scattering table number (iscat) for this particle type
        iscat = ptype(k,ix,iy,iz)
        if(Reff(k,ix,iy,iz) <= maxval(ReffTable(1:NReTab(iscat),iscat)) .and. &
           Reff(k,ix,iy,iz) >  minval(ReffTable(1:NReTab(iscat),iscat))) then 
           
          ! Binary search to find effective radius entry in table
          il = findIndex(Reff(k,ix,iy,iz), ReffTable(1:NReTab(iscat),iscat))

           ! Interpolate optical properties linearly in Reff for this component
          f = (Reff(k,ix,iy,iz)-ReffTable(il,iscat)) / (ReffTable(il+1,iscat)-ReffTable(il,iscat))
          extinct(ix,iy,iz,iscat) = MassCont(k,ix,iy,iz) * &  
                                    ((1-f)*ExtinctTable(il,iscat) + f*ExtinctTable(il+1,iscat))
          ssa(ix,iy,iz,iscat) = (1-f)*ssaTable(il,iscat) + f*ssaTable(il+1,iscat)
          ! Chose the closest phase function
          if (f < 0.5) then
            phaseFuncIndex(ix,iy,iz,iscat) = il
          else
            phaseFuncIndex(ix,iy,iz,iscat) = il+1
          endif
        
        else if(MassCont(k,ix,iy,iz) > 0.0) then 
          ! This effective radius isn't in the table
          print *, 'Warning: effective radius outside of table (ix,iy,iz,type,Reff):'
          print '(4(1x,i3),1x,f6.2)', ix, iy, iz, ptype(k,ix,iy,iz), Reff(k,ix,iy,iz)
          extinct(ix,iy,iz,iscat) = 0.
          ssa(ix,iy,iz,iscat) = 0.
        end if
      enddo
    enddo
   enddo
  enddo
  

  deallocate (MassCont, Reff, ptype, Ncomp)
  deallocate (Nretab, ReffTable, ExtinctTable, ssaTable)

  ! -----------------------------------------
  ! Package the optical properties in a domain object
  !
  !
  ! Create the domain
  !
  thisDomain = new_Domain (deltaX * (/ (i, i = 0, nX) /), &
                           deltaY * (/ (i, i = 0, nY) /), &
                           Zlevels(1:nZt+1), status)
  call printStatus(status)


   ! Add the optical properties for each particle component to the domain
  do i = 1, numScatTables
    print *, "Adding component ", i
    call addOpticalComponent (thisDomain, "Particle type " // trim(IntToChar(i)), &
                              extinct(:,:,:,i), &
                              ssa(:,:,:,i), phaseFuncIndex(:,:,:,i), &
                              phaseFuncTables(i), zLevelBase=izLevelBase, &
                              status = status)

    call printStatus(status)
    call finalize_PhaseFunctionTable (phaseFuncTables(i))
  end do
  deallocate (extinct, ssa, phaseFuncIndex, phaseFuncTables)

   ! Add the Rayleigh scattering optical properties
  allocate (ssaProf(nZt), phaseIndex(nZt))
  if (any(RaylExt(:) > 0.0)) then
    print *, "Adding Rayleigh scattering"
    ssaProf(:) = 1.0
    phaseIndex(:) = 1
    LegendreCoefs(1:2) = (/ 0.0, 0.5 /) / (/ 2.*1. + 1., 2.*2. + 1. /)
    PhaseFuncs(1) = new_PhaseFunction (LegendreCoefs(1:2), status=status)
    call printStatus(status)
    OnePhaseFuncTable = new_PhaseFunctionTable (PhaseFuncs(1:1), key=(/ 0.0 /),&
                                                tableDescription = "Rayleigh scattering", &
                                                status=status)
    call printStatus(status)
    call finalize_phaseFunction (phaseFuncObject)
    call addOpticalComponent (thisDomain, 'Rayleigh scattering', &
                              RaylExt(:), ssaProf(:), phaseIndex(:), &
                              OnePhaseFuncTable, status = status)
    call printStatus(status)
    call finalize_PhaseFunctionTable (OnePhaseFuncTable)
  endif

   ! Add the molecular absorption optical properties
  if (any(GasExt(:) > 0.0)) then
    print *, "Adding molecular absorption" 
    ssaProf(:) = 0.0
    phaseIndex(:) = 1
    LegendreCoefs(1:1) = 0.0
    PhaseFuncs(2) = new_PhaseFunction (LegendreCoefs(1:1), status=status)
    call printStatus(status)
    OnePhaseFuncTable = new_PhaseFunctionTable (PhaseFuncs(2:2), key=(/ 0.0 /),&
                                                tableDescription = "Molecular absorption", &
                                                status=status)
    call printStatus(status)
    call finalize_phaseFunction (phaseFuncObject)
    call addOpticalComponent (thisDomain, 'Molecular absorption', &
                              GasExt(:), ssaProf(:), phaseIndex(:), &
                              OnePhaseFuncTable, status = status)
    call printStatus(status)
    call finalize_PhaseFunctionTable (OnePhaseFuncTable)
  endif
  deallocate (ssaProf, phaseIndex)
  deallocate (RaylExt, GasExt, Zpar, TempPar)
  
  !
  ! Write the domain to a file
  !
  call write_Domain(thisDomain, outputFileName, status)
  call printStatus(status)
  call finalize_Domain(thisDomain)
end program ParticleFileToDomain
! --------------------------------------------------------------------------

subroutine read_particle_file_size (parfile, nx, ny, nzp)
  implicit none
  character(len=*), intent(in) :: parfile
  integer, intent(out) :: nx, ny, nzp
  
  open (unit=2, file=trim(parfile), status='old')
  read (2,*)
  read (2,*) nx, ny, nzp
  close (2)
end subroutine read_particle_file_size

! --------------------------------------------------------------------------

subroutine read_particle_file (parfile, nx, ny, nzp, nscattab, &
                               DropNumConc,  delx, dely, zpar, temppar,  &
                               ncomp, ptype, masscont, reff)
 ! Reads the particle physical properties file.  The file contains a header
 ! containing the array size (nx,ny,nzp), heights (zpar), and
 ! temperature profile (temppar).  Then each row has format: 
 !   iX iY iZ  Ncomp Ptype1 masscont1 Reff1 ... PtypeN masscontN ReffN
 ! where Ncomp is the number of components, Ptype is the particle type, 
 ! masscont is the mass content (g/m^3), and Reff is the particle 
 ! effective radius (micron).
 ! Also reads the 1 and 2 parameter LWC files with LWC and Reff at each cell.
 ! The droplet number concentration (cm^-3) is used to derive the 
 ! effective radius from LWC for the 1 parameter LWC file (assumes a
 ! gamma distribution with alpha=7).
  implicit none
  character(len=*), intent(in) :: parfile
  integer, intent(in) :: nx, ny, nzp, nscattab
  real,    intent(in) :: DropNumConc
  real,    intent(out) :: delx, dely, zpar(nzp+1), temppar(nzp+1)
  integer, intent(out) :: ncomp(nx,ny,nzp), ptype(nscattab,nx,ny,nzp)
  real,    intent(out) :: masscont(nscattab,nx,ny,nzp), reff(nscattab,nx,ny,nzp)
  
  ! Local variables
  integer :: i, filekind, ix, iy, iz, nc, pt(nscattab)
  real    :: mass(nscattab), re(nscattab)

   ! Initialize the output arrays (in case there are missing grid points)
  ncomp(:,:,:) = 0
  ptype(:,:,:,:) = 0
  masscont(:,:,:,:) = 0.0
  reff(:,:,:,:) = 0.0
  
   ! Read in the header
  open (unit=2, file=trim(parfile), status='old')
  read (2,*) filekind
  read (2,*) ! nx, ny, nzp
  read (2,*) delx, dely   
  read (2,*) (zpar(iz), iz=1, nzp+1)
  read (2,*) (temppar(iz), iz=1, nzp+1)
  
  select case(filekind)
    case(1, 2)
      if (nscattab /= 1) &
        stop 'read_particle_file: Must have only one scattering table to use 1 or 2 parameter LWC file.'
      ! Read in the data
      do while (.true.)  
        if (filekind == 1) then
          read (2,*,end=190) ix, iy, iz, mass(1)
          re(1) = 100* ( mass(1) *0.75*1.3889/(3.14159*DropNumConc) )**(1.0/3)
        else
          read (2,*,end=190) ix, iy, iz, mass(1), re(1)
        endif
        if (ix >= 1 .and. ix <= nx .and. iy >= 1 .and. iy <= ny .and. &
            iz >= 1 .and. iz <= nzp) then
          ncomp(ix,iy,iz) = 1
          ptype(1,ix,iy,iz) = 1
          masscont(1,ix,iy,iz) = mass(1) 
          reff(1,ix,iy,iz) = re(1)
        endif
      enddo  
    
    case(3)
      ! Read in the data
      do while (.true.)  
        read (2,*,end=190) ix, iy, iz, nc, &
                         (pt(i), mass(i), re(i), i=1,min(nc,nscattab))
        if (any(pt(:min(nc,nscattab)) > nscattab)) &
          stop 'read_particle_file: Particle type greater than number of scattering tables.'
        if (ix >= 1 .and. ix <= nx .and. iy >= 1 .and. iy <= ny .and. &
            iz >= 1 .and. iz <= nzp) then
          ncomp(ix,iy,iz) = nc
          ptype(1:nc,ix,iy,iz) = pt(1:nc)
          masscont(1:nc,ix,iy,iz) = mass(1:nc)
          reff(1:nc,ix,iy,iz) = re(1:nc)
        endif
      enddo  
    
  case default
    stop 'read_particle_file: Must be type 2 or 3 particle properties file.'
  end select 
  
190 continue
  close (2) 
end subroutine read_particle_file

! --------------------------------------------------------------------------

subroutine organize_levels (nZp, Zpar, TempPar, &
                            numOtherLevels, OtherHeights, OtherTemps, &
                            nZt, Zlevels, Temp, izLevelBase)
 ! Adds the "other" layers to the layers from the particle file.
 ! The other heights must be outside the range of heights in the particle file.
 ! The base level index of the particle levels is returned in izLevelBase.
  implicit none
  integer, intent(in) :: nZp, numOtherLevels, nZt
  real,    intent(in) :: Zpar(nZp+1), TempPar(nZp+1)
  real,    intent(in) :: OtherHeights(numOtherLevels), OtherTemps(numOtherLevels)
  real,    intent(out) :: Zlevels(nZt+1), Temp(nZt+1)
  integer, intent(out) :: izLevelBase
  integer :: i, j, k

  
  if(any(Zpar(2:nZp+1) - zpar(1:nzp) <= 0)) &
    stop 'organize_levels: Zpar must increase'
  if(any(OtherHeights(:) >= Zpar(1) .and. OtherHeights(:) <= Zpar(nZp+1))) &
    stop 'organize_levels: OtherHeights must be outside particle file height range'

  do j = 1, numOtherLevels-1
    if (OtherHeights(j) >= OtherHeights(j+1)) then
      stop 'organize_levels: OtherHeights must increase'
    endif
  enddo
  
  k = 1  ;  j = 1
  do while (j<=numOtherLevels)
    if(OtherHeights(j) > Zpar(1)) exit
    Zlevels(k) = OtherHeights(j)
    Temp(k) = OtherTemps(j)
    k = k + 1  ;  j = j + 1
  enddo
  izLevelBase = k
  do i = 1, nZp+1
    Zlevels(k) = Zpar(i)
    Temp(k) = TempPar(i)
    k = k + 1
  enddo
  do while (j<=numOtherLevels)
    Zlevels(k) = OtherHeights(j)
    Temp(k) = OtherTemps(j)
    k = k + 1  ;  j = j + 1
  enddo
end subroutine organize_levels


! --------------------------------------------------------------------------

subroutine read_molec_abs_file (MolecAbsFileName, nZt, Zlevels, GasExt)
 ! Reads the three line ascii molecular absorption file which contains
 ! the extinction profile in the format:
 !        nZ               [number of Z grid cells]
 !        Zlevels(1:nZ+1)  [increasing heights of cell boundaries in km]
 !        GasExt(1:nZ)     [molecular extinction in km^-1 for each cell]
 ! If MolecAbsFileName is 'NONE' or '' then no file is read and the
 ! GasExt profile is set to zeros.
  implicit none
  character(len=*), intent(in) :: MolecAbsFileName
  integer, intent(in) :: nZt
  real,    intent(in) :: Zlevels(1:nZt+1)
  real,    intent(out) :: GasExt(1:nZt+1)
  integer :: nZ
  real, allocatable :: Zlevin(:)

  GasExt(:) = 0.0
  if (trim(MolecAbsFileName) /= 'NONE' .and. len_trim(MolecAbsFIlename) > 0) then
    open (unit=2, file=trim(MolecAbsFileName), status='old')
    read (2,*) nZ
    allocate (Zlevin(nZ+1))
    read (2,*) Zlevin(1:nZ+1)
    if (nZ /= nZt .or. any(abs(Zlevin(:) - Zlevels(:)) > spacing(Zlevels))) then
      print *, 'read_molec_abs_file: input Z levels do not match'
      stop
    endif
    deallocate (Zlevin)
    read (2,*) GasExt(1:nZt)
  endif
end subroutine read_molec_abs_file


! --------------------------------------------------------------------------

subroutine rayleigh_extinct (nzt, Zlevels, Temp, wavelen, RaylExt)
 ! Computes the molecular Rayleigh extinction profile RaylExt [/km]
 ! from the temperature profile Temp [k] at Zlevels [km].  Assumes
 ! a linear lapse rate between levels to compute the pressure at
 ! each level.  The Rayleigh extinction is proportional to air
 ! density and depends on the wavelength [um].  The Rayleigh
 ! extinction is calculated at grid cell boundaries and then
 ! interpolated to get the correct average extinction assuming 
 ! an exponential behavior.
 ! The extinction profile is returned with zeros if wavelen<=0.
  implicit none
  integer, intent(in) :: nzt
  real,    intent(in) :: Zlevels(nzt+1), Temp(nzt+1), wavelen
  real,    intent(out) :: RaylExt(nzt)
  
  integer :: i
  real    :: raylcoef, pres, lapse, ts, dz, extlev(nzt+1)

  raylcoef = 2.97E-4*wavelen**(-4.15+0.2*wavelen)

   ! Find surface pressure by integrating hydrostatic relation
   !  for a dry atmosphere up to surface height.
  pres = 1013.
  ts = Temp(1)
  lapse = 6.5*0.001
  pres = pres*(ts/(ts+lapse*Zlevels(1)*1000.))**(9.8/(287.*lapse))

   ! Use layer mean temperature to compute fractional pressure change.
  do i = 1, nzt
    dz = 1000.*(Zlevels(i+1)-Zlevels(i))
    lapse = (Temp(i)-Temp(i+1))/dz
    if (abs(lapse) .gt. 0.0001) then
      pres = pres*(Temp(i+1)/Temp(i))**(9.8/(287.*lapse))
    else
      pres = pres*exp(-9.8*dz/(287.*Temp(i)))
    endif
  enddo
  extlev(:) = raylcoef*pres/Temp(:)
  RaylExt(1:nzt) = (extlev(1:nzt)-extlev(2:nzt+1)) & 
                  / log(extlev(1:nzt)/extlev(2:nzt+1))
end subroutine rayleigh_extinct

! --------------------------------------------------------------------------

