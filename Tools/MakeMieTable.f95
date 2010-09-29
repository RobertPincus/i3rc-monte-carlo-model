! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

PROGRAM MakeMiePhaseFunctionTable
! $Revision$
! $URL$
! 
! Does Mie computations to create an I3RC Monte Carlo phase function table 
! as a function of effective radius for gamma or lognormal size distributions
! of spherical particles.  The particles may be water or ice (in which case 
! the program provides the index of refraction depending on wavelength) or
! "aerosols" (in which case the index of refraction is user specified).
! For water or ice particles the scattering properties may be averaged 
! over the desired spectral range with Planck function weighting.  
! The phase functions in the output scattering table are represented 
! with Legendre series.
!
! Mie calculations are done at a large number of drop sizes, and then the 
!   single scattering parameters at each value of effetive radius are computed 
!   by integrating those properties across the drop size distribution
!
! Input parameters are specified with a namelist file.
!
!    Frank Evans    University of Colorado     December 2005
!      Revisions by Robert Pincus, Climate Diagnostics Center, Jan 2006
!

   ! Modules from the I3RC community Monte Carlo model framework
  use ErrorMessages
  use scatteringPhaseFunctions
  use UserInterface

  IMPLICIT NONE
  !
  ! Input parameters
  !
  INTEGER :: NRETAB = 0
  REAL    :: WAVELEN1 = 0., WAVELEN2 = 0., DELTAWAVE = 0., PARDENS = 0.
  REAL    :: SRETAB = 0, ERETAB = 0, ALPHA = 0, MAXRADIUS = 0
  COMPLEX :: RINDEX
  CHARACTER(LEN=1)   :: PARTYPE = "W", AVGFLAG = "C", DISTFLAG = "G"
  CHARACTER(LEN=256) :: phaseFunctionTableFile = "phaseFunctionTable.pft"
  NAMELIST /mie_table_input/ WAVELEN1, WAVELEN2, AVGFLAG, DELTAWAVE, &
                             PARTYPE, RINDEX, PARDENS, DISTFLAG, ALPHA, &   
                             NRETAB, SRETAB, ERETAB, MAXRADIUS,  phaseFunctionTableFile
  !
  ! Local variable
  !
  INTEGER :: NSIZE, MAXLEG, I, J, L, NL
  LOGICAL :: LOGSPACEDREFF
  REAL    :: PI, WAVELENCEN, XMAX, SCATTER
  CHARACTER(LEN=256) :: namelistFileName
  character(len=256) :: tableDescription
  
  INTEGER, ALLOCATABLE :: NLEG1(:), NLEG(:)
  REAL,    ALLOCATABLE :: RADII(:), ND(:)
  REAL,    ALLOCATABLE :: EXTINCT1(:), SCATTER1(:), LEGEN1(:,:)
  REAL,    ALLOCATABLE :: REFF(:), EXTINCT(:), SSALB(:), LEGCOEF(:,:)

  !
  ! Temperature at which to evaluate index of refraction for water and ice
  !
  real, parameter :: waterTemperature = 283., iceTemperature = 243. 

  type(ErrorMessage)               :: status
  type(phaseFunction), allocatable :: MiePhaseFuncs(:)
  type(phaseFunctionTable)         :: phaseFuncTable


   ! Get the input parameters
! ------------------------------------------------------------------------------  
  PI = ACOS(-1.0)
  RINDEX = CMPLX(0., 0.)

  namelistFileName = getOneArgument()
  OPEN (UNIT=1, FILE=trim(namelistFileName), STATUS='OLD')
  READ (1,NML=mie_table_input)
  CLOSE (1)

  !
  ! Check for valid input values and consistency among parameters
  !
  if(nReTab == 0) stop "MakeMieTable: Must specify at least one effective radius" 
  if(wavelen1 <= 0) &
    stop "MakeMieTable: Must specify a wavelength at which to calculate single scattering parameters."
  if(wavelen2 == 0.) wavelen2 = wavelen1  
  if(alpha <= 0.) stop "MakeMieTable: must specify parameter alpha for size distribution"
  if(sretab <= 0) stop "MakeMieTable: must specify a starting effective radius (sretab)"
  if(nretab == 0) stop "MakeMieTable: must specify at least one effective radius." 
  if(eretab <= 0.) then
    eretab = sretab
    if(nretab > 1) print *, "MakeMieTable: specified a range of effective radii but requested only 1. " // &
                            "Single scattering parameters will be computed at only one value."  
    nretab = 1
  end if
  if(maxradius <= 0) maxradius = 25 * max(sretab, eretab)
  LOGSPACEDREFF = NRETAB < 0
  NRETAB = ABS(NRETAB)

  ! 
  ! Default density for water and ice particles (g/cm^3)
  !
  select case(trim(partype))
    case('w', 'W')
      if(pardens /= 0.) print *, "Water particles are specified. Using default density." 
      pardens = 1. 
    case ('i', "I") 
      if(pardens /= 0.) print *, "Ice particles are specified. Using default density." 
      pardens = 0.916
    case default
      if(pardens <= 0.) stop "MakeMieTable: Must specify an particle density (g/cm^3) for aerosols"
      if(RINDEX == cmplx(0., 0.)) stop "MakeMieTable: must specify a refractive index for aerosols"
  end select
  
   ! Calculate the maximum size parameter and the max number of Legendre terms
  CALL GET_CENTER_WAVELEN (WAVELEN1, WAVELEN2, WAVELENCEN)
  IF (AVGFLAG == 'A' .or. AVGFLAG == 'a') THEN
    XMAX = 2*PI*MAXRADIUS/WAVELEN1
  ELSE
    XMAX = 2*PI*MAXRADIUS/WAVELENCEN
  ENDIF  
  !
  ! Maximum number of Legendre coefficients to compute (for largest
  !    size parameter) 
  ! Reference Wiscombe, W.J., 1980: Improved Mie scattering algorithms
  !   Appl. Opt. Vol 19, pg 1505-1509. 
  !
  MAXLEG = NINT(2*(XMAX + 4.0*XMAX**0.3334 + 2))

   ! Get the average index of refraction for water or ice
  select case(trim(partype))
    case("w", "W", "i", "I")
      if(rindex /= cmplx(0., 0.)) &
        print *, "Water or ice particles specified. Computing index of refraction " // &
                 "(ignoring the specified value). " 
      RINDEX =  GET_REFRACT_INDEX (PARTYPE, WAVELEN1, WAVELEN2)
  end select
 
! -----------------------------------------------------------
! Now we're ready to do the computation
! 
   ! Figure the number of radii there will be
  nsize = GET_NSIZE(SRETAB, MAXRADIUS, WAVELENCEN)

   ! Allocate all the arrays here
   ! Drop radius and number concentration (for all drops)
  ALLOCATE (RADII(NSIZE), ND(NSIZE))
   ! Extinction, single scattering albedo, number of Legendre coefficents, and 
   !   the value of those coefficients
  ALLOCATE (EXTINCT1(NSIZE), SCATTER1(NSIZE), NLEG1(NSIZE), LEGEN1(0:MAXLEG,NSIZE))
   ! Effective radius, extinction, and ssa at each effective radius
  ALLOCATE (REFF(NRETAB), EXTINCT(NRETAB), SSALB(NRETAB))
   ! Number of Legendre coefficients and their value at each tabulated effective radius
  ALLOCATE (NLEG(NRETAB), LEGCOEF(0:MAXLEG,NRETAB))


   ! Make up the discrete particle radii to use
  radii(:NSIZE) = GET_SIZES (SRETAB, MAXRADIUS, WAVELENCEN, NSIZE)
   ! Do the Mie computations for each radius, which may involve several
   !   Mie calculation over the wavelength integration
  CALL COMPUTE_MIE_ALL_SIZES (AVGFLAG, WAVELEN1, WAVELEN2, DELTAWAVE, PARTYPE, &
                              WAVELENCEN, RINDEX, NSIZE, RADII, MAXLEG, &
                              EXTINCT1, SCATTER1, NLEG1, LEGEN1)
  
  
  ! At which values of effective radius do we build the table? 
  IF (NRETAB == 1) THEN
    REFF(:nReTab) = SRETAB
  ELSE
    IF (LOGSPACEDREFF) THEN
      REFF(:) = SRETAB*(ERETAB/SRETAB)**(FLOAT( (/ (i, i = 0, nReTab - 1)/))/(NRETAB-1))
    ELSE
      REFF(:) = SRETAB + (ERETAB - SRETAB) * FLOAT((/ (i, i = 0, nReTab - 1) /) )/(NRETAB-1) 
    ENDIF
  ENDIF
  
  ! Loop over the number of output tabulated effective radii
  DO I = 1, NRETAB
    ! Set tabulated effective radius

    ! Calculate the discrete size number concentrations (ND), which vary
    !   according to a truncated gamma or lognormal distribution,
    !   that gives the desired effective radius (REFF) and LWC (1 g/m^3).
    CALL MAKE_SIZE_DIST (DISTFLAG, PARDENS, RADII, REFF(I), ALPHA, ND(:))

    ! Sum the scattering properties over the discrete size distribution
    extinct(i) = dot_product(nd(:nsize), extinct1(:nsize))
    scatter    = dot_product(nd(:nsize), scatter1(:nsize))
    LEGCOEF(:,I) = 0.0
    NL = 1
    do J = 1, NSIZE
      NL = max(NL, nLeg1(J))
      LEGCOEF(0:NL,I) = LEGCOEF(0:NL,I) + ND(J) * LEGEN1(0:NL,J)
    END do
    
    LEGCOEF(0:NL,I) = LEGCOEF(0:NL,I)/SCATTER
    DO L = 0, NL
      IF (LEGCOEF(L,I) .GT. 0.5E-5) NLEG(I) = L
    ENDDO
    
    !
    ! Sanity check - the first Legendre coefficient should be identically 1
    !
    IF (ABS(LEGCOEF(0,I)-1.0) > 0.0001) THEN
      PRINT *,'Phase function not normalized for Reff=',REFF(I),LEGCOEF(0,I)
      STOP
    ENDIF 
    IF (EXTINCT(I) > 0.0) THEN
      SSALB(I) = SCATTER/EXTINCT(I)
    ENDIF

  ENDDO  ! end of effective radius loop
  !
  ! Convert to units of extinction in g/m3
  !
  extinct(:nReTab) = 0.001 * extinct(:nReTab)

   ! Put the phase function, extinction, and single scattering albedo
   !  in a phase function object for each effective radius
  forall(i = 1:nReTab)  &
    LegCoef(0:Nleg(i), i) = LegCoef(0:Nleg(i), i) / (/ (2*l+1, l=0, Nleg(i)) /)
  !
  ! Trap values of single scattering albedo that are slightly bigger than 1
  !
  where(SSalb(1:Nretab) > 1. .and. SSalb(1:Nretab) <= 1. + 2. * SPACING(1.)) SSalb(1:Nretab) = 1
   
  ALLOCATE (MiePhaseFuncs(Nretab))
  DO i = 1, Nretab
    MiePhaseFuncs(i) = new_PhaseFunction( LegCoef(1:Nleg(i),i), &
                                          extinction=Extinct(i), &
                                          singleScatteringAlbedo=SSalb(i), &
                                          status=status)
  ENDDO
  CALL printStatus(status)
  
   ! Make the phase function table object, labeled with effective radius
   ! Provide a description of how the table was created. 
  tableDescription = "Mie phase function table for spheres made of" 
  select case (trim(partype))
   case('w', 'W')
     tableDescription = trim(tableDescription) // " water"
   case('i', 'I')
     tableDescription = trim(tableDescription) // " ice"
   case('a', 'A')
     tableDescription = trim(tableDescription) // " aerosol"
   case default
      tableDescription = trim(tableDescription) // " an unknown material "
  end select 
  tableDescription = trim(tableDescription) // " at a concentration of 1 g/m^3. Key is in microns. "  
  select case (trim(distflag))
   case('g', 'G')
     tableDescription = trim(tableDescription) // " Gamma size distribution. "
   case('l', 'L')
     tableDescription = trim(tableDescription) // " Lognormal size distribution. "
  end select
  
  phaseFuncTable = new_PhaseFunctionTable( &
                               MiePhaseFuncs(1:Nretab), key=reff(1:Nretab), &
                               tableDescription =trim (tableDescription),   &
                               status = status)
  CALL printStatus(status)
  
   ! Write the phase function table object to a new file
  CALL write_PhaseFunctionTable(phaseFuncTable, phaseFunctionTableFile, status)
  CALL printStatus(status)

   ! Clean up unused objects - not really necessary since the program is about 
   !   to end, but useful as an example
  DEALLOCATE (REFF, EXTINCT, SSALB, NLEG, LEGCOEF)
  DEALLOCATE (ND, EXTINCT1, SCATTER1, NLEG1, LEGEN1)
  DO i = 1, Nretab
    CALL finalize_PhaseFunction(MiePhaseFuncs(i))
  ENDDO
  DEALLOCATE (MiePhaseFuncs)
  call finalize_PhaseFunctionTable(phaseFuncTable)
contains 
! ------------------------------------------------------------------------------  

  elemental function planckRadiation(wavelength, temperature)
    implicit none
    real, intent(in) :: wavelength, temperature
    real             :: planckRadiation
  !
  ! Computes the Planck blackbody radiation (W/m2-str) as a function of 
  !   wavelength (microns) and black body temperature (K)
  !

    planckRadiation = (1.19E8/wavelength**5)/(EXP(1.439E4/(wavelength * temperature))-1.)

  end function planckRadiation

! ------------------------------------------------------------------------------  

   function effectiveBlackBodyTemp(wavelength1, wavelength2)
     implicit none
     real, intent(in) :: wavelength1, wavelength2 ! in microns
     real             :: effectiveBlackBodyTemp   ! K
     !
     ! Computes the black body temperature at the center
     !   of a wavelength interval. In the 
     !   shortwave (wavelengths < 3 microns) we use the sun's 
     !   temperature (5800K); in the longwave (> 5 microns) we use a temperaure
     !   of 270K, and in between we set the BB temperature = -1. 
     !
     real, parameter :: solarBBTemp = 5800., terrestrialBBTemp = 270., &
                        maxSolarBBWavelen = 3.0, minTerrestBBWavelen = 5.0
     real :: WAVECEN
     
     WAVECEN = 0.5*(wavelength1+wavelength2)
     IF (WAVECEN < maxSolarBBWavelen) THEN
       effectiveBlackBodyTemp = solarBBTemp
     ELSE IF (WAVECEN > minTerrestBBWavelen) THEN
       effectiveBlackBodyTemp = terrestrialBBTemp
     ELSE
       effectiveBlackBodyTemp = -1.0
     ENDIF 
     
   end function effectiveBlackBodyTemp

! ------------------------------------------------------------------------------  

  function computeNumPlanckWeightingWaves(wavelength1, wavelength2)
    implicit none
    real, intent(in) :: wavelength1, wavelength2
    integer          :: computeNumPlanckWeightingWaves
    !
    ! Compute how many indivdual wavelengths are required in order
    !   to weight a quantity by the Planck function across a wavelength 
    !   interval 
    !
    real :: wavecen, delwave
    
    IF (wavelength1 == wavelength2) THEN
      computeNumPlanckWeightingWaves = 1  
    ELSE
      WAVECEN = 0.5*(WAVELEN1+WAVELEN2)
      DELWAVE = MIN(WAVECEN/100., 0.1*ABS(WAVELEN2-WAVELEN1))
      DELWAVE = MAX(DELWAVE, WAVECEN*1.0E-5)
      computeNumPlanckWeightingWaves = int(ABS(WAVELEN2-WAVELEN1)/delwave)
    END IF 
      
  end function computeNumPlanckWeightingWaves

! ------------------------------------------------------------------------------  

  subroutine planckWeightingWavelengths(wavelength1, wavelength2, wavelengths)
    real, intent(in)  :: wavelength1, wavelength2
    real, dimension(0:), &
          intent(out) :: wavelengths
    !
    ! Returns the wavelengths requiredto weight a quantity by the Planck 
    !  function across a wavelength interval 
    ! 
    
    integer :: nWavelengths, i
    
    nWavelengths = size(wavelengths) - 1
    if(computeNumPlanckWeightingWaves(wavelength1, wavelength2) /= nWavelengths) &
      stop 'planckWeightingWavelengths: wrong number of wavelengths requested' 
      
    IF (wavelength1 == wavelength2) THEN
      wavelengths(0:) = wavelength1  
    ELSE 
      wavelengths(0:) = wavelength1 + (wavelength2 - wavelength1) * &
                                      float((/ (i, i = 0, nWavelengths) /))/ nWavelengths
    END IF
  end subroutine planckWeightingWavelengths
  
! ------------------------------------------------------------------------------  
  
  
  SUBROUTINE GET_CENTER_WAVELEN (WAVELEN1, WAVELEN2, WAVELENCEN)
  !  Returns the Planck weighted center wavelength averaged over the 
  ! wavelength interval (WAVELEN1 < WAVELEN2 [microns]).  A solar
  ! blackbody temperature of 5800 K is used for the Planck weighting
  ! if the average wavelength is less than 3 microns, no Planck weighting
  ! is done for an average wavelength between 3 and 5 microns, and a 
  ! blackbody temperature of 270 K is done for an average wavelength 
  ! greater than 5 microns.
    IMPLICIT NONE
    REAL, INTENT(IN)  :: WAVELEN1, WAVELEN2
    REAL, INTENT(OUT) :: WAVELENCEN
    
    ! Local variables
    REAL              :: BBTEMP
    integer           :: nSteps
    real, allocatable :: wavelengths(:), planckValues(:)
  
    IF (WAVELEN1 == WAVELEN2) THEN
      WAVELENCEN = WAVELEN1  
    ELSE
      BBTEMP = effectiveBlackBodyTemp(WAVELEN1, WAVELEN2)
      
      nSteps = computeNumPlanckWeightingWaves(WAVELEN1, WAVELEN2)
      allocate(wavelengths(0:nSteps), planckValues(0:nSteps))
      call planckWeightingWavelengths(WAVELEN1, WAVELEN2, wavelengths)

      if(bbtemp > 0) then
        planckValues(:) = planckRadiation(wavelengths(:), bbtemp)
      else
        planckValues(:) = 1. 
      end if 
      wavelencen = .001 * int(1000 * dot_product(planckValues(:), wavelengths(:)) / &
                                     sum(planckValues(:)) )
      deallocate(wavelengths, planckValues)
    ENDIF
  END SUBROUTINE GET_CENTER_WAVELEN
  
! ------------------------------------------------------------------------------    
  
  
  function GET_REFRACT_INDEX (PARTYPE, WAVELEN1, WAVELEN2) result(RINDEX)
   ! Returns the index of refraction for water or ice averaged over
   ! the wavelength interval (WAVELEN1 < WAVELEN2 [microns]).   The
   ! averaging is done at 0.05 micron intervals and is weighted by
   ! a Planck function at 5800 K temperature for central wavelengths
   ! less than 3 microns, a flat function between 3 and 5 microns, and
   ! 270 K Planck function beyond 5 microns.  
   ! The index of refraction is using parameters iceTemperature and 
   !   waterTemperature in the main program
   ! (the temperature dependence is important in the microwave).
    IMPLICIT NONE
    CHARACTER(LEN=1), INTENT(IN) :: PARTYPE
    REAL,             INTENT(IN) :: WAVELEN1, WAVELEN2
    COMPLEX                      :: RINDEX
    
    REAL    :: BBTEMP
    REAL    :: MRE, MIM, A
    integer :: i, nsteps
    real, allocatable :: M_real(:), M_Complex(:), wavelengths(:), planckValues(:)
  
    BBTEMP = effectiveBlackBodyTemp(WAVELEN1, WAVELEN2)

    nSteps = computeNumPlanckWeightingWaves(WAVELEN1, WAVELEN2)
    allocate(wavelengths(0:nSteps), planckValues(0:nSteps), &
             M_real(0:nSteps), M_Complex(0:nSteps))
    call planckWeightingWavelengths(WAVELEN1, WAVELEN2, wavelengths)
    
    if(BBTEMP > 0) then
      planckValues = planckRadiation(wavelengths, BBTEMP)
    else 
      planckValues = 1. 
    end if 

    IF (PARTYPE == 'I') THEN
      do i = 0, nSteps
        CALL REFICE (0, wavelengths(i),   iceTemperature, M_real(I), M_Complex(i), A, A)
      end do
    ELSE
      do i = 0, nSteps
        CALL REFWAT (0, wavelengths(i), waterTemperature, M_real(I), M_Complex(i), A, A)
      end do
    end if
    
    MRE = dot_product(planckValues(:), M_real(:)   ) / sum(planckValues(:))
    MIM = dot_product(planckValues(:), M_Complex(:)) / sum(planckValues(:))
    RINDEX = CMPLX(MRE, -MIM)
  END function GET_REFRACT_INDEX
  
! ------------------------------------------------------------------------------  

  function GET_NSIZE (SRETAB, MAXRADIUS, WAVELEN) result(nsize)
   ! Calculates the number of radii for which the Mie computation will be run.
   ! The formula and spacing in size parameter can be changed to trade
   ! off size distribution integration accuracy vs. computer time.
    IMPLICIT NONE
    REAL,    INTENT(IN)  :: SRETAB, MAXRADIUS, WAVELEN
    INTEGER              :: NSIZE
    
    ! Local variables
    REAL    :: TWOPI, RADMIN, RAD, X, DELX, DELRAD
  
    TWOPI = 2.0*ACOS(-1.0)
    RADMIN = 0.02*SRETAB
    RAD = RADMIN
    NSIZE = 1
    DO WHILE (RAD < MAXRADIUS)
      X = TWOPI*RAD/WAVELEN
      DELX = MAX(0.01,0.03*X**0.5)    ! coarser spacing at large size parameters
  !    DELX = 0.1                     ! One alternative method
      DELRAD = DELX*WAVELEN/TWOPI
      RAD = RAD + DELRAD
      NSIZE = NSIZE + 1
    ENDDO
  END function GET_NSIZE
  
! ------------------------------------------------------------------------------  
  
  function GET_SIZES (SRETAB, MAXRADIUS, WAVELEN, NSIZE) result(radii)
   ! Calculates the radii for which the Mie computation will be run and
   ! from which all the size distributions will be computed.
   ! The formula and spacing in size parameter can be changed to trade
   ! off size distribution integration accuracy vs. computer time.
    IMPLICIT NONE
    INTEGER, INTENT(IN ) :: NSIZE
    REAL,    INTENT(IN ) :: SRETAB, MAXRADIUS, WAVELEN
    REAL                 :: RADII(NSIZE)
    
    ! Local variables
    INTEGER :: N
    REAL    :: TWOPI, RAD, X, DELX, DELRAD
  
    TWOPI = 2.0*ACOS(-1.0)
    RAD = 0.02*SRETAB
    RADII(1) = RAD
    DO N = 2, NSIZE
      X = TWOPI*RAD/WAVELEN
      DELX = MAX(0.01,0.03*sqrt(X))    ! coarser spacing at large size parameters
  !    DELX = 0.1                     ! One alternative method
      DELRAD = DELX*WAVELEN/TWOPI
      RAD = RAD + DELRAD
      RADII(N) = RAD
    ENDDO
  END function GET_SIZES
  
! ------------------------------------------------------------------------------    
   
  
  SUBROUTINE COMPUTE_MIE_ALL_SIZES (AVGFLAG, WAVELEN1, WAVELEN2, DELTAWAVE, &
                                    PARTYPE, WAVELENCEN, RINDEX, NSIZE, RADII, &
                                    MAXLEG, EXTINCT1, SCATTER1, NLEG1, LEGEN1)
   ! Does a Mie computation for each particle radius in RADII and returns the
   ! optical properties in arrays EXTINCT1, SCATTER1, NLEG1, and LEGEN1.
   ! For AVGFLAG='C' the computation is done at a single wavelength (WAVELENCEN),
   ! using the input index of refraction (RINDEX).  For AVGFLAG='A' an
   ! integration of the Mie properties over wavelength is performed for
   ! each radius.  For each wavelength, with spacing DELTAWAVE, the water 
   ! or ice (depending on PARTYPE) index of refraction is obtained and
   ! used in the Mie computation for that wavelength, and the Mie optical
   ! properties are averaged with Planck function weighting (blackbody
   ! temperature depends on wavelength).  The Legendre coefficients are
   ! returned with the product of the phase function times the scattering
   ! coefficient.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NSIZE, MAXLEG
    REAL,    INTENT(IN) :: WAVELEN1, WAVELEN2, DELTAWAVE, WAVELENCEN
    REAL,    INTENT(IN) :: RADII(NSIZE)
    COMPLEX, INTENT(IN) :: RINDEX
    CHARACTER(LEN=1), &
             INTENT(IN) :: AVGFLAG, PARTYPE
    INTEGER, INTENT(OUT) :: NLEG1(NSIZE)
    REAL,    INTENT(OUT) :: EXTINCT1(NSIZE), SCATTER1(NSIZE)
    REAL,    INTENT(OUT) :: LEGEN1(0:MAXLEG,NSIZE)
    
    ! Local variables
    INTEGER :: I, NL
    REAL    :: WAVE, BBTEMP, PLANCK, SUMP, A
    REAL    :: MRE, MIM, EXT, SCAT, LEG(0:MAXLEG)
    COMPLEX :: REFIND
  
    IF (AVGFLAG == 'C') THEN
       ! For using one central wavelength: just call Mie routine for each radius
      DO I = 1, NSIZE
        CALL MIE_ONE (WAVELENCEN, RINDEX, RADII(I), MAXLEG, &
                      EXTINCT1(I), SCATTER1(I), NLEG1(I), LEGEN1(0,I) )
      ENDDO
  
    ELSE
      ! For averaging over wavelength range:
      BBTEMP = effectiveBlackBodyTemp(WAVELEN1, WAVELEN2)
      
      EXTINCT1(:) = 0.0
      SCATTER1(:) = 0.0
      NLEG1(:) = 1
      LEGEN1(:,:) = 0.0
      SUMP = 0.0
      WAVE = WAVELEN1
      DO WHILE (WAVE <= WAVELEN2)   ! Loop over the wavelengths
        IF (BBTEMP > 0) PLANCK = planckRadiation(WAVE, BBTEMP)
        SUMP = SUMP + PLANCK
        IF (PARTYPE == 'I') THEN   ! Get the index of refraction of water or ice
          CALL REFICE (0, WAVE, iceTemperature,   MRE, MIM, A, A)
        ELSE
          CALL REFWAT (0, WAVE, waterTemperature, MRE, MIM, A, A)
        ENDIF
        REFIND = CMPLX(MRE,-MIM)
        DO I = 1, NSIZE
          CALL MIE_ONE (WAVE, REFIND, RADII(I), MAXLEG, EXT, SCAT, NL, LEG)
          EXTINCT1(I) = EXTINCT1(I) + PLANCK*EXT
          SCATTER1(I) = SCATTER1(I) + PLANCK*SCAT
          NLEG1(I) = MAX(NLEG1(I),NL)
          LEGEN1(0:NL,I) = LEGEN1(0:NL,I) + PLANCK*LEG(0:NL)
        ENDDO
        WAVE = WAVE + DELTAWAVE
      ENDDO
      EXTINCT1(:) = EXTINCT1(:)/SUMP
      SCATTER1(:) = SCATTER1(:)/SUMP
      LEGEN1(:,:) = LEGEN1(:,:)/SUMP
    ENDIF
    
  END SUBROUTINE COMPUTE_MIE_ALL_SIZES
  
! ------------------------------------------------------------------------------  
   
  SUBROUTINE MAKE_SIZE_DIST (DISTFLAG, PARDENS, RADII, REFF, ALPHA, ND)
   ! Calculates the number concentrations (ND in cm^-3) at  
   ! discrete particle RADII (micron) of a gamma or lognormal size distribution
   ! with an effective radius of REFF (micron), gamma shape parameter or 
   ! lognormal standard deviation of ALPHA, and mass content of 1 g/m^3.
    IMPLICIT NONE
    CHARACTER(LEN=1), INTENT(IN ) :: DISTFLAG
    REAL,             INTENT(IN ) :: RADII(:), REFF, ALPHA, PARDENS
    REAL,             INTENT(OUT) :: ND(:)
    
    REAL,    PARAMETER :: TOL=0.001  ! fractional tolerance in achieving Reff
    integer, parameter :: maxIterations = 8
    INTEGER :: I, NSIZE
    REAL    :: TRUERE, F, REHI, RELO, REMID
    
    NSIZE = size(radii)
    if(size(nd) /= NSIZE) &
      stop 'MAKE_SIZE_DIST: vectors RADII and ND must be the same length'
  
     ! See if the true effective radius is already close enough
    CALL DO_SIZE_DIST (PARDENS, DISTFLAG, ALPHA, REFF, RADII, ND, TRUERE)
    IF (ABS(TRUERE-REFF) < TOL*REFF) RETURN
    F = REFF/TRUERE
  
    IF (TRUERE < REFF) THEN 
      ! Find Reff that gives true Reff above desired value
      RELO = REFF
      REHI = F*REFF
      I = 0
      TRUERE = REFF/F
      DO WHILE (TRUERE <= REFF .AND. I < maxIterations)
        REHI = F*REHI
        I = I + 1
        CALL DO_SIZE_DIST (PARDENS,DISTFLAG, ALPHA, REHI, RADII, ND, TRUERE)
      ENDDO
      IF (TRUERE <= REFF) THEN
        PRINT *, 'MAKE_SIZE_DIST: effective radius cannot be achieved',REFF,TRUERE
        STOP
      ENDIF
    ELSE 
      ! Find Reff that gives true Reff below desired value
      REHI = REFF
      RELO = F*REFF
      I = 0
      TRUERE = REFF/F
      DO WHILE (TRUERE >= REFF .AND. I < maxIterations)
        RELO = F*RELO
        I = I + 1
        CALL DO_SIZE_DIST (PARDENS,DISTFLAG, ALPHA, RELO, RADII, ND, TRUERE)
      ENDDO
      IF (TRUERE >= REFF) THEN
        PRINT *, 'MAKE_SIZE_DIST: effective radius cannot be achieved',REFF,TRUERE
        STOP
      ENDIF
    ENDIF
    ! Do bisection to get correct effective radius
    DO WHILE (ABS(TRUERE-REFF) > TOL*REFF)
      REMID = 0.5*(RELO+REHI)
      CALL DO_SIZE_DIST (PARDENS, DISTFLAG, ALPHA, REMID, RADII, ND, TRUERE)
      IF (TRUERE < REFF) THEN
        RELO = REMID
      ELSE
        REHI = REMID
      ENDIF
    ENDDO  
  END SUBROUTINE MAKE_SIZE_DIST
  
! ------------------------------------------------------------------------------  
  
  SUBROUTINE DO_SIZE_DIST (PARDENS, DISTFLAG, ALPHA, RE, RADII, ND, TRUERE)
   ! For the input effective radius (RE) [um], returns the number concentrations 
   ! ND [cm^-3] and the calculated effective radius TRUERE [um] for a 
   ! gamma or lognormal size distribution with mass content of 1 g/m^3.
    IMPLICIT NONE
    REAL,    INTENT(IN) :: PARDENS
    CHARACTER(LEN=*), &
             INTENT(IN)  :: DISTFLAG
    REAL,    INTENT(IN)  :: ALPHA, RE, RADII(:)
    REAL,    INTENT(OUT) :: ND(:), TRUERE
    
    INTEGER :: nSize
    REAL    :: PI, A, B, LWC, SUM2, SUM3
    real, dimension(:) &
            :: deltaR(size(RADII))
    ! External function
    real    :: GAMMLN
  
    nSize = size(radii) 
    if(size(nd) /= nSize) &
      stop 'DO_SIZE_DIST: vectors RADII and ND must be the same length'
    PI = ACOS(-1.0)
    
    deltaR(2:nSize-1) = sqrt(radii(2:nSize-1) * radii(3:nSize)) - &
                        sqrt(radii(2:nSize-1) * radii(1:nSize-2))
    deltaR(1        ) = sqrt(radii(2    )     * radii(3   )) - &
                        radii(1)
    deltaR(nSize    ) = radii(nsize) - &
                        sqrt(radii(nsize) * radii(nSize-1)) 
                        
    IF (DISTFLAG == 'G') THEN       ! Gamma distribution
      B = (ALPHA+3)/RE
      A = 1.E6/( (4*PI/3.)*PARDENS *B**(-ALPHA-4) *EXP(GAMMLN(ALPHA+4.)) )
      ND(:) = A * RADII(:)**ALPHA *EXP(-B*RADII(:)) * deltaR(:)
    ELSE                            ! Lognormal distribution
      B = RE*EXP(-2.5*ALPHA**2)
      A = 1.E6/( (4*PI/3.)*PARDENS *SQRT(2*PI)*ALPHA * B**3 *EXP(4.5*ALPHA**2) )
      ND(:) = (A/RADII(:))*EXP(-0.5*(LOG(RADII(:)/B))**2/ALPHA**2) * deltaR(:)
    ENDIF
    
    SUM2 = dot_product(nd(:), RADII(:)**2)
    SUM3 = dot_product(nd(:), RADII(:)**3)
    TRUERE = SUM3/SUM2
    LWC = 1.0E-6 * PARDENS * (4.*PI/3.) * SUM3
    ND(:) = (1.0/LWC) * ND(:)
    
  END SUBROUTINE DO_SIZE_DIST
! ------------------------------------------------------------------------------  

end PROGRAM MakeMiePhaseFunctionTable


