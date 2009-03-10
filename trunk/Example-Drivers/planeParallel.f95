! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

program computePlaneParallelRT
  ! $Revision: 1.7 $, $Date: 2009/03/09 19:17:22 $
  ! $Name: Cornish-Gilliflower $
  !
  ! Program to compute flux, absorption, and intensity for a 
  !   plane parallel slab  
  !   using the I3RC Community Monte Carlo Model. This might be used 
  !   to compare the results from new or modified Monte Carlo solvers 
  !   with the fluxes and intensities computer using other methods 
  !   (e.g. plane-parallel). It's also meant as a bare-bones
  !   example of how to use the I3RC Monte Carlo code. 
  !
  ! The phase function can be read from a table or one can use a 
  !   Henyey-Greenstein phase function. The latter can be defined using  
  !   moments or angle-value pairs according to the value of useMoments in the
  !   opticalParameters namelist. 
  !
  
  ! It would be cool if we could divide the optical thickness 
  !   unequally among a specified number of components. 

! ---------------------------------------------------
  !
  ! Modules from the I3RC community Monte Carlo model framework
  !
  use ErrorMessages
  use RandomNumbers
  use scatteringPhaseFunctions
  use opticalProperties
  use surfaceProperties
  use monteCarloIllumination
  use UserInterface
  !
  ! The Monte Carlo integrator
  !
  use monteCarloRadiativeTransfer
  implicit none
  
                               
  ! Radiative transfer parameters
  ! 
  real     :: solarMu = 0.5,  solarAzimuth = 0.
  real     :: surfaceAlbedo = 0.0
  integer  :: numIntensityAngles = 0
  integer, parameter &
           :: maxNumIntensityAngles = 20
  real, dimension(maxNumIntensityAngles) &
           :: intensityMus = 0., intensityPhis = 0. 
  logical  :: computeIntensity = .false. 
                        
  namelist /radiativeTransfer/ solarMu, solarAzimuth, surfaceAlbedo, &        
                               intensityMus, intensityPhis
  ! 
  ! Monte Carlo parameters
  !
  integer  :: numPhotonsPerBatch = 1E5, numBatches = 4
  integer  :: iseed = 10, nPhaseIntervals = 10000
  
  namelist /monteCarlo/ numPhotonsPerBatch, numBatches, iseed, nPhaseIntervals

  !
  ! Monte Carlo algorithmic choices 
  !
  logical              :: useRayTracing = .true., useRussianRoulette = .true. 
  logical              :: useHybridPhaseFunsForIntenCalcs = .false. 
  real                 :: hybridPhaseFunWidth = 7. 
  integer              :: numOrdersOrigPhaseFunIntenCalcs = 0
  logical              :: useRussianRouletteForIntensity = .true. 
  real                 :: zetaMin = 0. 

  namelist /algorithms/ useRayTracing, useRussianRoulette, &
                        useHybridPhaseFunsForIntenCalcs, hybridPhaseFunWidth, &
                        numOrdersOrigPhaseFunIntenCalcs,                      &
                        useRussianRouletteForIntensity, zetaMin
  ! 
  ! Problem specification
  !
  ! Optical parameters
  real     :: SSA = 1.0, opticalDepth = 1.
  !
  ! Phase function: Henyey-Greenstien
  !
  real     :: g = 0.85
  integer  :: nLegendreCoefficients = 64, nAngles = 5000
  logical  :: useMoments = .true.
  !
  ! Phase function: read from a table
  !
  character(len = 256) :: phaseFunctionTableFile = "" 
  integer  :: phaseFunctionTableIndex = 0
  
  namelist /problemOptics/ SSA, opticalDepth,                             &
                           g, nLegendreCoefficients, nAngles, useMoments, &
                           phaseFunctionTableFile, phaseFunctionTableIndex

  real     :: domainSize = 500., physicalThickness = 250.
  integer  :: nLayers = 1, nX = 1, nY = 1
  logical  :: useSurfaceProperties = .false. 

  namelist /problemDomain/ domainSize, nLayers, physicalThickness, nx, ny, &
                           useSurfaceProperties

  !
  ! The domain can be written to a file if desired
  !
  character(len=256) :: domainFileName = "" 
  namelist /filenames/ domainFileName
  
  ! -----------------------------------
  ! Local variables
  integer            :: batch, i
  
  real, dimension(:, :, :), allocatable &
                     :: fluxUp, fluxDown, fluxAbsorbed
  real               :: meanFluxUp, meanFluxDown, meanFluxAbsorbed, &
                        fluxUpStdDev, fluxDownStdDev, fluxAbsorbedStdDev
  real, dimension(:, :, :, :), allocatable &
                     :: intensity
  real               :: intensityMean

  type(randomNumberSequence) :: randoms
  type(photonStream)         :: incomingPhotons
  type(domain)               :: exampleDomain
  type(integrator)           :: mcIntegrator
  type(ErrorMessage)              :: status
  type(surfaceDescription)        :: surface
  ! ----------------------------------------------------------------
  call getUserInput
  !
  ! Create a domain object of specified size in which all cells 
  !   have the same specified optical properties - this can 
  !   be used to test MC against known plane-parallel solutions. 
  !
  call createDomain
  
  !
  ! Set up the integrator object - the integrator makes copies of the 3D distribution 
  !   of optical components, so we can release the resources used by the domain object
  !
  mcIntegrator = new_Integrator(exampleDomain, status = status)
  call finalize_Domain(exampleDomain)
  call printStatus(status)
  
  !
  ! Set the surface albedo. 
  !
  if(useSurfaceProperties) then 
    surface = new_surfaceDescription(surfaceParameters = (/ surfaceAlbedo /), status = status)
    call printStatus(status) 
    call specifyParameters(mcIntegrator, surfaceBDRF = surface, status = status)
    call finalize_surfaceDescription(surface)
  else 
    call specifyParameters(mcIntegrator, surfaceAlbedo = surfaceAlbedo, status = status)
  end if
  call printStatus(status) 
  
  if(computeIntensity) then
    call specifyParameters(mcIntegrator, &
                           intensityMus  = intensityMus( 1:numIntensityAngles), &
                           intensityPhis = intensityPhis(1:numIntensityAngles), &
                           status = status)
    call printStatus(status) 
  end if 

  !
  ! Make the algorithmic choices
  !
  call specifyParameters(mcIntegrator,                              &
                         useRayTracing      = useRayTracing,        &
                         useRussianRoulette = useRussianRoulette,   &
                         useHybridPhaseFunsForIntenCalcs =          &
                                 useHybridPhaseFunsForIntenCalcs,   &
                         hybridPhaseFunWidth = hybridPhaseFunWidth, &
                         useRussianRouletteForIntensity =           &
                                 useRussianRouletteForIntensity,    &
                         zetaMin = zetaMin,                         &
                         status = status)
  call printStatus(status) 

  if(numBatches > 0 .and. numPhotonsPerBatch > 0) then 
    !
    ! Allocate output arrays based on the size of the problem. 
    !
    allocate(      fluxUp(nX, nY, numBatches), &
                 fluxDown(nX, nY, numBatches), & 
             fluxAbsorbed(nX, nY, numBatches))
    if(computeIntensity) allocate(intensity(nX, nY, numIntensityAngles, numBatches))
  
    ! 
    ! The following loop can be used to estimate the uncertainty in the flux 
    !   estimates by looking at the variance between numBatches independent calculations. 
    !    
    !   The loop could be run as a set of parallel threads if the integrator and
    !   random number objects were copied and the flux arrays shared between 
    !   threads. 
    !
    do batch = 1, numBatches
      ! 
      ! Variable randoms holds the state of the random number generator. 
      !  It's seeded at the begining of each loop. 
      !
      randoms = new_RandomNumberSequence(seed = (/ batch, iseed /))
      !
      ! The initial direction and position of the photons are precomputed and 
      !   stored in an "illumination" object. 
      ! If solarAzimuth isn't specified the photons are introduced at random azimuths. 
      ! If solarMu isn't specified the solar flux on the horizontal is equally weighted
      !  in the cosine of the solar zenith angle 
      !
      incomingPhotons = new_PhotonStream(solarMu, solarAzimuth,                &
                                         numberOfPhotons = numPhotonsPerBatch, &
                                         randomNumbers = randoms, status = status)
      call printStatus(status)
      
      !
      ! Now we compute the radiative transfer for this batch of photons. 
      !
      call computeRadiativeTransfer(mcIntegrator, randoms, incomingPhotons, status)
      call printStatus(status)
  
      If(computeIntensity) then 
        call reportResults(mcIntegrator, intensity = intensity(:, :, :, batch), status = status)
      else 
        call reportResults(mcIntegrator, meanFluxUp, meanFluxDown, meanFluxAbsorbed,  &
                           fluxUp(:, :, batch), fluxDown(:, :, batch), fluxAbsorbed(:, :, batch), &
                           status = status)
      end if
      call printStatus(status)
      call finalize_RandomNumberSequence(randoms)
      call finalize_PhotonStream(incomingPhotons)
    end do
    
    ! 
    ! Output to screen
    !
    if(computeIntensity) then
      print *, "  tau  omega   g  theta0    mu   phi radiance    error"
      do i = 1, numIntensityAngles
        intensityMean = sum(intensity(:, :, i, :))/(numBatches * nX * nY)
        write(*, '(f6.2, 1x, 2(f5.3, 1x), 1x, f5.2, 1x, f7.5, 1x, i3, 1x, f8.6, 1x, f10.8)') &
          opticalDepth, ssa, g, acos(solarMu) * 180./acos(-1.),            &
          intensityMus(i), int(intensityPhis(i)), intensityMean,           & 
          sqrt(sum( (sum(sum(intensity(:, :, i, :), dim=1), dim=1)/(nX*nY) - &
                     intensityMean)**2)/numBatches)
      end do
    else
      ! 
      ! This integrator provides fluxes at the top and bottom 
      !   of the domain. Both the domain mean and the fluxes in each level 
      !   are available. 
      !
      meanFluxUp       = sum(fluxUp(:, :, :))/(numBatches * nX * nY)
      meanFluxDown     = sum(fluxDown(:, :, :))/(numBatches * nX * nY)
      meanFluxAbsorbed = sum(fluxAbsorbed(:, :, :))/(numBatches * nX * nY)
      if(numBatches > 1) then 
        fluxUpStdDev       = sqrt(sum((fluxUp(:, :, :) - meanFluxUp)**2)/((numBatches - 1) * nX * nY))
        fluxDownStdDev     = sqrt(sum((fluxDown(:, :, :) - meanFluxDown)**2)/((numBatches - 1) * nX * nY))
        fluxAbsorbedStdDev = sqrt(sum((fluxAbsorbed(:, :, :) - meanFluxAbsorbed)**2)/((numBatches - 1) * nX * nY))
      else
        fluxUpStdDev   = 0
        fluxDownStdDev = 0
      end if 
    
      print *, "  tau  omega   g  theta0   Fup      Fdn    FluxUpErr FluxDownErr FluxAbs FluxAbsErr"
      write(*, '(f6.2, 1x, 2(f5.3, 1x), 1x, f5.2, 1x, 6(f7.5, 3x))') &
        opticalDepth, ssa, g, acos(solarMu) * 180./acos(-1.),        &
        meanFluxUp, meanFluxDown, fluxUpStdDev, fluxDownStdDev, meanFluxAbsorbed, fluxAbsorbedStdDev
    end if
    
    deallocate(fluxUp, fluxDown, fluxAbsorbed)
    if(computeIntensity) deallocate(intensity)
  end if 
contains
 ! ---------------
   subroutine getUserInput
     integer, parameter   :: fileUnit = 10
     character(len = 256) :: namelistFileName
     
     namelistFileName = getOneArgument()
     open (unit=fileUnit, file = trim(namelistFileName), status='old')
     read (fileUnit, nml = radiativeTransfer); rewind(10)
     read (fileUnit, nml = monteCarlo);        rewind(10)
     read (fileUnit, nml = algorithms);        rewind(10)
     read (fileUnit, nml = filenames);         rewind(10)
     read (fileUnit, nml = problemOptics);     rewind(10)
     read (fileUnit, nml = problemDomain)
     close (fileUnit)

     numIntensityAngles = count(abs(intensityMus) > 0)
     if(numIntensityAngles > 0) computeIntensity = .true.
     
  end subroutine getUserInput
 ! ---------------
  subroutine createDomain
    integer :: i
    real,    dimension(:, :, :), allocatable :: extinction, singleScatteringAlbedo
    integer, dimension(:, :, :), allocatable :: phaseFunctionIndex
    real,    dimension(:),       allocatable :: scatteringAngle, value
    type(phaseFunction)             :: phase
    type(phaseFunctionTable)        :: table
  
    allocate(            extinction(nX, nY, nLayers), &
             singleScatteringAlbedo(nX, nY, nLayers), &
                 phaseFunctionIndex(nX, nY, nLayers))
    !
    ! Define the domain
    !
    exampleDomain =                                                                &
      new_Domain(xPosition = domainSize/nX * (/ 0., (real(i), i = 1, nX) /),       &
                 yPosition = domainSize/nY * (/ 0., (real(i), i = 1, nY) /),       &
                 zPosition = physicalThickness/real(nLayers) *                     &
                                             (/ 0., (real(i), i = 1, nLayers) /) , &
                 status = status)
    call printStatus(status)
    
    !
    ! Compute a Henyey-Greenstein phase function or read one from a file
    !
    if(len_trim(phaseFunctionTableFile) > 0) then
      !
      ! Read the phase function table from the file
      !
      call read_PhaseFunctionTable(trim(phaseFunctionTableFile), &
                                   table = table, status = status) 
      phaseFunctionIndex(:, :, :) = phaseFunctionTableIndex
    else 
      !
      ! Use a Henyey-Greenstein phase function
      !
      if(useMoments) then
        ! 
        ! Phase function using Legendre moments 
        !
        phase = &
          new_PhaseFunction(g**(/ (i, i = 1, nLegendreCoefficients )/), status = status)
         call printStatus(status)
      
        table = &
          new_PhaseFunctionTable((/ phase /), key = (/ 1. /), status = status)
      else
        !
        !  Phase function using angle-values pairs for Henyey-Greenstein phase function
        !
        allocate(scatteringAngle(nAngles), value(nAngles))
        scatteringAngle(:) = (/ (real(i), i = 0, nAngles - 1) /)/real(nAngles - 1) * acos(-1.)
        value(:) = (1 - g**2) / (1 + g**2 - 2 * g * cos(scatteringAngle))**(3./2.)
        table = &
          new_PhaseFunctionTable(scatteringAngle, spread(value, dim = 2, ncopies = 1), &
                                 key = (/ 1. /), status = status)
        deallocate(scatteringAngle, value)
      end if
      call printStatus(status) 
      phaseFunctionIndex(:, :, :) = 1
    end if 
    if(.not. stateIsFailure(status)) call setStateToCompleteSuccess(status)

    !
    ! Add the clouds to the domain
    !
    extinction(:, :, :)         = opticalDepth/physicalThickness
    singleScatteringAlbedo(:, :, :) = ssa
    
    call addOpticalComponent(exampleDomain, "cloud",  &
                             extinction, singleScatteringAlbedo, &
                             phaseFunctionIndex, table, status = status)
    call printStatus(status)
    deallocate(extinction, singleScatteringAlbedo, phaseFunctionIndex)
    
    if(len_trim(domainFileName) > 0) then
      call write_Domain(exampleDomain, trim(domainFileName), status)
      call printStatus(status)
      print *, "Wrote domain to file ", trim(domainFileName)
    end if 
  end subroutine createDomain
end program computePlaneParallelRT
