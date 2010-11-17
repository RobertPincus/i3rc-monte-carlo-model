! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

program monteCarloDriver
  ! $Revision$, $Date$
  ! $URL$

  ! Monte Carlo radiative transfer program that computes top and bottom
  ! of domain pixel level outgoing fluxes, domain average absorption
  ! profile, 3D absorption fields, and domain top radiances from an 
  ! input "domain" file describing the optical properties for the medium.  
  ! Does monochromatic solar radiative transfer only.
  ! 
  ! The input parameters are specified with a namelist file. 
  !
  ! Outputs four types of results: 
  !    Pixel level and domain average fluxes
  !    Pixel level radiances for the specified directions
  !    Domain average absorption profiles (flux/km)
  !    3D absorption field (flux/km)
  ! Each output has the mean value and standard error of the mean
  ! over numBatches of photon batches (to estimate the Monte Carlo noise).
  !
  ! Assumptions: 
  !   Monochromatic solar radiative transfer.
  !   Lambertian surface reflection.
  !   Periodic horizontal boundary conditions.
  !   Uniform optical properties in each grid cell.

  !    Frank Evans    University of Colorado     December 2005
  !      Modifications by Robert Pincus, Climate Diagnostics Center, January 2006
  !      Modifications by Robert Pincus and Frank Evans, June-July 2006
  !      Adapted for MPI by Robert Pincus, August 2007
  !      Merged with single processor version by Robert Pincus, January 2009

  ! Modules from the I3RC community Monte Carlo model framework
  use ErrorMessages,     only: ErrorMessage,   &
                               stateIsFailure, &
                               setStateToFailure, setStateToWarning, setStateToSuccess
  use MultipleProcesses, only: MasterProc, &
                               sumAcrossProcesses, initializeProcesses, finalizeProcesses, synchronizeProcesses
  use RandomNumbers,     only: RandomNumberSequence, &
                               new_RandomNumberSequence, finalize_RandomNumberSequence
  use opticalProperties, only: domain, &
                               read_Domain, getInfo_Domain, finalize_Domain
  use monteCarloIllumination, &
                         only: photonStream, &
                               new_PhotonStream, finalize_photonStream
  use monteCarloRadiativeTransfer, &
                         only: integrator, &
                               new_Integrator, specifyParameters, isReady_Integrator, finalize_Integrator,  &
                               computeRadiativeTransfer, reportResults
  use UserInterface,    only: printStatus, getOneArgument 

  implicit none

  ! Input parameters
  !   Radiative transfer 
  real                 :: solarFlux = 1., surfaceAlbedo = 0.
  real                 :: solarMu = 1., solarAzimuth = 0. 
  integer, parameter   :: maxNumRad = 20
  real                 :: intensityMus(maxNumRad)  = 0., &
                          intensityPhis(maxNumRad) = 0.
  !   Monte Carlo
  integer              :: numPhotonsPerBatch = 0, numBatches = 100, &
                          iseed = 10, nPhaseIntervals = 10001
                          
  !   Monte Carlo algorithmic choices 
  logical              :: useRayTracing = .true., useRussianRoulette = .true. 
  logical              :: useHybridPhaseFunsForIntenCalcs = .false. 
  real                 :: hybridPhaseFunWidth = 7. 
  integer              :: numOrdersOrigPhaseFunIntenCalcs = 0
  logical              :: useRussianRouletteForIntensity = .true. 
  real                 :: zetaMin = 0.3 
  logical              :: limitIntensityContributions = .false.
  real                 :: maxIntensityContribution    = 77.
  
  ! Control over output
  logical              :: reportVolumeAbsorption = .false., &
                          reportAbsorptionProfile = .false. 
  
  ! File names
  character(len=256)   :: domainFileName = ""
  character(len=256)   :: outputFluxFile = "", outputRadFile = "",  &
                          outputAbsProfFile = "", outputAbsVolumeFile = "", &
                          outputNetcdfFile = ""

  namelist /radiativeTransfer/ solarFlux, solarMu, solarAzimuth, surfaceAlbedo, &
                               intensityMus, intensityPhis
  namelist /monteCarlo/ numPhotonsPerBatch, numBatches, iseed, nPhaseIntervals
  namelist /algorithms/ useRayTracing, useRussianRoulette,                    &
                        useHybridPhaseFunsForIntenCalcs, hybridPhaseFunWidth, &
                        numOrdersOrigPhaseFunIntenCalcs,                      &
                        useRussianRouletteForIntensity, zetaMin,              &
                        limitIntensityContributions, maxIntensityContribution
                        
  namelist /output/ reportVolumeAbsorption, reportAbsorptionProfile
  
  namelist /fileNames/ domainFileName, &
                       outputRadFile, outputFluxFile, &
                       outputAbsProfFile, outputAbsVolumeFile, outputNetcdfFile


   ! Local variables
  character(len=256)   :: namelistFileName
  integer              :: nX, nY, nZ
  integer              :: batch
  integer              :: numRadDir
  logical              :: computeIntensity
  real                 :: cpuTime0, cpuTime1, cpuTime2, cpuTimeTotal, cpuTimeSetup
  real                 :: meanFluxUp, meanFluxDown, meanFluxAbsorbed
  real                 :: meanFluxUpStats(2), meanFluxDownStats(2), meanFluxAbsorbedStats(2)
  real, allocatable    :: xPosition(:), yPosition(:), zPosition(:)
  real, allocatable    :: fluxUp(:, :), fluxDown(:, :), fluxAbsorbed(:, :)
  real, allocatable    :: fluxUpStats(:, :, :), fluxDownStats(:, :, :), fluxAbsorbedStats(:, :, :)
  real, allocatable    :: absorbedProfile(:), absorbedProfileStats(:, :)
  real, allocatable    :: absorbedVolume(:, :, :), absorbedVolumeStats(:, :, :, :)
  real, allocatable    :: Radiance(:, :, :), RadianceStats(:, :, :, :)

   ! I3RC Monte Carlo code derived type variables
  type(domain)               :: thisDomain
  type(ErrorMessage)         :: status
  type(randomNumberSequence) :: randoms
  type(photonStream)         :: incomingPhotons
  type(integrator)           :: mcIntegrator

  ! Variables related to splitting up the job across processors
  integer            :: thisProc            ! ID OF CURRENT PROCESSOR; default 0
  integer            :: numProcs            ! TOTAL NUMBER OF PROCESSORS USED; default 1
  integer            :: batchesPerProcessor

  ! -----------------------------------------
  ! Start communications among multiple processes.
  ! 
  call initializeProcesses(numProcs, thisProc)

  ! -----------------------------------------
  ! Get the input variables from the namelist file
  !
  call cpu_time(cpuTime0)
  namelistFileName = getOneArgument()
  open (unit = 1, file = trim(namelistFileName), status='OLD')
  read (1, nml = radiativeTransfer); rewind(1)
  read (1, nml = monteCarlo);        rewind(1)
  read (1, nml = algorithms);        rewind(1)
  read (1, nml = output);            rewind(1)
  read (1, nml = fileNames)
  close (1)
  numRadDir = count(abs(intensityMus(:)) > 0.) 
  computeIntensity = numRadDir > 0 .and. &
                     (len_trim(outputRadFile) > 0 .or. len_trim(outputNetcdfFile) > 0)
  if(.not. computeIntensity) outputRadFile = ""

  ! -----------------------------------------
  !  Read the domain file
  !
  call read_Domain(domainFileName, thisDomain, status)
  call printStatus(status)
  call getInfo_Domain(thisDomain, numX = nx, numY = ny, numZ = nZ, status = status) 
  allocate(xPosition(nx+1), yPosition(ny+1), zPosition(nz+1))
  call getInfo_Domain(thisDomain,                                   &
                      xPosition = xPosition, yPosition = yPosition, &
                      zPosition = zPosition, status = status) 

  ! Set up the integrator object - the integrator makes copies of the 
  !   3D distribution of optical properties, so we can release the resources
  mcIntegrator = new_Integrator(thisDomain, status = status)
  call printStatus(status)
  call finalize_Domain(thisDomain)

   ! Set the surface albedo, table sizes, and maybe the radiance directions
  call specifyParameters (mcIntegrator,                          &
                          surfaceAlbedo = surfaceAlbedo,         &
                          minInverseTableSize = nPhaseIntervals, &
                          status = status)
  call printStatus(status) 

  if (computeIntensity) then
    call specifyParameters (mcIntegrator, &
                            minForwardTableSize=nPhaseIntervals, &
                            intensityMus=intensityMus(1:numRadDir), &
                            intensityPhis=intensityPhis(1:numRadDir), &
                            computeIntensity=computeIntensity, status=status)
    call printStatus(status) 
  endif

  !
  ! Make the algorithmic choices
  !
  call specifyParameters(mcIntegrator,                              &
                         useRayTracing      = useRayTracing,        &
                         useRussianRoulette = useRussianRoulette,   &
                         status = status)
  call printStatus(status) 
  
  !
  ! Algorithmic choices for intensity calculations
  !
  if (computeIntensity) then
    call specifyParameters(mcIntegrator,                            &
                         useHybridPhaseFunsForIntenCalcs =          &
                                 useHybridPhaseFunsForIntenCalcs,   &
                         hybridPhaseFunWidth = hybridPhaseFunWidth, &
                         numOrdersOrigPhaseFunIntenCalcs =          &
                                 numOrdersOrigPhaseFunIntenCalcs,   &
                         useRussianRouletteForIntensity =           &
                                 useRussianRouletteForIntensity,    &
                         zetaMin = zetaMin,                         &
                         limitIntensityContributions =              &
                                   limitIntensityContributions,     &
                         maxIntensityContribution =                 &
                                   maxIntensityContribution,        &
                         status = status)
    call printStatus(status) 
  end if

   ! Allocate and zero the arrays for radiative quantities and moments 
  allocate (fluxUp      (nX, nY), fluxUpStats      (nX, nY, 2))
  allocate (fluxDown    (nX, nY), fluxDownStats    (nX, nY, 2))
  allocate (fluxAbsorbed(nX, nY), fluxAbsorbedStats(nX, nY, 2))
  allocate (absorbedProfile(nZ), absorbedProfilestats(nZ, 2))
  allocate (absorbedVolume(nX, nY, nZ), absorbedVolumeStats(nX, nY, nZ, 2))
  meanFluxUpStats(:) = 0.0  ; meanFluxDownStats(:) = 0.0  ; meanFluxAbsorbedStats(:) = 0.0
  fluxUpStats(:, :, :) = 0.0  ; fluxDownStats(:, :, :) = 0.0  ; fluxAbsorbedStats(:, :, :) = 0.0
  absorbedProfilestats(:, :) = 0.0 ;  absorbedVolumeStats(:, :, :, :) = 0.0
  if (computeIntensity) then
    allocate (Radiance(nX, nY, numRadDir), RadianceStats(nX, nY, numRadDir, 2))
    RadianceStats(:, :, :, :) = 0.0
  endif  

  ! --------------------------------------------------------------------------
  ! Compute radiative transfer with a trivial number of photons. 
  !   This checks to see if the integrator is properly set up before 
  !   running all the batches, and also allows the integrator to 
  !   do any internal pre-computation. 

  ! Seed the random number generator.
  randoms = new_RandomNumberSequence(seed = (/ iseed, 0 /) )

   ! The initial direction and position of the photons are precomputed and 
   !   stored in an "illumination" object. 
  incomingPhotons = new_PhotonStream (solarMu, solarAzimuth, &
                                      numberOfPhotons = 1,   &
                                      randomNumbers = randoms, status=status)
  call printStatus(status)

  ! Now we compute the radiative transfer for a single photon 
  if(.not. isReady_Integrator (mcIntegrator)) stop 'Integrator is not ready.'
  call computeRadiativeTransfer (mcIntegrator, randoms, incomingPhotons, status)
  call printStatus(status) 
  call finalize_PhotonStream (incomingPhotons)

  call cpu_time(cpuTime1)
  call synchronizeProcesses
  cpuTimeSetup = sumAcrossProcesses(cpuTime1 - cpuTime0) 
  if (MasterProc) &
    print *, "Setup CPU time (secs, approx): ", int(cpuTimeSetup)
  ! --------------------------------------------------------------------------

  ! The  loop over batches is for estimating the uncertainty in the flux and
  !   radiance from the variance between numBatches independent calculations. 
  numBatches = max(numBatches,2)
  batchesPerProcessor = numBatches/numProcs
  ! If the number of batches doesn't divide among the processors evenly increase the 
  !   number until it does. 
  if(mod(numBatches, numProcs) /= 0) then 
    batchesPerProcessor = batchesPerProcessor + 1
    numBatches = batchesPerProcessor * numProcs
  end if 
  if (MasterProc) &
    print *, "Doing ", batchesPerProcessor, " batches on each of ", numProcs, " processors." 
  batches: do batch = thisProc*batchesPerProcessor + 1, thisProc*batchesPerProcessor + batchesPerProcessor
    ! Seed the random number generator.
    !   Variable randoms holds the state of the random number generator. 
    randoms = new_RandomNumberSequence(seed = (/ iseed, batch /) )

    ! The initial direction and position of the photons are precomputed and 
    !   stored in an "illumination" object. 
    incomingPhotons = new_PhotonStream (solarMu, solarAzimuth,                &
                                        numberOfPhotons = numPhotonsPerBatch, &
                                        randomNumbers = randoms, status = status)
    call printStatus(status)

    ! Now we compute the radiative transfer for this batch of photons. 
    call computeRadiativeTransfer (mcIntegrator, randoms, incomingPhotons, status)

     ! Get the radiative quantities:
     !   This particular integrator provides fluxes at the top and bottom 
     !   of the domain for both the domain mean and pixel level fluxes,
     !   the absorbed flux profile, 3D field of absorbed flux, and
     !   the pixel level radiances at top and/or bottom of domain.
    call reportResults (mcIntegrator, &
           meanFluxUp, meanFluxDown, meanFluxAbsorbed,       &
           fluxUp(:, :), fluxDown(:, :), fluxAbsorbed(:, :), &
           absorbedProfile(:), absorbedVolume(:, :, :), status = status)

     ! Accumulate the first and second moments of each quantity over the batches 
    meanFluxUpStats(1)       = meanFluxUpStats(1)       + meanFluxUp
    meanFluxUpStats(2)       = meanFluxUpStats(2)       + meanFluxUp**2
    meanFluxDownStats(1)     = meanFluxDownStats(1)     + meanFluxDown
    meanFluxDownStats(2)     = meanFluxDownStats(2)     + meanFluxDown**2
    meanFluxAbsorbedStats(1) = meanFluxAbsorbedStats(1) + meanFluxAbsorbed
    meanFluxAbsorbedStats(2) = meanFluxAbsorbedStats(2) + meanFluxAbsorbed**2
          fluxUpStats(:, :, 1) =       fluxUpStats(:, :, 1) + fluxUp(:, :)
          fluxUpStats(:, :, 2) =       fluxUpStats(:, :, 2) + fluxUp(:, :)**2
        fluxDownStats(:, :, 1) =     fluxDownStats(:, :, 1) + fluxDown(:, :)
        fluxDownStats(:, :, 2) =     fluxDownStats(:, :, 2) + fluxDown(:, :)**2
    fluxAbsorbedStats(:, :, 1) = fluxAbsorbedStats(:, :, 1) + fluxAbsorbed(:, :)
    fluxAbsorbedStats(:, :, 2) = fluxAbsorbedStats(:, :, 2) + fluxAbsorbed(:, :)**2
    absorbedProfileStats(:, 1) = absorbedProfileStats(:, 1) + absorbedProfile(:)
    absorbedProfileStats(:, 2) = absorbedProfileStats(:, 2) + absorbedProfile(:)**2
    absorbedVolumeStats(:, :, :, 1) = absorbedVolumeStats(:, :, :, 1) + absorbedVolume(:, :, :)
    absorbedVolumeStats(:, :, :, 2) = absorbedVolumeStats(:, :, :, 2) + absorbedVolume(:, :, :)**2

    if (computeIntensity) then
      call reportResults(mcIntegrator, intensity = Radiance(:, :, :), status = status)
      RadianceStats(:, :, :,1) = RadianceStats(:, :, :,1) + Radiance(:, :, :)
      RadianceStats(:, :, :,2) = RadianceStats(:, :, :,2) + Radiance(:, :, :)**2
    endif

     ! Release the photon "illumination" object memory
    call finalize_PhotonStream (incomingPhotons)
    call printStatus(status)
  end do batches
  
  !
  ! Accumulate statistics from across all the processors
  !
  ! Domain-mean fluxes
  
  meanFluxUpStats(:)       = sumAcrossProcesses(meanFluxUpStats)
  meanFluxDownStats(:)     = sumAcrossProcesses(meanFluxDownStats)
  meanFluxAbsorbedStats(:) = sumAcrossProcesses(meanFluxAbsorbedStats)

  ! Pixel-by-pixel fluxes
  fluxUpStats(:, :, :)       = sumAcrossProcesses(fluxUpStats)
  fluxDownStats(:, :, :)     = sumAcrossProcesses(fluxDownStats)
  fluxAbsorbedStats(:, :, :) = sumAcrossProcesses(fluxAbsorbedStats)
  
  ! Absorption (mean profile and cell-by-cell)
  absorbedProfileStats(:, :)      = sumAcrossProcesses(absorbedProfileStats)
  absorbedVolumeStats(:, :, :, :) = sumAcrossProcesses(absorbedVolumeStats)
  
  ! Radiance
  if (computeIntensity) &
    RadianceStats(:, :, :, :) = sumAcrossProcesses(RadianceStats)

  call synchronizeProcesses
  call cpu_time(cpuTime2)
  cpuTimeTotal = sumAcrossProcesses(cpuTime2 - cpuTime0)
  call finalizeProcesses

  if (MasterProc) print *, "Total CPU time (secs, approx): ", int(cpuTimeTotal)

   ! Calculate the mean and standard error of the radiative quantities from the two moments
  meanFluxUpStats(:)       = solarFlux * meanFluxUpStats(:)/numBatches
  meanFluxUpStats(2)       = sqrt( max(0.0, meanFluxUpStats(2) - meanFluxUpStats(1)**2) /(numBatches-1))
  meanFluxDownStats(:)     = solarFlux * meanFluxDownStats(:)/numBatches
  meanFluxDownStats(2)     = sqrt( max(0.0, meanFluxDownStats(2) - meanFluxDownStats(1)**2) /(numBatches-1))
  meanFluxAbsorbedStats(:) = solarFlux * meanFluxAbsorbedStats(:)/numBatches
  meanFluxAbsorbedStats(2) = sqrt( max(0.0, meanFluxAbsorbedStats(2) - meanFluxAbsorbedStats(1)**2) /(numBatches-1))
  fluxUpStats(:, :, :)       = solarFlux * fluxUpStats(:, :, :)/numBatches
  fluxUpStats(:, :, 2)       = sqrt( max(0.0, fluxUpStats(:, :,2) - fluxUpStats(:, :,1)**2) /(numBatches-1))
  fluxDownStats(:, :, :)     = solarFlux * fluxDownStats(:, :, :)/numBatches
  fluxDownStats(:, :, 2)     = sqrt( max(0.0, fluxDownStats(:, :,2) - fluxDownStats(:, :,1)**2) /(numBatches-1))
  fluxAbsorbedStats(:, :, :) = solarFlux * fluxAbsorbedStats(:, :, :)/numBatches
  fluxAbsorbedStats(:, :, 2) = sqrt( max(0.0, fluxAbsorbedStats(:, :,2) - fluxAbsorbedStats(:, :,1)**2) /(numBatches-1))
  absorbedProfileStats(:, :) = solarFlux * absorbedProfileStats(:, :)/numBatches
  absorbedProfileStats(:, 2) = sqrt( max(0.0, absorbedProfileStats(:,2) - absorbedProfileStats(:,1)**2) /(numBatches-1))
  absorbedVolumeStats(:, :, :, :) = solarFlux * absorbedVolumeStats(:, :, :, :)/numBatches
  absorbedVolumeStats(:, :, :, 2) = sqrt( max(0.0, absorbedVolumeStats(:, :, :, 2) - absorbedVolumeStats(:, :, :,1)**2) / &
                                    (numBatches-1))
  if (computeIntensity) then
    RadianceStats(:, :, :, :) = solarFlux * RadianceStats(:, :, :, :)/numBatches
    RadianceStats(:, :, :,2) = sqrt( max(0.0, RadianceStats(:, :, :,2) - RadianceStats(:, :, :,1)**2) /(numBatches-1))
  endif



  if(MasterProc) then ! Write a single output file. 
    if(any( (/ len_trim(outputFluxFile),      len_trim(outputAbsProfFile), &
               len_trim(outputAbsVolumeFile), len_trim(outputRadFile)      /) > 0)) then 
      call writeResults_ASCII(domainFileName,  numPhotonsPerBatch, numBatches,      &
                              useRayTracing, useRussianRoulette,                    &
                              useHybridPhaseFunsForIntenCalcs, hybridPhaseFunWidth, &
                              solarFlux, solarMu, solarAzimuth, surfaceAlbedo,      &
                              xPosition, yPosition, zPosition,                      &
                              outputFluxFile, meanFluxUpStats, meanFluxDownStats,   &
                              meanFluxAbsorbedStats,                                &
                              fluxUpStats, fluxDownStats, fluxAbsorbedStats,        &
                              outputAbsProfFile, absorbedProfileStats,              &
                              outputAbsVolumeFile, absorbedVolumeStats,             &
                              outputRadFile, intensityMus, intensityPhis, RadianceStats)
      print *, "Wrote ASCII results"  
    end if 
  
    if(len_trim(outputNetcdfFile) > 0) then
      if(computeIntensity) then 
        call writeResults_netcdf(domainFileName,  numPhotonsPerBatch, numBatches, &
                                 solarFlux, solarMu, solarAzimuth, surfaceAlbedo, &
                                 xPosition, yPosition, zPosition,                 &
                                 outputNetcdfFile,                                &
                                 fluxUpStats, fluxDownStats, fluxAbsorbedStats,   &
                                 absorbedProfileStats, absorbedVolumeStats,       &
                                 intensityMus, intensityPhis, RadianceStats)
      else
        call writeResults_netcdf(domainFileName,  numPhotonsPerBatch, numBatches, &
                                 solarFlux, solarMu, solarAzimuth, surfaceAlbedo, &
                                 xPosition, yPosition, zPosition,                 &
                                 outputNetcdfFile,                                &
                                 fluxUpStats, fluxDownStats, fluxAbsorbedStats,   &
                                 absorbedProfileStats, absorbedVolumeStats)
  
      end if 
      print *, "Wrote netcdf results"                         
    end if
  end if

  !
  ! Release all the memory. We should be able to finalize the derived types before we write 
  !   the output but this fails intermittently, so we want to be sure to get our results 
  !   before we take a chance on blowing up. 
  ! 
  deallocate (fluxUp, fluxUpStats, fluxDown, fluxDownStats, fluxAbsorbed, fluxAbsorbedStats)
  deallocate (absorbedProfile, absorbedProfileStats, absorbedVolume, absorbedVolumeStats)
  if (computeIntensity) deallocate (Radiance, RadianceStats)

  call finalize_RandomNumberSequence(randoms)
  call finalize_Integrator (mcIntegrator)

contains

! -------------------------------------------------------------------------------
  subroutine writeResults_ASCII(domainFileName,  numPhotonsPerBatch, numBatches,&
                          useRayTracing, useRussianRoulette,                    &
                          useHybridPhaseFunsForIntenCalcs, hybridPhaseFunWidth, &
                          solarFlux, solarMu, solarAzimuth, surfaceAlbedo,      &
                          xPosition, yPosition, zPosition,                      &
                          outputFluxFile, meanFluxUpStats, meanFluxDownStats,   &
                          meanFluxAbsorbedStats,                                &
                          fluxUpStats, fluxDownStats, fluxAbsorbedStats,        &
                          outputAbsProfFile, absorbedProfileStats,              &
                          outputAbsVolumeFile, absorbedVolumeStats,             &
                          outputRadFile, intensityMus, intensityPhis, RadianceStats)
    !
    ! Writes Monte Carlo results to ASCII files. 
    !   Fluxes, absorption profiles, 3D absorption, and intensities are written to separate files.
    !
    ! Variables describing the problem 
    character(len=*),   intent(in) :: domainFileName
    integer,            intent(in) :: numPhotonsPerBatch, numBatches
    logical,            intent(in) :: useRayTracing, useRussianRoulette, useHybridPhaseFunsForIntenCalcs
    real,               intent(in) :: hybridPhaseFunWidth
    real,               intent(in) :: solarFlux, solarMu, solarAzimuth, surfaceAlbedo
    real, dimension(:), intent(in) :: xPosition, yPosition, zPosition
    ! Flux variables
    character(len = *),       intent(in) :: outputFluxFile
    real, dimension(:),       intent(in) :: meanFluxUpStats, meanFluxDownStats, meanFluxAbsorbedStats
    real, dimension(:, :, :), intent(in) :: fluxUpStats, fluxDownStats, fluxAbsorbedStats
    ! Absorption variables
    character(len = *),       intent(in) :: outputAbsProfFile
    real, dimension(:, :),    intent(in) :: absorbedProfileStats
    character(len = *),       intent(in) :: outputAbsVolumeFile
    real, dimension(:, :, :, :), &
                              intent(in) :: absorbedVolumeStats
    ! Intensity variable
    character(len = *),       intent(in) :: outputRadFile
    real, dimension(:),       intent(in) :: intensityMus, intensityPhis
    real, dimension(:, :, :, :), &
                              intent(in) :: RadianceStats

    !
    ! Local variables
    !
    integer  :: i, j, k, nx, ny, nz, numRadDir

    nx = size(xPosition) - 1; ny = size(yPosition) - 1; nz = size(zPosition) - 1

    ! -----------------------------------------
    !  Output the flux results to an ASCII file
    !
    if(len_trim(outputFluxFile) > 0) then 
      open (unit=2, file=outputFluxFile, status='unknown')
      write (2,'(A)') '!   I3RC Monte Carlo 3D Solar Radiative Transfer: Flux'
      write (2,'(A,A60)') '!  Property_File=', domainFileName
      write (2,'(A,I10)')  '!  Num_Photons=', numPhotonsPerBatch * numBatches
      write (2,'(A,L1,A,L1)') '!  PhotonTracing=', useRayTracing, &
                              '    Russian_Roulette=',useRussianRoulette
      write (2,'(A,L1,A,F5.2)') '!  Hybrid_Phase_Func_for_Radiance=',useHybridPhaseFunsForIntenCalcs, &
                                '   Gaussian_Phase_Func_Width_deg=',hybridPhaseFunWidth
      write (2,'(A,E13.6,A,F10.7,A,F7.3)') '!  Solar_Flux=', SolarFlux, &
                          '   Solar_Mu=', SolarMu, '   Solar_Phi=', SolarAzimuth
      write (2,'(A,F7.4)') '!  Lambertian_Surface_Albedo=',surfaceAlbedo
      write (2,'(A)')  '!  Output_Type= Pixel Flux'
      write (2,'(A,F7.3,A,F7.3)') '!  Upwelling_Level=', zPosition(nZ+1), &
                                  '   Downwelling_level=', zPosition(1)
      write (2,'(A)') '!   X      Y           Flux_Up             Flux_Down            Flux_Absorbed '
      write (2,'(A)') '!                  Mean     StdErr       Mean     StdErr       Mean     StdErr'
      write (2,'(A14,3(1X,2(1X,F9.4)))') '!  Average:   ', meanFluxUpStats(1:2), &
                             meanFluxDownStats(1:2), meanFluxAbsorbedStats(1:2)
      do j = 1, nY
        do i = 1, nX  
          write (2,'(2(F7.3),3(1X,2(1X,F9.4)))') &
               sum(xPosition(i:i+1))/2., sum(yPosition(j:j+1))/2., &
               fluxUpStats(i,j,1:2), fluxDownStats(i,j,1:2), fluxAbsorbedStats(i,j,1:2)
        enddo
      enddo  
      close (2)
      end if 


    ! -----------------------------------------
    !  Output the absorption profile results to a file
    !
    if(len_trim(outputAbsProfFile) > 0) then 
      open (unit=2, file=outputAbsProfFile, status='unknown')
      write (2,'(A)') '!   I3RC Monte Carlo 3D Solar Radiative Transfer: Absorption Profile'
      write (2,'(A,A60)') '!  Property_File=', domainFileName
      write (2,'(A,I10)')  '!  Num_Photons=', numPhotonsPerBatch * numBatches
      write (2,'(A,L1,A,L1)') '!  PhotonTracing=', useRayTracing, &
                              '    Russian_Roulette=',useRussianRoulette
      write (2,'(A,L1,A,F5.2)') '!  Hybrid_Phase_Func_for_Radiance=',useHybridPhaseFunsForIntenCalcs, &
                                '   Gaussian_Phase_Func_Width_deg=',hybridPhaseFunWidth
      write (2,'(A,E13.6,A,F10.7,A,F7.3)') '!  Solar_Flux=', SolarFlux, &
                          '   Solar_Mu=', SolarMu, '   Solar_Phi=', SolarAzimuth
      write (2,'(A,F7.4)') '!  Lambertian_Surface_Albedo=',surfaceAlbedo
      write (2,'(A)')  '!  Output_Type= Absorption Profile'
      write (2,'(A)') '!   Z    Absorbed_Flux (flux/km) '
      write (2,'(A)') '!          Mean     StdErr '
      do k = 1, nZ
        write (2,'(F7.3,1X,2(1X,F9.4))') 0.5*(zPosition(k)+zPosition(k+1)), absorbedProfileStats(k,1:2) 
      enddo  
      close (2)
    end if


    ! -----------------------------------------
    !  Output the volume absorption results to a file
    !
    if(len_trim(outputAbsVolumeFile) > 0) then 
      open (unit=2, file=outputAbsVolumeFile, status='unknown')
      write (2,'(A)') '!   I3RC Monte Carlo 3D Solar Radiative Transfer: 3D Absorption Field'
      write (2,'(A,A60)') '!  Property_File=', domainFileName
      write (2,'(A,I10)')  '!  Num_Photons=', numPhotonsPerBatch * numBatches
      write (2,'(A,L1,A,L1)') '!  PhotonTracing=', useRayTracing, &
                              '    Russian_Roulette=',useRussianRoulette
      write (2,'(A,L1,A,F5.2)') '!  Hybrid_Phase_Func_for_Radiance=',useHybridPhaseFunsForIntenCalcs, &
                                '   Gaussian_Phase_Func_Width_deg=',hybridPhaseFunWidth
      write (2,'(A,E13.6,A,F10.7,A,F7.3)') '!  Solar_Flux=', SolarFlux, &
                          '   Solar_Mu=', SolarMu, '   Solar_Phi=', SolarAzimuth
      write (2,'(A,F7.4)') '!  Lambertian_Surface_Albedo=',surfaceAlbedo
      write (2,'(A)')  '!  Output_Type= Volume Absorption '
      write (2,'(A)') '!    X       Y        Z       Absorbed_Flux (flux/km)'
      write (2,'(A)') '!                               Mean     StdErr '
      do i = 1, nX  
        do j = 1, nY
          do k = 1, nZ
            write (2,'(3(F7.3,1X),2(1X,F9.4))') &
              sum(xPosition(i:i+1))/2., sum(yPosition(j:j+1))/2., sum(zPosition(k:k+1))/2., &
              absorbedVolumeStats(i,j,k,1:2) 
          enddo  
        enddo  
      enddo  
      close (2)
    end if

    ! -----------------------------------------
    !  Output the radiance results to a file
    !
    if (len_trim(outputRadFile) > 0) then
      numRadDir = count(abs(intensityMus(:)) > 0.)
      open (unit=2, file=outputRadFile, status='unknown')
      write (2,'(A)') '!   I3RC Monte Carlo 3D Solar Radiative Transfer: Radiance'
      write (2,'(A,A60)') '!  Property_File=', domainFileName
      write (2,'(A,I10)')  '!  Num_Photons=', numPhotonsPerBatch * numBatches
      write (2,'(A,L1,A,L1)') '!  PhotonTracing=', useRayTracing, &
                              '    Russian_Roulette=',useRussianRoulette
      write (2,'(A,L1,A,F5.2)') '!  Hybrid_Phase_Func_for_Radiance=',useHybridPhaseFunsForIntenCalcs, &
                                '   Gaussian_Phase_Func_Width_deg=',hybridPhaseFunWidth
      write (2,'(A,L1,A,F5.2)') '!  Intensity_uses_Russian_Roulette=', useRussianRouletteForIntensity, &
                                '   Intensity_Russian_Roulette_zeta_min=',zetaMin
      write (2,'(A,L1,A,F5.2)') '!  limited_intensity_contributions=', limitIntensityContributions, &
                                '   max_intensity_contribution=',maxIntensityContribution
      write (2,'(A,E13.6,A,F10.7,A,F7.3)') '!  Solar_Flux=', SolarFlux, &
                        '   Solar_Mu=', SolarMu, '   Solar_Phi=', SolarAzimuth
      write (2,'(A,F7.4)') '!  Lambertian_Surface_Albedo=',surfaceAlbedo
      write (2,'(A)')  '!  Output_Type= Pixel Radiance'
      write (2,'(A,F7.3,3(A,I4))') '!  RADIANCE AT Z=',  zPosition(nZ+1), &
              '   NXO=',nX, '   NYO=',nY, '   NDIR=',numRadDir
      write (2,'(A)') '!   X      Y         Radiance (Mean, StdErr)'
      do k = 1, numRadDir
        WRITE (2,'(A,1X,F8.5,1X,F6.2,2X,A)') &
            '! ', intensityMus(k), intensityPhis(k), '<- (mu,phi)'
        do j = 1, nY
          do i = 1, nX  
            write (2,'(2(F7.3),2(1X,F9.4))') sum(xPosition(i:i+1))/2., sum(yPosition(j:j+1))/2., &
               RadianceStats(i,j,k,1:2)
          enddo
        enddo  
      enddo  
      close (2)
    endif
  end subroutine writeResults_ASCII


! -------------------------------------------------------------------------------
  subroutine writeResults_netcdf(domainFileName,  numPhotonsPerBatch, numBatches, &
                                 solarFlux, solarMu, solarAzimuth, surfaceAlbedo, &
                                 xPosition, yPosition, zPosition,                 &
                                 outputFileName,                                  &
                                 fluxUpStats, fluxDownStats, fluxAbsorbedStats,   &
                                 absorbedProfileStats, absorbedVolumeStats,       &
                                 intensityMus, intensityPhis, RadianceStats)
    use netcdf
    !
    ! Writes Monte Carlo results to a netCDF file. 
    !   Writes top-of-dmain fluxes, absorption profiles, and intensities 
    !   Monte Carlo noise for each quantitiy is estimated as the std. dev. across the nBatches. 
    !
    ! Variables describing the problem 
    character(len = *), intent(in) :: domainFileName
    character(len = *), intent(in) :: outputFileName
    integer,            intent(in) :: numPhotonsPerBatch, numBatches
    real,               intent(in) :: solarFlux, solarMu, solarAzimuth, surfaceAlbedo
    ! The position variables mark the edges of the cells - we'll report the 
    !   results at the midpoints. 
    real, dimension(:), intent(in) :: xPosition, yPosition, zPosition

    ! Flux variables
    real, dimension(:, :, :),     intent(in) :: fluxUpStats, fluxDownStats, fluxAbsorbedStats
    ! Absorption variables
    real, dimension(:, :),        intent(in) :: absorbedProfileStats
    real, dimension(:, :, :, :),  intent(in) :: absorbedVolumeStats
    ! Intensity variable
    real, dimension(:), optional, intent(in) :: intensityMus, intensityPhis
    real, dimension(:, :, :, :), &
                        optional, intent(in) :: RadianceStats

    !
    ! Local variables
    !
    logical                :: computeIntensity
    integer, dimension(32) :: ncStatus = nf90_NoErr
    integer                :: ncFileId, xDimId, yDimId, zDimId, dirDimId, ncVarId, &
                              numX, numY, numZ, numIntensityDirs
    ! ---------------------------------------------------------------------------
    numX = size(xPosition) - 1; numY = size(yPosition) - 1; numZ = size(zPosition) - 1
    computeIntensity = present(intensityMus) .and. present(intensityPhis) .and. &
                       present(RadianceStats)

    ncStatus( 1) = nf90_create(trim(outputFileName), nf90_clobber, ncFileId)

    ! Store the simulation parameters as global attributes
    !
    ncStatus( 5) = nf90_put_att(ncFileId, NF90_Global, "description", &
                                                       "Output from I3RC Community Monte Carlo Model")
    ncStatus( 6) = nf90_put_att(ncFileId, NF90_Global, "Domain_filename", trim(domainFileName))
    ncStatus( 7) = nf90_put_att(ncFileId, NF90_Global, "Surface_albedo", surfaceAlbedo)
    ncStatus( 8) = nf90_put_att(ncFileId, NF90_Global, "Total_number_of_photons", &
                                                       numPhotonsPerBatch * numBatches)
    ncStatus( 9) = nf90_put_att(ncFileId, NF90_Global, "Number_of_batches", numBatches)
    ncStatus(10) = nf90_put_att(ncFileId, NF90_Global, "Solar_flux", SolarFlux)
    ncStatus(11) = nf90_put_att(ncFileId, NF90_Global, "Solar_mu", SolarMu)
    ncStatus(12) = nf90_put_att(ncFileId, NF90_Global, "Solar_phi", SolarAzimuth)

    ncStatus(13) = nf90_put_att(ncFileId, NF90_Global, "Random_number_seed", iseed)
    ncStatus(14) = nf90_put_att(ncFileId, NF90_Global, "Phase_function_table_sizes", &
                                                       nPhaseIntervals)
    if(useRayTracing) then
      ncStatus(15) = nf90_put_att(ncFileId, NF90_Global, "Algorithm", "Ray_tracing")
    else
      ncStatus(15) = nf90_put_att(ncFileId, NF90_Global, "Algorithm", "Max_cross_section")
    end if

    if(useHybridPhaseFunsForIntenCalcs) then 
      ncStatus(16) = nf90_put_att(ncFileId, NF90_Global, "Intensity_uses_hyrbid_phase_functions", 1)
      ncStatus(17) = nf90_put_att(ncFileId, NF90_Global, "Hybrid_phase_function_width", &
                                                          hybridPhaseFunWidth)
    else
      ncStatus(16) = nf90_put_att(ncFileId, NF90_Global, "Intensity_uses_hyrbid_phase_functions", 0)
      ncStatus(17) = nf90_put_att(ncFileId, NF90_Global, "Hybrid_phase_function_width", 0.)
    end if 

    if(useRussianRouletteForIntensity) then 
      ncStatus(18) = nf90_put_att(ncFileId, NF90_Global, "Intensity_uses_Russian_roulette", 1)
      ncStatus(19) = nf90_put_att(ncFileId, NF90_Global, "Intensity_Russian_roulette_zeta_min",  zetaMin)
    else
      ncStatus(18) = nf90_put_att(ncFileId, NF90_Global, "Intensity_uses_Russian_roulette", 0)
      ncStatus(19) = nf90_put_att(ncFileId, NF90_Global, "Intensity_Russian_roulette_zeta_min", 0.)
    end if 
    
    if(limitIntensityContributions) then 
      ncStatus(18) = nf90_put_att(ncFileId, NF90_Global, "limited_intensity_contributions", 1)
      ncStatus(19) = nf90_put_att(ncFileId, NF90_Global, "max_intensity_contribution",  maxIntensityContribution)
    else
      ncStatus(18) = nf90_put_att(ncFileId, NF90_Global, "limited_intensity_contributions", 0)
      ncStatus(19) = nf90_put_att(ncFileId, NF90_Global, "max_intensity_contribution", 0.)
    end if 
    ncStatus(20) = nf90_put_att(ncFileId, NF90_Global, "Cpu_time_total", CpuTimeTotal)
    ncStatus(21) = nf90_put_att(ncFileId, NF90_Global, "Cpu_time_setup", CpuTimeSetup)
    ncStatus(22) = nf90_put_att(ncFileId, NF90_Global, "Number_of_processors_used", numProcs)

    if(any(ncStatus(:) /= nf90_NoErr)) print *, "Attributes ", ncStatus(:19)

    !
    ! Dimensions
    !
    ncStatus( 1) = nf90_def_dim(ncFileId, "x", numX, xDimId)
    ncStatus( 2) = nf90_def_dim(ncFileId, "y", numY, yDimId)
    if(reportAbsorptionProfile .or. reportVolumeAbsorption) & 
      ncStatus( 3) = nf90_def_dim(ncFileId, "z", numZ, zDimId)
    !
    ! Dimension variables
    !
    ncStatus( 4) = nf90_def_var(ncFileId, "x", nf90_float, xDimId, ncVarId)
    ncStatus( 5) = nf90_def_var(ncFileId, "y", nf90_float, yDimId, ncVarId)
    if(reportAbsorptionProfile .or. reportVolumeAbsorption) & 
      ncStatus( 6) = nf90_def_var(ncFileId, "z", nf90_float, zDimId, ncVarId)
    !
    ! Flux variables
    !
    ncStatus( 7) = nf90_def_var(ncFileId, "fluxUp", &
                                nf90_float, (/ xDimId, yDimId /), ncVarId)
    ncStatus( 8) = nf90_def_var(ncFileId, "fluxDown", &
                                nf90_float, (/ xDimId, yDimId /), ncVarId)
    ncStatus( 9) = nf90_def_var(ncFileId, "fluxAbsorbed", &
                                nf90_float, (/ xDimId, yDimId /), ncVarId)
    ncStatus(10) = nf90_def_var(ncFileId, "fluxUp_StdErr", &
                                nf90_float, (/ xDimId, yDimId /), ncVarId)
    ncStatus(11) = nf90_def_var(ncFileId, "fluxDown_StdErr", &
                                nf90_float, (/ xDimId, yDimId /), ncVarId)
    ncStatus(12) = nf90_def_var(ncFileId, "fluxAbsorbed_StdErr", &
                                nf90_float, (/ xDimId, yDimId /), ncVarId)
    !
    ! Absorption profile
    !
    if(reportAbsorptionProfile) then 
      ncStatus(13) = nf90_def_var(ncFileId, "absorptionProfile", &
                                  nf90_float, zDimId, ncVarId)
      ncStatus(14) = nf90_def_var(ncFileId, "absorptionProfile_StdErr", &
                                  nf90_float, zDimId, ncVarId)
    end if 
    
    !
    ! Volume absorption
    !
    if(reportVolumeAbsorption) then 
      ncStatus(15) = nf90_def_var(ncFileId, "absorbedVolume", &
                                  nf90_float, (/ xDimId, yDimId, zDimID /), ncVarId)
      ncStatus(16) = nf90_def_var(ncFileId, "absorbedVolume_StdErr", &
                                  nf90_float, (/ xDimId, yDimId, zDimId /), ncVarId)
    end if 
    
    !
    ! Intensity
    !
    if(computeIntensity) then
      numIntensityDirs = size(RadianceStats, 3)
      ncStatus(17) = nf90_def_dim(ncFileId, "direction",     numIntensityDirs,      dirDimId)
      ncStatus(19) = nf90_def_var(ncFileId, "intensityMus",  nf90_float, dirDimId, ncVarId)
      ncStatus(19) = nf90_def_var(ncFileId, "intensityPhis", nf90_float, dirDimId, ncVarId)
      ncStatus(20) = nf90_def_var(ncFileId, "intensity", &
                                  nf90_float, (/ xDimId, yDimId, dirDimId /), ncVarId)
      ncStatus(21) = nf90_def_var(ncFileId, "intensity_StdErr", &
                                  nf90_float, (/ xDimId, yDimId, dirDimId /), ncVarId)
    end if

    ncStatus(22) = nf90_EndDef(ncFileId)
    if(any(ncStatus(:) /= nf90_NoErr)) print *, ncStatus(:22)
    ! ----------------------
    ! File is set up - now write each variable
    ! 
    ! Dimension variables
    !
    ncStatus( 1) = nf90_inq_varid(ncFileId, "x", ncVarId)
    ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, (xPosition(:numX) + xPosition(2:))/2)
    ncStatus( 3) = nf90_inq_varid(ncFileId, "y", ncVarId)
    ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, (yPosition(:numY) + yPosition(2:))/2)
    if(reportAbsorptionProfile .or. reportVolumeAbsorption) then 
      ncStatus( 5) = nf90_inq_varid(ncFileId, "z", ncVarId)
      ncStatus( 6) = nf90_put_var(ncFileId, ncVarId, (zPosition(:numZ) + ZPosition(2:))/2)
    end if 
    if(any(ncStatus(:) /= nf90_NoErr)) print *, "Position", ncStatus(:6)
    ncstatus(:) = nf90_NoErr
    !
    ! Upward flux
    !
    ncStatus( 1) = nf90_inq_varid(ncFileId, "fluxUp", ncVarId)
    ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, fluxUpStats(:, :,1))
    ncStatus( 3) = nf90_inq_varid(ncFileId, "fluxUp_StdErr", ncVarId)
    ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, fluxUpStats(:, :,2))
    if(any(ncStatus(:) /= nf90_NoErr)) print *, "Flux up", ncStatus(:4)

    !
    ! Downward flux
    !
    ncStatus( 1) = nf90_inq_varid(ncFileId, "fluxDown", ncVarId)
    ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, fluxDownStats(:, :,1))
    ncStatus( 3) = nf90_inq_varid(ncFileId, "fluxDown_StdErr", ncVarId)
    ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, fluxDownStats(:, :,2))
    if(any(ncStatus(:) /= nf90_NoErr)) print *, "Flux down", ncStatus(:4)

    !
    ! Absorbed flux
    !
    ncStatus( 1) = nf90_inq_varid(ncFileId, "fluxAbsorbed", ncVarId)
    ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, fluxAbsorbedStats(:, :,1))
    ncStatus( 3) = nf90_inq_varid(ncFileId, "fluxAbsorbed_StdErr", ncVarId)
    ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, fluxAbsorbedStats(:, :,2))
    if(any(ncStatus(:) /= nf90_NoErr)) print *, "Flux absorbed", ncStatus(:4)

    !
    ! Absorption profile
    !
    if(reportAbsorptionProfile) then 
      ncStatus( 1) = nf90_inq_varid(ncFileId, "absorptionProfile", ncVarId)
      ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, absorbedProfileStats(:,1))
      ncStatus( 3) = nf90_inq_varid(ncFileId, "absorptionProfile_StdErr", ncVarId)
      ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, absorbedProfileStats(:,2))
      if(any(ncStatus(:) /= nf90_NoErr)) print *, "Absorption profile", ncStatus(:4)
    end if 
    
    !
    ! Volume absorption
    !
    if(reportVolumeAbsorption) then 
      ncStatus( 1) = nf90_inq_varid(ncFileId, "absorbedVolume", ncVarId)
      ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, absorbedVolumeStats(:, :, :,1))
      ncStatus( 3) = nf90_inq_varid(ncFileId, "absorbedVolume_StdErr", ncVarId)
      ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, absorbedVolumeStats(:, :, :,2))
      if(any(ncStatus(:) /= nf90_NoErr)) print *, "Flux convergence", ncStatus(:4)
    end if 
    
    !
    ! Intensity
    !
    if(computeIntensity) then
      ncStatus( 1) = nf90_inq_varid(ncFileId, "intensityMus", ncVarId)
      ncStatus( 2) = nf90_put_var(ncFileId, ncVarId, intensityMus(:numIntensityDirs))
      ncStatus( 3) = nf90_inq_varid(ncFileId, "intensityPhis", ncVarId)
      ncStatus( 4) = nf90_put_var(ncFileId, ncVarId, intensityPhis(:numIntensityDirs))
      ncStatus( 5) = nf90_inq_varid(ncFileId, "intensity", ncVarId)
      ncStatus( 6) = nf90_put_var(ncFileId, ncVarId, RadianceStats(:, :, :,1))
      ncStatus( 7) = nf90_inq_varid(ncFileId, "intensity_StdErr", ncVarId)
      ncStatus( 8) = nf90_put_var(ncFileId, ncVarId, RadianceStats(:, :, :,2))
      if(any(ncStatus(:) /= nf90_NoErr)) print *, "Intensity", ncStatus(:8)
    end if

    ncStatus( 1) = nf90_close(ncFileId)
    if(any(ncStatus(:) /= nf90_NoErr)) print *, ncStatus(:1)

  end subroutine writeResults_netcdf
! -------------------------------------------------------------------------------
end program monteCarloDriver
