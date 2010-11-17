! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

program backwardsMonteCarloDriver
  ! $Revision$, $Date$
  ! $URL$

  ! Backwards Monte Carlo monochromatic radiative transfer program.
  !   Computes radiative transfer for a set of trajectories originating inside the domain
  !   These can originate in random (flux) or specific (intensity) directions
  !   We sum the local estimates (as we would for an intensity calculation in a 
  !   forward Monte Carlo calculation) to get the contribution for a solar position. 
  !  An input "domain" file describes the optical properties for the medium. 
  
  ! 
  ! The input parameters are specified with a namelist file. 
  !
  ! Each output has the mean value and standard error of the mean
  ! over numBatches of photon batches (to estimate the Monte Carlo noise).
  !
  ! Assumptions: 
  !   Monochromatic solar radiative transfer.
  !   Lambertian surface reflection.
  !   Periodic horizontal boundary conditions.
  !   Uniform optical properties in each grid cell.

  !    Robert Pincus   University of Colorado     October 2010
  !      Adapted from monteCarloDriver.f95 of the I3RC Community Monte Carlo model 
  
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
                               backwardsRadiativeTransfer, reportResults
  use UserInterface,     only: printStatus, getOneArgument 

  implicit none

  ! Input parameters
  !   Radiative transfer 
  real                 :: solarFlux = 1., surfaceAlbedo = 0.
  real                 :: detectorX = 0., detectorY = 0., detectorZ = 0., &
                          detDeltaX = 0., detDeltaY = 0.,                 &
                          detectorMu = 0., detectorPhi = 0., detDeltaTheta = 0. 
  logical              :: detectorPointsUp = .true. 
  integer, parameter   :: maxNumSuns = 20
  real                 :: solarMus(maxNumSuns)  = 0., &
                          solarPhis(maxNumSuns) = 0.
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
  
  ! File names
  character(len=256)   :: domainFileName = ""
  character(len=256)   :: outputFluxFile = "", outputRadFile = "",  &
                          outputAbsProfFile = "", outputAbsVolumeFile = "", &
                          outputNetcdfFile = ""

  namelist /radiativeTransfer/ solarFlux, surfaceAlbedo,                  &
                               detectorX, detectorY, detectorZ,           &
                               detectorMu, detectorPhi, detectorPointsUp, &
                               solarMus, solarPhis
  namelist /monteCarlo/ numPhotonsPerBatch, numBatches, iseed, nPhaseIntervals
  namelist /algorithms/ useRayTracing, useRussianRoulette,                    &
                        useHybridPhaseFunsForIntenCalcs, hybridPhaseFunWidth, &
                        numOrdersOrigPhaseFunIntenCalcs,                      &
                        useRussianRouletteForIntensity, zetaMin
                        
  namelist /fileNames/ domainFileName, &
                       outputRadFile, outputFluxFile, &
                       outputNetcdfFile


   ! Local variables
  logical              :: computeFlux
  character(len=256)   :: namelistFileName
  integer              :: nX, nY, nZ
  integer              :: batch
  integer              :: numSuns
  real                 :: cpuTime0, cpuTime1, cpuTime2, cpuTimeTotal, cpuTimeSetup
  real, allocatable    :: xPosition(:), yPosition(:), zPosition(:)
  real, allocatable    :: Outputs(:, :, :), AccumOutputs(:, :, :, :)

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
  read (1, nml = fileNames)
  close (1)
  numSuns = count(abs(solarMus(:)) > 0.) 

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

   ! Set the surface albedo, table sizes, and maybe the solar directions (radiance directions) 
  call specifyParameters (mcIntegrator,                          &
                          surfaceAlbedo = surfaceAlbedo,         &
                          minInverseTableSize = nPhaseIntervals, &
                          status = status)
  call printStatus(status) 

  call specifyParameters (mcIntegrator, &
                          minForwardTableSize=nPhaseIntervals, &
                          intensityMus  = solarMus(1:numSuns), &
                          intensityPhis = solarPhis(1:numSuns), &
                          computeIntensity=.true., status=status)
  call printStatus(status) 

  !
  ! Make the algorithmic choices
  !
  call specifyParameters(mcIntegrator,                              &
                         useRayTracing      = useRayTracing,        &
                         useRussianRoulette = useRussianRoulette,   &
                         useHybridPhaseFunsForIntenCalcs =          &
                                 useHybridPhaseFunsForIntenCalcs,   &
                         hybridPhaseFunWidth = hybridPhaseFunWidth, &
                         numOrdersOrigPhaseFunIntenCalcs =          &
                                 numOrdersOrigPhaseFunIntenCalcs,   &
                         useRussianRouletteForIntensity =           &
                                 useRussianRouletteForIntensity,    &
                         zetaMin = zetaMin,                         &
                         limitIntensityContributions = .false.,     &
                         status = status)
  call printStatus(status) 

  allocate (Outputs(nX, nY, numSuns), AccumOutputs(nX, nY, numSuns, 2))
  AccumOutputs(:, :, :, :) = 0.0

  !
  ! Are we computing flux or intensity?   
  !
  computeFlux = abs(detectorMu) <= tiny(detectorMu) 
  
  ! --------------------------------------------------------------------------
  ! Compute radiative transfer with a trivial number of photons. 
  !   This checks to see if the integrator is properly set up before 
  !   running all the batches, and also allows the integrator to 
  !   do any internal pre-computation. 

  ! Seed the random number generator.
  randoms = new_RandomNumberSequence(seed = (/ iseed, 0 /) )

  ! Create a single photon
  call getPhotons(1)
  call printStatus(status)

  ! Compute the radiative transfer for a single photon - this lets the 
  !   integrator make various initial calculations and checks that the 
  !   full calculation is likel to proceed. 
  if(.not. isReady_Integrator (mcIntegrator)) stop 'Integrator is not ready.'
  call backwardsRadiativeTransfer(mcIntegrator, randoms, incomingPhotons, status)
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

    call getPhotons(numPhotonsPerBatch) 
    call printStatus(status)

    ! Now we compute the radiative transfer for this batch of photons. 
    call backwardsRadiativeTransfer(mcIntegrator, randoms, incomingPhotons, status)

     ! Get the radiative quantities:
    call reportResults(mcIntegrator, intensity = Outputs(:, :, :), status = status)
    AccumOutputs(:, :, :, 1) = AccumOutputs(:, :, :, 1) + Outputs(:, :, :)
    AccumOutputs(:, :, :, 2) = AccumOutputs(:, :, :, 2) + Outputs(:, :, :)**2

     ! Release the photon "illumination" object memory
    call finalize_PhotonStream (incomingPhotons)
    call printStatus(status)
  end do batches
  
  !
  ! Accumulate statistics from across all the processors
  !
  AccumOutputs(:, :, :, :) = sumAcrossProcesses(AccumOutputs)

  call synchronizeProcesses
  call cpu_time(cpuTime2)
  cpuTimeTotal = sumAcrossProcesses(cpuTime2 - cpuTime0)
  call finalizeProcesses

  if (MasterProc) print *, "Total CPU time (secs, approx): ", int(cpuTimeTotal)

  ! Calculate the mean and standard error of the radiative quantities from the two moments
  AccumOutputs(:, :, :, :) = solarFlux * AccumOutputs(:, :, :, :)/numBatches
  AccumOutputs(:, :, :, 2) = sqrt( max(0., AccumOutputs(:, :, :, 2) - AccumOutputs(:, :, :,1)**2) /(numBatches-1))

  if(MasterProc) then ! Write a single output file. 
  end if

  !
  ! Release all the memory. We should be able to finalize the derived types before we write 
  !   the output but this fails intermittently, so we want to be sure to get our results 
  !   before we take a chance on blowing up. 
  ! 
  deallocate (outputs, AccumOutputs)
  call finalize_RandomNumberSequence(randoms)
  call finalize_Integrator (mcIntegrator)

contains 
  subroutine getPhotons(n) 
    integer, intent(in) :: n

    if(computeFlux) then 
      if(detDeltaX > 0. .or. detDeltaY > 0.) then 
        incomingPhotons = new_PhotonStream ((xPosition(nX + 1) - detectorX) / (xPosition(nX + 1) - xPosition(1)), &
                                            (yPosition(nY + 1) - detectorY) / (yPosition(nY + 1) - yPosition(1)), &
                                            (zPosition(nZ + 1) - detectorZ) / (zPosition(nZ + 1) - zPosition(1)), &
                                            detectorPointsUp,                                                     &
                                            (xPosition(nX + 1) - detDeltaX) / (xPosition(nX + 1) - xPosition(1)), &
                                            (yPosition(nY + 1) - detDeltaY) / (yPosition(nY + 1) - yPosition(1)), &
                                            numberOfPhotons = n, randomNumbers = randoms, status=status)
      else
        incomingPhotons = new_PhotonStream ((xPosition(nX + 1) - detectorX) / (xPosition(nX + 1) - xPosition(1)), &
                                            (yPosition(nY + 1) - detectorY) / (yPosition(nY + 1) - yPosition(1)), &
                                            (zPosition(nZ + 1) - detectorZ) / (zPosition(nZ + 1) - zPosition(1)), &
                                            detectorPointsUp,                                                     &
                                            numberOfPhotons = n, randomNumbers = randoms, status=status)
      end if
    else
      if(detDeltaX > 0. .or. detDeltaY > 0. .or. detDeltaTheta > 0.) then 
        incomingPhotons = new_PhotonStream ((xPosition(nX + 1) - detectorX) / (xPosition(nX + 1) - xPosition(1)), &
                                            (yPosition(nY + 1) - detectorY) / (yPosition(nY + 1) - yPosition(1)), &
                                            (zPosition(nZ + 1) - detectorZ) / (zPosition(nZ + 1) - zPosition(1)), &
                                            detectorMu, detectorPhi,                                              &
                                            (xPosition(nX + 1) - detDeltaX) / (xPosition(nX + 1) - xPosition(1)), &
                                            (yPosition(nY + 1) - detDeltaY) / (yPosition(nY + 1) - yPosition(1)), &
                                            detDeltaTheta,                                                        &
                                            numberOfPhotons = n, randomNumbers = randoms, status=status)
      else
        incomingPhotons = new_PhotonStream ((xPosition(nX + 1) - detectorX) / (xPosition(nX + 1) - xPosition(1)), &
                                            (yPosition(nY + 1) - detectorY) / (yPosition(nY + 1) - yPosition(1)), &
                                            (zPosition(nZ + 1) - detectorZ) / (zPosition(nZ + 1) - zPosition(1)), &
                                            detectorMu, detectorPhi,                                              &
                                            numberOfPhotons = n, randomNumbers = randoms, status=status)
      end if
    end if 
  end subroutine getPhotons
! -------------------------------------------------------------------------------

end program backwardsMonteCarloDriver
