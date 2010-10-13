! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision$, $Date$
! $URL$
module monteCarloRadiativeTransfer
  use CharacterUtils,    only: intToChar
  use ErrorMessages,     only: ErrorMessage,   &
                               stateIsFailure, &
                               setStateToFailure, setStateToWarning, setStateToSuccess, setStateToCompleteSuccess
  use RandomNumbers,     only: RandomNumberSequence, &
                               getRandomReal
  use numericUtilities,  only: findIndex
  use scatteringPhaseFunctions, &
                         only: phaseFunctionTable, &
                               getInfo_phaseFunctionTable, copy_PhaseFunctionTable, &
                               getPhaseFunctionValues, finalize_phaseFunctionTable
  use inversePhaseFunctions, &
                         only: computeInversePhaseFuncTable
  use opticalProperties, only: domain, &
                               read_Domain, getInfo_Domain, getOpticalPropertiesByComponent
  use surfaceProperties, only: surfaceDescription, &
                               copy_surfaceDescription, finalize_surfaceDescription, &
                               isReady_surfaceDescription, computeSurfaceReflectance
  use monteCarloIllumination, &
                       only: photonStream, &
                             getNextPhoton, morePhotonsExist
  implicit none
  private
  
  !------------------------------------------------------------------------------------------
  ! Constants
  !------------------------------------------------------------------------------------------
  integer, parameter :: defaultMinForwardTableSize = 9001, &
                        defaultMinInverseTableSize = 9001
  real,    parameter :: defaultHybridPhaseFunWidth =  7., & 
                        maxHybridPhaseFunWidth     = 30.
  integer, parameter :: defOrdersOrigPhaseFunIntenCalcs = 0
  real,    parameter :: defaultZetaMin             =  0.3
  real,    parameter :: defaultMaxIntensityContrib =  huge(1.) 
  real,    parameter :: Pi = 3.14159265358979312
  
  !------------------------------------------------------------------------------------------
  ! Type (object) definitions
  !------------------------------------------------------------------------------------------
  ! The public type
  !
  type integrator
    private
    ! Shortcuts - the status of the integrator could be obtained by looking to see if 
    !   arrays have been allocated, etc., but we also keep track with these flags
    logical                                 :: readyToCompute   = .false., &
                                               computeIntensity = .false. 
    ! Integrator specific parameters
    integer                                 :: minForwardTableSize = defaultMinForwardTableSize, &
                                               minInverseTableSize = defaultMinInverseTableSize
    ! -------------------------------------------------------------------------=                                           
    ! Algorithmic choices - these can be modified through specifyParameters()
    !
    ! Use ray-tracing? The alternative is max cross-section
    logical                                 :: useRayTracing = .true. 
    ! Use Russian roulette? 
    logical                                 :: useRussianRoulette = .true. 
    real                                    :: RussianRouletteW = 1. 
        
    ! -------------------------------------------------------------------------=                                           
    ! The atmosphere and surface 
    real                                    :: surfaceAlbedo = 0.  
    logical                                 :: xyRegularlySpaced = .false., &
                                                zRegularlySpaced = .false. 
    real                                    :: deltaX = 0., deltaY = 0. , deltaZ = 0., &
                                               x0 = 0., y0= 0., z0 = 0.
    real,    dimension(:),          pointer :: xPosition => null() 
    real,    dimension(:),          pointer :: yPosition => null()
    real,    dimension(:),          pointer :: zPosition => null()
    real,    dimension(:, :, :),    pointer :: totalExt           => null()
    real,    dimension(:, :, :, :), pointer :: cumulativeExt      => null()
    real,    dimension(:, :, :, :), pointer :: ssa                => null()
    integer, dimension(:, :, :, :), pointer :: phaseFunctionIndex => null()
    ! Surface reflection BDRF
    logical                          :: useSurfaceBDRF = .false.
    type(surfaceDescription)         :: surfaceBDRF

    !
    ! We store the original forward phase function tables even though calculations
    !   inside the module use the matrix representations below. The originals don't 
    !   take much rooms (a few Mb at most) and this allows us to recompute them to 
    !   arbitrary accuracy at any time. 
    !   
    type(phaseFunctionTable), &
             dimension(:),    pointer :: forwardTables => null()
    
    ! We store tabulated phase function and inverse (cumulative) phase functions in
    !   two-D arrays, but these arrays can be different sizes for different components. 
    !   We define a derived type to represent each 2D array and make a vector of these 
    !   (one matrix for each component). 
    !
    type(matrix), &
          dimension(:),       pointer :: tabulatedPhaseFunctions => null()
    type(matrix), &
          dimension(:),       pointer :: tabulatedOrigPhaseFunctions => null()
    type(matrix), &
          dimension(:),       pointer :: inversePhaseFunctions => null()
    ! -------------------------------------------------------------------------=                                           
    !
    ! Direction cosines at which to compute intensity
    !
    real, dimension(:, :),    pointer :: intensityDirections => null()
    
    ! -------------------------------------------------------------------------= 
    ! Variables for variance reduction for intensity calculations

    ! Build hybrid phase functions for intensity calculations for at least 
    !   one component

    logical :: useHybridPhaseFunsForIntenCalcs = .false. 
    real    :: hybridPhaseFunWidth             = defaultHybridPhaseFunWidth
    integer :: numOrdersOrigPhaseFunIntenCalcs = defOrdersOrigPhaseFunIntenCalcs

    ! Russian roulette for intensity 
    logical :: useRussianRouletteForIntensity = .False. 
    real    :: zetaMin                        = defaultZetaMin
    
    ! Redistribute large intenstity contributions in space
    logical                              :: limitIntensityContributions = .false. 
    real                                 :: maxIntensityContribution    = &
                                                        defaultMaxIntensityContrib
    real, dimension(:, :),       pointer :: intensityExcess => null()

    ! -------------------------------------------------------------------------=                                           
    ! Output arrays
    
    real, dimension(:, :),    pointer :: fluxUp => null(), fluxDown => null(), &
                                         fluxAbsorbed => null()
    real, dimension(:, :, :), pointer :: volumeAbsorption => null()
    real, dimension(:, :, :), pointer :: intensity => null()
    real, dimension(:, :, :, :), &
                              pointer :: intensityByComponent => null()
    
  end type integrator
  !------------------------------------------------------------------------------------------
  ! Matrix is a private type used only inside the module. 
  !   One of the components of the public type (integrator) is of type matrix. 
  
  type matrix
    integer                        :: numX = 0, numY  = 0
    real, dimension(:, :), pointer :: values
  end type matrix
  !------------------------------------------------------------------------------------------
  ! What is visible? 
  !------------------------------------------------------------------------------------------
  public :: integrator
  public :: new_Integrator, copy_Integrator, isReady_Integrator, finalize_Integrator, &
            specifyParameters, computeRadiativeTransfer, reportResults
contains
  !------------------------------------------------------------------------------------------
  ! Initialization: Routines to create new integrators, specifying 
  !   integrator-specific parameters, making copies. 
  !------------------------------------------------------------------------------------------
  function new_Integrator(atmosphere, status) result(new)
    type(domain),                       intent(in   ) :: atmosphere
    type(ErrorMessage),                 intent(inout) :: status
    type(integrator)                                  :: new
    ! 
    ! Create an object of type integrator. This holds all the information
    !   about the atmosphere and surface that's needed to do the 
    !   problem at hand. 
    
    ! Local variables
    integer :: numX, numY, numZ, numComponents
    real    :: deltaX, deltaY, deltaZ
      
    ! Sanity checks - be sure that the atmosphere is initialized. 
    !--------------------------------------------------------------
    ! Atmosphere
    !   How big is the atmospheric domain, and what are the cell boundaries? 
    !
    call getInfo_Domain(atmosphere, numX, numY, numZ, &
                        numberOfComponents = numComponents,  status = status)  
    if(.not. stateIsFailure(status)) then   
      allocate(new%xPosition(numX + 1), new%yPosition(numY + 1), &
               new%zPosition(numZ + 1)) 
      call getInfo_Domain(atmosphere,                       &
                          xPosition = new%xPosition, &
                          yPosition = new%yPosition, &
                          zPosition = new%zPosition, status = status)
      !
      ! Are the cells all the same size in the X and Y or in the Z directions? 
      !   This enables us to take some shortcuts 
      !
      new%x0 = new%xPosition(1)
      new%y0 = new%yPosition(1)
      new%z0 = new%zPosition(1)
      deltaX = new%xPosition(2)  - new%xPosition(1)
      deltaY = new%yPosition(2)  - new%yPosition(1)
      deltaZ = new%zPosition(2)  - new%zPosition(1)
      if(all(abs( (new%xPosition(2:) - new%xPosition(:numX)) -  deltaX ) <= &
             2. * spacing(new%xPosition(2:))) .and.                         &
         all(abs( (new%yPosition(2:) - new%yPosition(:numY)) -  deltaY ) <= &
             2. * spacing(new%yPosition(2:)))) then
        new%xyRegularlySpaced = .true. 
        new%deltaX            = deltaX
        new%deltaY            = deltaY
      end if
      if(all(abs( (new%zPosition(2:) - new%zPosition(:numZ)) -  deltaZ ) <= &
             spacing(new%zPosition(2:)))) then
        new%zRegularlySpaced = .true. 
        new%deltaZ           = deltaZ
      end if
         
      ! 
      ! Now we make room for the optical properties of the domain and copy 
      !   these fields from the atmophere object. 
      !
      allocate(new%totalExt          (numX, numY, numZ),                &
               new%cumulativeExt     (numX, numY, numZ, numComponents), &
               new%ssa               (numX, numY, numZ, numComponents), &
               new%phaseFunctionIndex(numX, numY, numZ, numComponents), &
               new%forwardTables(numComponents))
      call getOpticalPropertiesByComponent(atmosphere,       & 
                   new%totalExt, new%cumulativeExt, new%ssa, &
                   new%phaseFunctionIndex, new%forwardTables, status)
      !
      ! cumulativeExt(:, :, :, numComponents) should be 1. in every cell in which 
      !   there's extinction (and 0 in empty cells). We'll chose which component 
      !   does the scattering by comparing these elements to a random number r
      !   and finding l such that cumulativeExt(i, j, k, l-1) <= r < cumulativeExt(i, j, k, l)
      !   We increase cumulativeExt(:, :, :, numComponents) slightly to account for the 
      !   edge case of r == 1. 
      !
      where(abs(new%cumulativeExt(:, :, :, numComponents) - 1.) <= spacing(1.)) & 
        new%cumulativeExt(:, :, :, numComponents) = 1. + spacing(1.)

    end if
    if(stateIsFailure(status)) then
      call setStateToFailure(status, "new_Integrator: Problems reading domain.")
      new%readyToCompute = .false. 
    else
      !--------------------------------------------------------------
      ! Create and zero the output arrays
      !
      allocate(new%fluxUp(numX, numY), new%fluxDown(numX, numY), new%fluxAbsorbed(numX, numY))
      new%fluxUp(:, :) = 0.; new%fluxDown(:, :) = 0.; new%fluxAbsorbed(:, :) = 0.
      allocate(new%volumeAbsorption(numX, numY, numZ))
      new%volumeAbsorption(:, :, :) = 0.
      !--------------------------------------------------------------
      ! Ready to go? 
      !
      new%readyToCompute = .true. 
      call setStateToSuccess(status)
    end if
  end function new_Integrator
  ! ------------------------------------------------------- 

  !------------------------------------------------------------------------------------------
  ! Computation: Routines to trace the photons and compute the radiative
  !   fluxes. There actual work might be done in other subroutines (to reflect 
  !   algorithmic choices, say). 
  !------------------------------------------------------------------------------------------
  subroutine computeRadiativeTransfer(thisIntegrator, randomNumbers, incomingPhotons, status)
    type(integrator),           intent(inout) :: thisIntegrator
    type(randomNumberSequence), intent(inout) :: randomNumbers
    type(photonStream),         intent(inout) :: incomingPhotons  
    type(ErrorMessage),         intent(inout) :: status
    !
    ! Monte Carlo "integrator" to compute flux up at the top boundary, flux down at the
    !   bottom boundary, and colum absorption. (Absorption is calculated separately, 
    !   and should agree with the difference in boundary fluxes.) 
    ! Intensity at a set of zenith angle cosine and azimuth angles at the top and bottom boundaries
    !   may be computed using local estimation. 
    ! The bottom boundary is a Lambertian surface with default albedo 0. That 
    !   value may be changed with a call to specifyParameters(). 
    ! This routine calls one of two solvers: computeRT_MaxCrossSection for maximum
    !   cross-section and computeRT_PhotonTracing for simple ray tracing. 
    
    ! Local variables
    integer :: numPhotonsProcessed
    integer :: numIntensityDirections, numX, numY, numZ, numComponents, j, k, d
    real, dimension(:, :), allocatable &
            :: numPhotonsPerColumn

    ! Sanity checks
    if(.not. isReady_integrator(thisIntegrator)) then
      call setStateToFailure(status, "computeRadiativeTransfer: problem not completely specified.")
    else
      numX = size(thisIntegrator%xPosition) - 1
      numY = size(thisIntegrator%yPosition) - 1
      numZ = size(thisIntegrator%zPosition) - 1
      numComponents = size(thisIntegrator%cumulativeExt, 4) 
      if(thisIntegrator%computeIntensity) & 
        numIntensityDirections = size(thisIntegrator%intensityDirections, 2)
      
      ! Zero output arrays
      thisIntegrator%fluxUp(:, :)       = 0.
      thisIntegrator%fluxDown(:, :)     = 0.
      thisIntegrator%fluxAbsorbed(:, :) = 0.
      thisIntegrator%volumeAbsorption(:, :, :) &
                                        = 0.
      !
      ! We want the intensity arrays to be zero'd even if we're not computing new 
      !   intensity values
      !
      if(associated(thisIntegrator%intensity)) thisIntegrator%intensity(:, :, :) = 0. 
      if(associated(thisIntegrator%intensityByComponent)) &
                                               thisIntegrator%intensityByComponent(:, :, :, :) = 0. 
      if(associated(thisIntegrator%intensityExcess)) &
                                               thisIntegrator%intensityExcess(:, :) = 0. 
      
      !
      ! Compute tablulated forward and inverse phase functions
      !   Forward phase functions are only needed if we're going to compute intensity
      !
      call tabulateInversePhaseFunctions(thisIntegrator, status)
      if(thisIntegrator%computeIntensity .and. .not. stateIsFailure(status)) then
        call tabulateForwardPhaseFunctions(thisIntegrator, status)
      end if
            
      !------------------------------------------------------------------------------
      !
      ! Compute radiative transfer for this photon batch 
      !
      if(.not. stateIsFailure(status)) &
         call computeRT(thisIntegrator, randomNumbers, incomingPhotons, numPhotonsProcessed, status)

      if(thisIntegrator%computeIntensity .and. &
         thisIntegrator%limitIntensityContributions) then 
        !
        ! Redistribute large intensity contributions in space
        !
        do j  = 0, numComponents
          do d = 1, numIntensityDirections 
            if(thisIntegrator%intensityExcess(d, j) > 0.) then 
              
              thisIntegrator%intensity(:, :, d) =  thisIntegrator%intensity(:, :, d) + &
                (    thisIntegrator%intensityByComponent(:, :, d, j) /                 &
                 sum(thisIntegrator%intensityByComponent(:, :, d, j)) ) * thisIntegrator%intensityExcess(d, j)
              
              thisIntegrator%intensityByComponent(:, :, d, j) =   & 
                thisIntegrator%intensityByComponent(:, :, d, j) + &
                (    thisIntegrator%intensityByComponent(:, :, d, j) /                 &
                 sum(thisIntegrator%intensityByComponent(:, :, d, j)) ) * thisIntegrator%intensityExcess(d, j)
            end if
          end do 
        end do
      end if
      
      !------------------------------------------------------------------------------
      !
      ! Normalization - compute the average number of photons incident on each column
      !
      allocate(numPhotonsPerColumn(numX, numY))
      
      if(thisIntegrator%xyRegularlySpaced) then
        numPhotonsPerColumn(:, :) = numPhotonsProcessed / real(numX * numY)
      else
        forall(j = 1:numY)
          ! Relative area of each column 
          numPhotonsPerColumn(:, j) = ( (thisIntegrator%yPosition(j+1) - thisIntegrator%yPosition(j)) *       &
                                        (thisIntegrator%xPosition(2:) - thisIntegrator%xPosition(1:numX)) ) / & 
                                      ( (thisIntegrator%xPosition(numX+1) - thisIntegrator%xPosition(1)) *    &
                                        (thisIntegrator%yPosition(numY+1) - thisIntegrator%yPosition(1)) ) 
        end forall
        ! Now the number of photons incident per column
        numPhotonsPerColumn(:, :) = numPhotonsPerColumn(:, :) * numPhotonsProcessed  
     end if 
      
      ! 
      ! Normalize fluxes by average number of photons per column
      !
      thisIntegrator%fluxUp      (:, :) = thisIntegrator%fluxUp      (:, :) / numPhotonsPerColumn(:, :)
      thisIntegrator%fluxDown    (:, :) = thisIntegrator%fluxDown    (:, :) / numPhotonsPerColumn(:, :)
      thisIntegrator%fluxAbsorbed(:, :) = thisIntegrator%fluxAbsorbed(:, :) / numPhotonsPerColumn(:, :) 
      !
      ! Normalize absorption profile by cell depth
      !
      forall(k = 1:numZ)  &
        thisIntegrator%volumeAbsorption(:, :, k) =      &
          thisIntegrator%volumeAbsorption(:, :, k) /    &
               (numPhotonsPerColumn(:, :) * (thisIntegrator%zPosition(k+1) - thisIntegrator%zPosition(k)) )
  
      !
      ! Intensity is also normalized by the average number of photons per column. 
      !
      if(thisIntegrator%computeIntensity) then
        forall(d = 1:numIntensityDirections)
          thisIntegrator%intensity(:, :, d) =  thisIntegrator%intensity(:, :, d) / numPhotonsPerColumn(:, :)
        end forall
        forall(d = 1:numIntensityDirections, j = 1:numComponents)
          thisIntegrator%intensityByComponent(:, :, d, j) =  &
                                               thisIntegrator%intensityByComponent(:, :, d, j) / &
                                                                                   numPhotonsPerColumn(:, :)
        end forall
      end if
      deallocate(numPhotonsPerColumn)  
    end if
  end subroutine computeRadiativeTransfer
  !------------------------------------------------------------------------------------------
  subroutine computeRT(thisIntegrator, randomNumbers, incomingPhotons, &
                                numPhotonsProcessed, status)
    type(integrator),           intent(inout) :: thisIntegrator
    type(randomNumberSequence), intent(inout) :: randomNumbers
    type(photonStream),         intent(inout) :: incomingPhotons
    integer,                    intent(  out) :: numPhotonsProcessed
    type(ErrorMessage),         intent(inout) :: status
    !
    ! Implements a standard ray-tracing Monte Carlo algorthm or 
    !   the Marchuk (1980) maximum cross-section algorithm. The extinction
    !   along each path segment is scaled by the maximum extinction within the domain. 
    !   Scattering events are then classified as "mathemtical" or "physical" by comparing
    !   a random number to the ratio of local to maximum extinction. 

        
    ! Local variables
    real :: xPos, yPos, zPos, mu, phi
    real :: x0, y0, z0, xMax, yMax, zMax
    real :: tauToTravel, photonWeight, scatteringAngle, tauAccumulated, ssa, maxExtinction
    real :: initialMu, initialPhi
    logical :: useRayTracing, useMaxCrossSection, scatterThisEvent
    integer :: xIndex, yIndex, zIndex, &
               component, phaseFunctionIndex, nPhotons
    integer :: i, scatteringOrder 
    integer :: nBad
    real, dimension(3) :: directionCosines
    !
    ! Variables related to intensity calculations
    !
    integer            :: numIntensityDirections
    real, dimension(:), &
           allocatable ::  contributions
    integer, dimension(:), &
           allocatable :: xIndexF, yIndexF
    
    ! ---------------------------------------------------------------------------------------
    useRayTracing = thisIntegrator%useRayTracing; useMaxCrossSection = .not. useRayTracing
    scatterThisEvent = .true. 
    if(useMaxCrossSection) &
      maxExtinction = maxval(thisIntegrator%totalExt(:, :, :))
    x0 = thisIntegrator%x0; xMax = thisIntegrator%xPosition(size(thisIntegrator%xPosition))
    y0 = thisIntegrator%y0; yMax = thisIntegrator%yPosition(size(thisIntegrator%yPosition))
    z0 = thisIntegrator%z0; zMax = thisIntegrator%zPosition(size(thisIntegrator%zPosition))
    if(thisIntegrator%computeIntensity) then
      numIntensityDirections = size(thisIntegrator%intensityDirections, 2)
      allocate(contributions(numIntensityDirections), &
                xIndexF(numIntensityDirections), yIndexF(numIntensityDirections))
    end if 
    !
    ! Begin loop over photons
    !
    nPhotons = 0; nBad = 0
    photonLoop: do
      if(.not. morePhotonsExist(incomingPhotons)) exit photonLoop ! This means we've used all the photons
      call getNextPhoton(incomingPhotons, xPos, yPos, zPos, mu, phi, status)
      if(stateIsFailure(status)) exit photonLoop
      scatteringOrder = 0
      directionCosines(:) = makeDirectionCosines(mu, phi)
      photonWeight = 1. 
      nPhotons = nPhotons + 1
      !
      ! Incoming xPos, yPos are between 0 and 1. 
      !   Translate these positions to the local domain
      !
      xPos = x0 + xPos * (xMax - x0) 
      yPos = y0 + yPos * (yMax - y0)
      zPos = z0 + zPos * (zMax - z0)      

      XIndex = 1; yIndex = 1; zIndex = 1
      call findXYIndicies(thisIntegrator, xPos, yPos, xIndex, yIndex)
      call findZIndex(thisIntegrator, zPos, zIndex)
      !
      ! Loop over orders of scattering
      !
      scatteringLoop: do
        !
        ! The optical distance we need to travel. 
        !   It's possible for the random number generator to produce exactly 0; 
        !   we set a lower bound. 
        !
        tauToTravel = -log(max(tiny(tauToTravel), getRandomReal(randomNumbers)))
        if(useRayTracing) then 
          !
          ! Ray tracing  - travel until we have accumulated enough extinction
          !
          call accumulateExtinctionAlongPath(thisIntegrator, directionCosines, &
                                             xPos, yPos, zPos, xIndex, yIndex, zIndex, &
                                             tauAccumulated, tauToTravel)            
          if(tauAccumulated < 0.) nBad = nBad + 1 
          if(tauAccumulated < 0.) cycle photonLoop 
        else 
          !
          ! Max cross-section: move the photon to the new location according to the maximum extinction
          !
          xPos = makePeriodic(xPos + directionCosines(1) * tauToTravel/maxExtinction, x0, xMax)
          yPos = makePeriodic(yPos + directionCosines(2) * tauToTravel/maxExtinction, y0, yMax)
          zPos =              zPos + directionCosines(3) * tauToTravel/maxExtinction
        end if 
        
        if(zPos >= zMax) then
          !
          ! The photon has gone out the top of the domain.  
          !   Add to the upward flux at that point, then start on a new photon
          !
          if(useMaxCrossSection) then 
            !
            ! Trace backwards to domain top if using max cross-section 
            ! 
            xPos = makePeriodic(xPos - directionCosines(1) * abs((zPos - zMax)/directionCosines(3)), x0, xMax)
            yPos = makePeriodic(yPos - directionCosines(2) * abs((zPos - zMax)/directionCosines(3)), y0, yMax)
            call findXYIndicies(thisIntegrator, xPos, yPos, xIndex, yIndex)
          end if 
          
          thisIntegrator%fluxUp(xIndex, yIndex) = thisIntegrator%fluxUp(xIndex, yIndex) + photonWeight
          cycle photonLoop
        else if(zPos <= z0 + spacing(z0)) then
          scatteringOrder = scatteringOrder + 1
          !
          ! The photon is at the surface. Add it to the surface flux. 
          !   directionCosines(3) will always be non-zero since the photon was traveling vertically. 
          !
          if(useMaxCrossSection) then 
            !
            ! Trace backwards to domain base if using max cross-section 
            ! 
            xPos = makePeriodic(xPos - directionCosines(1) * abs((zPos - z0)/directionCosines(3)), x0, xMax)
            yPos = makePeriodic(yPos - directionCosines(2) * abs((zPos - z0)/directionCosines(3)), y0, yMax)
            call findXYIndicies(thisIntegrator, xPos, yPos, xIndex, yIndex)
          end if 
          zIndex = 1
          zPos = z0 + spacing(z0)
          thisIntegrator%fluxDown(xIndex, yIndex) = thisIntegrator%fluxDown(xIndex, yIndex) + photonWeight
          !
          ! Compute new photon weight and a new direction to travel. 
          !   Save the old directions in case we're using a BDRF
          !
          initialMu  = directionCosines(3)
          if(initialMu <= 1.) then 
            initialPhi = acos(directionCosines(1)/sqrt(1. - initialMu**2))
          else
            initialPhi = 0. 
          end if 
          do 
            !
            ! Ensure that new trajectory has at least some vertical component - otherwise 
            !   the trajectory can get stuck forever if the lowest layer has no extinction
            !
            mu  = sqrt(getRandomReal(randomNumbers))
            if(abs(mu) > 2 * tiny(mu)) exit
          end do 
          phi = 2 * Pi * getRandomReal(randomNumbers)
          !
          ! New weight from the surface reflectance. 
          !
          if(thisIntegrator%useSurfaceBDRF) then 
            photonWeight = photonWeight * &
                           computeSurfaceReflectance(thisIntegrator%surfaceBDRF, xPos, yPos, &
                                                     initialMu, mu, initialPhi, phi)
          else
            !   Special case: Lambertian surface
            photonWeight = photonWeight * thisIntegrator%surfaceAlbedo
          end if
          if(photonWeight <= tiny(photonWeight)) cycle photonLoop
          directionCosines(:) = makeDirectionCosines(mu, phi)
          !
          ! Add contribution of surface reflection to intensity
          !
          if(thisIntegrator%computeIntensity) then
            call computeIntensityContribution(thisIntegrator, photonWeight, &
                                              xPos,   yPos,   zPos,         & 
                                              xIndex, yIndex, zIndex,       &
                                              directionCosines, 0,          &
                                              randomNumbers, scatteringOrder, &
                                              contributions, xIndexF(:), yIndexF(:))              
            forall(i = 1:numIntensityDirections)
              thisIntegrator%intensity(xIndexF(i), yIndexF(i), i) =  &
                thisIntegrator%intensity(xIndexF(i), yIndexF(i), i) + contributions(i)
              thisIntegrator%intensityByComponent(xIndexF(i), yIndexF(i), i, 0) = &
                thisIntegrator%intensityByComponent(xIndexF(i), yIndexF(i), i, 0) + contributions(i)
            end forall 
          end if 
        else 
          !
          ! Scattering event. 
          !
          
          ! Max cross-section - test for "Physical scattering event" 
          if(useMaxCrossSection) scatterThisEvent = & 
            getRandomReal(randomNumbers) < thisIntegrator%totalExt(xIndex, yIndex, zIndex)/maxExtinction

          if(useRayTracing .or. scatterThisEvent ) then
            scatteringOrder = scatteringOrder + 1
            !
            ! Time for the next scattering event.
            ! The photon might have accumulated the right amount of extinction 
            !   just on the boundary of a cell. If it's traveling in the positive
            !   direction, the upper boundary is considered in the next cell, and 
            !   that cell might not have any extinction in it. 
            ! In this (very rare) case we take the smallest possible step backwards, 
            !   since the optical properties aren't defined in cells with no extinction. 
            !  Don't need to worry about max cross-section, since we can't have a physical 
            !    scattering event if extinction is 0. 
            ! 
            ! We need to enforce periodicity here 
            !
      
            if(thisIntegrator%totalExt(xIndex, yIndex, zIndex) <= 0.) then
              if(xPos - thisIntegrator%xPosition(xIndex) <= 0. .and. directionCosines(1) > 0. ) then
                xPos = xPos - spacing(xPos) 
                xIndex = xIndex - 1
                if (xIndex <= 0) then
                  xIndex = size(thisIntegrator%xPosition) - 1
                  xPos =  thisIntegrator%xPosition(xIndex) 
                  xPos = xpos - 2. * spacing(xPos)
                end if
              end if 
              if(yPos - thisIntegrator%yPosition(yIndex) <= 0. .and. directionCosines(2) > 0. ) then
                yPos = yPos - spacing(yPos) 
                yIndex = yIndex - 1
                if (yIndex <= 0) then
                  yIndex = size(thisIntegrator%yPosition) - 1
                  yPos =  thisIntegrator%xPosition(yIndex) 
                  yPos = xpos - 2. * spacing(yPos)
                end if 
              end if 
              if(zPos - thisIntegrator%zPosition(zIndex) <= 0. .and. directionCosines(3) > 0. ) then
                zPos = zPos - spacing(zPos) 
                zIndex = zIndex - 1
              end if 
              !
              ! No need to worry about edge case - photons won't have come from below the domain
              !
            end if 
            !
            !   Figure out which component does the extinction, 
            !   and compute the new direction and weight of the photon. 
            !
            component = findIndex(getRandomReal(randomNumbers), &
                                  (/ 0., thisIntegrator%cumulativeExt(xIndex, yIndex, zIndex, :) /))
            !
            ! Absorption 
            !
            ssa = thisIntegrator%ssa(xIndex, yIndex, zIndex, component)
            if(ssa < 1.) then 
              thisIntegrator%fluxAbsorbed(xIndex, yIndex) =   &
                thisIntegrator%fluxAbsorbed(xIndex, yIndex)             + photonWeight * (1. - ssa)
              thisIntegrator%volumeAbsorption(xIndex, yIndex, zIndex) =   &
                thisIntegrator%volumeAbsorption(xIndex, yIndex, zIndex) + photonWeight * (1. - ssa)
              photonWeight = photonWeight * ssa
            end if 
  
            !
            ! Compute the contribution to intensity from this scattering event if need be
            !
            if(thisIntegrator%computeIntensity) then
              call computeIntensityContribution(thisIntegrator, photonWeight, &
                                                xPos,   yPos,   zPos,         & 
                                                xIndex, yIndex, zIndex,       &
                                                directionCosines, component,  &
                                                randomNumbers, scatteringOrder, &
                                                contributions, xIndexF(:), yIndexf(:))              
    
              do i = 1, numIntensityDirections
                thisIntegrator%intensity(xIndexF(i), yIndexf(i), i) =  &
                  thisIntegrator%intensity(xIndexF(i), yIndexf(i), i) + contributions(i)
                thisIntegrator%intensityByComponent(xIndexF(i), yIndexf(i), i, component) = &
                  thisIntegrator%intensityByComponent(xIndexF(i), yIndexf(i), i, component) + contributions(i)
              end do
            end if  ! end of local estimation for radiance contribution
    
            !
            ! "Russian roulette" 
            !
            if(thisIntegrator%useRussianRoulette .and. photonWeight < thisIntegrator%RussianRouletteW/2. ) then 
              if(getRandomReal(randomNumbers) >= photonWeight/thisIntegrator%RussianRouletteW) then 
                photonWeight = 0. 
              else
                photonWeight = thisIntegrator%RussianRouletteW
              end if 
            end if
            if(photonWeight <= tiny(photonWeight)) cycle photonLoop
            !
            ! Scattering - look up the scattering angle
            !
            phaseFunctionIndex = thisIntegrator%phaseFunctionIndex(xIndex, yIndex, zIndex, component)
            scatteringAngle = computeScatteringAngle(getRandomReal(randomNumbers), &
                                      thisIntegrator%inversePhaseFunctions(component)%values(:, phaseFunctionIndex)) 
            call next_direct(randomNumbers, cos(scatteringAngle), directionCosines)
          end if
        end if
      end do scatteringLoop
    end do photonLoop
    
    if(thisIntegrator%computeIntensity) deallocate(contributions, xIndexF, yIndexF)

    !
    ! Status is only set when getting photons, so if that hasn't failed we know we're ok. 
    !
    if(.not. stateIsFailure(status) .and. nPhotons > 0) then
      call setStateToCompleteSuccess(status, "computeRadiativeTransfer: finished with photons")
      numPhotonsProcessed = nPhotons
    else if(nPhotons == 0) then 
        call setStateToFailure(status, "computeRadiativeTransfer: Didn't process any photons.")
    else 
        call setStateToFailure(status, "computeRadiativeTransfer: Error.")
    end if
    
  end subroutine computeRT
  !------------------------------------------------------------------------------------------
  ! Reporting 
  !------------------------------------------------------------------------------------------
  subroutine reportResults(thisIntegrator,                             &
                           meanFluxUp, meanFluxDown, meanFluxAbsorbed, &
                               fluxUp,     fluxDown,     fluxAbsorbed, & 
                           absorbedProfile, volumeAbsorption,          &
                           meanIntensity, intensity, status)
    type(integrator),                   intent(in   ) :: thisIntegrator
    real,                     optional, intent(  out) :: meanFluxUp, meanFluxDown, meanFluxAbsorbed
    real, dimension(:, :),    optional, intent(  out) ::     fluxUp,     fluxDown,     fluxAbsorbed
    real, dimension(:),       optional, intent(  out) :: absorbedProfile
    real, dimension(:, :, :), optional, intent(  out) :: volumeAbsorption
    real, dimension(:),       optional, intent(  out) :: meanIntensity
    real, dimension(:, :, :), optional, intent(  out) ::     intensity
    type(ErrorMessage),                 intent(inout) :: status
    
    !
    ! Local variables
    !
    integer :: direction, numDirections, numColumns
    
    !
    ! This integrator computes fluxes at the top and bottom boundaries, the total 
    !   column absorption, and the absorption profile.  Users can ask for any of 
    !   the pixel level fluxes as a domain mean and/or on a column-by-column basis. 
    !
    numColumns = size(thisIntegrator%fluxUp) 
    
    ! Domain averaged fluxes
    !
    if(present(meanFluxUp))   meanFluxUp       = sum(thisIntegrator%fluxUp)       / numColumns
    if(present(meanFluxDown)) meanFluxDown     = sum(thisIntegrator%fluxDown)     / numColumns
    if(present(meanFluxAbsorbed)) &
                              meanFluxAbsorbed = sum(thisIntegrator%fluxAbsorbed) / numColumns
    !
    ! Pixel-by-pixel fluxes
    !
    if(present(fluxUp)) then
      if(any((/ size(fluxUp, 1), size(fluxUp, 2) /) /= &
             (/ size(thisIntegrator%fluxUp, 1), size(thisIntegrator%fluxUp, 2) /))) then
        call setStateToFailure(status, "reportResults: fluxUp array is the wrong size")
      else
        fluxUp(:, :) = thisIntegrator%fluxUp(:, :)
      end if 
    end if 
    
    if(present(fluxDown)) then
      if(any((/ size(fluxDown, 1), size(fluxDown, 2) /) /= &
             (/ size(thisIntegrator%fluxDown, 1), size(thisIntegrator%fluxDown, 2) /))) then
        call setStateToFailure(status, "reportResults: fluxDown array is the wrong size")
      else
        fluxDown(:, :) = thisIntegrator%fluxDown(:, :)
      end if 
    end if 
    
    if(present(fluxAbsorbed)) then
      if(any((/ size(fluxAbsorbed, 1), size(fluxAbsorbed, 2) /) /= &
             (/ size(thisIntegrator%fluxAbsorbed, 1), size(thisIntegrator%fluxAbsorbed, 2) /))) then
        call setStateToFailure(status, "reportResults: fluxAbsorbed array is the wrong size")
      else
        fluxAbsorbed(:, :) = thisIntegrator%fluxAbsorbed(:, :)
      end if 
    end if 
    
    !
    ! Absorption - heating rate profile and volume absorption
    !
    if(present(absorbedProfile)) then
      if (size(absorbedProfile) /= size(thisIntegrator%volumeAbsorption, 3)) then
        call setStateToFailure(status, "reportResults: absorbedProfile array is the wrong size")
      else
        absorbedProfile(:) = sum(sum(thisIntegrator%volumeAbsorption(:, :, :), dim = 1), dim = 1) / numColumns
      end if 
    end if 
     
    if(present(volumeAbsorption)) then
      if(any((/ size(volumeAbsorption, 1), size(volumeAbsorption, 2), size(volumeAbsorption, 3) /) /= &
             (/ size(thisIntegrator%volumeAbsorption, 1), size(thisIntegrator%volumeAbsorption, 2),   &
                size(thisIntegrator%volumeAbsorption, 3)  /))) then
        call setStateToFailure(status, "reportResults: volumeAbsorption array is the wrong size")
      else
        volumeAbsorption(:, :, :) = thisIntegrator%volumeAbsorption(:, :, :)
      end if 
    end if 
    !
    ! Domain-averaged intensity
    !
    if(present(meanIntensity)) then
      if(.not. associated(thisIntegrator%intensity)) then 
        call setStateToFailure(status, "reportResults: intensity information not available") 
      else if (size(thisIntegrator%intensity, 3) /= size(meanIntensity)) then 
        call setStateToFailure(status, "reportResults: requesting mean intensity in the wrong number of directions.") 
      else
        numDirections = size(thisIntegrator%intensity, 3)
        forall(direction = 1:numDirections)
          meanIntensity(direction) = sum(thisIntegrator%intensity(:, :, direction)) / numColumns
        end forall
      end if 
    end if
    
    !
    ! Pixel-by-pixel intensity
    !
    if(present(intensity)) then
      if(.not. associated(thisIntegrator%intensity)) then 
        call setStateToFailure(status, "reportResults: intensity information not available") 
      else if (any( (/ size(               intensity, 1), size(               intensity, 2), &
                       size(               intensity, 3) /) /=                               &
                    (/ size(thisIntegrator%intensity, 1), size(thisIntegrator%intensity, 2), &
                       size(thisIntegrator%intensity, 3) /) ))  then 
        call setStateToFailure(status, "reportResults: intensity array has wrong dimensions.") 
      else
        intensity(:, :, :) = thisIntegrator%intensity(:, :, :)
      end if
    end if
    
    if(.not. stateIsFailure(status)) call setStateToSuccess(status)
  end subroutine reportResults 
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  subroutine specifyParameters(thisIntegrator, surfaceAlbedo, surfaceBDRF,    &
                               minForwardTableSize, minInverseTableSize,      &
                               intensityMus, intensityPhis, computeIntensity, &
                               useRayTracing,  useRussianRoulette,                   &
                               useRussianRouletteForIntensity, zetaMin,              &
                               useHybridPhaseFunsForIntenCalcs, hybridPhaseFunWidth, &
                               numOrdersOrigPhaseFunIntenCalcs,                      &
                               limitIntensityContributions, maxIntensityContribution, &
                               status)
    type(integrator),   intent(inout) :: thisIntegrator
    real,     optional, intent(in   ) :: surfaceAlbedo
    type(surfaceDescription), &
              optional, intent(in   ) :: surfaceBDRF
    integer,  optional, intent(in   ) :: minForwardTableSize, minInverseTableSize
    real, dimension(:), &
              optional, intent(in   ) :: intensityMus, intensityPhis
    logical,  optional, intent(in   ) :: computeIntensity, useRayTracing, useRussianRoulette
    logical,  optional, intent(in   ) :: useRussianRouletteForIntensity
    real,     optional, intent(in   ) :: zetaMin
    logical,  optional, intent(in   ) :: useHybridPhaseFunsForIntenCalcs
    real,     optional, intent(in   ) :: hybridPhaseFunWidth
    integer,  optional, intent(in   ) :: numOrdersOrigPhaseFunIntenCalcs
    logical,  optional, intent(in   ) :: limitIntensityContributions
    real,     optional, intent(in   ) :: maxIntensityContribution
    type(ErrorMessage), intent(inout) :: status
    
    !
    !  Set integrator-specific parameter. The parameters might or might not need to be persistent. 
    !    If the medium is known to be highly variable within some vertical subset
    !    of the domain, for example, we might store the upper and lower boundaries 
    !    of this region, then do maximum cross section within it and photon tracing
    !    outside. 
    !  This might also be useful outside the core integrators, i.e. those that work 
    !    on a single batch of photons.
    
    ! Local variables
    integer :: i

    ! ------------------------------------------------------------------------          
    ! 
    ! Sanity checks for input variables
    !
    if(present(surfaceBDRF) .and. present(surfaceAlbedo)) &
      call setStateToFailure(status, "specifyParameters: only one surface specification can be provided") 
    !
    !   If surface albedo is supplied it should be between 0 and 1. 
    !
    if(present(surfaceAlbedo)) then
      if(surfaceAlbedo > 1. .or. surfaceAlbedo < 0.) &
        call setStateToFailure(status, "specifyParameters: surface albedo out of range.") 
    end if
    if(present(surfaceBDRF)) then
      if(.not. isReady_surfaceDescription(surfaceBDRF)) & 
        call setStateToFailure(status, "specifyParameters: surface description isn't valid.") 
    end if 
    
    !
    ! Table sizes must be bigger than 0; could put a minimum bound
    !
    if( present(minForwardTableSize)) then
      if(minForwardTableSize < defaultMinForwardTableSize) &
        call setStateToWarning(status, "specifyParameters: minForwardTableSize less than default. Value ignored.") 
    end if
    if( present(minInverseTableSize)) then
      if(minInverseTableSize < defaultMinInverseTableSize) &
        call setStateToWarning(status, "specifyParameters: minInverseTableSize less than default. Value ignored.") 
    end if
    
    if(present(hybridPhaseFunWidth)) then
      if(hybridPhaseFunWidth > maxhybridPhaseFunWidth .or. hybridPhaseFunWidth < 0.)        &
        call setStateToWarning(status,                                                         &
                               "specifyParameters: hybridPhaseFunWidth out of range (0 to " // &
                               trim(intToChar(int(maxhybridPhaseFunWidth))) // "degrees)."  // &
                               "Using default (" // trim(intToChar(int(defaultHybridPhaseFunWidth))) // ")")
    end if
    if(present(numOrdersOrigPhaseFunIntenCalcs)) then
      if(numOrdersOrigPhaseFunIntenCalcs < 0)                                                        &
        call setStateToWarning(status,                                                               &
                               "specifyParameters: numOrdersExactPhaseFunIntenCalcs less than 0." // &
                               "Using default (" // trim(intToChar(defOrdersOrigPhaseFunIntenCalcs)) // ")")
    end if
    
    if(present(maxIntensityContribution)) then
      if(maxIntensityContribution <= 0.)                                                      &
        call setStateToWarning(status,                                                       &
                               "specifyParameters: maxIntensityContribution <= 0. Value is unchanged.")
        
    end if 


    ! intensity direction arrays should be the same length
    !  both or neither supplied
    ! intensityMus must be between -1 and 1; can't be identically 0
    ! intensityPhis must be between 0 and 360 
    !
    if(present(intensityMus) .neqv. present(intensityPhis)) &
      call setStateToFailure(status, "specifyParameters: Both or neither of intensityMus and intensityPhis must be supplied") 
    if(present(intensityMus)) then
      if(size(intensityMus) /= size(intensityPhis)) &
        call setStateToFailure(status, "specifyParameters: intensityMus, intensityPhis must be the same length.") 
      if(any(intensityMus < -1.) .or. any(intensityMus > 1.)) &
        call setStateToFailure(status, "specifyParameters: intensityMus must be between -1 and 1") 
      if(any(abs(intensityMus) < tiny(intensityMus))) &
        call setStateToFailure(status, "specifyParameters: intensityMus can't be 0 (directly sideways)") 
      if(any(intensityPhis < 0.) .or. any(intensityPhis > 360.)) &
        call setStateToFailure(status, "specifyParameters: intensityPhis must be between 0 and 360") 
    end if
    
    !
    ! If someone specifies intensityDirections but trys to set computeIntensity to false we ignore them
    !
    if(present(computeIntensity)) then
      if(.not. computeIntensity .and. present(intensityMus))                                                             &
        call setStateToWarning(status, "specifyParameters: intensity directions *and* computeIntensity set to false." // &
                                       "Will compute intensity at given angles.") 
      if(computeIntensity .and. .not. present(intensityMus) .and. .not. associated(thisIntegrator%intensityDirections)) &
        call setStateToFailure(status, "specifyParameters: Can't compute intensity without specifying directions.") 
    end if
    
    ! ------------------------------------------------------------------------          
    if( .not. StateIsFailure(status)) then 
      if(present(surfaceAlbedo)) then 
        thisIntegrator%surfaceAlbedo = surfaceAlbedo
        thisIntegrator%useSurfaceBDRF = .false. 
      else if(present(surfaceBDRF)) then
        thisIntegrator%surfaceBDRF = copy_SurfaceDescription(surfaceBDRF)
        thisIntegrator%useSurfaceBDRF = .true. 
      end if 
      !
      ! Algorithmic choices
      !
      if(present(useRayTracing)) thisIntegrator%useRayTracing = useRayTracing
      if(present(minForwardTableSize)) &
        thisIntegrator%minForwardTableSize    = max(minForwardTableSize, defaultMinForwardTableSize)
      if(present(minInverseTableSize)) &
        thisIntegrator%minInverseTableSize    = max(minInverseTableSize, defaultMinInverseTableSize)

      if(present(useRussianRoulette)) &
        thisIntegrator%useRussianRoulette = useRussianRoulette
        
      !
      ! Russian roulette for intensity
      !
      if(present(useRussianRouletteForIntensity)) &
         thisIntegrator%useRussianRouletteForIntensity = useRussianRouletteForIntensity
      if(present(zetaMin)) then 
        if(zetaMin < 0.) then 
          call setStateToWarning(status, "specifyParameters: zetaMin must be >= 0. Value is unchanged.")
        else 
          thisIntegrator%zetaMin = zetaMin
          if(zetaMin > 1.) &
            call setStateToWarning(status, "specifyParameters: zetaMin > 1. That's kind of large.")
        end if 
      end if 

      !
      ! Hyrbid phase function for local estimation
      !
      if(present(useHybridPhaseFunsForIntenCalcs)) &
        thisIntegrator%useHybridPhaseFunsForIntenCalcs = useHybridPhaseFunsForIntenCalcs
      if(present(hybridPhaseFunWidth)) then
        if(hybridPhaseFunWidth > 0 .and. hybridPhaseFunWidth < maxhybridPhaseFunWidth) then
          thisIntegrator%hybridPhaseFunWidth = hybridPhaseFunWidth
        else
           thisIntegrator%hybridPhaseFunWidth = defaultHybridPhaseFunWidth
        end if
        !
        ! Need to re-tabulate the phase functions
        !
        if(associated(thisIntegrator%tabulatedPhaseFunctions)) then 
          do i = 1, size(thisIntegrator%tabulatedPhaseFunctions)
            call finalize_Matrix(thisIntegrator%tabulatedPhaseFunctions(i))
            deallocate(thisIntegrator%tabulatedPhaseFunctions)
          end do
        end if
      end if 
      if(present(numOrdersOrigPhaseFunIntenCalcs)) then
        if(numOrdersOrigPhaseFunIntenCalcs >= 0 ) then
          thisIntegrator%numOrdersOrigPhaseFunIntenCalcs = numOrdersOrigPhaseFunIntenCalcs
        else
          thisIntegrator%numOrdersOrigPhaseFunIntenCalcs = defOrdersOrigPhaseFunIntenCalcs
        end if
      end if 

      !
      ! Limited maximum local estimate
      !
      if(present(limitIntensityContributions)) & 
        thisIntegrator%limitIntensityContributions = limitIntensityContributions
      if(present(maxIntensityContribution)) then 
        if(maxIntensityContribution > 0.) &
          thisIntegrator%maxIntensityContribution = maxIntensityContribution
      end if
      !
      ! Intensity 
      !
      if(present(intensityMus)) then 
        if(associated(thisIntegrator%intensityDirections))  deallocate(thisIntegrator%intensityDirections)
        allocate(thisIntegrator%intensityDirections(3, size(intensityMus)))

        if(associated(thisIntegrator%intensity))            deallocate(thisIntegrator%intensity)
        allocate(thisIntegrator%intensity(size(thisIntegrator%totalExt, 1), &
                                          size(thisIntegrator%totalExt, 2), &
                                          size(intensityMus)))

        if(associated(thisIntegrator%intensityByComponent)) deallocate(thisIntegrator%intensityByComponent) 
        allocate(thisIntegrator%intensityByComponent(size(thisIntegrator%totalExt, 1), &
                                                     size(thisIntegrator%totalExt, 2), &
                                                     size(intensityMus),               &
                                                     0:size(thisIntegrator%cumulativeExt, 4)))

        forall(i = 1:size(intensityMus)) &
          thisIntegrator%intensitydirections(:, i) = makeDirectionCosines(intensityMus(i), &
                                                                          intensityPhis(i) * Pi/180.)
        thisIntegrator%computeIntensity = .true.
      end if 
      
      if(present(computeIntensity)) then 
        !
        ! If computeintensity is true we've already assured that the directions have been supplied at some point, 
        !    so there's nothing to do. 
        !
        if(.not. computeIntensity .and. .not. present(intensityMus)) then
          if(associated(thisIntegrator%intensityDirections))  deallocate(thisIntegrator%intensityDirections)
          if(associated(thisIntegrator%intensity))            deallocate(thisIntegrator%intensity)
          if(associated(thisIntegrator%intensityByComponent)) deallocate(thisIntegrator%intensityByComponent)
          thisIntegrator%computeIntensity = .false.
        end if 
      end if 
            
      if(thisIntegrator%computeIntensity .and. thisIntegrator%limitIntensityContributions) then
        if(associated(thisIntegrator%intensityExcess))      deallocate(thisIntegrator%intensityExcess)
        allocate(thisIntegrator%intensityExcess(size(thisIntegrator%intensityDirections, 2), &
                                                0:size(thisIntegrator%cumulativeExt, 4)))
      end if

      call setStateToSuccess(status)
    end if 
    
  end subroutine specifyParameters 
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  function isReady_Integrator(thisIntegrator)
    type(integrator), intent( in) :: thisIntegrator
    logical                       :: isReady_Integrator
    
    isReady_Integrator = thisIntegrator%readyToCompute
  end function isReady_Integrator
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  function copy_Integrator(original) result(copy)
    type(integrator), intent(in) :: original
    type(integrator)             :: copy
    !
    ! Copy the state of one integrator to another
    !
    
    ! Local variables
    integer :: i
    
    ! Scalars
    copy%readyToCompute    = original%readyToCompute
    copy%computeIntensity  = original%computeIntensity
    
    copy%surfaceAlbedo     = original%surfaceAlbedo
    if(isReady_surfaceDescription(original%surfaceBDRF)) &
      copy%surfaceBDRF     = copy_surfaceDescription(original%surfaceBDRF)
    copy%useSurfaceBDRF    = original%useSurfaceBDRF
    
    copy%xyRegularlySpaced = original%xyRegularlySpaced
    copy%zRegularlySpaced  = original%zRegularlySpaced
    
    copy%deltaX = original%deltaX; copy%deltaY = original%deltaY; copy%deltaZ = original%deltaZ
    copy%x0     = original%x0;     copy%y0     = original%y0;     copy%z0     = original%z0 
    
    copy%minForwardTableSize = original%minForwardTableSize
    copy%minInverseTableSize = original%minInverseTableSize
    !
    ! Algorithmic choices
    !
    copy%useRayTracing = original%useRayTracing
    
    copy%useRussianRoulette = original%useRussianRoulette
    copy%RussianRouletteW   = original%RussianRouletteW

    copy%useRussianRouletteForIntensity   = original%useRussianRouletteForIntensity
    copy%zetaMin                          = original%zetaMin
    
    !
    ! Evans phase function truncation
    !
    copy%useHybridPhaseFunsForIntenCalcs = original%useHybridPhaseFunsForIntenCalcs
    copy%hybridPhaseFunWidth             = original%hybridPhaseFunWidth
    copy%numOrdersOrigPhaseFunIntenCalcs = original%numOrdersOrigPhaseFunIntenCalcs

    !
    ! Location vectors
    !
    if(associated(original%xPosition)) then
      allocate(copy%xPosition(size(original%xPosition)))
      copy%xPosition(:) = original%xPosition(:)
    end if

    if(associated(original%yPosition)) then
      allocate(copy%yPosition(size(original%yPosition)))
      copy%yPosition(:) = original%yPosition(:)
    end if

    if(associated(original%zPosition)) then
      allocate(copy%zPosition(size(original%zPosition)))
      copy%zPosition(:) = original%zPosition(:)
    end if
    
    !
    ! Extinction, ssa, phase function index (three-D and four-D arrays)
    !
    if(associated(original%totalExt)) then
      allocate(copy%totalExt(size(original%totalExt, 1), &
                             size(original%totalExt, 2), &
                             size(original%totalExt, 3)))
      copy%totalExt(:, :, :) = original%totalExt(:, :, :)
    end if

    if(associated(original%cumulativeExt)) then
      allocate(copy%cumulativeExt(size(original%cumulativeExt, 1), &
                                  size(original%cumulativeExt, 2), &
                                  size(original%cumulativeExt, 3), &
                                  size(original%cumulativeExt, 4)))
      copy%cumulativeExt(:, :, :, :) = original%cumulativeExt(:, :, :, :)
    end if

    if(associated(original%ssa)) then
      allocate(copy%ssa(size(original%ssa, 1), size(original%ssa, 2), &
                        size(original%ssa, 3), size(original%ssa, 4)))
      copy%ssa(:, :, :, :) = original%ssa(:, :, :, :)
    end if

    if(associated(original%phaseFunctionIndex)) then
      allocate(copy%phaseFunctionIndex(size(original%phaseFunctionIndex, 1), &
                                  size(original%phaseFunctionIndex, 2), &
                                  size(original%phaseFunctionIndex, 3), &
                                  size(original%phaseFunctionIndex, 4)))
      copy%phaseFunctionIndex(:, :, :, :) = original%phaseFunctionIndex(:, :, :, :)
    end if
    
    !
    ! Intensity directions
    !
    if(associated(original%intensityDirections)) then
      allocate(copy%intensityDirections(size(original%intensityDirections, 1), size(original%intensityDirections, 2)))
      copy%intensityDirections(:, :)  = original%intensityDirections(:, :)
    end if 
    !
    ! Forward  phase function tables
    !   Have to loop over the tables for each component
    !
    if(associated(original%forwardTables)) then 
      allocate(copy%forwardTables(size(original%forwardTables)))
      do i = 1, size(original%forwardTables)
        copy%forwardTables(i) = copy_PhaseFunctionTable(original%forwardTables(i))
      end do
    end if
    
    !
    ! Tabulated forward and inverse phase function tables
    !
    if(associated(original%tabulatedPhaseFunctions)) then 
      allocate(copy%tabulatedPhaseFunctions(size(original%tabulatedPhaseFunctions)))
      do i = 1, size(original%tabulatedPhaseFunctions)
        copy%tabulatedPhaseFunctions(i) = new_Matrix(original%tabulatedPhaseFunctions(i)%values)
      end do
    end if
    if(associated(original%tabulatedOrigPhaseFunctions)) then 
      allocate(copy%tabulatedOrigPhaseFunctions(size(original%tabulatedOrigPhaseFunctions)))
      do i = 1, size(original%tabulatedOrigPhaseFunctions)
        copy%tabulatedOrigPhaseFunctions(i) = new_Matrix(original%tabulatedOrigPhaseFunctions(i)%values)
      end do
    end if
    if(associated(original%inversePhaseFunctions)) then 
      allocate(copy%inversePhaseFunctions(size(original%inversePhaseFunctions)))
      do i = 1, size(original%inversePhaseFunctions)
        copy%inversePhaseFunctions(i) = new_Matrix(original%inversePhaseFunctions(i)%values)
      end do
    end if
           
    
    !
    ! Output arrays
    !
    if(associated(original%fluxUp)) then
      allocate(copy%fluxUp(size(original%fluxUp, 1), &
                           size(original%fluxUp, 2)))
      copy%fluxUp(:, :) = original%fluxUp(:, :)
    end if

    if(associated(original%fluxDown)) then
      allocate(copy%fluxDown(size(original%fluxDown, 1), &
                             size(original%fluxDown, 2)))
      copy%fluxDown(:, :) = original%fluxDown(:, :)
    end if

    if(associated(original%fluxAbsorbed)) then
      allocate(copy%fluxAbsorbed(size(original%fluxAbsorbed, 1), &
                                 size(original%fluxAbsorbed, 2)))
      copy%fluxAbsorbed(:, :) = original%fluxAbsorbed(:, :)
    end if
    
    if(associated(original%volumeAbsorption)) then
      allocate(copy%volumeAbsorption(size(original%volumeAbsorption, 1), &
                                     size(original%volumeAbsorption, 2), &
                                     size(original%volumeAbsorption, 3)))
      copy%volumeAbsorption(:, :, :) = original%volumeAbsorption(:, :, :)
    end if
    
    if(associated(original%intensity)) then
      allocate(copy%intensity(size(original%intensity, 1), &
                              size(original%intensity, 2), &
                              size(original%intensity, 3)))
      copy%intensity(:, :, :) = original%intensity(:, :, :)
    end if

  end function copy_Integrator

  !------------------------------------------------------------------------------------------
  ! Finalization 
  !------------------------------------------------------------------------------------------
  subroutine finalize_Integrator(thisIntegrator)
    type(integrator), intent(out) :: thisIntegrator
    
    !
    ! Finalize by copying component values from a variable that's never been used
    !
    type(integrator)               :: pristineI
    
    ! Local variable
    integer :: i
    
    thisIntegrator%readyToCompute    =  pristineI%readyToCompute
    thisIntegrator%computeIntensity  =  pristineI%computeIntensity 
    
    thisIntegrator%xyRegularlySpaced = pristineI%xyRegularlySpaced 
    thisIntegrator%zRegularlySpaced  = pristineI%zRegularlySpaced
    
    thisIntegrator%useRayTracing              = pristineI%useRayTracing     
    thisIntegrator%useRussianRoulette         = pristineI%useRussianRoulette 
    thisIntegrator%RussianRouletteW           = pristineI%RussianRouletteW
    
    thisIntegrator%useRussianRouletteForIntensity = pristineI%useRussianRouletteForIntensity 
    thisIntegrator%zetaMin                        = pristineI%zetaMin 

    thisIntegrator%useHybridPhaseFunsForIntenCalcs = pristineI%useHybridPhaseFunsForIntenCalcs 
    thisIntegrator%hybridPhaseFunWidth             = pristineI%hybridPhaseFunWidth
    thisIntegrator%numOrdersOrigPhaseFunIntenCalcs = pristineI%numOrdersOrigPhaseFunIntenCalcs 

    thisIntegrator%deltaX = 0.; thisIntegrator%deltaY = 0.; thisIntegrator%deltaZ = 0. 
    thisIntegrator%x0     = 0.; thisIntegrator%y0     = 0.; thisIntegrator%z0     = 0. 
    
    thisIntegrator%minForwardTableSize = pristineI%minForwardTableSize
    thisIntegrator%minInverseTableSize = pristineI%minInverseTableSize
    
    thisIntegrator%surfaceAlbedo =  pristineI%surfaceAlbedo
    thisIntegrator%useSurfaceBDRF = pristineI%useSurfaceBDRF 
    call finalize_surfaceDescription(thisIntegrator%surfaceBDRF)
    
    if(associated(thisIntegrator%xPosition))          deallocate(thisIntegrator%xPosition)
    if(associated(thisIntegrator%yPosition))          deallocate(thisIntegrator%yPosition)
    if(associated(thisIntegrator%zPosition))          deallocate(thisIntegrator%zPosition)
    if(associated(thisIntegrator%totalExt))           deallocate(thisIntegrator%totalExt)
    if(associated(thisIntegrator%cumulativeExt))      deallocate(thisIntegrator%cumulativeExt)
    if(associated(thisIntegrator%ssa))                deallocate(thisIntegrator%ssa)
    if(associated(thisIntegrator%phaseFunctionIndex)) deallocate(thisIntegrator%phaseFunctionIndex)

    !
    ! Forward and inverse phase functions - finalize each element in the array (i.e. 
    !   the table for each component) to free the underlying memory, 
    !   then deallocate the array that holds the tables
    !
    if(associated(thisIntegrator%forwardTables)) then 
      do i = 1, size(thisIntegrator%forwardTables)
        call finalize_PhaseFunctionTable(thisIntegrator%forwardTables(i))
      end do
      deallocate(thisIntegrator%forwardTables)
    end if 
    
    if(associated(thisIntegrator%intensityDirections))    deallocate(thisIntegrator%intensityDirections)
    
    !
    ! Tabulated forward and inverse phase functions - stored as one matrix per
    !   component. Finalize each matrix, then the array that holds them. 
    !
    if(associated(thisIntegrator%tabulatedPhaseFunctions)) then 
      do i = 1, size(thisIntegrator%tabulatedPhaseFunctions)
        call finalize_Matrix(thisIntegrator%tabulatedPhaseFunctions(i))
      end do 
      deallocate(thisIntegrator%tabulatedPhaseFunctions)
    end if
    if(associated(thisIntegrator%tabulatedOrigPhaseFunctions)) then 
      do i = 1, size(thisIntegrator%tabulatedOrigPhaseFunctions)
        call finalize_Matrix(thisIntegrator%tabulatedOrigPhaseFunctions(i))
      end do 
      deallocate(thisIntegrator%tabulatedOrigPhaseFunctions)
    end if
    if(associated(thisIntegrator%inversePhaseFunctions)) then 
      do i = 1, size(thisIntegrator%inversePhaseFunctions)
        call finalize_Matrix(thisIntegrator%inversePhaseFunctions(i))
      end do 
      deallocate(thisIntegrator%inversePhaseFunctions)
    end if
    
    !
    ! Output arrays
    !
    if(associated(thisIntegrator%fluxUp))           deallocate(thisIntegrator%fluxUp)
    if(associated(thisIntegrator%fluxDown))         deallocate(thisIntegrator%fluxDown)
    if(associated(thisIntegrator%fluxAbsorbed))     deallocate(thisIntegrator%fluxAbsorbed)
    if(associated(thisIntegrator%volumeAbsorption)) deallocate(thisIntegrator%volumeAbsorption)
    if(associated(thisIntegrator%intensity))        deallocate(thisIntegrator%intensity)
  end subroutine finalize_Integrator
  !------------------------------------------------------------------------------------------
  ! Functions for use inside the module 
  !------------------------------------------------------------------------------------------
  pure subroutine findXYIndicies(thisIntegrator, xPos, yPos, xIndex, yIndex)
    type(integrator), intent(in ) :: thisIntegrator
    real,             intent(in ) :: xPos, yPos
    integer,          intent(inout) :: xIndex, yIndex
    
    
    if(thisIntegrator%xyRegularlySpaced) then
      xIndex = min(int((xPos - thisIntegrator%x0)/thisIntegrator%deltaX) + 1, & 
                   size(thisIntegrator%totalExt, 1))
      yIndex = min(int((yPos - thisIntegrator%y0)/thisIntegrator%deltaY) + 1, &
                   size(thisIntegrator%totalExt, 2))
      !
      ! Insure against rounding errors
      if(abs(thisIntegrator%xPosition(xIndex+1) - xPos) < spacing(xPos)) xIndex = xIndex + 1
      if(abs(thisIntegrator%yPosition(yIndex+1) - yPos) < spacing(yPos)) yIndex = yIndex + 1
      if(xIndex == size(thisIntegrator%xPosition)) xIndex = 1
      if(yIndex == size(thisIntegrator%yPosition)) yIndex = 1
    else
      xIndex = findIndex(xPos, thisIntegrator%xPosition, xIndex)
      yIndex = findIndex(yPos, thisIntegrator%yPosition, yIndex)
    end if
  end subroutine findXYIndicies
  !------------------------------------------------------------------------------------------
  pure subroutine findZIndex(thisIntegrator, zPos, zIndex)
    type(integrator), intent(in ) :: thisIntegrator
    real,             intent(in ) :: zPos
    integer,          intent(out) :: zIndex
    if(thisIntegrator%zRegularlySpaced) then
      zIndex = min(int((zPos - thisIntegrator%z0)/thisIntegrator%deltaZ) + 1, &
                   size(thisIntegrator%totalExt, 3))
      ! Insure against rounding errors
      if(abs(thisIntegrator%zPosition(zIndex+1) - zPos) < spacing(zPos)) zIndex = zIndex + 1
    else
      zIndex = findIndex(zPos, thisIntegrator%zPosition, zIndex)
    end if
  end subroutine findZIndex
  !------------------------------------------------------------------------------------------
  pure function computeScatteringAngle(randomDeviate, inversePhaseFunctionTable) result(scatteringAngle)
    !
    ! Linearly interpolate the scattering angle from a vector containing the 
    !   the angle as a function of the cumulative distribution (inverse phase function). 
    !   Recall that the first entry in the table is for CDF = 0 (hence angleIndex - 1) 
    ! 
    real,              intent (in ) :: randomDeviate
    real, dimension(:), intent(in ) :: inversePhaseFunctionTable
    real                            :: scatteringAngle
  
    ! Local variables
    integer :: angleIndex, numIntervals
    real    :: leftOver
    
    numIntervals = size(inversePhaseFunctionTable)
    angleIndex = int(randomDeviate * numIntervals) + 1
    if(angleIndex < numIntervals) then
      !
      ! Interpolate between entries in the inverse phase function table. 
      !   The first entry in the table is for CDF = 0 (hence angleIndex - 1) 
      !
      leftOver           = randomDeviate - real(angleIndex - 1)/real(numIntervals)
      scatteringAngle = (1. - leftOver) * inversePhaseFunctionTable(angleIndex) + &
                              leftOver  * inversePhaseFunctionTable(angleIndex + 1)
    else
      scatteringAngle = inversePhaseFunctionTable(numIntervals) 
    end if 
  end function computeScatteringAngle
  !------------------------------------------------------------------------------------------
  subroutine computeIntensityContribution(thisIntegrator, photonWeight,    &
                                          xPos,   yPos,   zPos,         & 
                                          xIndex, yIndex, zIndex,       &
                                          directionCosines, component,  &
                                          randomNumbers, scatteringOrder, &
                                          contributions, xIndexF, yIndexF)
    !
    ! Compute the contribution to the intensity in each direction 
    !   from a scattering event at xPos, xPos, zPos from a specified
    !   component. 
    !
    type(integrator),      intent(inout) :: thisIntegrator
    real,                  intent(in   ) :: photonWeight, xPos,   yPos,   zPos
    integer,               intent(in   ) ::               xIndex, yIndex, zIndex
    real,    dimension(:), intent(in   ) :: directionCosines
    integer,               intent(in   ) :: component
    type(randomNumberSequence), &
                           intent(inout) :: randomNumbers ! Needed for Russian roulette, if used
    integer,               intent(in   ) :: scatteringOrder
    real,    dimension(:), intent(  out) :: contributions
    integer, dimension(:), intent(  out) :: xIndexF, yIndexF
    
    ! Local variables
    integer :: phaseFunctionIndex, numIntensityDirections, i
    real, dimension(size(thisIntegrator%intensityDirections, 2)) &
            :: projections, scatteringAngles, phaseFunctionVals, tausToBoundary
    integer,  dimension(size(thisIntegrator%intensityDirections, 2)) &
            :: zIndexF
            
    !
    ! Variables for Russian Roulette as applied to intensity calculations
    !   Notation follows H. Iwabuchi, JAS 2006, Vol 63, pp 2324-2339
    !   contributions(:) corresponds to w_n zeta_n/Pi in his Eq 9
    !   We omit the cancelling factors of Pi that appear in his Eqs 9 and 10
    !
    real    :: tauFree, tauMax
    integer :: zIndexMax
    real    :: xTemp, yTemp, zTemp, xPosI, yPosI, zPosI
    real, dimension(size(thisIntegrator%intensityDirections, 2)) &
            :: normalizedPhaseFunc
    
    ! -------------------------------------
    numIntensityDirections = size(thisIntegrator%intensityDirections, 2)
    zIndexMax = size(thisIntegrator%zPosition)
    xPosI = xPos; yPosI = yPos; zPosI = zPos
    
    xIndexF(:) = xIndex; yIndexF(:) = yIndex; zIndexF(:) = zIndex
    
    ! The photon weight as supplied has already been adjusted to account for the 
    !    albedo of the surface or scattering component
    !
    ! We rely on the calling routines setting component to 0 to indicate surface reflection. 
    !   There might be better ways to trap for that. 
    !
    if(component < 1) then ! We're reflecting from the surface
      !
      ! This normalization may only work for Lambertian surfaces
      !   It's supposed to be the ratio of the BRDF to the albedo at 
      !   the incident solar zenith angle
      !
      normalizedPhaseFunc(:) = 1/Pi
    else
      !
      ! Determine the scattering angle (in radians) from the photon direction to each intensity direction
      !   as the acos of dot product of the photon's direction cosine with the intensity direction's cosine
      !
      ! Numerics are such that once in a great while the projection can be slightly larger than 1. 
      !
      projections(:) = matmul(directionCosines, thisIntegrator%intensityDirections(:, :))
      where(abs(projections(:)) > 1.) projections(:) = sign(1., projections(:))
      scatteringAngles(:) = acos(projections)
      
      !
      ! Look up the phase function values in each of the scattering directions 
      !   from the tabulated phase functions
      !
      phaseFunctionIndex = thisIntegrator%phaseFunctionIndex(xIndex, yIndex, zIndex, component)
      if(thisIntegrator%useHybridPhaseFunsForIntenCalcs .and. & 
         scatteringOrder <= thisIntegrator%numOrdersOrigPhaseFunIntenCalcs) then 
        phaseFunctionVals(:) =                                                                                          &
          lookUpPhaseFuncValsFromTable(thisIntegrator%tabulatedOrigPhaseFunctions(component)%values(:, phaseFunctionIndex), &
                                       scatteringAngles)
      else 
        !
        ! If useHybridPhaseFunsForIntenCalcs is false then the tabulated phase functions are the orginals
        ! 
        phaseFunctionVals(:) =                                                                                          &
          lookUpPhaseFuncValsFromTable(thisIntegrator%tabulatedPhaseFunctions(component)%values(:, phaseFunctionIndex), &
                                       scatteringAngles)
      end if 
      normalizedPhaseFunc(:) = phaseFunctionVals(:) / (4 * Pi * abs(thisIntegrator%intensityDirections(3, :)))
    end if 
    
    if(.not. thisIntegrator%useRussianRouletteForIntensity) then 
      !
      ! Find the integrated extinction from the current location to the boundary along each
      !   direction at which intensity is desired. 
      !
      do i = 1, numIntensityDirections
        xTemp = xPosI; yTemp = yPosI; zTemp = zPosI
        call accumulateExtinctionAlongPath(thisIntegrator, thisIntegrator%intensityDirections(:, i), &
                                           xTemp, yTemp, zTemp, xIndexF(i), yIndexF(i), zIndexF(i),  &
                                           tausToBoundary(i))            
      end do
      !
      ! The contribution to intensity in each direction is the product of 
      !   photonWeight * phase function value * 1/(4 pi) * 1/abs(mu) * transmission
      !   The photon weight already contains the single scattering albedo for 
      !   this collision
      ! This contribution will be added to the intensity at the x-y location
      where(tausToBoundary(:) >= 0.) 
        contributions(:) = photonWeight * normalizedPhaseFunc * exp(-tausToBoundary(:)) 
      elsewhere
        ! The extinction subroutine signals errors by producing extinction values < 0 
        !   We could trap this but for the moment we'll just safeguard against disaster
        contributions(:) = 0. 
      end where 
    else
      ! 
      ! Russian roulette for intensity calculations (H. Iwabuchi, JAS 2006, Vol 63, pp 2324-2339) 
      !
      do i = 1, numIntensityDirections
        xTemp = xPosI; yTemp = yPosI; zTemp = zPosI
        tauFree = -log(max(tiny(tauFree), getRandomReal(randomNumbers)))
        !
        !   Small phase function contributions (Iwabuchi Eq 13)
        ! 
        if(Pi * normalizedPhaseFunc(i) <= thisIntegrator%zetaMin) then 
          call accumulateExtinctionAlongPath(thisIntegrator, thisIntegrator%intensityDirections(:, i), &
                                             xTemp, yTemp, zTemp, xIndexF(i), yIndexF(i), zIndexF(i),  &
                                             tausToBoundary(i), tauFree)            
          !
          ! accumulateExtinctionAlongPath stops tracing once tausToBoundary(i) == tauFree
          !    so we use zF to see if the photon has escaped (i.e. if tau <= tauFree)
          !
          if(getRandomReal(randomNumbers) <= Pi * normalizedPhaseFunc(i)/thisIntegrator%zetaMin .and. &
             zIndexF(i) >= zIndexMax) then
            contributions(i) = photonWeight * thisIntegrator%zetaMin/Pi
          else
            contributions(i) = 0. 
          end if 
        else 
          !
          ! Use the full contribution if the optical depth is small and 
          !   play Russian roulette with contributions that undergo large extinction. 
          !   (Iwabuchi Eq 14). 
          !
          tauMax = -log(thisIntegrator%zetaMin/max(tiny(normalizedPhaseFunc), Pi * normalizedPhaseFunc(i)))
          call accumulateExtinctionAlongPath(thisIntegrator, thisIntegrator%intensityDirections(:, i), &
                                             xTemp, yTemp, zTemp, xIndexF(i), yIndexF(i), zIndexF(i),  &
                                             tausToBoundary(i), tauMax)            
          if(zIndexF(i) >= zIndexMax .and. tausToBoundary(i) >= 0.) then 
            !
            !  This means tau <= tauMax
            ! 
            contributions(i) =  photonWeight * normalizedPhaseFunc(i) * exp(-tausToBoundary(i)) 
          else if (tausToBoundary(i) >= 0.) then 
            call accumulateExtinctionAlongPath(thisIntegrator, thisIntegrator%intensityDirections(:, i), &
                                               xTemp, yTemp, zTemp, xIndexF(i), yIndexF(i), zIndexF(i),  &
                                               tausToBoundary(i), tauFree)            
            !
            ! accumulateExtinctionAlongPath stops tracing once tausToBoundary(i) == tauFree
            !    so we use zF to see if the photon has escaped (i.e. if tau <= tauMax + tauFree)
            !
            if(zIndexF(i) >= zIndexMax) then 
              contributions(i) = photonWeight * thisIntegrator%zetaMin/Pi
            else 
              contributions(i) = 0. 
            end if
          else
            !
            ! This branch means the ray tracing to tauMax failed
            !
            contributions(i) = 0.
          end if
        end if
      end do 
    end if 
    
    if(thisIntegrator%limitIntensityContributions) then 
      !
      ! Limit local estimate contribution to intensity; keep track of excess 
      !   so it can be redistributed after the fact
      !
      where(contributions(:) > thisIntegrator%maxIntensityContribution) 
        thisIntegrator%intensityExcess(:, component) = &
          thisIntegrator%intensityExcess(:, component) + & 
          contributions(:) - thisIntegrator%maxIntensityContribution
        contributions(:) = thisIntegrator%maxIntensityContribution
      end where
    end if 
    
  end subroutine computeIntensityContribution
  !------------------------------------------------------------------------------------------
  pure function lookUpPhaseFuncValsFromTable(tablulatedPhaseFunction, scatteringAngles) &
                    result(phaseFunctionVals)
    real, dimension(:),         intent(in ) :: tablulatedPhaseFunction, scatteringAngles
    real, dimension(size(scatteringAngles)) :: phaseFunctionVals
    !
    ! Find the value of phase function for the component doing the scattering at the 
    !   angle from the photon direction to each intensity direction.  
    ! Interpolate linearly in angle between the two closest tabulated values - 
    !   here we have passed in only the single relevant tabulated phase function
    !   (i.e. we have picked the correct element from the table corresponding to the 
    !    right component.) 
    ! The phase functions are tabulated equally in angle, with the first element corresponding 
    !   to scattering angle 0 and the last to scattering angle of pi. 
    !
    
    ! Local variables
    integer                                    :: nAngleSteps
    real                                       :: deltaTheta
    integer, dimension(size(scatteringAngles)) :: angleIndicies
    real,    dimension(size(scatteringAngles)) :: weights

    nAngleSteps = size(tablulatedPhaseFunction)
    deltaTheta  = Pi / (nAngleSteps - 1) 
    angleIndicies(:) = int(scatteringAngles(:) / deltaTheta) + 1

    !
    ! Since angleIndices could be >= nAngleSteps...
    !
    where(angleIndicies(:) < nAngleSteps)  
      ! Angle at corresponding index 
      weights(:) = 1. - (scatteringAngles(:) - (angleIndicies(:) - 1) * deltaTheta)/ deltaTheta
      phaseFunctionVals(:) =                                            &
              weights(:)  * tablulatedPhaseFunction(angleIndicies(:)) + &
        (1. - weights(:)) * tablulatedPhaseFunction(angleIndicies(:) + 1)
    elsewhere
      phaseFunctionVals(:) = tablulatedPhaseFunction(nAngleSteps) 
    end where
    
    
  end function lookUpPhaseFuncValsFromTable
  !------------------------------------------------------------------------------------------
  pure subroutine accumulateExtinctionAlongPath(thisIntegrator, directionCosines,         &
                                                xPos, yPos, zPos, xIndex, yIndex, zIndex, &
                                                extAccumulated, extToAccumulate)
    !
    ! Trace through the medium in a given direction accumulating extinction 
    !   along the path. Tracing stops either when the boundary is reached or
    !   when the accumulated extinction reaches the (optional) value extToAccumulate.
    !   Reports the final position and optical path accumulated. 
    !
    type(integrator),   intent(in   ) :: thisIntegrator
    real, dimension(3), intent(in   ) :: directionCosines
    real,               intent(inout) :: xPos, yPos, zPos
    integer,            intent(inout) :: xIndex, yIndex, zIndex
    real,               intent(  out) :: extAccumulated
    real, optional,     intent(in   ) :: extToAccumulate
    
    ! Local variables
    integer               :: nXcells, nYcells, nZcells
    real                  :: thisStep, thisCellExt, totalPath, z0, zMax
    real,    dimension(3) :: step
    integer, dimension(3) :: SideIncrement, CellIncrement
    
    extAccumulated = 0.; totalPath = 0.

     ! Make some useful parameters
    nXcells = size(thisIntegrator%xPosition) - 1
    nYcells = size(thisIntegrator%yPosition) - 1
    nZcells = size(thisIntegrator%zPosition) - 1
    ! Side increment is 1 where directionCosines is > 0; 0 otherwise
    sideIncrement(:) = merge(1, 0, directionCosines(:) >= 0.) 
    ! Cell increment is +1 where direction cosines > 0; -1 otherwise
    CellIncrement(:) = merge(1, -1, directionCosines(:) >= 0.) 
    
    z0   = thisIntegrator%zPosition(1)
    zMax = thisIntegrator%zPosition(nZCells + 1) 

    accumulationLoop: do
      !
      ! step - how far away is the closest cell boundary along the direction of travel
      !     
      ! How big a step must we take to reach the next edge?  Go the right direction (+ or -) 
      !    and find out for each dimension (but don't divide by zero).
      !
      where(abs(directionCosines) >= 2. * tiny(directionCosines))
        step(:) = (/ thisIntegrator%xPosition(xIndex + SideIncrement(1)) - xPos,     &
                     thisIntegrator%yPosition(yIndex + SideIncrement(2)) - yPos,     &
                     thisIntegrator%zPosition(zIndex + SideIncrement(3)) - zPos /) / &
                  directionCosines(:)
      elsewhere
        step(:) = huge(step)
      end where

       ! The step length across the cell is the smallest of the three directions
       !   We guard against step being negative, which can happen if the 
       !   direction cosine or the distance to the boundary is very small 
       !
      thisStep = minval(step(:))
      if (thisStep <= 0.0) then
        extAccumulated = -2.0  ! Error flag
        exit accumulationLoop
      end if

       ! If this cell pushes the optical path past the desired extToAccumulate,
       !  then find how large a step is needed, take it, and exit
       
      thisCellExt = thisIntegrator%totalExt(xIndex, yIndex, zIndex)

      if (present(extToAccumulate)) then 
        if(extAccumulated + thisStep * thisCellExt > extToAccumulate) then
          thisStep = (extToAccumulate - extAccumulated) / thisCellExt
          xPos = xPos + thisStep * directionCosines(1)
          yPos = yPos + thisStep * directionCosines(2)
          zPos = zPos + thisStep * directionCosines(3)
          totalPath = totalPath + thisStep
          extAccumulated = extToAccumulate
          exit accumulationLoop
        end if
      end if

       ! Add this cell crossing to the accumulated optical path and distance
       
      extAccumulated = extAccumulated + thisStep*thisCellExt
      totalPath = totalPath + thisStep

       ! Determine which side of the cell we're going to hit first (the smallest step).
       !   Set the position to that side of the cell and increment the cell indices to 
       !   the next cell, which means extra code for periodic boundaries in X and Y.
       ! If we wind up within spacing() of the coming boundary just say the position is 
       !   in the next cell. (This slows things down but protects againt rounding errors.)
       
      if(step(1) <= thisStep) then
        xPos = thisIntegrator%xPosition(xIndex + SideIncrement(1)) 
        xIndex = xIndex + CellIncrement(1)
      else
        xPos = xPos + thisStep * directionCosines(1)
        if(abs(thisIntegrator%xPosition(xIndex + sideIncrement(1)) - xPos) <= 2 * spacing(xPos)) &
             xIndex = xIndex + cellIncrement(1)      !
      end if
      
      if(step(2) <= thisStep) then
        yPos = thisIntegrator%yPosition(yIndex+SideIncrement(2))
        yIndex = yIndex + CellIncrement(2)
      else
        yPos = yPos + thisStep * directionCosines(2)
        if(abs(thisIntegrator%yPosition(yIndex + sideIncrement(2)) - yPos) <= 2 * spacing(yPos)) &
                yIndex = yIndex + cellIncrement(2)
      end if

      if(step(3) <= thisStep) then
        zPos = thisIntegrator%zPosition(zIndex+SideIncrement(3))
        zIndex = zIndex + CellIncrement(3)
      else
        zPos = zPos + thisStep * directionCosines(3)
        if(abs(thisIntegrator%zPosition(zIndex + sideIncrement(3)) - zPos) <= 2 * spacing(zPos)) &
          zIndex = zIndex + cellIncrement(3)
      end if

      !
      ! Enforce periodicity
      !
      if (xIndex <= 0) then
        xIndex = nXcells
        xPos = thisIntegrator%xPosition(xIndex+1) + cellIncrement(1) * 2 * spacing(xPos)
      else if (xIndex >= nXcells+1) then
        xIndex = 1
        xPos = thisIntegrator%xPosition(xIndex) + cellIncrement(1) * 2 * spacing(xPos)
      end if
      
      if (yIndex <= 0) then
        yIndex = nYcells
        yPos = thisIntegrator%yPosition(yIndex+1) + cellIncrement(1) * 2 * spacing(yPos)
      else if (yIndex >= nYcells+1) then
        yIndex = 1
        yPos = thisIntegrator%yPosition(yIndex) + cellIncrement(1) * 2 * spacing(yPos)
      end if

      !
      ! Out the top? 
      !
      if(zIndex > nZCells) then
        zPos = zMax + 2 * spacing(zMax) 
        exit accumulationLoop
      end if 
      
      !
      ! Hit the bottom? 
      !
      if(zIndex < 1) then
        zPos = z0
        exit accumulationLoop
      end if 

    end do accumulationLoop
  end subroutine accumulateExtinctionAlongPath  
  !------------------------------------------------------------------------------------------
  subroutine tabulateInversePhaseFunctions(thisIntegrator, status)
    !
    ! Tabulate the inverse (cummulative) phase functions (i.e. scattering angle
    !   as a function of the position in the cumulative distribution). 
    !
    type(integrator),  intent(inout) :: thisIntegrator
    type(ErrorMessage),intent(inout) :: status
    
    ! Local variables
    integer                            :: i, numComponents, nEntries, nSteps
    logical                            :: computeThisTable
    real, dimension(:, :), allocatable :: tempMatrix

    numComponents = size(thisIntegrator%forwardTables)
    if(.not. associated(thisIntegrator%inversePhaseFunctions)) &
        allocate(thisIntegrator%inversePhaseFunctions(numComponents))
        
    componentLoop: do i = 1, numComponents
      !
      ! Does the table already exist at high enough resolution? 
      !
      computeThisTable = .true. 
      if(associated(thisIntegrator%inversePhaseFunctions(i)%values)) then
        if(thisIntegrator%inversePhaseFunctions(i)%numX >= thisIntegrator%minInverseTableSize) computeThisTable = .false.
      end if
        
      if(computeThisTable) then
        !
        ! Compute at the minimum desired "resolution" (number of intervals bet. 0 and 1)
        !   computeInversePhaseFuncTable expects a simple 2D array, dimensioned nSteps, nEntries
        !   We need to copy this into our "matrix" type
        !
        call getInfo_phaseFunctionTable(thisIntegrator%forwardTables(i), nEntries = nEntries, status = status)
        if(stateIsFailure(status)) exit componentLoop
        nSteps = thisIntegrator%minInverseTableSize
        allocate(tempMatrix(nSteps, nEntries))
        !
        ! Use the code in the inversePhaseFunctions module to compute the inverse table 
        !
        call computeInversePhaseFuncTable(thisIntegrator%forwardTables(i), tempMatrix, status)
        if(stateIsFailure(status)) exit componentLoop
        thisIntegrator%inversePhaseFunctions(i) = new_Matrix(tempMatrix)
        deallocate(tempMatrix)
      end if
    end do componentLoop
    
    if(stateIsFailure(status)) then
      call setStateToFailure(status, "tabulateInversePhaseFunctions: failed on component" // trim(intToChar(i)) )
    else
      call setStateToSuccess(status)
    end if 
  
  end subroutine tabulateInversePhaseFunctions
  !------------------------------------------------------------------------------------------
  subroutine tabulateForwardPhaseFunctions(thisIntegrator, status)
    type(integrator),  intent(inout) :: thisIntegrator
    type(ErrorMessage),intent(inout) :: status
    
    ! Local variables
    integer                            :: i, j, numComponents, nEntries, nSteps
    logical                            :: computeThisTable
    real, dimension(:),    allocatable :: angles
    real, dimension(:, :), allocatable :: tempMatrix

    !
    ! We store the forward tables as matrices evenly spaced in angle to make interpolation simple.
    !
    numComponents = size(thisIntegrator%forwardTables)
    if(.not. associated(thisIntegrator%tabulatedPhaseFunctions)) &
        allocate(thisIntegrator%tabulatedPhaseFunctions(numComponents))
    if(.not. associated(thisIntegrator%tabulatedOrigPhaseFunctions)) &
        allocate(thisIntegrator%tabulatedOrigPhaseFunctions(numComponents))

    componentLoop: do i = 1, numComponents
      !
      ! Does the matrix form of the table already exist high enough resolution? 
      !   The tabulated phase functions and the original versions are always computed at the same resolution 
      !
      computeThisTable = .true. 
      if(associated(thisIntegrator%tabulatedPhaseFunctions(i)%values)) then
        if(thisIntegrator%tabulatedPhaseFunctions(i)%numX >= thisIntegrator%minForwardTableSize) computeThisTable = .false.
      end if 
      
      if(computeThisTable) then 
        !
        ! Compute the phase function values at nSteps points equally spaced in angle from 0 to pi radians
        !
        nSteps = thisIntegrator%minForwardTableSize
        call getInfo_PhaseFunctionTable(thisIntegrator%forwardTables(i), nEntries = nEntries, status = status)
        if(stateIsFailure(status)) exit componentLoop
        allocate(angles(nSteps), tempMatrix(nSteps, nEntries))
        angles(:) = (/ (j, j = 0, nSteps - 1) /) / real(nSteps - 1) * Pi
        call getPhaseFunctionValues(thisIntegrator%forwardTables(i), angles(:), tempMatrix(:, :), status)
        
        thisIntegrator%tabulatedOrigPhaseFunctions(i) = new_Matrix(tempMatrix)
        
        if(thisIntegrator%useHybridPhaseFunsForIntenCalcs) then 
          if(thisIntegrator%hybridPhaseFunWidth > 0.)                                   &
            tempMatrix(:, :) = computeHydridPhaseFunctions(angles(:), tempMatrix(:, :), &
                                                           thisIntegrator%hybridPhaseFunWidth)
        end if
        !
        ! Copy the tabulated phase functions into a matrix
        !
        thisIntegrator%tabulatedPhaseFunctions(i) = new_Matrix(tempMatrix)
        deallocate(angles, tempMatrix)
      end if
    end do componentLoop
  
    if(stateIsFailure(status)) then
      call setStateToFailure(status, "tabulatePhaseFunctions: failed on component" // trim(intToChar(i)) )
    else
      call setStateToSuccess(status)
    end if 
  end subroutine tabulateForwardPhaseFunctions
  !------------------------------------------------------------------------------------------
  pure function computeHydridPhaseFunctions(angles, values,  GaussianWidth) result(newValues)
    real, dimension(:),    intent( in) :: angles         ! phase function angles in radians
    real, dimension(:, :), intent( in) :: values         ! Phase  function, nEntries by nAngles
    real,                  intent( in) :: GaussianWidth  ! In degrees
    real, dimension(size(values, 1), size(values, 2) ) :: newValues
    !
    ! Creates a phase function that's a hybrid of the original and a Gaussian of
    !   specified width. The Gaussian part replaces the forward  peak, and is
    !   continuous with the original phase function

    ! Local variables.
    real, dimension(size(angles)) :: gaussianValues, angleCosines
    real                          :: P0, lowDiff, upDiff, midDiff
    integer                       :: nAngles, transitionIndex
    integer                       :: i, lowerBound, upperBound,  midPoint, increment

    nAngles = size(angles)
    angleCosines(:) = cos(angles(:))
    !
    ! Gaussian phase function values - we won't need most of these
    !   but it's easier to compute them all than keep adding
    !
    gaussianValues(:) = exp(-( angles(:)/(GaussianWidth * Pi/180) )**2)

    ! First set the output phase function in input one in case there is no root
    newValues(:, :) = values(:, :)
    entryLoop: do i = 1, size(values, 2)

      ! Set the lower transition bound according to width of the Gaussian
      lowerBound = findIndex(GaussianWidth * Pi/180., angles(:)) + 1
      if(lowerBound >= nAngles - 2) exit entryLoop 
      ! We haven't found the position of the Gaussian width in the table

      ! We want the transition angle at which the two phase functions are the same
      !   (i.e. the difference between them is 0).
      !
      ! First we "hunt", expanding the range in which we're trying to bracket the value
      !
      lowDiff = phaseFuncDiff(angleCosines(:), values(:, i),  gaussianValues(:), lowerBound)
      increment = 1
      huntingLoop: do
        upperBound = min(lowerBound + increment, nAngles - 1)
        upDiff = phaseFuncDiff(angleCosines(:), values(:, i),  gaussianValues(:), upperBound)

        if (lowerBound == nAngles - 1) cycle entryLoop  ! There's no root, so use the original phase function
        if (lowDiff * upDiff < 0) exit huntingLoop

        lowerBound = upperBound
        lowDiff = upDiff
        increment = increment * 2
      end do huntingLoop

      ! Bisection: figure out which half of the remaining interval holds the
      !   root, discard the other half, and repeat
      bisectionLoop: do
        if (upperBound <= lowerBound + 1) exit bisectionLoop
        midPoint = (lowerBound + upperBound)/2
        midDiff = phaseFuncDiff(angleCosines(:), values(:, i),  gaussianValues(:), midPoint)
        if (midDiff * upDiff < 0) then
          lowerBound = midPoint
          lowDiff = midDiff
        else
          upperBound = midPoint
          upDiff = midDiff
        end if
      end do bisectionLoop

      transitionIndex = lowerBound
      P0 = computeNormalization(angleCosines(:), values(:, i),  gaussianValues(:), transitionIndex)
      newValues(:transitionIndex,   i) = P0 * gaussianValues(:transitionIndex)
      newValues(transitionIndex+1:, i) = values(transitionIndex+1:, i)
    end do entryLoop

  end function computeHydridPhaseFunctions
  !------------------------------------------------------------------------------------------
  pure function computeNormalization(angleCosines, values,  gaussianValues, transitionIndex) result(P0)
    real, dimension(:), intent(in) :: angleCosines, values,  gaussianValues
    integer,            intent(in) :: transitionIndex
    real                           :: P0

    integer :: nAngles
    real    :: IntegralGaus, IntegralOrig

    ! Normalization for the Gaussian part of the phase function, computed by
    !   forcing the complete phase function to be normalized.
    !
    nAngles = size(angleCosines)
    IntegralGaus = dot_product( &
                   0.5*(gaussianValues(1:transitionIndex-1) + gaussianValues(2:transitionIndex)), &
                   angleCosines(1:transitionIndex-1) - angleCosines(2:transitionIndex) )
    IntegralOrig = dot_product( &
                   0.5*(values(transitionIndex:nAngles-1) + values(transitionIndex+1:nAngles)), &
                   angleCosines(transitionIndex:nAngles-1) - angleCosines(transitionIndex+1:nAngles) )
    if (IntegralOrig >= 2.0) then
      P0 = 1.0/IntegralGaus
    else
      P0 = (2. - IntegralOrig) / IntegralGaus
    end if
  end function computeNormalization
  !------------------------------------------------------------------------------------------
  pure function phaseFuncDiff(angleCosines, values, gaussianValues,  transitionIndex) &
                                        result(d)
    !
    ! Compute the difference between the normalized Gaussian phase function and the
    !   original phase function at the transition index. 
    !
    real, dimension(:), intent(in) :: angleCosines, values,  gaussianValues
    integer,            intent(in) :: transitionIndex
    real                           :: d

    real :: P0

    P0 = computeNormalization(angleCosines, values, gaussianValues,  transitionIndex)
    d = P0 * gaussianValues(transitionIndex) - values(transitionIndex)
  end function phaseFuncDiff
  !------------------------------------------------------------------------------------------
  pure function makeDirectionCosines(mu, phi)
    real,   intent(in) :: mu, phi
    real, dimension(3) :: makeDirectionCosines
    !
    ! Direction cosines S are defined so that
    !  S(1) = sin(theta) * cos(phi), projection in X direction 
    !  S(2) = sin(theta) * sin(phi), projection in Y direction
    !  S(3) = cos(theta),            projection in Z direction
    ! 
    ! Input is mu = cos(theta) and phi
    !
    
    real :: sinTheta, cosPhi, sinPhi
    
    sinTheta = sqrt(1. - mu**2)
    cosPhi   = cos(Phi) 
    sinPhi   = sin(Phi) ! sqrt(1 - cosPhi) is ambiguous at 90, 270 degrees. 
    makeDirectionCosines(:) = (/ sinTheta * cosPhi, sinTheta * sinPhi, mu /)
  end function makeDirectionCosines
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  elemental function makePeriodic(a, aMin, aMax)
    real, intent(in) :: a, aMin, aMax
    real             :: makePeriodic
    !
    ! Ensure that a position is within domain when the boundary conditions are periodic
    !
    ! makePeriodic = aMin + mod(a - aMin, aMax - aMin)
    ! if(makePeriodic < aMin) makePeriodic = aMax - abs(makePeriodic - aMin)
    makePeriodic = a
    do 
      if(makePeriodic <= aMax .and. makePeriodic > aMin) exit
      if(makePeriodic > aMax) then
        makePeriodic = makePeriodic - (aMax - aMin)
      else if (makePeriodic == aMin) then 
        makePeriodic = aMax 
      else
        makePeriodic = makePeriodic + (aMax - aMin)
      end if 
    end do
  end function makePeriodic
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  SUBROUTINE NEXT_DIRECT (randomNumbers, scatteringCosine, S)
    !
    !   Finds new set of direction cosines S from the original set and 
    !   the cosine of the scattering angle. 
    !   Algorithm from Marchuk et al (1980); implementation by Frank Evans
    !
    type(randomNumberSequence), intent(inout) :: randomNumbers
    real,                       intent(in   ) :: scatteringCosine
    real, dimension(:),         intent(inout) :: S
     
    REAL :: D, AX, AY, B
  
    D = 2.0
    DO WHILE (D .GT. 1.0)
      AX = 1.0 - 2.0*getRandomReal(randomNumbers)
      AY = 1.0 - 2.0*getRandomReal(randomNumbers)
      D = AX**2 + AY**2
    ENDDO
    B = SQRT((1.0 - scatteringCosine**2)/D)
    AX = AX*B
    AY = AY*B
    B = S(1)*AX-S(2)*AY
    D = scatteringCosine - B/(1.0+ABS(S(3)))
    
    S(1) = S(1)*D + AX
    S(2) = S(2)*D - AY
    S(3) = S(3)*scatteringCosine - sign(b, s(3) * b)
  END SUBROUTINE NEXT_DIRECT
  !------------------------------------------------------------------------------------------
  function new_Matrix(array) 
    real, dimension(:, :) :: array
    type(matrix)          :: new_Matrix
    
    new_Matrix%numX         = size(array, 1)
    new_Matrix%numY         = size(array, 2)
    allocate(new_Matrix%values(size(array, 1), size(array, 2)))
    new_Matrix%values(:, :) = array(:, :)
  end function new_Matrix
  !------------------------------------------------------------------------------------------
  subroutine finalize_Matrix(thisMatrix)
    type(matrix), intent(out) :: thisMatrix
    
    thisMatrix%numX = 0; thisMatrix%numY = 0
    if(associated(thisMatrix%values)) deallocate(thisMatrix%values)
  end subroutine finalize_Matrix
  !------------------------------------------------------------------------------------------
end module monteCarloRadiativeTransfer
