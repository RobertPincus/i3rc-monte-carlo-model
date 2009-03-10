! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision$, $Date$
! $URL$
program i3rcStepCloud
  use ErrorMessages
  use scatteringPhaseFunctions
  use inversePhaseFunctions
  use opticalProperties
  use userInterface
  implicit none
  !
  ! Write the domains corresponding to the I3RC step cloud. 
  !   The domain is .5 km wide and 32 columns across. The first 16 
  !   columns have optical depth 2, the second 16 are optical depth 18. 
  !   We use a Henyey-Greenstein phase function with g=0.85. 
  !   We write two domains, one with single scattering albedo = 0.99 
  !   and the other entirely non-absorbing. 
  ! Two parameters weed need that aren't specified are the number of 
  !    Legendre moments and the physical thickness of the cloud. 
  !
  
  ! I3RC parameters
  real,    parameter :: domainSize = 500, g = 0.85
  integer, parameter :: nColumns   = 32, nSSAs = 2
  real, dimension(nSSAs), &
           parameter :: SSAs = (/ 1., 0.99 /) 
  

  ! Other parameters
  integer, parameter :: nLegendreCoefficients = 64, nLayers = 32
  real,    parameter :: physicalThickness = 250
  real               :: deltaX, deltaZ 
  character(len = 31), dimension(2), &
           parameter :: fileNames = (/ "Data/StepCloud_NonAbsorbing.opt", &
                                       "Data/StepCloud_Absorbing.opt   " /)
  
  ! Variables
  integer                           :: i
    
  real,    dimension(nColumns, 1, nLayers) :: extinction, singleScatteringAlbedo
  integer, dimension(nColumns, 1, nLayers) :: phaseFunctionIndex
  
  type(ErrorMessage)              :: status
  type(phaseFunction)             :: phase
  type(phaseFunctionTable)        :: table
  ! type(InversePhaseFunctionTable) :: inverseTable
  type(domain)                    :: stepCloud

  ! ------------------------------------------
  phase = &
    new_PhaseFunction(g**(/ (i, i = 1, nLegendreCoefficients )/), status = status)
  table = &
    new_PhaseFunctionTable((/ phase /), key = (/ 1. /), status = status)
  if(.not. stateIsFailure(status)) call setStateToCompleteSuccess(status)
  
  deltaX = domainSize/real(nColumns)
  deltaZ = physicalThickness/real(nLayers)
  extinction(:, 1, :) = spread((/ (2, i = 1, nColumns/2), (18, i = 1, nColumns/2) /), &
                               dim = 2, nCopies = nLayers) / physicalThickness
  phaseFunctionIndex(:, :, :) = 1
  
  singleScatteringAlbedo(:, :, :) = SSAs(1)
  
  !
  ! Define the domain
  !
  stepCloud =                                                            &
    new_Domain(xPosition = deltaX * (/ 0., (real(i), i = 1, nColumns) /), &
               yPosition = (/ 0., 500.0 /),                               &
               zPosition = deltaZ * (/ 0., (real(i), i = 1, nLayers) /) , &
               status = status)
  !
  ! Add the clouds
  !
  call addOpticalComponent(stepCloud, "cloud: non-absorbing",  &
                           extinction, singleScatteringAlbedo, &
                           phaseFunctionIndex, table, status = status)
  call printStatus(status)
  !
  ! Write it to the file
  !
  call write_Domain(stepCloud, trim(fileNames(1)), status = status)
  print *, "Writing non-absorbing domain"; call printStatus(status)
  if(stateIsSuccess(status)) call setStateToCompleteSuccess(status)
  
  !
  ! Now write out the same domain but with SSA = =.99
  !
  singleScatteringAlbedo(:, :, :) = SSAs(2)
  call replaceOpticalComponent(stepCloud, 1, "cloud: absorbing", &
                               extinction, singleScatteringAlbedo,   &
                               phaseFunctionIndex, table, status = status)
  call write_Domain(stepCloud, trim(fileNames(2)), status = status)
  print *, "Writing absorbing domain"; call printStatus(status)
  if(stateIsSuccess(status)) call setStateToCompleteSuccess(status)
end program i3rcStepCloud
