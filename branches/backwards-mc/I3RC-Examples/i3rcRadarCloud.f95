! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision$, $Date$
! $URL$

program i3rcRadarCloud
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
  integer,                parameter :: nColumns   = 640, nLayers = 54, nSSAs = 2
  real,                   parameter :: deltaX = 50., deltaZ = 45., g = 0.85
  real, dimension(nSSAs), parameter :: SSAs = (/ 1., 0.99 /) 
  integer,                parameter :: nLegendreCoefficients = 299, nScatteringAngles = 1801
  

  ! Other parameters
  integer,             parameter :: fileUnit = 10
  character(len =  5), parameter :: i3rcDir = "Data/"
  character(len = 20), parameter ::           tauFile = "mmcr_tau_32km_020898"
  character(len = 12), parameter :: phaseFunctionFile = "C.1_PF      ", &
                                         legendreFile = "C.1_leg_coef" 
  character(len =  3), dimension(2), &
                       parameter :: phaseFunctionNames = (/ "HG_", "C1_" /)
  character(len = 12), dimension(2), &
                       parameter :: ssaNames = (/ "NonAbsorbing", "Absorbing   " /)
  character(len = 11), parameter :: outputFilePrefix =  "RadarCloud_" 
  
  ! Variables
  integer                           :: i, j
  real,    dimension(nLegendreCoefficients) :: coefficients
  real,    dimension(nScatteringAngles)     :: scatteringAngle, value
    
  real,    dimension(nColumns, 1, nLayers)  :: extinction, singleScatteringAlbedo
  integer, dimension(nColumns, 1, nLayers)  :: phaseFunctionIndex
  
  type(ErrorMessage)              :: status
  type(phaseFunction)             :: HG, C1_tabulated, C1_expanded
  type(phaseFunctionTable)        :: HGTable, C1Table_Tabulated, C1Table_Expanded
  type(InversePhaseFunctionTable) :: inverseTable
  type(domain)                    :: radarCloud

  ! ------------------------------------------
  call initializeState(status)
  ! 
  ! Phase functions
  !
  !   Henyey-Greenstein
  HG = new_PhaseFunction(g**(/ (i, i = 1, nLegendreCoefficients )/), status = status)
  call printStatus(status)
    
  ! Tabulated version of C1, read in from file
  open(unit = 10, file = trim(i3rcDir) // trim(phaseFunctionFile), status = 'old')
  do i = 1, nScatteringAngles
    read(fileUnit, *) scatteringAngle(i), value(i)
  end do
  close(fileUnit)
  C1_tabulated = new_PhaseFunction(scatteringAngle * acos(-1.)/180., value, status = status)
  call printStatus(status)
  
  ! Legendre expansion of C1 - skip the first value and read from file
  open(unit = 10, file = trim(i3rcDir) // trim(legendreFile), status = 'old')
  read(10, *) coefficients(1)
  do i = 1, nLegendreCoefficients
    read(fileUnit, *) coefficients(i)
  end do
  close(fileUnit)
  ! Frank supplies his coefficients as (2 l + 1) * mine
  C1_expanded = new_PhaseFunction(coefficients / (/ (2*i + 1, i = 1, nLegendreCoefficients) /), &
                                  status = status)
  call printStatus(status)
  !
  ! Phase function tables
  !
  HGTable           = new_PhaseFunctionTable((/ HG /), (/ 1. /), &
                                             tableDescription = "Henyey-Greenstein with g = 0.85", &
                                             status =  status)
  call printStatus(status)
  C1Table_Tabulated = new_PhaseFunctionTable(scatteringAngle * acos(-1.)/180.,     &
                                             spread(value, 2, nCopies = 1),        &
                                             key = (/ 1. /),                       &
                                             tableDescription = "Dermeindjian C1", &
                                             status = status)
  call printStatus(status)
  C1Table_Expanded  = new_PhaseFunctionTable((/ C1_expanded /), (/ 1. /), &
                                             tableDescription = "Dermeindjian C1",  &
                                             status = status)
  call printStatus(status)
  
  open(unit = FileUnit, file = trim(i3rcDir) // trim(tauFile), status = 'old')
  do j = nLayers, 1, -1
    read(fileUnit, '(640f8.3)') extinction(:, 1, j)
  end do
  close(fileUnit)
  ! This I3RC case has defined optical depth in each cell, and we need to convert it 
  !   to extinction
  extinction(:, :, :) = extinction(:, :, :) / deltaZ
  phaseFunctionIndex(:, :, :) = 1
  singleScatteringAlbedo(:, :, :) = 1. 
  
  !
  ! Define the domain
  !
  radarCloud =                                                             &
    new_Domain(xPosition = deltaX * (/ 0., (real(i), i = 1, nColumns) /), &
               yPosition = (/ 0., deltaX * nColumns /),                   &
               zPosition = deltaZ * (/ 0., (real(i), i = 1, nLayers) /) , &
               status = status)
  !
  ! Add the clouds. We're just going to replace them below...
  !
  call addOpticalComponent(radarCloud, "cloud: non-absorbing",  &
                           extinction, singleScatteringAlbedo, &
                           phaseFunctionIndex, HGTable, status = status)
  call printStatus(status)
  
  !
  ! Loop over each combination of phase function and single scattering albedo 
  !   and write each domain
  !
  do i = 1, size(SSAs) 
    singleScatteringAlbedo = SSAs(i)
    call replaceOpticalComponent(radarCloud, 1, "cloud: HG",  &
                                 extinction, singleScatteringAlbedo, &
                                 phaseFunctionIndex, HGTable, status = status)
    call printStatus(status)
    call write_Domain(radarCloud, trim(i3rcDir) // trim(outputFilePrefix) // &
                            trim(phaseFunctionNames(1)) // trim(ssaNames(i)) // ".opt", status)
    call printStatus(status)

    call replaceOpticalComponent(radarCloud, 1, "cloud: C1",  &
                                 extinction, singleScatteringAlbedo, &
                                 phaseFunctionIndex, C1Table_tabulated, status = status)
    call printStatus(status)
    call write_Domain(radarCloud, trim(i3rcDir) // trim(outputFilePrefix) // &
                            trim(phaseFunctionNames(2)) // trim(ssaNames(i)) // ".opt", status)
    call printStatus(status)
  end do
  
end program i3rcRadarCloud
