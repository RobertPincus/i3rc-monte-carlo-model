! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision: 1.7 $, $Date: 2009/03/09 19:17:22 $
! $Name: Cornish-Gilliflower $
program i3rcLandsatCloud
  use ErrorMessages
  use scatteringPhaseFunctions
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
  integer,                parameter :: nX= 128, nY = 128, nSSAs = 2
  real,                   parameter :: deltaXY = 30., g = 0.85
  real, dimension(nSSAs), parameter :: SSAs = (/ 1., 0.99 /) 
  

  ! Other parameters
  integer,             parameter :: nLegendreCoefficients = 299
  integer,             parameter :: deltaZ = 20, maxThickness = 2380 ! max thickness in file
  integer,             parameter :: nLayers = (maxThickness + deltaZ/2)/deltaZ
  
  integer,             parameter :: fileUnit = 10
  character(len =  5), parameter :: i3rcDir = "Data/"
  character(len = 19), parameter :: tauFile = "scene43.tau.128x128", &
                                     dzFile = "scene43.dz.128x128 "
  character(len = 13), parameter :: outputFilePrefix =  "LandsatCloud_" 
  
  ! Variables
  integer                           :: i, j
  
  real,    dimension(nX, nY)           :: opticalDepth, thickness
  real,    dimension(nX, nY, nLayers)  :: extinction = 0., singleScatteringAlbedo = 0.
  integer, dimension(nX, nY, nLayers)  :: phaseFunctionIndex = 0
  
  type(ErrorMessage)       :: status
  type(phaseFunction)      :: HG
  type(phaseFunctionTable) :: HGTable
  type(domain)             :: landsatCloud

  ! ------------------------------------------
  call initializeState(status)
  ! 
  ! Phase function and phsae function table
  !
  !   Henyey-Greenstein
  HG = new_PhaseFunction(g**(/ (i, i = 1, nLegendreCoefficients )/), status = status)
  call printStatus(status)
  HGTable           = new_PhaseFunctionTable((/ HG /), (/ 1. /), &
                                             tableDescription = "Henyey-Greenstein with g = 0.85", status = status)
  call printStatus(status)

  !
  !  Read in optical depth and geometrical thickness arrays
  !
  open(unit = FileUnit, file = trim(i3rcDir) // trim(tauFile), status = 'old')
  do i = 1, nY
    read(fileUnit, '(128f7.2)') opticalDepth(:, i)
  end do
  close(fileUnit)
  
  open(unit = FileUnit, file = trim(i3rcDir) // trim(dzFile),  status = 'old')
  do i = 1, nY
    read(fileUnit, '(128f7.2)') thickness(:, i)
  end do
  close(fileUnit)
  thickness(:, :) = thickness(:, :) * 1000. ! km to m conversion

  !
  ! Optical depth and thickness should be greater then 0 in all the same places
  !
  if(any(thickness(:, :) > 0. .neqv. opticalDepth(:, :) > 0.)) &
    print *, "Thickness and optical depth not zero in all the same places." 
  
  !
  ! Define the domain
  !
  landsatCloud =                                                            &
    new_Domain(xPosition = deltaXY * (/ (real(i), i = 0, nX)      /),       &
               yPosition = deltaXY * (/ (real(i), i = 0, nY)      /),       &
               zPosition = deltaZ  * (/ (real(i), i = 0, nLayers) /) + 200, &
               status = status)
  !
  ! Cloud extinction, single scattering albedo, phase function
  !
  forall(i = 1:nX, j = 1:nY)
    extinction(i, j, :nint(thickness(i, j) / deltaZ)) =                      &
        merge(opticalDepth(i, j)/(nint(thickness(i, j) / deltaZ) * deltaZ), &
              0., opticalDepth(i, j) > tiny(opticalDepth))
  end forall
  where(extinction(:, :, :) > 0.) 
    singleScatteringAlbedo(:, :, :) = SSAs(1)
    phaseFunctionIndex(:, :, :) = 1
  end where
  
  !
  ! Let's just be sure we have the same optical depth, shall we ?
  !
  if(any(abs(opticalDepth - count(extinction > 0., dim = 3) * deltaZ * extinction(:, :, 1)) >             &
         spacing(opticalDepth)))                                                                          &
     print *, minval(abs(opticalDepth - count(extinction > 0., dim = 3) * deltaZ * extinction(:, :, 1))), &
              maxval(abs(opticalDepth - count(extinction > 0., dim = 3) * deltaZ * extinction(:, :, 1)))

  !
  ! Add the cloud to the domain, write it out
  !
  call addOpticalComponent(landsatCloud, "cloud: non-absorbing", &
                           extinction, singleScatteringAlbedo, &
                           phaseFunctionIndex, HGTable, status = status)
  call printStatus(status)
  call write_Domain(landsatCloud, &
                          trim(i3rcDir) // trim(outputFilePrefix) // "NonAbsorbing.opt", status)
  call printStatus(status)
  
  !
  ! Now the absorbing cloud. 
  !
  where(extinction(:, :, :) > 0.) &
    singleScatteringAlbedo(:, :, :) = SSAs(2)
  call replaceOpticalComponent(landsatCloud, 1, "cloud: absorbing",  &
                               extinction, singleScatteringAlbedo, &
                               phaseFunctionIndex, HGTable, status = status)
  call printStatus(status)
  call write_Domain(landsatCloud, &
                          trim(i3rcDir) // trim(outputFilePrefix) // "Absorbing.opt", status)
  call printStatus(status)  
end program i3rcLandsatCloud
