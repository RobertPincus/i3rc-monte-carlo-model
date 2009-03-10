! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision$, $Date$
! $URL$
  ! This module represent surface  models in which the reflectance 
  !   (the "bi-directional reflectance distribution function," or BRDF) can be modeled with
  !    just a few parameters. These parameters can vary with horizontal position.
  ! The BRDF is a function of the incoming and outgoing polar cosines and azimuthal 
  !   angles.
  ! This particular module implements Lambertian surface reflectance but 
  !    can provide a template for more complicated surface BRDFs. 
  ! The horizontal coordinate system is local to the surface reflection object. 
  !   In particular, it's up to the user to guarantee that the coordinate system
  !   used to define the surface model is the same as the one used to define 
  !   the atmosphere. 
! --------------------------------------------
module surfaceProperties
  use ErrorMessages
  use numericUtilities
  implicit none
  private
  !------------------------------------------------------------------------------------------
  ! Module constants
  !------------------------------------------------------------------------------------------
  integer, parameter :: numberOfParameters = 1
  !------------------------------------------------------------------------------------------
  ! Type (object) definitions
  !------------------------------------------------------------------------------------------
  type surfaceDescription
    private
    real, dimension(:),       pointer :: xPosition, yPosition
    real, dimension(:, :, :), pointer :: BRDFParameters
  end type surfaceDescription
  
  !------------------------------------------------------------------------------------------
  ! Overloading
  !------------------------------------------------------------------------------------------
  !  Surface properties can vary with location or be constant. 
  ! 
  interface new_SurfaceDescription
    module procedure newSurfaceDescriptionXY, newSurfaceUniform
  end interface ! new_SurfaceDescription
  
  !------------------------------------------------------------------------------------------
  ! What is visible? 
  !------------------------------------------------------------------------------------------
  public :: surfaceDescription
  public :: new_SurfaceDescription, copy_surfaceDescription, finalize_SurfaceDescription, &
            isReady_surfaceDescription, &
            computeSurfaceReflectance
contains  
  !------------------------------------------------------------------------------------------
  ! Initialization: Routines to create new surface description objects
  !------------------------------------------------------------------------------------------
  function newSurfaceDescriptionXY(surfaceParameters, xPosition, yPosition, status) &
    result(thisSurfaceDescription)
    real, dimension(:, :, :), intent(in   ) :: surfaceParameters
    real, dimension(:),       intent(in   ) :: xPosition, yPosition
    type(ErrorMessage),       intent(inout) :: status
    type(surfaceDescription)                :: thisSurfaceDescription 

    integer :: numX, numY
    
    ! Checks : array sizes 
    numX = size(xPosition); numY = size(yPosition)
    if(size(surfaceParameters, 1) /= numberOfParameters) &
      call setStateToFailure(status, "new_SurfaceDescription: Wrong number of parameters supplied for surface BRDF.") 
    if(size(surfaceParameters, 2) /= numX - 1 .or. &
       size(surfaceParameters, 3) /= numY - 1)     &
       call setStateToFailure(status, "new_SurfaceDescription: position vector(s) are incorrect length.")
    if(any(xPosition(2:) - xPosition(:numX-1) <= 0.) .or. &
       any(yPosition(2:) - yPosition(:numY-1) <= 0.))     &
      call setStateToFailure(status, "new_SurfaceDescription: positions must be unique, increasing.")

    ! Check to ensure that surface property parameter values make sense
    !   For a Lambertian surface (our example) the reflectance must be between 0 and 1, incl. 
    if(any(surfaceParameters(1, :, :) < 0.) .or. &
       any(surfaceParameters(1, :, :) > 1.)) &
      call setStateToFailure(status, "new_SurfaceDescription: surface reflectance must be between 0 and 1")
       
       
    if(.not. stateIsFailure(status)) then
      allocate(thisSurfaceDescription%xPosition(numX), &
               thisSurfaceDescription%yPosition(numY), &
               thisSurfaceDescription%BRDFParameters(numberOfParameters, numX - 1, numY - 1)) 
      thisSurfaceDescription%xPosition(:) = xPosition(:)
      thisSurfaceDescription%yPosition(:) = yPosition(:)  
      thisSurfaceDescription%BRDFParameters(:, :, :) = surfaceParameters(:, :, :)
      call setStateToSuccess(status) 
    end if 
  end function newSurfaceDescriptionXY
  !------------------------------------------------------------------------------------------
  function newSurfaceUniform(surfaceParameters, status) result(thisSurfaceDescription)
    real, dimension(:), intent(in   ) :: surfaceParameters
    type(ErrorMessage), intent(inout) :: status
    type(surfaceDescription)          :: thisSurfaceDescription
    !
    ! Define a horizontally uniform surface
    !
    real, dimension(2), parameter :: xPosition = (/ 0., huge(surfaceParameters) /), &
                                     yPosition = (/ 0., huge(surfaceParameters) /)
    real, dimension(numberOfParameters, 1, 1) :: surfaceParams
 
    if(size(surfaceParameters) /= numberOfParameters) then
      call setStateToFailure(status, "new_SurfaceDescription: Wrong number of parameters supplied for surface BRDF.") 
    else 
      surfaceParams(:, 1, 1) = surfaceParameters(:)
      thisSurfaceDescription = &
        newSurfaceDescriptionXY(surfaceParams, xPosition, yPosition, status)
    end if 
    
  end function newSurfaceUniform
  !------------------------------------------------------------------------------------------
  ! Compute surface reflectance at a given x-y location
  !------------------------------------------------------------------------------------------
  pure function computeSurfaceReflectance(thisSurfaceDescription, xPos, yPos, &
                                     incomingMu, outgoingMu, incomingPhi, outgoingPhi) result(surfaceReflectance)
    type(surfaceDescription), intent(in) :: thisSurfaceDescription
    real,                     intent(in) :: xPos, yPos, incomingMu, outgoingMu, incomingPhi, outgoingPhi
    real                                 :: surfaceReflectance

    real    :: x0, xMax, y0, yMax
    integer :: xIndex, yIndex
    
    !
    ! Find the square on the surface that's doing the reflecting 
    !
    x0 = thisSurfaceDescription%xPosition(1)
    y0 = thisSurfaceDescription%yPosition(1)
    xMax = thisSurfaceDescription%xPosition(size(thisSurfaceDescription%xPosition))
    yMax = thisSurfaceDescription%yPosition(size(thisSurfaceDescription%yPosition))
    xIndex = findIndex(makePeriodic(xPos, x0, xMax), thisSurfaceDescription%xPosition)
    yIndex = findIndex(makePeriodic(yPos, y0, yMax), thisSurfaceDescription%yPosition)
    
    !
    ! Compute the reflectance from the BRDF parameters
    !   Developers should make no assumptions about the sign of value 
    !   of incomingMu and outgoingMu and should be prepared to treat 
    !   values of phi between 360 and -360. 
    !
    surfaceReflectance = R(thisSurfaceDescription%BRDFParameters(:, xIndex, yIndex))
    
  end function computeSurfaceReflectance
  !------------------------------------------------------------------------------------------
  !   Compute reflectance given a set of BRDF parameters 
  !     This is where the work is done and it's the main section of the code
  !     that needs changing when implementing a new BRDF. 
  !------------------------------------------------------------------------------------------
  pure function R(surfaceParameters) 
    real, dimension(numberOfParameters), intent(in) :: surfaceParameters
    real                                            :: R
    !
    ! This example implements a Lambertian surface albedo
    !
    R = surfaceParameters(1)
    
  end function R 
  !------------------------------------------------------------------------------------------
  !   Readiness
  !------------------------------------------------------------------------------------------
  elemental function isReady_surfaceDescription(thisSurface) 
    type(surfaceDescription), intent(in) :: thisSurface
    logical                              :: isReady_surfaceDescription
    
    isReady_surfaceDescription = associated(thisSurface%xPosition) .and. &
                                 associated(thisSurface%yPosition) .and. &
                                 associated(thisSurface%BRDFParameters)

  end function isReady_surfaceDescription
  !------------------------------------------------------------------------------------------
  !   Copy
  !------------------------------------------------------------------------------------------
  function copy_surfaceDescription(original) result(copy)
    type(surfaceDescription), intent( in) :: original
    type(surfaceDescription)              :: copy
    
    integer :: numX, numY 
    
    if(isReady_surfaceDescription(original)) then 
      numX = size(original%xPosition); numY = size(original%yPosition)
      allocate(copy%xPosition(numX), &
               copy%yPosition(numY), &
               copy%BRDFParameters(numberOfParameters, numX - 1, numY - 1)) 
      copy%xPosition(:) = original%xPosition(:)
      copy%yPosition(:) = original%yPosition(:)  
      copy%BRDFParameters(:, :, :) = original%BRDFParameters(:, :, :)
    end if
  end function copy_surfaceDescription

  !------------------------------------------------------------------------------------------
  !   Finalization
  !------------------------------------------------------------------------------------------
  subroutine finalize_surfaceDescription(thisSurface)
    type(surfaceDescription), intent(out) :: thisSurface
    
    if(associated(thisSurface%xPosition)) deallocate(thisSurface%xPosition)
    if(associated(thisSurface%yPosition)) deallocate(thisSurface%yPosition)
    if(associated(thisSurface%BRDFParameters)) &
                                          deallocate(thisSurface%BRDFParameters)

  end subroutine finalize_surfaceDescription
  !------------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------------
  ! Utility procedures private to the module 
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
end  module surfaceProperties
