! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision$, $Date$
! $URL$
! --------------------------------------------
module opticalProperties
  ! Provide a representation of the three dimensional optical properties 
  !   of the atmosphere with each component (e.g. cloud, aerosol) 
  !   represented separately. 
  ! There are two objects: the domain in which the optical properties are specified, 
  !   and a set of one or more optical components. 
  use characterUtils
  use ErrorMessages
  use scatteringPhaseFunctions
  use netcdf
  implicit none
  private
  
  integer, parameter :: maxNameLength = 256
  
  !------------------------------------------------------------------------------------------
  ! Type (object) definitions
  !------------------------------------------------------------------------------------------
  type opticalComponent
    ! Each component has a name. The optical properties may be defined for a subset of the 
    !   z levels in the domain, and may be uniform or variable  horiontally. 
    !   the values of the optical properties (extinction, single scattering 
    !   albedo, phase function) of each component at each point. 
    ! The vectors defining the domain's cell boundaries should be one element
    !   larger than the corresponding array dimension. TIn other words, a field (extinction, say)
    !   defined everywhere is an an array of size nx, ny, nz, while the xPosition, xPosition, 
    !   and zPostion vectors are of length nx+1, ny+1, and nz+1.
    !
    private
    character (len = maxNameLength)      :: name       = ""
    integer                              :: zLevelBase = 0
    logical                              :: horizontallyUniform = .false.
    real,    dimension(:, :, :), pointer :: extinction             => null()
    real,    dimension(:, :, :), pointer :: singleScatteringAlbedo => null()
    integer, dimension(:, :, :), pointer :: phaseFunctionIndex     => null()
    type(phaseFunctionTable)             :: table
  end type opticalComponent

  type domain
    ! The domain contains the x, y, and z cell boundaries, flags to indicate whether the x and y 
    !   and/or z grids are evenly spaced, and pointers to the optical components. 
    !
    private
    real, pointer, dimension(:)          :: xPosition => null()
    real, pointer, dimension(:)          :: yPosition => null()
    real, pointer, dimension(:)          :: zPosition => null()
    logical                              :: xyRegularlySpaced = .false., zRegularlySpaced = .false. 
    type(opticalComponent), &
                   dimension(:), pointer :: components  => null()
  end type domain
  !------------------------------------------------------------------------------------------
  ! Overloading
  !------------------------------------------------------------------------------------------
  
  interface addOpticalComponent
    module procedure addOpticalComponent3D, addOpticalComponent1D
  end interface ! addOpticalComponent

  interface replaceOpticalComponent
    module procedure replaceOpticalComponent3D, replaceOpticalComponent1D
  end interface ! replaceOpticalComponent
  !------------------------------------------------------------------------------------------
  ! What is visible? 
  !------------------------------------------------------------------------------------------

  ! The types...
  public :: domain 
  ! ... and these procedures
  public :: new_Domain, getInfo_Domain, write_Domain, read_Domain, finalize_Domain, &
            addOpticalComponent, deleteOpticalComponent, replaceOpticalComponent,   &
            getOpticalPropertiesByComponent !,  getAverageOpticalProperties
            

contains
  !------------------------------------------------------------------------------------------
  ! Initialization: Routine to create new domains 
  !------------------------------------------------------------------------------------------
  function new_Domain(xPosition, yPosition, zPosition, status)
    real,    dimension(:), intent(in   ) :: xPosition, yPosition, zPosition
    type(ErrorMessage),    intent(inout) :: status
    type(domain)                         :: new_Domain
    
    ! Local variables
    integer :: numX, numY, numZ
     
    ! -------------------------
    ! Checks: always increasing, within limits; 
    !   other checks performed by addOpticalComponent
    numX = size(xPosition); numY = size(yPosition); numZ = size(zPosition)
    if(any(xPosition(2:) - xPosition(:numX-1) <= 0.) .or. &
       any(yPosition(2:) - yPosition(:numY-1) <= 0.) .or. &
       any(zPosition(2:) - zPosition(:numZ-1) <= 0.))     &
      call setStateToFailure(status, "new_Domain: Positions must be increasing, unique.")
     
    ! -------------------------
     if(.not. stateIsFailure(status)) then
       allocate(new_Domain%xPosition(numX), new_Domain%yPosition(numY), new_Domain%zPosition(numZ))
       new_Domain%xPosition(:) = xPosition(:)
       new_Domain%yPosition(:) = yPosition(:)
       new_Domain%zPosition(:) = zPosition(:)
       
       ! Are the grids regularly spaced? Compare the distance between each pair 
       !   of array elements to the distance between the first two. 
       ! The default value is false. 
       !
       if(all(abs( (xPosition(2:) - xPosition(:numX-1)) -                                &
                   (xPosition(2)  - xPosition(1)) ) <= 2 * spacing(xPosition(2:))) .and. &
          all(abs( (yPosition(2:) - yPosition(:numY-1)) -                                &
                   (yPosition(2)  - yPosition(1)) ) <= 2 * spacing(yPosition(2:))))      &
         new_Domain%xyRegularlySpaced = .true.
       if(all(abs( (zPosition(2:) - zPosition(:numZ-1)) -                           &
                   (zPosition(2)  - zPosition(1)) ) <= 2 * spacing(zPosition(2:)))) &
         new_Domain%zRegularlySpaced = .true.
       call setStateToSuccess(status)
     end if
  end function new_Domain
  !------------------------------------------------------------------------------------------
  subroutine addOpticalComponent3D(thisDomain, componentName,          &
                                   extinction, singleScatteringAlbedo, &
                                   phaseFunctionIndex, phaseFunctions, &
                                   zLevelBase, status)
    type(domain),                intent(inout) :: thisDomain
    character (len = *),         intent(in   ) :: componentName
    real,    dimension(:, :, :), intent(in   ) :: extinction, singleScatteringAlbedo
    integer, dimension(:, :, :), intent(in   ) :: phaseFunctionIndex
    type(phaseFunctionTable),    intent(in   ) :: phaseFunctions
    integer, optional,           intent(in   ) :: zLevelBase
    type(ErrorMessage),          intent(inout) :: status
    ! 
    ! Add a new optical component to the domain. We check the optical 
    !   properties to make sure they make sense, and also that they can be 
    !   associated with the domain. 
    ! Implementation note: the domain type has an array containing the optical components; 
    !   each new component means allocating a slightly bigger array and copying the 
    !   optical component objects.  

    ! Local variables
    integer                               :: nComponents
    type(opticalComponent), &
                    dimension(:), pointer :: tempComponents
    integer                               :: baseLevel 

    ! -----------------
    ! Checks are done in validateOpticalComponent; zBaseLevel is assumed to be 1
    !   if it isn't supplied, if this doesn't make sense it'll get caught in validateOpticalComponent
    if(present(zLevelBase)) then
      baseLevel = zLevelBase
    else
      baseLevel = 1
    end if
    call validateOpticalComponent(thisDomain, componentName,          &
                                  extinction, singleScatteringAlbedo, &
                                  phaseFunctionIndex, phaseFunctions, &
                                  baseLevel, status)
    ! -----------------
    if(.not. stateIsFailure(status)) then 
      ! How many components will we have when we're done adding this one? 
      !   Allocate new memory, then copy the data from the pre-existing array. 
      ! 
      if(containsComponents(thisDomain)) then
        nComponents = size(thisDomain%components)
        allocate(tempComponents(nComponents + 1))
        tempComponents(:nComponents) = thisDomain%components(:)
      else
        nComponents = 0
        allocate(tempComponents(1))
      end if 
      
      tempComponents(nComponents + 1) =                                          &
          newOpticalComponent(componentName, extinction, singleScatteringAlbedo, &
                              phaseFunctionIndex, baseLevel, phaseFunctions)
      
      ! Free the memory associated with the old array
      if(associated(thisDomain%components)) deallocate(thisDomain%components)       
      thisDomain%components => tempComponents !   and point to the new array
      
      call setStateToSuccess(status)
    else
      call setStateToFailure(status, "addOpticalComponent: optical properties aren't valid.")
    end if
  end subroutine addOpticalComponent3D
  !------------------------------------------------------------------------------------------
  subroutine addOpticalComponent1D(thisDomain, componentName,          &
                                   extinction, singleScatteringAlbedo, &
                                   phaseFunctionIndex, phaseFunctions, &
                                   zLevelBase, status)
    type(domain),                intent(inout) :: thisDomain
    character (len = *),         intent(in   ) :: componentName
    real,    dimension(:),       intent(in   ) :: extinction, singleScatteringAlbedo
    integer, dimension(:),       intent(in   ) :: phaseFunctionIndex
    type(phaseFunctionTable),    intent(in   ) :: phaseFunctions
    integer,           optional, intent(in   ) :: zLevelBase
    type(ErrorMessage),          intent(inout) :: status
    !
    ! Add a new one-dimensional component to the domain. 
    !   We do this by pretending it's a 3D component with dimension size 1 in the 
    !   x and y directions.

    ! -----------------
    ! Local variables
    integer :: baseLevel

    if(present(zLevelBase)) then
      baseLevel = zLevelBase
    else
      baseLevel = 1
    end if

    call addOpticalComponent3D(thisDomain, componentName,                                                 &
                               reshape(extinction,             (/ 1, 1, size(extinction) /)),             &
                               reshape(singleScatteringAlbedo, (/ 1, 1, size(singleScatteringAlbedo) /)), &
                               reshape(phaseFunctionIndex,     (/ 1, 1, size(phaseFunctionIndex) /)),     &
                               phaseFunctions, zLevelBase = baseLevel, status = status)

  end subroutine addOpticalComponent1D
  !------------------------------------------------------------------------------------------
  subroutine replaceOpticalComponent3D(thisDomain, componentNumber, componentName,                &
                                       extinction, singleScatteringAlbedo,                        &
                                       phaseFunctionIndex, phaseFunctions, &
                                       zLevelBase, status)
    type(domain),                intent(inout) :: thisDomain
    integer,                     intent(in   ) :: componentNumber
    character (len = *),         intent(in   ) :: componentName
    real,    dimension(:, :, :), intent(in   ) :: extinction, singleScatteringAlbedo
    integer, dimension(:, :, :), intent(in   ) :: phaseFunctionIndex
    type(phaseFunctionTable),    intent(in   ) :: phaseFunctions
    integer, optional,           intent(in   ) :: zLevelBase
    type(ErrorMessage),          intent(inout) :: status
    ! 
    ! Replace one of the domains optical components. This is like adding except we don't have 
    !   to allocate memory.  
    integer                               :: baseLevel 

    ! -----------------
    if(.not. containsComponents(thisDomain)) then
      call setStateToFailure(status, "replaceOpticalComponent: no components to replace" )
    else if(componentNumber < 1 .or. componentNumber > size(thisDomain%components)) then
      call setStateToFailure(status, "replaceOpticalComponent: no components to replace" )
    end if
    ! Checks are done in validateOpticalComponent; zBaseLevel is assumed to be 1
    !   if it isn't supplied, if this doesn't make sense it'll get caught in validateOpticalComponent
    if(present(zLevelBase)) then
      baseLevel = zLevelBase
    else
      baseLevel = 1
    end if
    call validateOpticalComponent(thisDomain, componentName,          &
                                  extinction, singleScatteringAlbedo, &
                                  phaseFunctionIndex, phaseFunctions, &
                                  baseLevel, status)

    ! -----------------
    if(.not. stateIsFailure(status)) then 
      call finalizeComponent(thisDomain%components(componentNumber))
      thisDomain%components(componentNumber) =                                   &
          newOpticalComponent(componentName, extinction, singleScatteringAlbedo, &
                             phaseFunctionIndex, baseLevel, phaseFunctions)
      call setStateToSuccess(status)
    end if
  end subroutine replaceOpticalComponent3D
  !------------------------------------------------------------------------------------------
  subroutine replaceOpticalComponent1D(thisDomain, componentNumber, componentName, &
                                       extinction, singleScatteringAlbedo,         &
                                       phaseFunctionIndex, phaseFunctions,         &
                                       zLevelBase, status)
    type(domain),                intent(inout) :: thisDomain
    integer,                     intent(in   ) :: componentNumber
    character (len = *),         intent(in   ) :: componentName
    real,    dimension(:),       intent(in   ) :: extinction, singleScatteringAlbedo
    integer, dimension(:),       intent(in   ) :: phaseFunctionIndex
    type(phaseFunctionTable),    intent(in   ) :: phaseFunctions
    integer,           optional, intent(in   ) :: zLevelBase
    type(ErrorMessage),          intent(inout) :: status
    !
    ! Replace a new component in the domain with something one-dimensional. It's like adding but with no 
    !   memory headaches. 

    ! -----------------
    ! Local variables
    integer :: baseLevel

    if(present(zLevelBase)) then
      baseLevel = zLevelBase
    else
      baseLevel = 1
    end if

    call replaceOpticalComponent3D(thisDomain, componentNumber, componentName,                            &
                               reshape(extinction,             (/ 1, 1, size(extinction) /)),             &
                               reshape(singleScatteringAlbedo, (/ 1, 1, size(singleScatteringAlbedo) /)), &
                               reshape(phaseFunctionIndex,     (/ 1, 1, size(phaseFunctionIndex) /)),     &
                               phaseFunctions, zLevelBase = baseLevel, status = status)

  end subroutine replaceOpticalComponent1D
  !------------------------------------------------------------------------------------------
  subroutine deleteOpticalComponent(thisDomain, componentNumber, status)
    ! Delete a component from the domain
    type(domain),                intent(inout) :: thisDomain
    integer,                     intent(in   ) :: componentNumber
    type(ErrorMessage),          intent(inout) :: status
    !
    ! Delete a component from the domain. By analogy to the process when 
    !   we add a new component, we create a shorter array in the domain 
    !   object, then copy the remaining optical components to the new array 
    !   and loose the old one. 
    ! 

    ! Local variables
    type(opticalComponent), &
      dimension(size(thisDomain%components) - 1), target :: tempComponents
    
    ! --------------------
    if(.not. isValid(thisDomain)) &
      call setStateToFailure(status, "deleteOpticalComponent: domain hasn't been initialized.")
    if(containsComponents(thisDomain)) then
      if(componentNumber < 1 .or. componentNumber > size(thisDomain%components)) &
        call setStateToFailure(status, "deleteOpticalComponent: non-existent component.")
    else
      call setStateToFailure(status, "deleteOpticalComponent: no components to delete.")
    end if 
    
    ! --------------------
    if(.not. stateIsFailure(status)) then
      if(size(thisDomain%components) == 1) then 
        ! There is only one component, so we deallocate (and nullify) the 
        !   array that holds the components. 
        !
        call finalizeComponent(thisDomain%components(1))
        deallocate(thisDomain%components)
      else 
        tempComponents(:) = (/ thisDomain%components(:(componentNumber - 1)), &
                               thisDomain%components((componentNumber + 1):) /)
        ! This deallocates the pointer to the array that holds the optical components
        !   but not the underlying data. 
        call finalizeComponent(thisDomain%components(componentNumber))
        deallocate(thisDomain%components); allocate(thisDomain%components(size(tempComponents)))
        thisDomain%components(:) = tempComponents(:)
      end if 
      
      call setStateToSuccess(status)
    end if
  end subroutine deleteOpticalComponent
  !------------------------------------------------------------------------------------------
  !  Getting information back from the object
  !------------------------------------------------------------------------------------------
  subroutine getInfo_Domain(thisDomain, numX, numY, numZ,    &
                            xPosition, yPosition, zPosition, &
                            numberOfComponents, componentNames, status) 
    type(domain),                    intent(in   ) :: thisDomain
    integer,               optional, intent(  out) :: numX, numY, numZ
    real,    dimension(:), optional, intent(  out) :: xPosition, yPosition, zPosition
    integer,               optional, intent(  out) :: numberOfComponents
    character(len = *), &
             dimension(:), optional, intent(  out) :: componentNames
    type(ErrorMessage),              intent(inout) :: status
    
    ! What can you get back from the domain? The number of cells in the arrays, 
    !   and the locations of the x, y, z, boundaries (which is one bigger). 
    ! Also the number and names of the various optical components. 
    
    ! --------------------
    ! Checks: domain is valid, position arrays are big enough (the right size?)
    if(.not. isValid(thisDomain)) then
      call setStateToFailure(status, "getInfo_Domain: domain hasn't been initialized.")
    else
      ! Number of positions in each dimension for the property arrays. It's one shorter 
      !   than the position arrays, which describe the boundaries. 
      if(present(numX)) numX = size(thisDomain%xPosition) - 1
      if(present(numY)) numY = size(thisDomain%yPosition) - 1 
      if(present(numZ)) numZ = size(thisDomain%zPosition) - 1
      
      ! Location of boundaries in each dimension
      if(present(xPosition)) then 
        if(size(xPosition) /= size(thisDomain%xPosition)) then
          call setStateToFailure(status, "getInfo_Domain: vector for x positions is wrong length.")
        else
          xPosition(:) = thisDomain%xPosition(:)
        end if
      end if
      if(present(yPosition)) then 
        if(size(yPosition) /= size(thisDomain%yPosition)) then
          call setStateToFailure(status, "getInfo_Domain: vector for y positions is wrong length.")
        else
          yPosition(:) = thisDomain%yPosition(:)
        end if
      end if
      if(present(zPosition)) then 
        if(size(zPosition) /= size(thisDomain%zPosition)) then
          call setStateToFailure(status, "getInfo_Domain: vector for z positions is wrong length.")
        else
          zPosition(:) = thisDomain%zPosition(:)
        end if
      end if
      
      if(containsComponents(thisDomain)) then 
        if(present(numberOfComponents)) numberOfComponents = size(thisDomain%components)
        if(present(componentNames)) then 
          if(size(componentNames) < size(thisDomain%components)) then
            call setStateToFailure(status, "getInfo_Domain: component names array is wrong length")
          else 
            componentNames(:size(thisDomain%components)) = thisDomain%components(:)%name
          end if
        end if
      else 
        if(present(numberOfComponents)) numberOfComponents = 0
        if(present(componentNames))     componentNames(:) = ""
      end if
      if(.not. stateIsFailure(status)) call setStateToSuccess(status)
    end if 
  end subroutine getInfo_Domain
  !------------------------------------------------------------------------------------------
  ! Getting the optical properties of the domain
  !------------------------------------------------------------------------------------------
  subroutine getOpticalPropertiesByComponent(thisDomain,                        &
                 totalExtinction, cumulativeExtinction, singleScatteringAlbedo, &
                 phaseFunctionIndex, phaseFunctions, status)
    type(domain),                    intent(in   ) :: thisDomain
    real,    dimension(:, :, :),     intent(  out) :: totalExtinction
    real,    dimension(:, :, :, :),  intent(  out) :: cumulativeExtinction, singleScatteringAlbedo
    integer, dimension(:, :, :, :),  intent(  out) :: phaseFunctionIndex
    type(phaseFunctionTable), &
             dimension(:), optional, intent(  out) :: phaseFunctions
    type(ErrorMessage),              intent(inout) :: status
    ! 
    ! Return the optical properties of the domain, component by component. 
    !   The properties returned are defined at all points in the domain. 
    !   If the component is horizontally homogeneous the properties are expanded into 
    !   3D; if it exists only on some subset of vertical layers the properties at 
    !   other layers are set to 0. 
    ! Arrays are dimensioned x, y, z, component. 
    ! Cumulative extinction for the i'th component is the proportion of the total extinction due 
    !   to components 1 to i. It's intended to be used to pick the component doing the extinction
    !   at each scattering/absorption event in a Monte Carlo scheme.  
    ! Phase function tables are returned in arrays. 
    !   
    
    ! Local variable
    integer :: numX, numY, numZ, numComponents
    integer :: i, minZ, maxZ
    
    !------------
    ! Sanity checks
    !   Domain is in use and contains components
    !   x, y, z, dimensions, number of components correct for all arrays
    !   number of entries in each table is correct
    !
    if(.not. isValid(thisDomain)) then
      call setStateToFailure(status, "getOpticalPropertiesByComponent: domain is not initialized.")
    else if (.not. containsComponents(thisDomain)) then
      call setStateToFailure(status, "getOpticalPropertiesByComponent: domain contains no optical components.")
    else 
      !
      ! Checks for x, y, z, sizes
      numX = size(thisDomain%xPosition) - 1; numY = size(thisDomain%yPosition) - 1
      numZ = size(thisDomain%zPosition) - 1; numComponents = size(thisDomain%components)
      if(any( (/ size(totalExtinction,        1), size(cumulativeExtinction, 1),          &
                 size(singleScatteringAlbedo, 1), size(phaseFunctionIndex,   1) /) /= numX)) &
        call setStateToFailure(status, "getOpticalPropertiesByComponent: extent of one or more arrays incorrect in x dimension.")
      if(any( (/ size(totalExtinction,        2), size(cumulativeExtinction, 2),          &
                 size(singleScatteringAlbedo, 2), size(phaseFunctionIndex,   2) /) /= numY)) &
        call setStateToFailure(status, "getOpticalPropertiesByComponent: extent of one or more arrays incorrect in y dimension.")
      if(any( (/ size(totalExtinction,        3), size(cumulativeExtinction, 3),  &
                 size(singleScatteringAlbedo, 3), size(phaseFunctionIndex,   3) /) /= numZ)) &
        call setStateToFailure(status, "getOpticalPropertiesByComponent: extent of one or more arrays incorrect in z dimension.")
      !
      ! Checks for number of components
      if(any( (/ size(cumulativeExtinction, 4),  size(singleScatteringAlbedo, 4), &
                 size(phaseFunctionIndex,   4) /)  /= numComponents))             &
        call setStateToFailure(status, "getOpticalPropertiesByComponent: number of components in one or more arrays incorrect.")
      if(present(phaseFunctions)) then
        if(size(phaseFunctions) /= numComponents) &
          call setStateToFailure(status, "getOpticalPropertiesByComponent: number of components in phaseFunctions array.")
      end  if 
    end if 

    ! -----------------------
    !   Everything looks ok for now. 
    !
    if(.not. stateIsFailure(status)) then 
      totalExtinction       (:, :, :)    = 0.
      cumulativeExtinction  (:, :, :, :) = 0.
      singleScatteringAlbedo(:, :, :, :) = 0.
      phaseFunctionIndex    (:, :, :, :) = 0
      
      do i = 1, numComponents
        minZ = thisDomain%components(i)%zLevelBase
        maxZ = minZ + size(thisDomain%components(i)%extinction, 3) - 1
        
        !
        ! In this loop cumulativeExtinction is the extinction by component 
        !
        if(.not. thisDomain%components(i)%horizontallyUniform) then 
          cumulativeExtinction  (:, :, minZ:maxZ, i) = thisDomain%components(i)%extinction(:, :, :)
          singleScatteringAlbedo(:, :, minZ:maxZ, i) = thisDomain%components(i)%singleScatteringAlbedo(:, :, :)
          phaseFunctionIndex    (:, :, minZ:maxZ, i) = thisDomain%components(i)%phaseFunctionIndex(:, :, :)
        else
          cumulativeExtinction(:, :, minZ:maxZ, i) =                                                      &
                                             spread(spread(thisDomain%components(i)%extinction(1, 1, :),  &
                                                           1, nCopies = numX), 2, nCopies = numY)
          singleScatteringAlbedo(:, :, minZ:maxZ, i) =                                                               &
                                             spread(spread(thisDomain%components(i)%singleScatteringAlbedo(1, 1, :), &
                                                           1, nCopies = numX), 2, nCopies = numY)
          phaseFunctionIndex    (:, :, minZ:maxZ, i) =                                                           &
                                             spread(spread(thisDomain%components(i)%phaseFunctionIndex(1, 1, :), &
                                                           1, nCopies = numX), 2, nCopies = numY)
        end if
        ! We don't finalize the array element phaseFunctions(i) in case something else is pointing to the 
        !   underlying memory. Users should finalize before they pass the array in. 
        ! 
        if(present(phaseFunctions)) phaseFunctions(i) = copy_PhaseFunctionTable(thisDomain%components(i)%table)
        
      end do 
      !
      ! Compute total extinction and convert cumulative extinction to fractional cumulative extinction. 
      !   Every element of cumulativeExtinction(:, :, :, numComponents) should equal 1 by definition
      !
      totalExtinction(:, :, minZ:maxZ) = sum(cumulativeExtinction  (:, :, minZ:maxZ, :), dim = 4)
      where(spread(totalExtinction, 4, nCopies = numComponents) > tiny(totalExtinction)) &
        cumulativeExtinction(:, :, :, :) = cumulativeExtinction(:, :, :, :)  / spread(totalExtinction, 4, nCopies = numComponents)
    end if
  end subroutine getOpticalPropertiesByComponent
  !------------------------------------------------------------------------------------------
!   subroutine getAverageOpticalProperties(thisDomain, extinction, singleScatteringAlbedo, &
!                 phaseFunctionIndex, phaseFunctions, status)
!     type(domain),                 intent( in) :: thisDomain
!     real,    dimension(:, :, :),  intent(out) :: extinction, singleScatteringAlbedo
!     integer, dimension(:, :, :),  intent(out) :: phaseFunctionIndex
!     type(phaseFunctionTable),     intent(out) :: phaseFunctions
!     type(ErrorMessage),           intent(inout) :: status
!     
!     
!   end subroutine getAverageOpticalProperties
  !------------------------------------------------------------------------------------------
  ! Storage and retrieval
  !------------------------------------------------------------------------------------------
  subroutine write_Domain(thisDomain, fileName, status)
    type(domain),       intent(in   ) :: thisDomain
    character(len = *), intent(in   ) :: fileName
    type(ErrorMessage), intent(inout) :: status
    
    ! Local variables
    integer                        :: i, j
    logical                        :: fillsDomainInVertical
    
    ! Netcdf-related local variables
    integer, dimension(16) :: ncStatus
    integer                :: ncFileId, xEdgeDimID, yEdgeDimId, zEdgeDimId,  &
                                        xGridDimId, yGridDimId, zGridDimId,  &
                                        extinctionVarId, ssaVarId, indexVarId, ncVarId
    
    ! Checks: domain is valid, contains components(?)
    ncStatus(:) = nf90_NoErr
    if(.not. isValid(thisDomain)) then
      call setStateToFailure(status, "write_Domain: domain hasn't been initialized.") 
    else
      ncStatus( 1) = nf90_create(trim(fileName), nf90_Clobber, ncFileId)
      !
      ! Domain dimensions 
      !
      ncStatus( 2) = nf90_def_dim(ncFileId, "x-Edges", size(thisDomain%xPosition),    xEdgeDimId) 
      ncStatus( 3) = nf90_def_dim(ncFileId, "y-Edges", size(thisDomain%yPosition),    yEdgeDimId) 
      ncStatus( 4) = nf90_def_dim(ncFileId, "z-Edges", size(thisDomain%zPosition),    zEdgeDimId) 
      ncStatus( 5) = nf90_def_dim(ncFileId, "x-Grid",  size(thisDomain%xPosition) - 1, xGridDimId) 
      ncStatus( 6) = nf90_def_dim(ncFileId, "y-Grid",  size(thisDomain%yPosition) - 1, yGridDimId) 
      ncStatus( 7) = nf90_def_dim(ncFileId, "z-Grid",  size(thisDomain%zPosition) - 1, zGridDimId)
      !
      ! Domain variables
      ! 
      ncStatus( 8) = nf90_def_var(ncFileId, "x-Edges", nf90_float, xEdgeDimId, ncVarId) 
      ncStatus( 9) = nf90_def_var(ncFileId, "y-Edges", nf90_float, yEdgeDimId, ncVarId) 
      ncStatus(10) = nf90_def_var(ncFileId, "z-Edges", nf90_float, zEdgeDimId, ncVarId) 
      !
      ! Domain attributes
      !
      ncStatus(11) = nf90_put_att(ncFileID, nf90_Global, "xyRegularlySpaced", asInt(thisDomain%xyRegularlySpaced)) 
      ncStatus(12) = nf90_put_att(ncFileID, nf90_Global,  "zRegularlySpaced", asInt(thisDomain%zRegularlySpaced))
      if(any(ncStatus(:) /= nf90_NoErr)) &
        call setStateToFailure(status, "write_Domain: error writing domain information") 
        
      if(.not. stateIsFailure(status) .and. containsComponents(thisDomain)) then
        ncStatus( 1) = nf90_put_att(ncFileID, nf90_Global, "numberOfComponents", size(thisDomain%components))
        do i = 1, size(thisDomain%components)
          ! 
          ! For each component: Attributes 
          !
          ncStatus( 1) = nf90_put_att(ncFileId, nf90_global, trim(makePrefix(i)) // "Name",       &
                                      trim(thisDomain%components(i)%name))
          ncStatus( 2) = nf90_put_att(ncFileId, nf90_global, trim(makePrefix(i)) // "zLevelBase", &
                                      thisDomain%components(i)%zLevelBase)
          !
          ! Dimension definition -  the length of the z dimension can be different than the domains
          ! 
          fillsDomainInVertical = thisDomain%components(i)%zLevelBase == 1 .and. &
                                  size(thisDomain%components(i)%extinction, 3) == (size(thisDomain%zPosition) -1)
          if(fillsDomainInVertical) then
            ncStatus( 6) = nf90_inq_dimid(ncFileId, "z-Grid", zGridDimId)
          else
            ncStatus( 6) = nf90_def_dim(ncFileId, trim(makePrefix(i)) // "z-Grid", &
                                        size(thisDomain%components(i)%extinction, 3), zGridDimId)
          end if 
          !
          ! Variables
          !   Doesn't seem like there will ever be more than 2^15 possible phase functions, 
          !   so we can use a two-byte integer to store the phaseFunctionIndex. 
          ! 
          if(thisDomain%components(i)%horizontallyUniform) then
            ! Variables are 1D
            !
            ncStatus( 7) = nf90_def_var(ncFileId, trim(makePrefix(i)) // "Extinction",             nf90_float, &
                                                       zGridDimId, extinctionVarId)
            ncStatus( 8) = nf90_def_var(ncFileId, trim(makePrefix(i)) // "SingleScatteringAlbedo", nf90_float, &
                                                       zGridDimId, ssaVarId)
            ncStatus( 9) = nf90_def_var(ncFileId, trim(makePrefix(i)) // "PhaseFunctionIndex",     nf90_short, &
                                                       zGridDimId, indexVarId)
          else
            ! Variables are 3D
            !
            ncStatus( 7) = nf90_def_var(ncFileId, trim(makePrefix(i)) // "Extinction",             nf90_float, &
                                        (/ xGridDimId, yGridDimId, zGridDimId /) , extinctionVarId)
            ncStatus( 8) = nf90_def_var(ncFileId, trim(makePrefix(i)) // "SingleScatteringAlbedo", nf90_float, &
                                        (/ xGridDimId, yGridDimId, zGridDimId /) , ssaVarId)
            ! Doesn't seem like there will ever be more than 2^15 possible phase functions...
            ncStatus( 9) = nf90_def_var(ncFileId, trim(makePrefix(i)) // "PhaseFunctionIndex",     nf90_short, &
                                        (/ xGridDimId, yGridDimId, zGridDimId /) , indexVarId)
          end if
          if(any(ncStatus(:) /= nf90_NoErr)) &
            call setStateToFailure(status,   &
                                   "write_Domain: Error creating definitions for component" // trim(IntToChar(i)))
        end do 
      end if
      !
      ! All the attributes have been written, and the dimensions and variables defined. Now we'll take the
      !   file out of define mode and write the domain data, then the data for each component.
      !
      ncStatus( 1) = nf90_EndDef(ncFileId)
      ncStatus( 2) = nf90_inq_varid(ncFileId, "x-Edges", ncVarID)
      ncStatus( 3) = nf90_put_var(ncFileId, ncVarId, thisDomain%xPosition)
      ncStatus( 4) = nf90_inq_varid(ncFileId, "y-Edges", ncVarID)
      ncStatus( 5) = nf90_put_var(ncFileId, ncVarId, thisDomain%yPosition)
      ncStatus( 6) = nf90_inq_varid(ncFileId, "z-Edges", ncVarID)
      ncStatus( 7) = nf90_put_var(ncFileId, ncVarId, thisDomain%zPosition)
      if(any(ncStatus(:) /= nf90_NoErr)) &
        call setStateToFailure(status, "write_Domain: error writing domain data") 
        
      if(.not. stateIsFailure(status) .and. containsComponents(thisDomain)) then
        do i = 1, size(thisDomain%components)
          ncstatus( 1) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "Extinction",  extinctionVarId)
          ncstatus( 2) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "SingleScatteringAlbedo",  ssaVarId)
          ncstatus( 3) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "PhaseFunctionIndex",  indexVarId)
          if(thisDomain%components(i)%horizontallyUniform) then
            ncStatus( 4) = nf90_put_var(ncFileId, extinctionVarId, &
                                        thisDomain%components(i)%extinction(1, 1, :))
            ncStatus( 5) = nf90_put_var(ncFileId, ssaVarId,        &
                                        thisDomain%components(i)%singleScatteringAlbedo(1, 1, :))
            ncStatus( 6) = nf90_put_var(ncFileId, indexVarId,      &
                                        thisDomain%components(i)%phaseFunctionIndex(1, 1, :))
          else
            ncStatus( 4) = nf90_put_var(ncFileId, extinctionVarId, &
                                        thisDomain%components(i)%extinction(:, :, :))
            ncStatus( 5) = nf90_put_var(ncFileId, ssaVarId,        &
                                        thisDomain%components(i)%singleScatteringAlbedo(:, :, :))
            ncStatus( 6) = nf90_put_var(ncFileId, indexVarId,      &
                                        thisDomain%components(i)%phaseFunctionIndex(:, :, :))
          end if
          call add_PhaseFunctionTable(thisDomain%components(i)%table,            &
                                      fileId = ncFileId,                         &
                                      prefix = "Component" // trim(IntToChar(i)) // "_", &
                                      status = status)
          do j = 1, 6
            if(ncStatus(j) /= nf90_NoErr) &
              call setStateToFailure(status, "write_Domain: " // trim(nf90_StrError(ncStatus(j))))
          end do
        end do
      end if
      ncStatus(1) = nf90_close(ncFileId)
      if(.not. stateIsFailure(status)) call setStateToSuccess(status)
    end if 
    
  end subroutine write_Domain
  !------------------------------------------------------------------------------------------
  subroutine read_Domain(fileName, thisDomain, status)
    character(len = *), intent(in   ) :: fileName
    type(domain),       intent(  out) :: thisDomain
    type(ErrorMessage), intent(inout) :: status
    
    ! Local variables
    integer(kind = selected_int_kind(2)) &
                       :: oneByte
    integer            :: i
    integer            :: nXEdges, nYEdges, nZEdges, nComponents, nZGrid
    real, dimension(:), &
           allocatable :: xEdges, yEdges, zEdges
    
    ! Variable for each component
    real,    dimension(:, :, :), allocatable :: extinction, singleScatteringAlbedo
    integer, dimension(:, :, :), allocatable :: phaseFunctionIndex
    type(phaseFunctionTable)                 :: table
    character(len = maxNameLength)           :: name
    integer                                  :: zLevelBase
    logical                                  :: fillsVerticalDomain, horizontallyUniform
    
    ! Netcdf-related local variables
    integer, dimension(16) :: ncStatus
    integer                :: ncFileId, ncDimId, ncVarId, zGridDimId, nDims
    integer, dimension(3)  :: dimIds
   
    ! There is nothing to check a priori
    ncStatus(:) = nf90_NoErr
    if(nf90_open(trim(fileName), nf90_NoWrite, ncFileID) /= nf90_NoErr) then
      call setStateToFailure(status, "read_Domain: Can't open file " // trim(fileName)) 
    end if
    
    if(.not. stateIsFailure(status)) then 
      ncStatus( 1) = nf90_inq_dimid(ncFileId, "x-Edges", ncDimId) 
      ncStatus( 2) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nXEdges)
      ncStatus( 3) = nf90_inq_dimid(ncFileId, "y-Edges", ncDimId) 
      ncStatus( 4) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nYEdges)
      ncStatus( 5) = nf90_inq_dimid(ncFileId, "z-Edges", ncDimId) 
      ncStatus( 6) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nZEdges)
      allocate(xEdges(nXEdges), yEdges(nYEdges), zEdges(nZEdges))
      ncStatus( 7) = nf90_inq_varid(ncFileId, "x-Edges", ncVarId)
      ncStatus( 8) = nf90_get_var(ncFileId, ncVarId, xEdges)
      ncStatus( 9) = nf90_inq_varid(ncFileId, "y-Edges", ncVarId)
      ncStatus(10) = nf90_get_var(ncFileId, ncVarId, yEdges)
      ncStatus(11) = nf90_inq_varid(ncFileId, "z-Edges", ncVarId)
      ncStatus(12) = nf90_get_var(ncFileId, ncVarId, zEdges)
      ncStatus(13) = nf90_inq_dimId(ncFileId, "z-Grid", zGridDimId) 
      if(any(ncStatus(:) /= nf90_NoErr)) &
        call setStateToFailure(status, "read_Domain: " // trim(fileName) // &
                               " doesn't look an optical properties file.") 
    end if 
    
    if(.not. stateIsFailure(status)) then 
      !
      ! Create a new domain using the initialization procedure. 
      !   The domain may have been regularly spaced when it was written to the file, but 
      !   this isn't guaranteed any more. 
      !
      call finalize_Domain(thisDomain)
      thisDomain = new_Domain(xEdges, yEdges, zEdges, status)
      ncStatus( 1) = nf90_get_att(ncFileID, nf90_Global, "xyRegularlySpaced", oneByte) 
      if(asLogical(oneByte) .neqv. thisDomain%xyRegularlySpaced) &
        call setStateToWarning(status, "read_Domain: file and new domain don't agree on regularity of x-y spacing.") 
      ncStatus( 2) = nf90_get_att(ncFileID, nf90_Global,  "zRegularlySpaced", oneByte)
      if(asLogical(oneByte) .neqv. thisDomain%zRegularlySpaced) &
        call setStateToWarning(status, "read_Domain: file and new domain don't agree on regularity of z spacing.") 
      deallocate(xEdges, yEdges, zEdges)
    end if 

    if(.not. stateIsFailure(status)) then 
      !
      ! Read in the data for each component, and then add the component to the domain using the 
      !   same function we'd use from outside the module. 
      !
      ncStatus( 1) = nf90_get_att(ncFileID, nf90_Global, "numberOfComponents", nComponents)
      do i = 1, nComponents
        ncStatus( 1) = nf90_get_att(ncFileId, nf90_global, trim(makePrefix(i)) // "Name", name)
        ncStatus( 2) = nf90_get_att(ncFileId, nf90_global, trim(makePrefix(i)) // "zLevelBase", zLevelBase)
        !
        ! Is the component horizontally homogeneous? Does it fill the entire domain vertically? 
        ! 
        ncStatus( 4) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "Extinction", ncVarId)
        ncStatus( 5) = nf90_Inquire_Variable(ncFileId, ncVarId, ndims = nDims, dimids = dimIds)
        horizontallyUniform = (ndims == 1)
        fillsVerticalDomain = (dimIds(ndims) == zGridDimId)
        if(fillsVerticalDomain) then
          nZGrid = nZEdges - 1
        else
          ncStatus( 6) = nf90_inq_dimId(ncFileId, trim(makePrefix(i)) // "z-Grid", ncDimId)
          ncStatus( 7) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nZGrid)
        end if
        !
        ! Read in the scalar 3D or 1D fields
        !
        if(horizontallyUniform) then
          allocate(extinction(1, 1, nZGrid), singleScatteringAlbedo(1, 1, nZGrid), phaseFunctionIndex(1, 1, nZGrid)) 
          ncStatus( 8) = nf90_get_var(ncFileId, ncVarId, extinction(1, 1, :))
          ncStatus( 9) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "SingleScatteringAlbedo", ncVarId)
          ncStatus(10) = nf90_get_var(ncFileId, ncVarId, singleScatteringAlbedo(1, 1, :))
          ncStatus(11) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "PhaseFunctionIndex", ncVarId)
          ncStatus(12) = nf90_get_var(ncFileId, ncVarId, phaseFunctionIndex(1, 1, :))
        else 
          allocate(            extinction(nXEdges - 1, nYEdges - 1, nZGrid), &
                   singleScatteringAlbedo(nXEdges - 1, nYEdges - 1, nZGrid), &
                       phaseFunctionIndex(nXEdges - 1, nYEdges - 1, nZGrid))
          ncStatus( 8) = nf90_get_var(ncFileId, ncVarId, extinction(:, :, :))
          ncStatus( 9) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "SingleScatteringAlbedo", ncVarId)
          ncStatus(10) = nf90_get_var(ncFileId, ncVarId, singleScatteringAlbedo(:, :, :))
          ncStatus(11) = nf90_inq_varid(ncFileId, trim(makePrefix(i)) // "PhaseFunctionIndex", ncVarId)
          ncStatus(12) = nf90_get_var(ncFileId, ncVarId, phaseFunctionIndex(:, :, :))
        end if
        if(any(ncStatus(:) /= nf90_NoErr))                                                                &
          call setStateToFailure(status, "read_Domain: Error reading scalar fields from file " // &
                                         trim(fileName))
        !
        ! Read in the phase function table(s) 
        !
        if(.not. stateIsFailure(status)) &
          call read_PhaseFunctionTable(fileId = ncFileId, table = table,                  &
                                       prefix = "Component" // trim(IntToChar(i)) // "_", &
                                       status = status)
          
        !
        ! Add the new component to the domain. 
        !
        if(.not. stateIsFailure(status)) then
          call addOpticalComponent(thisDomain, name, extinction, singleScatteringAlbedo, &
                                   phaseFunctionIndex, table, zLevelBase = zLevelBase,   &
                                   status = status)
        else
          call setStateToFailure(status, "read_Domain: Error reading phase function table.") 
        end if
        deallocate(extinction, singleScatteringAlbedo, phaseFunctionIndex)
      end do
      if(.not. stateIsFailure(status)) call setStateToSuccess(status)
    end if
  end subroutine read_Domain
  !------------------------------------------------------------------------------------------
  ! Finalization
  !------------------------------------------------------------------------------------------
  subroutine finalize_Domain(thisDomain) 
    ! Return the variable to it uninitialized state
    type(domain), intent(out) :: thisDomain
    
    ! Loca variable
    integer :: i
    
    thisDomain%xyRegularlySpaced = .false.
    thisDomain%zRegularlySpaced  = .false.
    
    ! Position vectors
    if(associated(thisDomain%xPosition)) deallocate(thisDomain%xPosition)
    if(associated(thisDomain%yPosition)) deallocate(thisDomain%yPosition)
    if(associated(thisDomain%zPosition)) deallocate(thisDomain%zPosition)
    
    !  Optical components
    if(containsComponents(thisDomain)) then
      do i = 1, size(thisDomain%components)
        call finalizeComponent(thisDomain%components(i))
      end do
      deallocate(thisDomain%components)
    end if
    
  end subroutine finalize_Domain
  !------------------------------------------------------------------------------------------
  subroutine finalizeComponent(component) 
    ! Return the variable to it uninitialized state
    type(opticalComponent), intent(out) :: component
    
    if(associated(component%extinction))             deallocate(component%extinction)
    if(associated(component%singleScatteringAlbedo)) deallocate(component%singleScatteringAlbedo)
    if(associated(component%phaseFunctionIndex))     deallocate(component%phaseFunctionIndex)
    
    call finalize_PhaseFunctionTable(component%table)

    component%horizontallyUniform     = .false. 
    component%zLevelBase              = 0
    component%name                    = ""
    
  end subroutine finalizeComponent
  !------------------------------------------------------------------------------------------
  ! Utilities, used inside the module
  !------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------
  ! Optical component objects aren't visible outside the module. so this creation routine isn't
  !   visible either.
  !------------------------------------------------------------------------------------------
  function newOpticalComponent(componentName, extinction, singleScatteringAlbedo, &
                               phaseFunctionIndex, zLevelBase, phaseFunctions) result(newOne)
    character (len = *),         intent( in) :: componentName
    real,    dimension(:, :, :), intent( in) :: extinction, singleScatteringAlbedo
    integer, dimension(:, :, :), intent( in) :: phaseFunctionIndex
    integer,           optional, intent( in) :: zLevelBase
    type(phaseFunctionTable),    intent( in) :: phaseFunctions
    type(opticalComponent)                   :: newOne
        
    ! Local variables
    integer :: nx, ny, nz
    
    ! Optical components are always associated with a domain. 
    !   So all error checking (grid sizes, etc.) that might be done here is instead done
    !   in validateOpticalComponent, which knows about the grid). 
    nx = size(extinction, 1); ny = size(extinction, 2); nz = size(extinction, 3)
    newOne%name = trim(componentName)
    newOne%zLevelBase = zLevelBase
    allocate(newOne%extinction(nx, ny, nz), newOne%singleScatteringAlbedo(nx, ny, nz), &
             newOne%phaseFunctionIndex(nx, ny, nz))
    newOne%extinction(:, :, :)             = extinction(:, :, :)
    newOne%singleScatteringAlbedo(:, :, :) = singleScatteringAlbedo(:, :, :)
    newOne%phaseFunctionIndex(:, :, :)     = phaseFunctionIndex(:, :, :)
    newOne%table = copy_PhaseFunctionTable(phaseFunctions)
    if(nx == 1 .and. ny == 1) newOne%horizontallyUniform = .true.
    
    if(present(zLevelBase)) then
      newOne%zLevelBase = zLevelBase
    else 
      newOne%zLevelBase = 1
    end if  
    
  end function newOpticalComponent
  !------------------------------------------------------------------------------------------
  subroutine validateOpticalComponent(thisDomain, componentName,          &
                                      extinction, singleScatteringAlbedo, &
                                      phaseFunctionIndex, table,          &
                                      zLevelBase, status)
    type(domain),                intent(inout) :: thisDomain
    character (len = *),         intent(in   ) :: componentName
    real,    dimension(:, :, :), intent(in   ) :: extinction, singleScatteringAlbedo
    integer, dimension(:, :, :), intent(in   ) :: phaseFunctionIndex
    type(phaseFunctionTable),    intent(in   ) :: table
    integer,                     intent(in   ) :: zLevelBase
    type(ErrorMessage),          intent(inout) :: status
    !
    ! Check to be sure that data being used to add or replace components is compatible 
    !   with the existing domain. This is folded into a separate subroutine because the logic
    !   is complicated enough that we don't want to replicate it everywhere. 
    

    ! Local variables
    integer :: numX, numY, numZ, numPhaseFunctions
        
    if(isValid(thisDomain)) then 
      numX = size(extinction, 1); numY = size(extinction, 2); numZ = size(extinction, 3)
      
      ! Are the arrays all the same size? 
      if(any((/ size(singleScatteringAlbedo, 1), size(phaseFunctionIndex, 1) /) /= numX) .or. &
         any((/ size(singleScatteringAlbedo, 2), size(phaseFunctionIndex, 2) /) /= numY) .or. &
         any((/ size(singleScatteringAlbedo, 3), size(phaseFunctionIndex, 3) /) /= numZ)) &
        call setStateToFailure(status, "validateOpticalComponent: optical property grids must be the same size.") 
        
      ! Do the arrays conform to the grid in the domain? 
      if(.not. any(numX == (/ 1, size(thisDomain%xPosition) - 1 /)) .or. &
         .not. any(numY == (/ 1, size(thisDomain%yPosition) - 1/)))      &
        call setStateToFailure(status, "validateOpticalComponent: arrays don't conform to horizontal extent of domain.")
      if(zLevelBase + numZ - 1 > size(thisDomain%zPosition) .or. zLevelBase < 1) &
        call setStateToFailure(status, "validateOpticalComponent: arrays don't conform to vertical extent of domain.")
     
      ! Resonable values for the properties
      if(any(extinction(:, :, :) < 0.)) &
        call setStateToFailure(status, "validateOpticalComponent: extinction must be >= 0.")
      if(any(singleScatteringAlbedo(:, :, :) < 0.) .or. any(singleScatteringAlbedo(:, :, :) > 1. )) &
        call setStateToFailure(status, "validateOpticalComponent: singleScatteringAlbedo must be between 0 and 1")
      
      ! Check the phase function table
      if(.not. stateIsFailure(status)) then 
        call getInfo_PhaseFunctionTable(table, nEntries = numPhaseFunctions, status = status)
        if(any(phaseFunctionIndex(:, :, :) < 0 .or. phaseFunctionIndex(:, :, :) > numPhaseFunctions)) &
          call setStateToFailure(status, "validateOpticalComponent: phase function index is out of bounds")
        ! Are the phase functions ready to go ? 
        if(.not. isReady_PhaseFunctionTable(table)) &
          call setStateToFailure(status, "validateOpticalComponent: phase function table is not ready.")
      end if 
      
      ! We could check to see if the component names overlap. 
    else
      call setStateToFailure(status, "validateOpticalComponent: domain hasn't been initialized.")
    end if

    if(.not. stateIsFailure(status)) call setStateToSuccess(status)
  end subroutine validateOpticalComponent
  !------------------------------------------------------------------------------------------
  logical function isValid(thisDomain) 
    type(domain), intent(in) :: thisDomain
    ! Checks to see if domain is associated with an initialized
    !   object. 
    
    isValid = associated(thisDomain%xPosition) .and. &
              associated(thisDomain%yPosition) .and. &
              associated(thisDomain%zPosition)
  end function isValid
  !------------------------------------------------------------------------------------------
  logical function containsComponents(thisDomain) 
    type(domain), intent(in) :: thisDomain
    ! Checks to see if this domain contains any optical components. 
    
    containsComponents = associated(thisDomain%components)
  end function containsComponents
  !------------------------------------------------------------------------------------------
  function makePrefix(i) 
    integer, intent( in) :: i
    character(len = 32)  :: makePrefix
    ! Constructs a unique prefix for each component
    !   We put this in a function so reading and writing subroutine are consistent
    !   Uses utility IntToChar from module characterUtils 
    
    character(len = 16), parameter :: prefixBase = "Component"
    
    makePrefix = trim(prefixBase) // trim(IntToChar(i)) // "_"
  end function makePrefix
  !------------------------------------------------------------------------------------------
  function makePhaseFunctionFileName(baseFileName, i)
    character(len = *), intent( in) :: baseFileName 
    integer,            intent( in) :: i
    character(len = nf90_max_name)  :: makePhaseFunctionFileName
    
    makePhaseFunctionFileName  = trim(baseFileName) // "_C" // trim(IntToChar(i)) // "_pft" 
  end function makePhaseFunctionFileName
  !------------------------------------------------------------------------------------------
  function makeInvPhaseFunctionFileName(baseFileName, i)
    character(len = *), intent( in) :: baseFileName 
    integer,            intent( in) :: i
    character(len = nf90_max_name)  :: makeInvPhaseFunctionFileName
    
    makeInvPhaseFunctionFileName = trim(baseFileName) // "_C" // trim(IntToChar(i)) // "_ipf" 
  end function makeInvPhaseFunctionFileName
  !------------------------------------------------------------------------------------------
  elemental function asInt(inValue)
    logical,                 intent( in) :: inValue
    integer(kind = selected_Int_Kind(2)) :: asInt
    
    asInt = 0; if(inValue) asInt = 1
    
  end function asInt
  !------------------------------------------------------------------------------------------
  elemental function asLogical(inValue)
    integer(kind = selected_Int_Kind(2)), intent( in) :: inValue
    logical :: asLogical
    
    asLogical = .not. (inValue == 0)
    
  end function asLogical
  !------------------------------------------------------------------------------------------
end module opticalProperties   
