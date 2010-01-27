! $Revision$, $Date$
! $URL$
module kDistributions
  !
  ! This module encapsulates a k-distribution applied to a specified atmosphere to create
  !    a weighted set of profiles of gaseous absorption within some portion of the spectrum. 
  !    Some (complicated) external program takes a profile of atmospheric temperature, 
  !    pressure, and gas concentrations, and produces a set of absorption profiles and the 
  !    fraction of the band ("weight") corresponding to each profile of absorption. 
  !    The k-distribution applies to some portion of the spectrum (a "band"); typically this is 
  !    a small enough region that the optical properties of cloud drops and ice crystals 
  !    can be considered constant. A k-distribution could also be generated for 
  !    an instrument channel. 
  ! This module provides a way to store the sets of absorption profiles produced by 
  !   by an external program. One k-distribution must be created for each band. 
  !
  use ErrorMessages
  implicit none 
  private
  
  !
  ! Module parameters
  !
  integer, parameter :: maxDescriptionLength = 256
  real,    parameter :: defaultValue = 0. 
  
  !
  ! Module type (object)
  !
  type kDistribution
    private
    real, dimension(:, :), pointer        :: absorptionProfiles => null() ! Dimensions numZs, numKs
    real, dimension(:),    pointer        :: absorptionWeights  => null(), & ! numKs
                                             zPosition          => null()    ! numZs
    character(len = maxDescriptionLength) :: description = ""
    real, dimension(2)                    :: wavelengthLimits = (/ defaultValue, defaultValue /) 
    real                                  :: spectralFraction = defaultValue 
  end type kDistribution
  
  !
  ! What is visible? The type defined by the module...
  !
  public :: kDistribution 
  !
  ! ... and functions to operate on that type.  
  !
  public :: new_kDistribution, getInfo_kDistribution, copy_kDistribution,   &
            read_kDistribution, write_kDistribution, isReady_kDistribution, &
            getAbsorptionProfile, getAbsorptionProfiles
contains
! ------------------------------------------
  function new_kDistribution(zPosition, absorptionProfiles, absorptionWeights, &
                             description, wavelengthLimits, spectralFraction, status) result(kDist)
    real, dimension(:),    intent(in   ) :: zPosition
    real, dimension(:, :), intent(in   ) :: absorptionProfiles
    real, dimension(:),    intent(in   ) :: absorptionWeights
    character(len = *),    intent(in   ) :: description
    real, dimension(2), optional, & 
                           intent(in   ) :: wavelengthLimits
    real,               optional, &
                           intent(in   ) :: spectralFraction
    type(ErrorMessage),    intent(inout) :: status
    type(kDistribution)                :: kDist
    
    ! Local variables
    integer :: numZs, numKs 
    
    numZs = size(zPosition); numKs = size(absorptionWeights)
    !
    ! Input validation - zPosition, absorptionProfiles, absorptionWeights size must match; 
    !   zPosition must increase (negative values are allowed)
    !   wavelengthLimits have to be positive and increasing
    !   spectralFraction must be between 0 and 1
    !
    if(size(absorptionProfiles, 1) /= numZs .or. size(absorptionProfiles, 2) /= numKs) &
      call setStateToFailure(status, "new_kDistribution: absorptionProfiles must be dimensioned (nLevels, nKs).")
    if(any(zPosition(2:) - zPosition(:numZs - 1) <= 0.)) &
      call setStateToFailure(status, "new_kDistribution: zPosition must be increasing, unique.")
    if(present(wavelengthLimits)) then 
      if(any(wavelengthLimits(:) <= 0) .or. wavelengthLimits(2) < wavelengthLimits(1)) &
        call setStateToFailure(status, "new_kDistribution: wavelengthLimits must be positive, increasing.")
    end if
    if(present(spectralFraction)) then 
      if(spectralFraction <= 0 .or. spectralFraction > 1) &
        call setStateToFailure(status, "new_kDistribution: spectralFraction must be between 0 and 1")
    end if
    
    if(.not. stateIsFailure(status)) then 
      allocate(kDist%zPosition(numZs),      &
               kDist%absorptionWeights(numKs), &
               kDist%absorptionProfiles(numZs, numKs))
      kDist%zPosition(:)             = zPosition(:)
      kDist%absorptionProfiles(:, :) = absorptionProfiles(:, :) 
      kDist%absorptionWeights(:)     = absorptionWeights(:)
      kDist%description              = trim(description)
      if(present(wavelengthLimits)) kDist%wavelengthLimits = wavelengthLimits
      if(present(spectralFraction)) kDist%spectralFraction = spectralFraction
      call setStateToSuccess(status)
    end if 
  end function new_kDistribution
! ------------------------------------------                  
  subroutine getInfo_kDistribution(kdist, nLevels, zPosition,                        &
                                   nProfiles, absorptionProfiles, absorptionWeights, &
                                   description, wavelengthLimits, spectralFraction, status)
    type(kDistribution),             intent(in  ) :: kDist
    integer,               optional, intent(  out) :: nLevels
    real, dimension(:),    optional, intent(  out) :: zPosition
    integer,               optional, intent(  out) :: nProfiles
    real, dimension(:, :), optional, intent(  out) :: absorptionProfiles
    real, dimension(:),    optional, intent(  out) :: absorptionWeights
    character (len = *),   optional, intent(  out) :: description
    real, dimension(:),    optional, intent(  out) :: wavelengthLimits
    real,                  optional, intent(  out) :: spectralFraction
    type(ErrorMessage),              intent(inout) :: status
    !
    !  Provides information about the k-distribution object, including the native 
    !    grid on which the extinction profiles are defined. 
    !
    
    integer :: numKs, numZs
    !
    ! Input validation - kDistribution must be ready 
    !   zPosition, absorptionProfiles, absorptionWeights sizes must match; 
    !   wavelengthLimits should be 2 elements long
    !   description string should be long enough
    !
    if(.not. isReady_kDistribution(kDist)) then
      call setStateToFailure(status, "getInfo_kDistribution: kDistribution has not been initialized.") 
    else
      numKs = size(kDist%absorptionWeights); numZs = size(kDist%zPosition)
      
      if(present(zPosition)) then
        if(size(zPosition) /= numZs) &
          call setStateToFailure(status, "getInfo_kDistribution: wrong number of heights requested.") 
      end if
      if(present(absorptionWeights)) then
        if(size(absorptionWeights) /= numKs) &
          call setStateToFailure(status, "getInfo_kDistribution: wrong number of weights requested.") 
      end if
      if(present(absorptionProfiles)) then
        if(any( (/ size(absorptionProfiles, 1), size(absorptionProfiles, 2) /)  /= (/ numZs, numKs /) )) &
          call setStateToFailure(status, "getInfo_kDistribution: absorptionProfiles has the wrong dimensions.") 
      end if
      if(present(wavelengthLimits)) then
        if(size(wavelengthLimits) < 2) &
          call setStateToFailure(status, "getInfo_kDistribution: wavelengthLimits must be at least two elements long.") 
        if(size(wavelengthLimits) > 2) &
          call setStateToWarning(status, "getInfo_kDistribution: wavelengthLimits should only be two elements long.") 
      end if
      if(present(description)) then
        if(len(description) < len_trim(kDist%description)) &
          call setStateToWarning(status, "getInfo_kDistribution: space allocated for description is too small.") 
      end if
    end if 
    
    if(.not. stateIsFailure(status)) then 
      if(present(nLevels))            nLevels                  = numZs
      if(present(zPosition))          zPosition(:)             = kDist%zPosition(:)
      if(present(nProfiles))          nProfiles                = numKs
      if(present(absorptionProfiles)) absorptionProfiles(:, :) = kDist%absorptionProfiles(:, :)
      if(present(absorptionWeights))  absorptionWeights(:)     = kDist%absorptionWeights(:)
      if(present(description))        description              = trim(kdist%description)
      if(present(wavelengthLimits))   wavelengthLimits(:2)     = kDist%wavelengthLimits(:)
      if(present(spectralFraction))   spectralFraction         = kDist%spectralFraction
      call setStateToSuccess(status)
    end if
  end subroutine getInfo_kDistribution
! ------------------------------------------
  subroutine getAbsorptionProfile(kdist, profileNumber, zPosition, &
                                  absorptionProfile, weight, status)
    !
    !  Provide profile of gasesous extinction for the i'th k value. The profile
    !    can be requested on an arbitrary grid specified by zPosition; the original 
    !    profiles will be interpolated/integrated onto this grid. 
    !
    type(kDistribution), intent(in   ) :: kDist
    integer,             intent(in   ) :: profileNumber
    real, dimension(:),  intent(in   ) :: zPosition
    real, dimension(:),  intent(  out) :: absorptionProfile
    real,                intent(  out) :: weight
    type(ErrorMessage),  intent(inout) :: status
    
    ! Warning if: 
    ! Fails if: 
    !  kDist isn't ready
    !  user asks for a profile number < 1 or > numKs
    !  zPosition and absorptionProfile aren't the same length
    ! What happens if 
    !   some of the requested levels are outside the region where the k-distribution is defined? (warning)
    !   all  of the requested levels are outside the region where the k-distribution is defined? (failure)
    !   all the zPositions are within one layer of the region where the k-distribution is defined? (warning? nothing?) 
    
  end subroutine getAbsorptionProfile
! ------------------------------------------
  subroutine getAbsorptionProfiles(kdist, zPosition, &
                                   absorptionProfiles, weights, status)
    type(kDistribution),   intent(in   ) :: kDist
    real, dimension(:),    intent(in   ) :: zPosition
    real, dimension(:, :), intent(  out) :: absorptionProfiles
    real, dimension(:),    intent(  out) :: weights
    type(ErrorMessage),    intent(inout) :: status
    
    ! Local variables
    integer :: numKs, i
    
    ! Warning if: 
    ! Fails if: 
    
    do i = 1, numKs
      call getAbsorptionProfile(kdist, i, zPosition, &
                                absorptionProfiles(i), weights(i), status)
    end do 
  end subroutine getAbsorptionProfiles
! ------------------------------------------
  function isReady_kDistribution(kdist)
    type(kDistribution), intent(in) :: kDist
    logical                         :: isReady_kDistribution
    !
    ! Has this variable been initialized, so we can use it? 
    !

    isReady_kDistribution = associated(kdist%absorptionProfiles)    
    
  end function isReady_kDistribution
! ------------------------------------------
  function copy_kDistribution(kdist) result(new)
    type(kDistribution), intent(in) :: kDist
    type(kDistribution)             :: new
    !
    ! Create a new copy of the data stored in the kDistribution object. 
    !   If the copy fails for some reason the user won't hear about it. 
    
    ! Local variables
    type(ErrorMessage) :: status

    if(isReady_kDistribution(kDist)) then 
      new = new_kDistribution(kDist%zPosition,                                &
                              kDist%absorptionProfiles, kDist%absorptionWeights, &
                              kDist%description, kDist%wavelengthLimits,      &
                              kDist%spectralFraction, status)
    else
      call setStateToFailure(status, "copy: source kDistribution hasn't been initialized.") 
    end if 
  end function copy_kDistribution
! ------------------------------------------
  subroutine finalize_kDistribution(kdist)
    type(kDistribution), intent(inout) :: kDist
    !
    ! Leave the variable in a pristine state - free all memory and 
    !   reset values to initial 
    !
    
    if(associated(kdist%absorptionProfiles)) deallocate(kdist%absorptionProfiles)
    if(associated(kdist%absorptionWeights    )) deallocate(kdist%absorptionWeights)
    if(associated(kdist%zPosition         )) deallocate(kdist%zPosition)
    kDist%description = ""
    kDist%wavelengthLimits = (/ defaultValue, defaultValue /) 
    kDist%spectralFraction = defaultValue 
  end subroutine finalize_kDistribution
! ------------------------------------------
  subroutine read_kDistribution(fileName, kDist, status)
    character (len = *), intent(in   ) :: fileName
    type(kDistribution), intent(  out) :: kDist
    type(ErrorMessage),  intent(inout) :: status
  
    ! Fails if: kdist isn't ready
    !          file can't be written 
end subroutine read_kDistribution
! ------------------------------------------
  subroutine write_kDistribution(kDist, fileName, status)
    type(kDistribution), intent(in   ) :: kDist
    character (len = *), intent(in   ) :: fileName
    type(ErrorMessage),  intent(inout) :: status
    
    ! Fails if: file can't be read 
    !  end module k-Distribution
  end subroutine write_kDistribution  
end module kDistributions