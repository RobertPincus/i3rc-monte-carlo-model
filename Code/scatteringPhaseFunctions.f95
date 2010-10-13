! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision$, $Date$
! $URL$
module scatteringPhaseFunctions
  ! Provides a representation of the scattering phase function for 
  !   use in radiative transfer calculations. Users may supply 
  !   either values of the phase function versus scattering angle, or
  !   the coefficients in the Legendre expansion of the phase function.
  !   They may request a stored phase function in either representation. 
  ! Users may optionally provide single scattering albedo and extinction in m^2/particle. 
  
  use ErrorMessages,    only: ErrorMessage,   &
                              stateIsFailure, &
                              setStateToFailure, setStateToWarning, setStateToSuccess
  use CharacterUtils,   only: charToInt, intToChar
  use numericUtilities, only: findIndex, computeLobattoTerms, computeLegendrePolynomials
  implicit none
  private 
  
  ! Scattering angle is in radians
  ! Is there a way to make this transparent?
  real,    parameter :: Pi = 3.141592654, radiansToDegrees = 180. / (Pi/2.)
  real,    parameter :: minScatteringAngle = 0., maxScatteringAngle = Pi ! radians
  real,    parameter :: defaultExtinction = 0.,  defaultSSA = 0.
  integer, parameter :: maxSingleDescriptorLength = 64, maxTableDescriptorLength = 1024

  !------------------------------------------------------------------------------------------
  ! Type (object) definitions
  !------------------------------------------------------------------------------------------
  type phaseFunction
    ! This holds the phase function as supplied. 
    !   It has members for both types of representations
    !   but only the one in use gets allocated. 
    private
    real, dimension(:), pointer :: scatteringAngle  => null()
    real, dimension(:), pointer :: value            => null()
    real, dimension(:), pointer :: legendreCoefficients  => null()
    real                        :: extinction             = defaultExtinction, &
                                   singleScatteringAlbedo = defaultSSA
    character(len = maxSingleDescriptorLength) &
                                :: description = ""
  end type phaseFunction
  
  type phaseFunctionTable
    ! This holds a series of phase functions, a real valued key into the table, 
    !   and an optional description
    private 
    type(phaseFunction), dimension(:), pointer :: phaseFunctions  => null()
    real,                dimension(:), pointer :: key             => null()
    character(len = maxSingleDescriptorLength), &
                         dimension(:), pointer :: phaseFunctionDescriptions  => null()
    character(len = maxTableDescriptorLength)  :: description = ""
    logical                                    :: oneAngleSet = .false.
  end type phaseFunctionTable
  
  !------------------------------------------------------------------------------------------
  ! Overloading
  !------------------------------------------------------------------------------------------
  ! The initialization procedure for individual phase functions ("newPhaseFunction") 
  !   can be called with scatteringAngle and phase function value, or with 
  !   the Legendre expansion coefficients, but not with both. 
  interface new_PhaseFunction
    module procedure newphaseFunctionTabulated, newPhaseFunctionExpansion
  end interface ! new_PhaseFunction
  
  ! The initialization procedure for phase function tables ("newPhaseFunctionTable") 
  !   can be called with a vector of N scatteringAngles and an NxM matrix of M phase 
  !   function values, or with a vector of phase functions. 
  interface new_PhaseFunctionTable
    module procedure newPhaseFunctionTableTabulated, newPhaseFunctionTableGeneral
  end interface ! newPhaseFunctionTable
  
  interface getPhaseFunctionValues
    module procedure getPhaseFunctionValues_one, getPhaseFunctionValues_table
  end interface ! getPhaseFunctionValues
  !------------------------------------------------------------------------------------------
  ! What is visible? 
  !------------------------------------------------------------------------------------------
  ! The types...
  public :: phaseFunction, phaseFunctionTable
  ! ... and these procedures
  public :: new_PhaseFunction,       new_PhaseFunctionTable,      &
            copy_PhaseFunction,      copy_PhaseFunctionTable,     &
            getInfo_PhaseFunction,   getInfo_PhaseFunctionTable,  &
            read_PhaseFunctionTable, add_PhaseFunctionTable,      &
            write_PhaseFunctionTable,                             &
            finalize_PhaseFunction,  finalize_PhaseFunctionTable, &
            isReady_PhaseFunction,   isReady_PhaseFunctionTable,  &
            getElement, getExtinction, getSingleScatteringAlbedo, &
            getPhaseFunctionValues, getPhaseFunctionCoefficients
            
contains
  !------------------------------------------------------------------------------------------
  ! Code
  !------------------------------------------------------------------------------------------
  ! Initialization: Routines to create new phase function objects
  !------------------------------------------------------------------------------------------
  function newPhaseFunctionTabulated(scatteringAngle, value, &
                                     extinction, singleScatteringAlbedo, description, status) result(newPhaseFunction)
    real, dimension(:),           intent(in   ) :: value, scatteringAngle
    real,               optional, intent(in   ) :: extinction, singleScatteringAlbedo
    character(len = *), optional, intent(in   ) :: description
    type(ErrorMessage),           intent(inout) :: status
    type(phaseFunction)                         :: newPhaseFunction
    ! ------------------------
    ! Initialization - scattering angle and phase function supplied
    
    ! Local variables
    integer :: nValues
    
    ! ------------------------
    ! Sanity checks
    !   angles in bounds, ascending, unique, have correct limits; 
    !   arrays same length; value always > 0. 
    nValues = size(scatteringAngle)
    if(any(scatteringAngle(:) < minScatteringAngle) .or. &
       any(scatteringAngle(:) > maxScatteringAngle))     &
      call setStateToFailure(status, "newPhaseFunction: ScatteringAngle out of bounds.")
    if(abs(scatteringAngle(1      ) - minScatteringAngle) > spacing(minScatteringAngle)) &
      call setStateToFailure(status, "newPhaseFunction: First scattering angle must be min value")
    if(abs(scatteringAngle(nValues) - maxScatteringAngle) > spacing(maxScatteringAngle)) &
      call setStateToFailure(status, "newPhaseFunction: Last scattering angle must be max value")
    if(any(scatteringAngle(2:) - scatteringAngle(:nValues-1) <= 0.)) &
      call setStateToFailure(status, "newPhaseFunction: Scattering angle must be increasing, unique.")
    if(any(value(:) < 0.)) &
      call setStateToFailure(status, "newPhaseFunction: Negative phase function values supplied.")
    if(size(scatteringAngle) /= size(value)) &
      call setStateToFailure(status,         &
            "newPhaseFunction: Number of scattering angles and phase function values must match.")
    if(present(extinction)) then
      if(extinction             < 0.) call setStateToFailure(status, "newPhaseFunction: negative extinction supplied.") 
    end if
    if(present(singleScatteringAlbedo)) then
      if(singleScatteringAlbedo < 0. .or. singleScatteringAlbedo > 1.) &
        call setStateToFailure(status, "newPhaseFunction: singleScatteringAlbedo out of bounds.") 
    end if
    if(present(description)) then
      if(len_trim(description) > maxSingleDescriptorLength) &
        call setStateToWarning(status, "newPhaseFunction: description will be trunctated.")
    end if
            
    ! ------------------------
    !  All the values passed are reasonable. Allocate arrays for the angles and the phase 
    !    function value at those angles. 
    !
    if(.not. stateIsFailure(status)) then
      allocate(newPhaseFunction%scatteringAngle(nValues))
      allocate(newPhaseFunction%value(          nValues))
      newPhaseFunction%scatteringAngle(:) = scatteringAngle(:)
      newPhaseFunction%value(:)   = normalizePhaseFunction(scatteringAngle, value)
      
      if(present(extinction))             newPhaseFunction%extinction             = extinction 
      if(present(singleScatteringAlbedo)) newPhaseFunction%singleScatteringAlbedo = singleScatteringAlbedo
      if(present(description))            newPhaseFunction%description            = description
      
      call setStateToSuccess(status)
    end if
  end function newPhaseFunctionTabulated
  !------------------------------------------------------------------------------------------
  function newPhaseFunctionExpansion(legendreCoefficients, extinction, singleScatteringAlbedo, &
                                     description, status) result(newPhaseFunction)
    real, dimension(:),           intent(in   ) :: legendreCoefficients
    real,               optional, intent(in   ) :: extinction, singleScatteringAlbedo
    character(len = *), optional, intent(in   ) :: description
    type(ErrorMessage),           intent(inout) :: status
    type(phaseFunction)                         :: newPhaseFunction
    
    integer, parameter :: numTestSteps = 1801  
    integer            :: i 
    real, dimension(numTestSteps) &
                       :: testAngles, testValues
    ! ------------------------
    ! Initialization - Legendre moments supplied
    
    ! Sanity checks     
    ! P0 is identically 1, and the supplied coefficients should start with P1. 
    if(size(legendreCoefficients) > 1) then
      if(legendreCoefficients(1) > 1. .or. legendreCoefficients(1) < -1.) &
        call setStateToFailure(status, "newPhaseFunction: Asymmetery parameter out of bounds.") 

      if(abs(legendreCoefficients(1) - 1.0) < spacing(1.)) &
        call setStateToWarning(status, "newPhaseFunction: Should the first legendre moment (P1) really be 1?")
    end if
    if(present(extinction)) then
      if(extinction             < 0.) call setStateToFailure(status, "newPhaseFunction: negative extinction supplied.") 
    end if
    if(present(singleScatteringAlbedo)) then
      if(singleScatteringAlbedo < 0. .or. singleScatteringAlbedo > 1.) &
        call setStateToFailure(status, "newPhaseFunction: singleScatteringAlbedo out of bounds.") 
    end if
    if(present(description)) then
      if(len_trim(description) > maxSingleDescriptorLength) &
        call setStateToWarning(status, "newPhaseFunction: description will be trunctated.")
    end if
    ! ------------------------
    ! All the values passed in look reasonable
    ! 
    if(.not. stateIsFailure(status)) then
      allocate(newPhaseFunction%legendreCoefficients(size(legendreCoefficients)))
      newPhaseFunction%legendreCoefficients(:) = legendreCoefficients(:)

      if(present(extinction))             newPhaseFunction%extinction             = extinction 
      if(present(singleScatteringAlbedo)) newPhaseFunction%singleScatteringAlbedo = singleScatteringAlbedo
      if(present(description))            newPhaseFunction%description            = description
      
      !
      ! Sanity check - expand the phase function using numTestSteps and see if there are any 
      !   negative values
      !
      testAngles(:) = (/ (i, i = 0, numTestSteps - 1) /) / real(numTestSteps - 1) * Pi
      call getPhaseFunctionValues(newPhaseFunction, testAngles, testValues, status)
      if(any(testValues(:) < 0.)) then 
        call setStateToWarning(status, &
                               "newPhaseFunction: Phase function coefficients give "           // &
                               trim(IntToChar((100 * count(testValues(:) < 0.))/numTestSteps)) // &
                               "% negative phase function values" )
      else
        call setStateToSuccess(status)
      end if
    end if
  end function newPhaseFunctionExpansion
  !------------------------------------------------------------------------------------------
  function newPhaseFunctionTableTabulated(scatteringAngle, values, key,       &
                                          extinction, singleScatteringAlbedo, &
                                          phaseFunctionDescriptions, tableDescription, status) result(table)
    real, dimension(:),            intent(in   ) :: scatteringAngle
    real, dimension(:, :),         intent(in   ) :: values
    real, dimension(:),            intent(in   ) :: key
    real, dimension(:), optional,  intent(in   ) :: extinction, singleScatteringAlbedo
    character(len = *), &
          dimension(:), optional,  intent(in   ) :: phaseFunctionDescriptions
    character(len = *), optional,  intent(in   ) :: tableDescription
    type(ErrorMessage),            intent(inout) :: status
    type(phaseFunctionTable)                     :: table
    ! ------------------------
    ! Initializing a table of phase functions given one set of scattering angles and a set
    !   of phase functions
    ! Array values is dimensioned(nAngles, nEntries)
    
    ! Local variables
    integer :: nAngles, nEntries, i
    
    ! ------------------------
    ! Sanity checks
    nAngles = size(scatteringAngle); nEntries = size(values, 2)
    ! Sanity checks
    !   angles in bounds, ascending, unique, have correct limits; 
    !   arrays same length; value always > 0. 
    if(any(scatteringAngle(:) < minScatteringAngle) .or. &
       any(scatteringAngle(:) > maxScatteringAngle))     &
      call setStateToFailure(status, "newPhaseFunctionTable: ScatteringAngle out of bounds.")
    if(abs(scatteringAngle(1      ) - minScatteringAngle) > spacing(minScatteringAngle)) &
      call setStateToFailure(status, "newPhaseFunctionTable: First scattering angle must be min value")
    if(abs(scatteringAngle(nAngles) - maxScatteringAngle) > spacing(maxScatteringAngle)) &
      call setStateToFailure(status, "newPhaseFunctionTable: Last scattering angle must be max value")
    if(any(scatteringAngle(2:) - scatteringAngle(:nAngles-1) <= 0.)) &
      call setStateToFailure(status, "newPhaseFunctionTable: Scattering angle must be increasing, unique.")
    if(any(values(:, :) < 0.)) &
      call setStateToFailure(status, "newPhaseFunctionTable: Negative phase function values supplied.")
    if(size(values, 1) /= nAngles) &
      call setStateToFailure(status,         &
            "newPhaseFunctionTable: Number of scattering angles and phase function values must match.")
    if(size(key) /= nEntries) &
      call setStateToFailure(status,       &
          "newPhaseFunctionTable: Number of phase functions and key values must match.")
    if(any(key(2:) - key(1:nEntries-1) <= 0.)) &
      call setStateToFailure(status,           &
          "newPhaseFunctionTable: Key values must be unique, increasing.")
    if(present(extinction)) then
      if(size(extinction) /= nEntries ) &
        call setStateToFailure(status, "newPhaseFunctionTable: extinction must be provided for each phase function.")
      if(any(extinction(:)             < 0.)) call setStateToFailure(status, "newPhaseFunction: negative extinction supplied.")
    end if 
    if(present(singleScatteringAlbedo)) then 
      if(size(singleScatteringAlbedo) /= nEntries) &
        call setStateToFailure(status, "newPhaseFunctionTable: single scattering albedo must be provided for each phase function.")
      if(any(singleScatteringAlbedo(:) < 0. .or. singleScatteringAlbedo(:) > 1.)) &
                                      call setStateToFailure(status, "newPhaseFunction: singleScatteringAlbedo must be > 0, <= 1.") 
    end if
    if(present(phaseFunctionDescriptions)) then
      if(size(phaseFunctionDescriptions) /= nEntries) &
        call setStateToFailure(status, "newPhaseFunctionTable: number of descriptions must match number of phase functions.")
      if(len_trim(phaseFunctionDescriptions(1)) > maxSingleDescriptorLength) &
        call setStateToWarning(status, "newPhaseFunctionTable: phaseFunctionDescriptions will be trunctated.")
    end if
    if(present(tableDescription)) then
      if(len_trim(tableDescription) > maxTableDescriptorLength) &
        call setStateToWarning(status, "newPhaseFunction: tableDescription will be trunctated.")
    end if
            
    ! ------------------------
    ! Everything looks OK to this point. 
    ! 
    if(.not. stateIsFailure(status)) then
      ! 
      ! Allocate the vector that holds the individual phase functions
      !
      allocate(table%phaseFunctions(nEntries), table%key(nEntries))
      !
      ! Now allocate each phase function in turn
      !   Because the scattering angles are the same for all phase functions we allocate 
      !   memory for the set only once, then have the scatteringAngle members of each 
      !   element in the array point to the first
      ! 
      allocate(table%phaseFunctions(1)%scatteringAngle(nAngles))
      do i = 1, nEntries
        allocate(table%phaseFunctions(i)%value(nAngles))
      end do
      
      table%phaseFunctions(1)%scatteringAngle(:) = scatteringAngle(:)
      table%phaseFunctions(1)%value(:)           = values(:, 1)
      forall(i = 2:nEntries)
        ! Phase function is normalized as in routine newPhaseFunctionTablulated
        table%phaseFunctions(i)%value(:) = values(:, i)
        table%phaseFunctions(i)%scatteringAngle => table%phaseFunctions(1)%scatteringAngle(:)
      end forall
      ! Normalize the phase functions
      do i = 1, nEntries
        table%phaseFunctions(i)%value(:) = normalizePhaseFunction(scatteringAngle(:), &
                                                                  table%phaseFunctions(i)%value(:))
      end do
      
      table%key(:) = key(:)
      ! Optional components
      if(present(extinction))                table%phaseFunctions(:)%extinction             = extinction(:)
      if(present(singleScatteringAlbedo))    table%phaseFunctions(:)%singleScatteringAlbedo = singleScatteringAlbedo(:)
      if(present(phaseFunctionDescriptions)) table%phaseFunctions(:)%description            = phaseFunctionDescriptions(:)
      if(present(tableDescription)) table%description = tableDescription
      
      table%oneAngleSet = .true. 
      call setStateToSuccess(status)
    end if
  end function newPhaseFunctionTableTabulated
  !------------------------------------------------------------------------------------------
  function newPhaseFunctionTableGeneral(phaseFunctions, key, phaseFunctionDescriptions, tableDescription, status) &
      result(table)
    type(phaseFunction), dimension(:), intent(in   ) :: phaseFunctions
    real, dimension(:),                intent(in   ) :: key
    character(len = *), &
          dimension(:), optional,      intent(in   ) :: phaseFunctionDescriptions
    character(len = *), optional,      intent(in   ) :: tableDescription
    type(ErrorMessage),                intent(inout) :: status
    type(phaseFunctionTable)                         :: table
    ! ------------------------
    ! Initializing a table of phase functions given a vector of type phaseFunction
    ! Array values is dimensioned(nAngles, nEntries)
    
    ! Local variables
    integer :: nEntries, i
    
    ! ------------------------
    ! Sanity checks on the input data
    !
    nEntries = size(phaseFunctions)
    ! Sanity checks - few are here; most were done when setting elements of the vector
    if(size(key) /= nEntries) &
      call setStateToFailure(status,       &
          "newPhaseFunctionTable: Number of phase functions and key values must match.")
    if(any(key(2:) - key(1:nEntries-1) <= 0.)) &
      call setStateToFailure(status,           &
          "newPhaseFunctionTable: Key values must be unique, increasing.")
    if(present(tableDescription)) then
      if(len_trim(tableDescription) > maxTableDescriptorLength) &
        call setStateToWarning(status, "newPhaseFunctionTable: tableDescription will be trunctated.")
    end if
    
    ! ------------------------
    if(.not. stateIsFailure(status)) then
      allocate(table%phaseFunctions(nEntries), table%key(nEntries))
     
      do i = 1, nEntries
        table%phaseFunctions(i) = copy_PhaseFunction(phaseFunctions(i))
      end do
      table%key(:) = key(:)
      if(present(phaseFunctionDescriptions)) table%phaseFunctions(:)%description = phaseFunctionDescriptions(:)
      if(present(tableDescription))          table%description                   = tableDescription
      table%oneAngleSet = .false. 
      call setStateToSuccess(status)
    end if
  
  end function newPhaseFunctionTableGeneral 
  !------------------------------------------------------------------------------------------
  function copy_PhaseFunction(original) result(thisCopy)
    ! Makes a copy of a scattering phase function. 
    !   Needed because simple assignment of the derived type doesn't duplicate
    !   the underlying data; this function does. 
    type(phaseFunction), intent( in) :: original
    type(phaseFunction)              :: thisCopy

    ! Local variables
    type(ErrorMessage) :: status
    
    ! ------------------------
    if(storedAsLegendre(original)) then
      thisCopy = new_PhaseFunction(original%legendreCoefficients,                        &
                                   original%extinction, original%singleScatteringAlbedo, &
                                   original%description, status)
    else if(storedAsTabulated(original)) then
      thisCopy = new_PhaseFunction(original%scatteringAngle, original%value,             &
                                   original%extinction, original%singleScatteringAlbedo, &
                                   original%description, status)
    else 
      call setStateToFailure(status, "copy: Can't determine how phase function is stored")
    end if 
  end function copy_PhaseFunction
  !------------------------------------------------------------------------------------------
  function copy_PhaseFunctionTable(original) result(thisCopy)
    ! Makes a copy of a scattering phase function. 
    !   Needed because simple assignment of the derived type doesn't duplicate
    !   the underlying data; this function does. 
    type(phaseFunctionTable), intent( in) :: original
    type(phaseFunctionTable)              :: thisCopy

    ! Local variables
    type(ErrorMessage)                 :: status
    integer                            :: nAngles, nEntries, i
    real, dimension(:, :), allocatable :: values
    real, dimension(:),    allocatable :: extinction, singleScatteringAlbedo
    
    ! ------------------------
    if(original%oneAngleSet) then
      nAngles = size(original%phaseFunctions(1)%scatteringAngle); nEntries = size(original%phaseFunctions)
      allocate(extinction(nEntries), singleScatteringAlbedo(nEntries), values(nAngles, nEntries))
      forall (i = 1:nEntries)
        extinction(i) = original%phaseFunctions(i)%extinction
        singleScatteringAlbedo(i) = original%phaseFunctions(i)%singleScatteringAlbedo
        values(:, i) = original%phaseFunctions(i)%value(:)
      end forall
      thisCopy = new_PhaseFunctionTable(original%phaseFunctions(1)%scatteringAngle, values(:, :), &
                                        original%key, extinction, singleScatteringAlbedo,         &
                                        tableDescription = original%description, status = status)
      deallocate(extinction, singleScatteringAlbedo, values)
    else 
      thisCopy = new_PhaseFunctionTable(original%phaseFunctions, original%key, &
      !                                 phaseFunctionDescriptions = original%phaseFunctionDescriptions,    &
                                        tableDescription = original%description, status = status)
    end if 
  end function copy_PhaseFunctionTable
  !------------------------------------------------------------------------------------------
  ! Getting the phase function back
  !------------------------------------------------------------------------------------------
  subroutine getPhaseFunctionValues_one(phaseFunctionVar, scatteringAngle, value, status)
    type(phaseFunction), intent(in   ) :: phaseFunctionVar
    real, dimension(:),  intent(in   ) :: scatteringAngle
    real, dimension(:),  intent(  out) :: value
    type(ErrorMessage),  intent(inout) :: status
    ! ------------------------
    ! Return the value of the phase function given a set of scattering angles. 
    !
    
    ! Local variables
    integer                               :: l, maxL, nValues, nStored
    integer, dimension(:),    allocatable :: tableIndicies, indiciesPlus1
    real,    dimension(:),    allocatable :: weights, dMu
    real,    dimension(:, :), allocatable :: legendreP
    
    ! ------------------------
    ! Sanity checks: 
    !   phase function contains values; angles in bounds; angles in order; arrays same length
    
    nValues = size(scatteringAngle)
    if(.not. isValid(phaseFunctionVar)) &
      call setStateToFailure(status, &
                             "getPhaseFunctionValues: Phase function variable has not been initialized.") 
    if(any(scatteringAngle(:) < minScatteringAngle) .or. &
       any(scatteringAngle(:) > maxScatteringAngle))     &
      call setStateToFailure(status, "getPhaseFunctionValues: ScatteringAngle out of bounds.")
    if(size(scatteringAngle) /= size(value)) &
      call setStateToFailure(status,         &
            "getPhaseFunctionValues: Number of scattering angles and phase function values must match.")
    
    ! ------------------------
    if(.not. stateIsFailure(status)) then 
      if(storedAsLegendre(phaseFunctionVar)) then
        !
        ! Phase function is stored as Legendre moments        
        !    We expand the values using the Legendre coefficients
        !
        maxL = size(phaseFunctionVar%legendreCoefficients)
        if(maxL == 0) then 
          ! 
          ! This is the phase function for isotropic scattering, which has only 
          !   P0 = 1. 
          !
          value(:) = 1/2. 
        else
          allocate(legendreP(0:maxL, nValues))
          legendreP(:, :) = computeLegendrePolynomials(maxL, cos(scatteringAngle(:)))
          value(:) = matmul((/ 1., phaseFunctionVar%legendreCoefficients(:) /) * (/ (2*l + 1, l = 0, maxL) /), &
                            legendreP(:,:))
          deallocate(legendreP)
        end if
      else
        !
        ! Phase function is stored as table     
        !    Return values interpolated in cosine of the scattering angle from the points in the original table
        !
        nStored = size(phaseFunctionVar%scatteringAngle)
        allocate(tableIndicies(nValues), weights(nValues), dMu(nValues), indiciesPlus1(nValues))
        do l = 1, nValues
          tableIndicies(l) = findIndex(scatteringAngle(l), phaseFunctionVar%scatteringAngle)
        end do
        indiciesPlus1(:) = tableIndicies(:) + 1

        ! There is the possibility that the some entry requested in scatteringAngle
        !   will be exactly Pi radians, which would make that entry in tableIndicies
        !   be nStored, and we don't want to refer to array elements beyond the table. 
        !
        where(tableIndicies(:) < nStored)
          dMu(:) = cos(phaseFunctionVar%scatteringAngle(indiciesPlus1(:))) - &
                   cos(phaseFunctionVar%scatteringAngle(tableIndicies(:)))
        elsewhere ! Weight is 0 when dAngle is large (tablesIndices >= nStored)
          dMu(:) = huge(dMu)
          indiciesPlus1(:) = tableIndicies(:) ! Don't want to try to access beyond the end of the table
        end where
        weights(:) = 1. - (cos(scatteringAngle(:)) -                                  &
                           cos(phaseFunctionVar%scatteringAngle(tableIndicies(:)))) / &
                          dMu(:)
        value(:) =       weights(:) *  phaseFunctionVar%value(tableIndicies(:)) + &
                   (1. - weights(:)) * phaseFunctionVar%value(indiciesPlus1(:)) 
        deallocate(tableIndicies, weights, dMu, indiciesPlus1)
      end if
      call setStateToSuccess(status)
    end if
  end subroutine getPhaseFunctionValues_one
  !------------------------------------------------------------------------------------------
  subroutine getPhaseFunctionValues_table(table, scatteringAngle, values, status)
    type(phaseFunctionTable), intent(in   ) :: table
    real, dimension(:),       intent(in   ) :: scatteringAngle
    real, dimension(:, :),    intent(  out) :: values
    type(ErrorMessage),       intent(inout) :: status
    ! ------------------------
    ! Return the values of the phase functions in a table given a set of scattering angles.
    !   The output array (values) has dimensions numScatteringAngles by numPhaseFunctions
    !
    
    ! Local variables
    integer                               :: i, l, maxL, nValues, nStored
    integer, dimension(:),    allocatable :: tableIndicies, indiciesPlus1
    real,    dimension(:),    allocatable :: weights, dMu
    real,    dimension(:, :), allocatable :: legendreP
    
    ! ------------------------
    ! Sanity checks: 
    !   phase function contains values; angles in bounds; angles in order; arrays same length
    nValues = size(scatteringAngle)
    if(.not. isReady_phaseFunctionTable(table)) &
      call setStateToFailure(status, &
                             "getPhaseFunctionValues: Phase function table has not been initialized.") 
    if(any(scatteringAngle(:) < minScatteringAngle) .or. &
       any(scatteringAngle(:) > maxScatteringAngle))     &
      call setStateToFailure(status, "getPhaseFunctionValues: ScatteringAngle out of bounds.")
    if(any(scatteringAngle(2:) - scatteringAngle(:nValues-1) <= 0.)) &
      call setStateToFailure(status, "getPhaseFunctionValues: Scattering angle must be increasing, unique.")
    if(size(scatteringAngle)                   /= size(values, 1) .or. &
       size(table%phaseFunctions) /= size(values, 2)) &
      call setStateToFailure(status,         &
            "getPhaseFunctionValues: Number of scattering angles and phase function values must match.")
    
    ! ------------------------
    if(.not. stateIsFailure(status)) then
    
      if(any(storedAsLegendre(table%phaseFunctions))) then
        !
        ! What's the largest number of Legendre coefficients in any of the phase functions? 
        !
        maxL = 0
        do i = 1, size(table%phaseFunctions)
          if(storedAsLegendre(table%phaseFunctions(i))) &
            maxL = max(maxL, size(table%phaseFunctions(i)%legendreCoefficients))
        end do 
        !
        ! Compute all the Legendre terms we'll ever need 
        !   Include the 2 L + 1 term so we don't recalculate it for every phase function
        !
        if(maxL > 0) then
          allocate(legendreP(0:maxL, nValues))
          legendreP(:, :) = spread( (/ (2*l + 1, l = 0, maxL) /), dim = 2, ncopies = nValues) * &
                            computeLegendrePolynomials(maxL, cos(scatteringAngle(:)))
        end if
      else if(table%oneAngleSet) then
        !
        ! Precompute mapping from native phase function angles to the angles requested
        !
        nStored = size(table%phaseFunctions(1)%scatteringAngle(:))
        allocate(tableIndicies(nValues), weights(nValues), dMu(nValues), indiciesPlus1(nValues))
        tableIndicies(1) = findIndex(scatteringAngle(1), &
                                     table%phaseFunctions(1)%scatteringAngle(:))
        do l = 2, nValues
          tableIndicies(l) = findIndex(scatteringAngle(l),                         &
                                       table%phaseFunctions(1)%scatteringAngle(:), &
                                       firstGuess = tableIndicies(l-1))
        end do
        indiciesPlus1(:) = tableIndicies(:) + 1

        ! There is the possibility that the last entry requested in scatteringAngle
        !   will be exactly Pi radians, which would make that entry in tableIndicies
        !   be nStored, and we don't want to refer to array elements beyond the table. 
        !
        where(tableIndicies(:) < nStored)
          dMu(:) = cos(table%phaseFunctions(1)%scatteringAngle(indiciesPlus1(:))) - &
                   cos(table%phaseFunctions(1)%scatteringAngle(tableIndicies(:)))
        elsewhere ! Weight is 0 when dAngle is large (tablesIndices >= nStored)
          dMu(:) = huge(dMu)
          indiciesPlus1(:) = tableIndicies(:) ! Don't want to try to access beyond the end of the table
        end where
        weights(:) = 1. - (cos(scatteringAngle(:)) -                                         &
                           cos(table%phaseFunctions(1)%scatteringAngle(tableIndicies(:)))) / &
                          dMu(:)
      end if
      
      do i = 1, size(table%phaseFunctions)
        if(storedAsLegendre(table%phaseFunctions(i))) then
          !
          ! Use previously computed Legendre polynomials
          !
          maxL = size(table%phaseFunctions(i)%legendreCoefficients)
          if(maxL == 0) then 
            values(:, i) = 1/2. 
          else
            values(:, i) = matmul((/ 1., table%phaseFunctions(i)%legendreCoefficients(:) /), &
                                  legendreP(0:maxL,:))
          end if
        else if(table%oneAngleSet) then 
          !
          ! Use previously computed map from native angles to requested angles
          !
          values(:, i) =       weights(:)  * table%phaseFunctions(i)%value(tableIndicies(:)) + &
                         (1. - weights(:)) * table%phaseFunctions(i)%value(indiciesPlus1(:)) 
        else 
          !
          ! General case (phase functions stored as values at different angles sets) 
          !
          call getPhaseFunctionValues_one(table%phaseFunctions(i), &
                                          scatteringAngle(:), values(:, i), status)
          if(stateIsFailure(status)) exit
        end if 
      end do 

      if(allocated(legendreP))     deallocate(legendreP)
      if(allocated(tableIndicies)) deallocate(tableIndicies, indiciesPlus1, weights, dMu) 
      if(.not. stateIsFailure(status)) call setStateToSuccess(status)
    end if
  end subroutine getPhaseFunctionValues_table
  !------------------------------------------------------------------------------------------
  subroutine getPhaseFunctionCoefficients(phaseFunctionVar, legendreCoefficients, status)
    type(phaseFunction), intent(in   ) :: phaseFunctionVar
    real, dimension(:),  intent(  out) :: legendreCoefficients
    type(ErrorMessage),  intent(inout) :: status
    ! ------------------------
    ! Return the legendre moments of the phase function, starting with P1 
    !   (because P0 is identically 1).  
    
    ! Local variables
    real, dimension(:, :), allocatable :: legendreP
    real, dimension(:),    allocatable :: weights, mus, values
    integer                            :: maxL, nAngles
    ! How many quadrature points should we use when computing Legendre
    !   series from tabulated phase functions? This number is completely ad-hoc
    integer, parameter                 :: resolutionMuliplier = 2
    
    ! ------------------------
    ! Sanity checks: 
    !   phase function contains values; array for coefficients is big enough
    call setStateToSuccess(status)
    if(.not. isValid(phaseFunctionVar)) &
      call setStateToFailure(status, &
                             "getPhaseFunctionCoefficients: Phase function variable has not been initialized.") 
      
    ! ------------------------
    if(.not. stateIsFailure(status)) then 
      if(storedAsLegendre(phaseFunctionVar)) then
        if(size(legendreCoefficients) < size(phaseFunctionVar%legendreCoefficients)) then
          legendreCoefficients(:) = phaseFunctionVar%legendreCoefficients(:size(legendreCoefficients))
        else
          legendreCoefficients(:size(phaseFunctionVar%legendreCoefficients)) = &
            phaseFunctionVar%legendreCoefficients(:)
          legendreCoefficients(size(phaseFunctionVar%legendreCoefficients) + 1:) = 0
        end if
      else ! phase function is stored as angle/value pairs
        ! The Legendre moments are computed as the projection of the phase function on Legendre polynomials. 
        !   An alternative would be to interpolate the phase function to Gaussian quadrature points
        maxL    = size(legendreCoefficients)
        nAngles = size(phaseFunctionVar%value)
        allocate(legendreP(0:maxL, resolutionMuliplier * size(phaseFunctionVar%scatteringAngle)), & 
                 weights(resolutionMuliplier * size(phaseFunctionVar%scatteringAngle)),           &
                 mus    (resolutionMuliplier * size(phaseFunctionVar%scatteringAngle)),           &
                 values (resolutionMuliplier * size(phaseFunctionVar%scatteringAngle)) )
        call computeLobattoTerms(mus, weights); mus = mus(size(mus):1:-1) 
        legendreP(:, :) = computeLegendrePolynomials(maxL, mus)
        call getPhaseFunctionValues(phaseFunctionVar, acos(mus), values, status)
        
        legendreCoefficients(:) = 0.5 * matmul(legendreP(1:, :), values(:))
        !
        ! We're integrating the product of each Legendre polynomial and the phase function over mu 
        ! The first factor of 1/2 is for the projection (see, eg., Thomas and Stamnes, eq 6.29),
        ! the second is for the numerical integration; the minus sign because we're doing the integral 
        ! ovr mu from -1 to 1
        ! We ignore legendreP(0:, :) because don't need to compute P0 
        !
!           legendreCoefficients(:) =                                                                       &
!               -0.5 * matmul(0.5 * (legendreP(1:, 2:)         * spread(phaseFunctionVar%value(2:),         &
!                                                                       dim = 1, ncopies = maxL)   +        &
!                                    legendreP(1:, :nAngles-1) * spread(phaseFunctionVar%value(:nAngles-1), &
!                                                                       dim = 1, ncopies = maxL) ),         &
!                             cos(phaseFunctionVar%scatteringAngle(2:)) - cos(phaseFunctionVar%scatteringAngle(:nAngles-1)))


        deallocate(legendreP, weights, mus, values)
      end if
    call setStateToSuccess(status)
    end if
  end subroutine getPhaseFunctionCoefficients
  !------------------------------------------------------------------------------------------
  elemental function getExtinction(phaseFunctionVar)
    type(phaseFunction), intent( in) :: phaseFunctionVar
    real                             :: getExtinction
    
    getExtinction = phaseFunctionVar%extinction
  end function getExtinction
  !------------------------------------------------------------------------------------------
  elemental function getSingleScatteringAlbedo(phaseFunctionVar)
    type(phaseFunction), intent( in) :: phaseFunctionVar
    real                             :: getSingleScatteringAlbedo
    
    getSingleScatteringAlbedo = phaseFunctionVar%singleScatteringAlbedo
  end function getSingleScatteringAlbedo  
  !------------------------------------------------------------------------------------------
  elemental function isReady_phaseFunction(testPhaseFunction)
    type(phaseFunction), intent( in) :: testPhaseFunction
    logical                          :: isReady_phaseFunction
    
    ! Has this phase function been filled in? 
    !   All the other checks should have been done when the 
    !   variable is initialized
    isReady_phaseFunction = associated(testPhaseFunction%value) .or. &
                            associated(testPhaseFunction%legendreCoefficients)

  end function isReady_phaseFunction
  !------------------------------------------------------------------------------------------
  elemental function isReady_phaseFunctionTable(table)
    type(phaseFunctionTable), intent( in) :: table
    logical                               :: isReady_phaseFunctionTable
    !
    ! A phase function table is ready to use if there is room for one or more phase functions
    !   and if each of those phase functions has been initialized. 
    !
    integer :: i
    
    isReady_phaseFunctionTable = associated(table%phaseFunctions)
    if(isReady_phaseFunctionTable) then
      do i = 1, size(table%phaseFunctions)
        isReady_phaseFunctionTable = isReady_phaseFunctionTable .and. &
                                     isReady_phaseFunction(table%phaseFunctions(i))
      end do
    end if

  end function isReady_phaseFunctionTable
  !------------------------------------------------------------------------------------------
  subroutine getInfo_PhaseFunction(phaseFunctionVar, nCoefficients, nAngles, nativeAngles, status)
    type(phaseFunction),          intent(in   ) :: phaseFunctionVar
    integer,            optional, intent(  out) :: nCoefficients, nAngles
    real, dimension(:), optional, intent(  out) :: nativeAngles
    type(ErrorMessage), optional, intent(inout) :: status
    
    integer :: nNativeAngles, nScatteringAngles

    ! ------------------------
    ! Sanity checks
    if(.not. isValid(phaseFunctionVar)) &
      call setStateToFailure(status, "getInfo_PhaseFunction: phase function hasn't been initialized.")
    
    ! ------------------------
    if(.not. stateIsFailure(status)) then
      if(storedAsLegendre(phaseFunctionVar)) then
        if(present(nCoefficients)) nCoefficients   = size(phaseFunctionVar%legendreCoefficients) 
        if(present(nAngles))       nAngles         = 0
        if(present(nativeAngles))  nativeAngles(:) = 0.
        call setStateToSuccess(status)
      else ! must be stored as angle-value pairs
        if(present(nCoefficients))  nCoefficients  = 0
        if(present(nAngles))        nAngles        = size(phaseFunctionVar%scatteringAngle)
        if(present(nativeAngles)) then
          nNativeAngles     = size(nativeAngles)
          nScatteringAngles = size(phaseFunctionVar%scatteringAngle)
          if(nNativeAngles == nScatteringAngles) then
            nativeAngles(:) = phaseFunctionVar%scatteringAngle(:)
            call setStateToSuccess(status)
          else
            nativeAngles(:min(nNativeAngles, nScatteringAngles)) = &
              phaseFunctionVar%scatteringAngle(:min(nNativeAngles, nScatteringAngles))
            if(nNativeAngles > nScatteringAngles) nativeAngles(nNativeAngles + 1: ) = 0. 
            call setStateToWarning(status, &
                                   "getInfo_PhaseFunction: wrong number of native angles requested")
          end if
        else 
          call setStateToSuccess(status)
        end if
        
      end if
    end if
        
  end subroutine getInfo_PhaseFunction  
  !------------------------------------------------------------------------------------------
  subroutine getInfo_PhaseFunctionTable(table, nEntries, key, extinction, singleScatteringAlbedo, &
                                        phaseFunctionDescriptions, tableDescription, status)
    type(phaseFunctionTable),           intent(in   ) :: table
    integer,                  optional, intent(  out) :: nEntries
    real,    dimension(:),    optional, intent(  out) :: key, extinction, singleScatteringAlbedo
    character(len = *), &
              dimension(:),   optional, intent(  out) :: phaseFunctionDescriptions
    character(len = *),       optional, intent(  out) :: tableDescription
    type(ErrorMessage),                 intent(inout) :: status
    
    !
    ! Local variables
    !
    integer             :: i 
    ! ------------------------
    if(.not. associated(table%phaseFunctions)) &
      call setStateToFailure(status, "getKeyInfo: phase function table hasn't been initialized.")
    
    ! ------------------------
    if(.not. stateIsFailure(status)) then 
      if(present(nEntries)) nEntries = size(table%key)
      
      if(present(key)) then
        if(size(key) /= size(table%key)) &
          call setStateToWarning(status, "getKeyInfo: supplied key length not the same as table size.") 
        key(:) = 0
        key(:min(size(key), size(table%key))) = table%key(:min(size(key), size(table%key)))
      end if
      
      if(present(extinction)) then
        if(size(extinction) /= size(table%phaseFunctions)) &
          call setStateToWarning(status, "getKeyInfo: vector extinction isn't the same size as the table.") 
        forall(i = 1:max(size(extinction), size(table%phaseFunctions))) &
          extinction(i) = getExtinction(table%phaseFunctions(i))
        extinction(size(table%phaseFunctions)+1:) = defaultExtinction 
      end if
      
      if(present(singleScatteringAlbedo)) then
        if(size(singleScatteringAlbedo) /= size(table%phaseFunctions)) &
          call setStateToWarning(status, "getKeyInfo: vector singleScatteringAlbedo isn't the same size as the table.") 
        forall(i = 1:max(size(singleScatteringAlbedo), size(table%phaseFunctions))) &
          singleScatteringAlbedo(i) = getSingleScatteringAlbedo(table%phaseFunctions(i))
        singleScatteringAlbedo(size(table%phaseFunctions)+1:) = defaultSSA 
      end if
      
      if(present(phaseFunctionDescriptions)) then
        if(associated(table%phaseFunctionDescriptions)) then 
          if(size(phaseFunctionDescriptions) /= size(table%phaseFunctionDescriptions)) &
            call setStateToWarning(status, "getKeyInfo: number of phaseFunctionDescriptions not the same as table size.") 
          phaseFunctionDescriptions(:) = ""
          phaseFunctionDescriptions(:min(size(phaseFunctionDescriptions), size(table%phaseFunctionDescriptions))) = &
            table%phaseFunctionDescriptions(:min(size(phaseFunctionDescriptions), size(table%phaseFunctionDescriptions)))
        else 
          phaseFunctionDescriptions(:) = ""
        end if
      end if
      
      if(present(tableDescription)) tableDescription = table%description 
      
      call setStateToSuccess(status)
    end if
     
  end subroutine getInfo_PhaseFunctionTable
  !------------------------------------------------------------------------------------------
  function getElement(n, table, status) 
    integer,                  intent(in   ) :: n
    type(phaseFunctionTable), intent(in   ) :: table
    type(ErrorMessage),       intent(inout) :: status
    type(phaseFunction)                     :: getElement
    ! Returns element n from a table of scattering phase functions
    !   Users should ensure that the variable is properly finalized 
    
    ! ------------------------
    if(.not. associated(table%phaseFunctions)) &
      call setStateToFailure(status, "getElement: phase function table hasn't been initialized.")
    if(.not. stateIsFailure(status)) then
      if(n < 1 .or. n > size(table%phaseFunctions)) &
        call setStateToFailure(status, "getElement: asking for non-existent element.")
    end if 

    ! ------------------------
    if(.not. stateIsFailure(status)) then
      getElement = copy_PhaseFunction(table%phaseFunctions(n))
      call setStateToSuccess(status)
    end if 

  end function getElement
  !------------------------------------------------------------------------------------------
  ! Storage and retrieval
  !------------------------------------------------------------------------------------------
  subroutine write_PhaseFunctionTable(table, fileName, status)
    use netcdf
    type(phaseFunctionTable), intent(in   ) :: table
    character(len = *),       intent(in   ) :: fileName
    type(ErrorMessage),       intent(inout) :: status
    ! Write a phase function table to an external file. 
    !   Here we use netCDF files, and only make provisions for a 
    !   set of tablulated phase functions that share a common angle set, or
    !   a set of phase functions stored as Legendre polynomials
    
    integer :: ncStatus, ncFileId
    
    ncStatus = nf90_create(trim(fileName), nf90_clobber, ncFileID)
    if(ncStatus == nf90_NoErr) then 
      ncStatus = nf90_EndDef(ncFileId)
      call add_PhaseFunctionTable(table, ncFileId, status = status)
      ncStatus = nf90_close(ncFileId)
      if(ncStatus /= nf90_NoErr) then
        call setStateToFailure(status, "write_PhaseFunctionTable: " // trim(nf90_StrError(ncStatus)))
        ncstatus = nf90_close(ncFileId) 
        open(20, file = trim(fileName))
        close(20, status = "delete")
      end if   
    else
      call setStateToFailure(status, "write_PhaseFunctionTable: can't create file " // trim(fileName)) 
    end if 
    if(.not. stateIsFailure(status)) call setStateToSuccess(status)
  end subroutine write_PhaseFunctionTable
  !------------------------------------------------------------------------------------------
  subroutine add_PhaseFunctionTable(table, fileId, prefix, status)
    use netcdf
    type(phaseFunctionTable),   intent(in   ) :: table
    integer,                    intent(in   ) :: fileId
    character(len=*), optional, intent(in   ) :: prefix
    type(ErrorMessage),         intent(inout) :: status
    ! Add a phase function table to an external file.
    !   Here we use netCDF files, and only make provisions for a 
    !   set of tablulated phase functions that share a common angle set, or
    !   a set of phase functions stored as Legendre polynomials
    
    ! fileId should point to an open file; in the netcdf version this file should
    !   not be in define mode
    
    
    ! Local variables
    integer, dimension(16)                :: ncStatus
    integer                               :: ncFileId, entryDimID, keyVarId, &
                                             extinctionVarId, ssaVarId
    integer                               :: nEntries, i
    real,    dimension(:),    allocatable :: extinction, singleScatteringAlbedo
    ! For angle-value pairs
    integer                               :: angleDimID, angleVarId, phaseFunctionVarID 
    integer                               :: nAngles
    real,    dimension(:, :), allocatable :: phaseFunctionArray
    ! For sets of Legendre coefficients
    integer                               :: coefficientDimID, coefficientVarID, &
                                             startVarId, lengthVarId
    integer, dimension(:),    allocatable :: start, length
    real,    dimension(:),    allocatable :: legendreCoefficients
    
    ! To learn about the existing file
    integer                               :: numDimensions
    character(len = 16)                   :: thisPrefix
    

    ! ------------------------
    ncStatus(:) = nf90_NoErr
    thisPrefix = ""
    if(present(prefix)) thisPrefix = prefix
    ncFileId = fileId
    if(.not. isReady_phaseFunctionTable(table)) then 
      call setStateToFailure(status, "add_PhaseFunctionTable: phase function table hasn't been initialized.") 
    else if (.not.(table%oneAngleSet .or. &
                   all(storedAsLegendre(table%phaseFunctions(:))))) then
      !
      ! We're not including a general way to write a table of phase functions. 
      !   (It seems like the two supplied here should work for most people.) 
      !
      call setStateToFailure(status, "add_PhaseFunctionTable: Can't write general phase function tables to files.")
    else
      ncStatus( :) = nf90_NoErr
      !
      ! Is the file being pointed to a reasonable netcdf file?  
      !
      ncstatus( 1) = nf90_inquire(ncFileId, numDimensions)
      if(ncStatus(1) /= nf90_noErr) &
        call setStateToFailure(status, "add_PhaseFunctionTable: trying to add to the wrong kind of file")
      !
      ! Are the variables we want to define already in place? That would be bad.
      !   We're only going to check one dimension, and assume that all the other
      !   dimensions and variables are also fine
      !
      ncStatus( 1) = nf90_inq_dimid(ncFileId, trim(thisPrefix) // "phaseFunctionNumber", entryDimID)
      if(ncStatus(1) == nf90_noErr) &
        call setStateToFailure(status, "add_PhaseFunctionTable: trying to add to the wrong kind of file") 
    end if
    
    ! ------------------------
    if(.not. stateIsFailure(status)) then 
      nEntries = size(table%phaseFunctions)
      
      ncStatus( :) = nf90_NoErr
      ncStatus( 1) = nf90_redef(ncFileID)
      ncStatus( 2) = nf90_def_dim(ncFileId, trim(thisPrefix) // "phaseFunctionNumber", nEntries, entryDimId)
      ncStatus( 3) = nf90_def_var(ncFileId, trim(thisPrefix) // "phaseFunctionKeyT",       &
                                  nf90_float, entryDimId, keyVarID)
      ncStatus( 4) = nf90_def_var(ncFileId, trim(thisPrefix) // "extinctionT",             &
                                  nf90_float, entryDimId, extinctionVarId)
      ncStatus( 5) = nf90_def_var(ncFileId, trim(thisPrefix) // "singleScatteringAlbedoT", &
                                  nf90_float, entryDimId, ssaVarId)
      if(len_trim(table%description) > 0) &
        ncStatus( 6) = nf90_put_att(ncFileID, nf90_Global, &
                                              trim(thisPrefix) // "description", trim(table%description))
     
      if(table%oneAngleSet) then
        ! Table is stored as one set of angles and phase function values at each angle.
        !   We'll write the angles and keys as vectors, and the values as a matrix
        nAngles = size(table%phaseFunctions(1)%scatteringAngle)
        
        ! Specific variable definitions
        ncStatus( 7) = nf90_def_dim(ncFileId, trim(thisPrefix) // "scatteringAngle",     &
                                    nAngles,  angleDimID)
        ncStatus( 8) = nf90_def_var(ncFileId, trim(thisPrefix) // "scatteringAngle",     &
                                    nf90_float, angleDimId, angleVarID)
        ncStatus( 9) = nf90_def_var(ncFileId, trim(thisPrefix) // "phaseFunctionValues", &
                                    nf90_float, (/ angleDimID, entryDimID /), phaseFunctionVarID)
        ncStatus(10) = nf90_put_att(ncFileID, nf90_Global, trim(thisPrefix) // "phaseFunctionStorageType", &
                                    "Angle-Value")
      else 
        ! All the phase functions are stored as Legendre polynomials. We'll store the 
        !   coefficients in a single vector, keeping track of the starting position and length
        !   of each set of coefficients in two separate vectors (yes, this is redundent). 
        allocate(start(nEntries), length(nEntries))

        ! What is the length of each Legendre series? How many terms are there in total? 
        do i = 1, nEntries
          length(i) = size(table%phaseFunctions(i)%legendreCoefficients)
        end do
        start(1) = 1
        do i = 2, nEntries
          start(i) = start(i-1) + length(i-1)
        end do
        
        ! Specific variable definitions
        ncStatus( 6) = nf90_def_dim(ncFileId, trim(thisPrefix) // "coefficents",          &       
                                    start(nEntries) + length(nEntries) - 1, coefficientDimID)
        ncStatus( 7) = nf90_def_var(ncFileId, trim(thisPrefix) // "start",                &
                                    nf90_int,   entryDimId,  startVarId)
        ncStatus( 8) = nf90_def_var(ncFileId, trim(thisPrefix) // "length",               &
                                    nf90_int,   entryDimId, lengthVarId)
        ncStatus( 9) = nf90_def_var(ncFileId, trim(thisPrefix) // "legendreCoefficients", &
                                    nf90_float, coefficientDimID, coefficientVarID)
        ncStatus(10) = nf90_put_att(ncFileID, nf90_Global, &
                                              trim(thisPrefix) // "phaseFunctionStorageType", "LegendreCoefficients")
      end if
      ncStatus(11) = nf90_EndDef(ncFileID)
      
      if(any(ncStatus(:) /= nf90_noErr)) then
        do i = 1, size(ncStatus)
          if(ncStatus(i) /= nf90_NoErr) &
            call setStateToFailure(status, "add_PhaseFunctionTable: " // trim(nf90_StrError(ncStatus(i))))
        end do
      end if
    end if ! End of file definitions
      
    ! Write the variables
    if(.not. stateIsFailure(status)) then ! File writing 
      ! These variables are common to both formats
      allocate(extinction(nEntries), singleScatteringAlbedo(nEntries))
      extinction             = table%phaseFunctions(:)%extinction
      singleScatteringAlbedo = table%phaseFunctions(:)%singleScatteringAlbedo
      ncStatus( 1) = nf90_put_var(ncFileId,        keyVarId, table%key) 
      ncStatus( 2) = nf90_put_var(ncFileID, extinctionVarId, extinction)
      ncStatus( 3) = nf90_put_var(ncFileID,        ssaVarId, singleScatteringAlbedo)
      deallocate(extinction, singleScatteringAlbedo)
      
      if(table%oneAngleSet) then
        !
        ! Copying the phase function values into temporary, local memory is lots 
        !   faster than writing each phase function in turn. 
        !
        allocate(phaseFunctionArray(nAngles, nEntries))
        do i = 1, nEntries
          phaseFunctionArray(:, i) = table%phaseFunctions(i)%value(:)
        end do
        
        ncStatus( 4) = nf90_put_var(ncFileId,         angleVarId, table%phaseFunctions(1)%scatteringAngle) 
        ncStatus( 5) = nf90_put_var(ncFileId, phaseFunctionVarId, phaseFunctionArray) 
        deallocate(phaseFunctionArray)
      else 
        !
        ! Copy all the Legendre coefficients into one long array
        !
        allocate(legendreCoefficients(start(nEntries) + length(nEntries) - 1))
        forall( i = 1:nEntries)
          legendreCoefficients(start(i):(start(i) + length(i) - 1)) = &
            table%phaseFunctions(i)%legendreCoefficients(:)
        end forall
        
        ncStatus( 4) = nf90_put_var(ncFileId,       startVarId, start) 
        ncStatus( 5) = nf90_put_var(ncFileId,      lengthVarId, length) 
        ncStatus( 6) = nf90_put_var(ncFileId, coefficientVarID, legendreCoefficients) 
        deallocate(start, length, legendreCoefficients)
      end if
    end if ! File writing is done
    
    ! Did everything go OK? 
    if(any(ncStatus(:) /= nf90_noErr)) then
      do i = 1, size(ncStatus)
        if(ncStatus(i) /= nf90_NoErr) &
          call setStateToFailure(status, "add_PhaseFunctionTable: " // trim(nf90_StrError(ncStatus(i))))
      end do
    end if
    if(.not. stateIsFailure(status)) call setStateToSuccess(status)

  end subroutine add_PhaseFunctionTable
  !------------------------------------------------------------------------------------------
  subroutine read_PhaseFunctionTable(fileName, fileId, table, prefix, status)
    use netcdf
    character(len = *), optional, intent(in   ) :: fileName
    integer,            optional, intent(in   ) :: fileId
    type(phaseFunctionTable),     intent(  out) :: table
    character(len = *), optional, intent(in   ) :: prefix
    type(ErrorMessage),           intent(inout) :: status
    ! Read a phase function table from an external file. 
    !   Here we use netCDF files. 
    
    ! Local variables
    integer, dimension(24)                     :: ncStatus
    integer                                    :: ncFileId, ncDimId, ncVarId
    integer                                    :: nEntries, i
    character(len = 32)                        :: storageType
    real,    dimension(:),    allocatable      :: key, extinction, singleScatteringAlbedo
    character(len = maxTableDescriptorLength)  :: description
    ! For angle-value pairs
    integer                               :: nAngles
    real,    dimension(:),    allocatable :: scatteringAngle
    real,    dimension(:, :), allocatable :: phaseFunctionArray
    ! For sets of Legendre coefficients
    integer                               :: nTotalCoefficients
    integer, dimension(:),    allocatable :: start, length
    real,    dimension(:),    allocatable :: legendreCoefficients
    type(phaseFunction), &
             dimension(:),    allocatable :: phaseFunctions
    
    character(len = 16)                   :: thisPrefix
             
    ! ------------------------
    description = "" 
    thisPrefix = ""
    if(present(prefix)) thisPrefix = trim(prefix)
    if(present(fileName) .eqv. present(fileId)) then 
      call setStateToFailure(status, "read_PhaseFunctionTable: must supply either fileName or fileId")
    else if (present(fileName)) then   
      ! Open the file, see if it is of the correct format
      ncStatus( :) = nf90_NoErr
      ncStatus( 1) = nf90_open(trim(fileName), nf90_NoWrite, ncFileID)
      ncStatus( 2) = nf90_get_att(ncFileId, nf90_Global, trim(thisPrefix) // "phaseFunctionStorageType", storageType)
      if(any(ncStatus(1:2) /= nf90_noErr)) & 
        call setStateToFailure(status, "read_PhaseFunctionTable: " // trim(fileName) // " is not a phase function table file.")
    else
      ncFileId = fileId
      ncStatus( 1) = nf90_get_att(ncFileId, nf90_Global, trim(thisPrefix) // "phaseFunctionStorageType", storageType)
      if(ncStatus(  1) /= nf90_noErr) & 
        call setStateToFailure(status, "read_PhaseFunctionTable: " // trim(fileName) // " doesn't contain this phase function.")
    end if
    
    ! Read the data from the file
    ncStatus(:) = nf90_NoErr
    if(.not. stateIsFailure(status)) then
      ! What's common to all the possible file formats? 
      ncStatus( 1) = nf90_inq_dimid(ncFileId, trim(thisPrefix) // "phaseFunctionNumber", ncDimId)
      ncStatus( 2) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nEntries)
      allocate(key(nEntries), extinction(nEntries), singleScatteringAlbedo(nEntries))
      ncStatus( 3) = nf90_inq_varid(ncFileId, trim(thisPrefix) // "phaseFunctionKeyT", ncVarId)
      ncStatus( 4) = nf90_get_var(ncFileId, ncVarId, key)
      ncStatus( 5) = nf90_inq_varid(ncFileId, trim(thisPrefix) // "extinctionT", ncVarId)
      ncStatus( 6) = nf90_get_var(ncFileId, ncVarId, extinction)
      ncStatus( 7) = nf90_inq_varid(ncFileId, trim(thisPrefix) // "singleScatteringAlbedoT", ncVarId)
      ncStatus( 8) = nf90_get_var(ncFileId, ncVarId, singleScatteringAlbedo)
      ncStatus( 9) = nf90_get_att(ncFileId, nf90_Global, trim(thisPrefix) // "description", description)
      if (ncStatus( 9) /= nf90_NoErr) then ! If a zero-length description was provided the attribute wasn't even written
        description = ""
        ncStatus( 9) = nf90_NoErr
      end if
      
      if(index(trim(storageType), "Angle-Value") == 1) then
        ncStatus(10) = nf90_inq_dimid(ncFileId, trim(thisPrefix) // "scatteringAngle", ncDimId)
        ncStatus(11) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nAngles)
        allocate(scatteringAngle(nAngles), phaseFunctionArray(nAngles, nEntries))
       
        ncStatus(12) = nf90_inq_varid(ncFileId, trim(thisPrefix) //  "scatteringAngle", ncVarId)
        ncStatus(13) = nf90_get_var(ncFileId, ncVarId, scatteringAngle)
        ncStatus(14) = nf90_inq_varid(ncFileId, trim(thisPrefix) //  "phaseFunctionValues", ncVarId)
        ncStatus(15) = nf90_get_var(ncFileId, ncVarId, phaseFunctionArray)
      else if (index(trim(storageType), "LegendreCoefficients") == 1) then
        ncStatus(10) = nf90_inq_dimid(ncFileId, trim(thisPrefix) //  "coefficents", ncDimId)
        ncStatus(11) = nf90_Inquire_Dimension(ncFileId, ncDimId, len = nTotalCoefficients)
        allocate(start(nEntries), length(nEntries), legendreCoefficients(nTotalCoefficients), &
                 phaseFunctions(nEntries))
       
        ncStatus(12) = nf90_inq_varid(ncFileId, trim(thisPrefix) // "start", ncVarId)
        ncStatus(13) = nf90_get_var(ncFileId, ncVarId, start)
        ncStatus(14) = nf90_inq_varid(ncFileId, trim(thisPrefix) // "length", ncVarId)
        ncStatus(15) = nf90_get_var(ncFileId, ncVarId, length)
        ncStatus(16) = nf90_inq_varid(ncFileId, trim(thisPrefix) // "legendreCoefficients", ncVarId)
        ncStatus(17) = nf90_get_var(ncFileId, ncVarId, legendreCoefficients)
      else
        call setStateToFailure(status, "read_PhaseFunctionTable: " // trim(fileName) // " is of unknown format.")
      end if

      !
      ! If the name was supplied we're writing a stand-alone file, so it's time to close it
      !
      if(present(fileName)) ncStatus(18) = nf90_close(ncfileId)

      !
      ! Report any netcdf errors in the status variable
      !
      if(any(ncStatus(:) /= nf90_noErr)) then
        do i = 1, size(ncStatus)
          if(ncStatus(i) /= nf90_NoErr) &
            call setStateToFailure(status, "read_PhaseFunctionTable: " // trim(nf90_StrError(ncStatus(i))))
        end do
      end if 
      
      !
      ! Put the data into the phase function table, depending on how it was stored
      !
      if(.not. stateIsFailure(status)) then
        if(index(trim(storageType), "Angle-Value") == 1) then
          table = new_PhaseFunctionTable(scatteringAngle, phaseFunctionArray,     &
                                        key, extinction, singleScatteringAlbedo, &
                                        tableDescription = description, status = status)
          deallocate(scatteringAngle, phaseFunctionArray)
        else 
          do i = 1, nEntries
            phaseFunctions(i) = new_PhaseFunction(legendreCoefficients(start(i):(start(i) + length(i) - 1)), &
                                                  extinction(i), singleScatteringAlbedo(i), status = status)
          end do
          if(.not. StateIsFailure(status)) &
            table = new_PhaseFunctionTable(phaseFunctions, key, tableDescription = description, status = status)
         
          ! Each element in the array of phase function needs to give its memory up
          do i = 1, nEntries
            call finalize_PhaseFunction(phaseFunctions(i))
          end do
          deallocate(start, length, legendreCoefficients, phaseFunctions)
        end if  ! Which kind of phase function table
      end if
      deallocate(key, extinction, singleScatteringAlbedo)
      if(.not. stateIsFailure(status)) call setStateToSuccess(status)
    end if
  end subroutine read_PhaseFunctionTable
  !------------------------------------------------------------------------------------------
  ! Finalization
  !------------------------------------------------------------------------------------------
  subroutine finalize_PhaseFunction(phaseFunctionVar)
    type(phaseFunction), intent(inout) :: phaseFunctionVar
    ! Free memory and nullify pointers. This leaves the variable in 
    !   a pristine state
    
    ! ------------------------
    if(associated(phaseFunctionVar%legendreCoefficients)) deallocate(phaseFunctionVar%legendreCoefficients)
    if(associated(phaseFunctionVar%scatteringAngle))      deallocate(phaseFunctionVar%scatteringAngle)
    if(associated(phaseFunctionVar%value))                deallocate(phaseFunctionVar%value)
    phaseFunctionVar%extinction             = defaultExtinction
    phaseFunctionVar%singleScatteringAlbedo = defaultSSA
    phaseFunctionVar%description            = ""
  end subroutine finalize_PhaseFunction 
  !------------------------------------------------------------------------------------------
  subroutine finalize_PhaseFunctionTable(table)
    type(phaseFunctionTable), intent(inout) :: table
    ! Free memory and nullify pointers. This leaves the variable in 
    !   a pristine state
    
    ! Local variables
    integer :: i
    
    ! ------------------------
    if(isReady_phaseFunctionTable(table)) then
      ! First we finalize each of the phase functions we've stored
      if(table%oneAngleSet) then
        do i = 2, size(table%phaseFunctions)
          if(associated(table%phaseFunctions(i)%value))           deallocate(table%phaseFunctions(i)%value)
          if(associated(table%phaseFunctions(i)%scatteringAngle)) deallocate(table%phaseFunctions(i)%scatteringAngle)
        end do
        call finalize_PhaseFunction(table%phaseFunctions(1))
      else 
        do i = 1, size(table%phaseFunctions) 
          call finalize_PhaseFunction(table%phaseFunctions(i))
        end do
      end if
      
      ! Next we delete the members of the table itself. 
      if(associated(table%phaseFunctions)) deallocate(table%phaseFunctions)
      if(associated(table%key))            deallocate(table%key)
      if(associated(table%phaseFunctionDescriptions)) &
                                           deallocate(table%phaseFunctionDescriptions)
      table%description = ""
      table%oneAngleSet = .false.
    end if
  end subroutine finalize_PhaseFunctionTable
  !------------------------------------------------------------------------------------------
  ! Utilities, used inside the module
  !------------------------------------------------------------------------------------------
  elemental function isValid(testPhaseFunction)
    type(phaseFunction), intent(in) :: testPhaseFunction
     logical                        :: isValid
     ! Has this phase function been filled in? 
     !   All the other checks should have been done when the 
     !   variable is initialized
     isValid = associated(testPhaseFunction%value) .or. &
               associated(testPhaseFunction%legendreCoefficients)
  end function isValid
  !------------------------------------------------------------------------------------------
  elemental function storedAsLegendre(testPhaseFunction)
    type(phaseFunction), intent(in) :: testPhaseFunction
     logical                                  :: storedAsLegendre
     
     storedAsLegendre = associated(testPhaseFunction%legendreCoefficients)
  end function storedAsLegendre
  !------------------------------------------------------------------------------------------
  elemental function storedAsTabulated(testPhaseFunction)
    type(phaseFunction), intent(in) :: testPhaseFunction
     logical                                  :: storedAsTabulated
     
     storedAsTabulated = associated(testPhaseFunction%value)
  end function storedAsTabulated
  !------------------------------------------------------------------------------------------
  pure function normalizePhaseFunction(scatteringAngle, values) result(normalized)
    real, dimension(:), intent( in) :: scatteringAngle, values
    real, dimension(size(values))   :: normalized

    ! Local variables
    integer :: nAngles
    
    ! ------------------------
    nAngles = size(scatteringAngle)
    ! Normalization with respect to scattering angle
!     normalized(:) = values(:) * 2. / dot_product((scatteringAngle(2:) - scatteringAngle(:nAngles-1)), &
!                                                  0.5 * (values(2:) + values(:nAngles-1)))
    ! Normalization with respect to cosine of the scattering angle
    ! Minus sign because we're doing the integral backwards (from mu=1 to mu=-1)
    normalized(:) = -values(:) * 2. / dot_product((cos(scatteringAngle(2:)) - cos(scatteringAngle(:nAngles-1))), &
                                                  0.5 * (values(2:) + values(:nAngles-1)))
  end function normalizePhaseFunction
  !------------------------------------------------------------------------------------------
end module scatteringPhaseFunctions
