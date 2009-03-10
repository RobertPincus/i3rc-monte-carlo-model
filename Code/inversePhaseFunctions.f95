! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision: 1.7 $, $Date: 2009/03/09 19:17:22 $
! $Name: Cornish-Gilliflower $
! --------------------------------------------
module inversePhaseFunctions
  ! Provides a way to compute, store, and retrieve inverse (cummulative) phase functions 
  !   from objects of type scatteringPhaseFunction. The inverse phase function is a
  !   table of scattering angle as a 
  use ErrorMessages
  use scatteringPhaseFunctions, only: phaseFunction, phaseFunctionTable,                 &
                                      isReady_PhaseFunctionTable,                        &
                                      getInfo_PhaseFunction, getInfo_PhaseFunctionTable, &
                                      isReady_PhaseFunctionTable,                        &
                                      getElement, getPhaseFunctionValues
  use numericUtilities,         only: findIndex, computeLobattoTerms
  implicit none
  
  private
  public :: computeInversePhaseFuncTable
contains
  !------------------------------------------------------------------------------------------
  subroutine computeInversePhaseFuncTable(forwardTable, inverseTable, status)
    !
    ! Compute the scattering angle as a function of the nSteps steps between 0 and 1
    !   for every entry in a phaseFunctionTable
    !
    type(phaseFunctionTable), intent(in   ) :: forwardTable
    real, dimension(:, :),    intent(  out) :: inverseTable
    type(ErrorMessage),       intent(inout) :: status
    
    ! Local variable
    integer            :: nTableEntries, nSteps, i
    type(ErrorMessage) :: localStatus
    
    ! --------------------------------
    nSteps = size(inverseTable, 1) 
    if(.not. isReady_phaseFunctionTable(forwardTable)) then
      call setStateToFailure(status, "computeInversePhaseFunctionTable: Forward table isn't ready.")
    else   
      call getInfo_phaseFunctionTable(forwardTable, nEntries = nTableEntries, status = status)
      if(size(inverseTable, 2) /= nTableEntries) &
        call setStateToFailure(status, "computeInversePhaseFunctionTable: Array for inverse table has the wrong number of entries")
    end if 
    
    if(.not. stateIsFailure(status)) then 
      entryLoop: do i = 1, nTableEntries
        !
        ! We're safe not checking the result from localStatus because we've protected against 
        !   the two possible error conditions (table not being ready or asking for an out-of-bound element)
        !
        call computeInversePhaseFunction(getElement(i, forwardTable, localStatus), inverseTable(:, i), status)
        if(stateIsFailure(status)) then 
          call setStateToFailure(status, "computeInversePhaseFunctionTable: Can't compute inverse tables.")
          exit entryLoop
        end if 
      end do entryLoop
    end if
    if(.not. stateIsFailure(status)) call setStateToSuccess(status)
    
  end subroutine computeInversePhaseFuncTable
  !------------------------------------------------------------------------------------------
  subroutine computeInversePhaseFunction(thisPhaseFunction, inverseTable, status) 
    !
    ! Compute the scattering angle as a function of the nSteps steps between 0 and 1. 
    !   We do this by finding the angle at which the fraction of the integrated phase function
    !   is equal to the desired probability using an analytic relationship I worked out. 
    !
    type(phaseFunction), intent(in   ) :: thisPhaseFunction
    real, dimension(:),  intent(  out) :: inverseTable
    type(ErrorMessage),  intent(inout) :: status

    ! Local variables
    integer                            :: i, index, nAngles, nMoments, nSteps
    real,    dimension(:), allocatable :: angles, values, mus, cdf
    integer, dimension(size(inverseTable)) &
                                       :: indicies
    real                               :: probabilityValue
    
    ! ------------------------
    nSteps = size(inverseTable)
    call getInfo_PhaseFunction(thisPhaseFunction, nMoments, nAngles, status = status)
    
    if(.not. stateIsFailure(status)) then
      if(nAngles > 0) then
        !
        ! The phase function is stored as angle-value pairs.  We ask for the phase function 
        !   at the native resolution. We'll calculate the CDF at each of those angles, then 
        !   find the CDF for the requested points by interpolation. 
        !
        allocate(angles(nAngles), values(nAngles), mus(nAngles), cdf(nAngles))
        call getInfo_PhaseFunction(thisPhaseFunction, nativeAngles = angles, status = status)
        if(.not. stateIsFailure(status)) &
          call getPhaseFunctionValues(thisPhaseFunction, angles, values, status)
        mus(:) = cos(angles(nAngles:1:-1)); values = values(nAngles:1:-1)
      else 
        !
        ! The phase function is stored as Legendre polynomials. We'll ask for it at the 
        !   zeros of the L-1 polynomial. 
        !
        ! The phase function for isotropic scattering has no moments beyond P0, so 
        !    it's possible that nMoments = 0. 
        !
        nAngles = max(nMoments, 2)
        allocate(angles(nAngles), values(nAngles), mus(nAngles), cdf(nAngles))
        ! We're going to ignore the weights (which we store for now in the array angles(:))
        call computeLobattoTerms(mus = mus(:), weights = angles(:))
        call getPhaseFunctionValues(thisPhaseFunction, acos(mus(nAngles:1:-1)), values, status)
        values = values(nAngles:1:-1)
      end if 
      if(.not. stateIsFailure(status)) then 
        !
        ! Compute the cumulative phase function (CDF) at each value of the phase function. 
        !   Integrate using trapezoidal integration in cosine of the
        !   scattering angle. Reorder the angles and values to increase in mu. 
        !
        cdf(1) = 0. 
        do i = 2, nAngles
          cdf(i) = cdf(i - 1) + (mus(i) - mus(i-1)) * 0.5 * (values(i) + values(i-1))
        end do
        ! The CDF should be normalized to 2 (i.e. cdf(nAngles) should equal 2) but the 
        !   important thing is to have cdf(:) run from 0 to 1. 
        !
        cdf(:) = cdf(:)/cdf(nAngles)

        ! FInd the location of each probability value in the CDF
        !
        indicies(1) = findIndex(0., cdf)
        do i = 2, nSteps
          probabilityValue = float(i - 1)/float(nSteps - 1)
          indicies(i) = findIndex(probabilityValue, cdf, firstGuess = indicies(i-1))
        end do
        
        do i = 1, nSteps - 1
          probabilityValue = float(i - 1)/float(nSteps - 1)
          index = indicies(i)
          !
          ! Following is an analytic inversion for scattering angle given an arbitrary value of the CDF from 
          !   a tabulated CDF, consistent with a trapezoidal integration of the phase function with repect to 
          !   cosine of the scattering angle. 
          !   I'm not aware of a reference. 
          !
          ! Special cases: the phase function cloud be identically 0 and the CDF locally constant...
          !
          if(cdf(index+1) - cdf(index) <= spacing(cdf(index))) then 
            inverseTable(i) = acos(mus(index))
          !
          ! ... or the phase function might be constant with angle locally...
          !
          else if(abs(values(index) - values(index + 1)) <= spacing(values(index))) then
            inverseTable(i) = acos(mus(index) +                  &
                                   (mus(index+1) - mus(index)) * &
                                   (probabilityValue - cdf(index))/(cdf(index+1) - cdf(index)))
          !
          !... Otherwise we use the generic relationship.  
          !
          else
            inverseTable(i) = acos(mus(index) +                                                                       &
                                   (mus(index+1) - mus(index))/(values(index) - values(index + 1)) *                  &
                                   (values(index) - sqrt(((cdf(index+1) - probabilityValue) * values(index    )**2 +  &
                                                          (probabilityValue - cdf(index)  ) * values(index + 1)**2) / &
                                                         (cdf(index+1) - cdf(index)))))
          end if
        end do
        inverseTable(nSteps) = 0. 
        call setStateToSuccess(status)
        deallocate(angles, values, mus, cdf)
      end if 
    end if 

  end subroutine computeInversePhaseFunction
  ! ------------------------------
end module inversePhaseFunctions