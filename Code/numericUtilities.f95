! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision$, $Date$
! $URL$
module numericUtilities
  implicit none
  private 
  public :: computeLobattoTerms, computeGaussLegendreTerms, &
            computeLegendrePolynomials, findIndex
contains
  !------------------------------------------------------------------------------------------
  pure subroutine computeLobattoTerms(mus, weights)
    real, dimension(:), intent(out) :: mus, weights
    ! Find the weights and abscissas for Lobatto on the interval [-1, 1]. 
    !   For n point quadrature the weights are the 0s of the derivative of the 
    !   (n-1)th degree Legendre polynomial. We find the roots using Newton's method. 
    
    ! Local variables
    real,    parameter                 :: relativeAccuracy = 3.
    integer, parameter                 :: maxIterations = 25
    real, dimension(:, :), allocatable :: legendreP
    real, dimension(:),    allocatable :: trialMus, derivatives, secondDerivs, lastGuess
    integer                            :: nTerms, midPoint, i
    real                               :: pi, c1

    pi = acos(-1.)
    nTerms = min(size(mus), size(weights))
    midPoint = (nTerms + 1)/2
    
    allocate(trialMus(midPoint-1),    lastGuess(midPoint-1),    &
             derivatives(midPoint-1), secondDerivs(midPoint-1), &
             legendreP(0:nTerms-1, midPoint-1))
             
    ! First guess, depending on whether the number of points is even or odd. 
    if(mod(nTerms, 2) == 1) then
      c1 = 1. 
    else
      c1 = 0.5
    end if
    ! c1 = merge(1., .5, mod(nTerms, 2) == 1)
    trialMus(:) = sin(pi * ((/ (real(i), i = 1, midPoint-1) /) - c1)/(nTerms - 1. + .5))
    
    ! One iteration of Newton's method
    legendreP(:,:) = computeLegendrePolynomials(nTerms-1, trialMus(:))
    ! First and second derivatives, computed analytically 
    derivatives(:) = (nTerms - 1) *                                                       &
                    (trialMus(:) * legendreP(nTerms - 1, :) - legendreP(nTerms - 2, :)) / &
                    (trialMus(:)**2 - 1.)
    
    secondDerivs(:) = (2. * trialMus * derivatives(:) -                      &
                       (nTerms * (nTerms - 1) * legendreP(nTerms - 1, :))) / &
                      (1. - trialMus(:)**2)
    lastGuess(:) = trialMus(:)
    trialMus(:) = trialMus(:) - derivatives(:)/secondDerivs(:)
    
    ! Second and following iterations of Newton's method
    i = 0
    findRoots: do
      if(all(abs(trialMus(:) - lastGuess(:)) <= relativeAccuracy * spacing(trialMus))) &
        exit findRoots
      ! I can't figure out how to do compute the Legendre coefficients within the where
      !   construct; this involves some wasted computation computing Legendre
      !   polynomials for mus that have already converged. 
      legendreP(:,:) = computeLegendrePolynomials(nTerms-1, trialMus(:))
      
      where(abs(trialMus(:) - lastGuess(:)) > relativeAccuracy * spacing(trialMus))
        ! Compute the derivative "using a well-known relation" 
        derivatives(:) = (nTerms - 1) *                                                       &
                        (trialMus(:) * legendreP(nTerms - 1, :) - legendreP(nTerms - 2, :)) / &
                        (trialMus(:)**2 - 1.)
        
        secondDerivs(:) = (2. * trialMus * derivatives(:) -                      &
                           (nTerms * (nTerms - 1) * legendreP(nTerms - 1, :))) / &
                          (1. - trialMus(:)**2)
        lastGuess(:) = trialMus(:)
        trialMus(:) = trialMus(:) - derivatives(:)/secondDerivs(:)
      end where
      i = i + 1
      if (i > maxIterations) exit findRoots
    end do findRoots
  
    
    mus(1)     = -1
    weights(1) = 2./(nTerms * (nTerms - 1))
    mus(midPoint:2:-1) = - trialMus(:)
    weights(midPoint:2:-1) = 2./(nTerms * (nTerms - 1) * legendreP(nTerms - 1, :)**2)

    ! Calculation is symmetric
    if(mod(nTerms, 2) == 0) then
      mus(midPoint+1:nTerms)     = -mus(midPoint:1:-1)
      weights(midPoint+1:nTerms) = weights(midPoint:1:-1)
    else
      mus(midPoint:nTerms)     = -mus(midPoint:1:-1)
      weights(midPoint:nTerms) = weights(midPoint:1:-1)
    end if
    mus(nTerms+1:) = 0.; weights(nTerms+1:) = 0.
 
    deallocate(trialMus, derivatives, secondDerivs, lastguess, legendreP)
  end subroutine computeLobattoTerms
  !------------------------------------------------------------------------------------------
  pure subroutine computeGaussLegendreTerms(mus, weights)
    real, dimension(:), intent(out) :: mus, weights
    ! Find the weights and abscissas for Gauss-Legendre integration on the interval (-1, 1). 
    !   For n point quadrature the weights are the 0s of the nth degree Legendre polynomial. 
    !   We find the roots using Newton's method. 
        
    ! Local variables
    real,    parameter                 :: relativeAccuracy = 2.0
    integer, parameter                 :: maxIterations = 25
    real, dimension(:, :), allocatable :: legendreP
    real, dimension(:),    allocatable :: trialMus, derivatives, lastGuess
    integer                            :: nTerms, midPoint, i
    real                               :: pi
    
    pi = acos(-1.)
    nTerms = min(size(mus), size(weights))
    midPoint = (nTerms + 1)/2
    
    allocate(trialMus(midPoint), derivatives(midPoint), lastGuess(midPoint), &
             legendreP(0:nTerms, midPoint))
             
    ! First guess
    trialMus(:) = cos(pi * ((/ (real(i), i = 1, midPoint) /) - .25)/(nTerms + .5))
    
    ! One iteration of Newton's method
    legendreP(:,:) = computeLegendrePolynomials(nTerms, trialMus(:))
    ! The derivative is analytic 
    derivatives(:) = nTerms *                                                         &
                    (trialMus(:) * legendreP(nTerms, :) - legendreP(nTerms - 1, :)) / &
                    (trialMus(:)**2 - 1.)
    lastGuess(:) = trialMus(:)
    trialMus(:) = trialMus(:) - legendreP(nTerms, :)/derivatives(:)
    
    ! Second and following iterations of Newton's method
    i = 0
    findRoots: do
      if(all(abs(trialMus(:) - lastGuess(:)) <= relativeAccuracy * spacing(trialMus))) &
        exit findRoots
      ! I can't figure out how to do compute the Legendre coefficients within the where
      !   construct; this involves some wasted computation computing Legendre
      !   polynomials for mus that have already converged. 
      legendreP(:,:) = computeLegendrePolynomials(nTerms, trialMus(:))
      
      where(abs(trialMus(:) - lastGuess(:)) > relativeAccuracy * spacing(trialMus))
        ! Compute the derivative "using a well-known relation" 
        derivatives(:) = nTerms *                                                         &
                        (trialMus(:) * legendreP(nTerms, :) - legendreP(nTerms - 1, :)) / &
                        (trialMus(:)**2 - 1.)
        lastGuess(:) = trialMus(:)
        trialMus(:) = trialMus(:) - legendreP(nTerms, :)/derivatives(:)
      end where
    i = i+ 1
    if (i > maxIterations) exit findRoots
    end do findRoots
  
    mus(    :midPoint) = -trialMus(:)
    weights(:midPoint) = 2./((1. - trialMus(:)**2) * derivatives(:)**2)
    
    ! Calculation is symmetric
    if(mod(nTerms, 2) == 0) then
      mus(midPoint+1:nTerms)     = -mus(midPoint:1:-1)
      weights(midPoint+1:nTerms) = weights(midPoint:1:-1)
    else
      mus(midPoint:nTerms)     = -mus(midPoint:1:-1)
      weights(midPoint:nTerms) = weights(midPoint:1:-1)
    end if
    mus(nTerms+1:) = 0.; weights(nTerms+1:) = 0.
    
    deallocate(trialMus, derivatives, lastguess, legendreP)
  end subroutine computeGaussLegendreTerms
  !------------------------------------------------------------------------------------------
  pure function computeLegendrePolynomials(maxL, mus) result(legendreP)
    integer,            intent( in) :: maxL
    real, dimension(:), intent( in) :: mus
    real, dimension(0:maxL, size(mus)) :: legendreP
    ! Compute the value of the Legendre polynomials from order 0 to order maxL  
    !   at the values of mu provided. 
    
    ! Local variable
    integer :: l
    
    legendreP(0, :) = 1
    legendreP(1, :) = mus(:)
    ! Legendre polynomials from recursion relation 
    !   See, for example, Numerical Recipies in Fortran 2nd Ed., page 172
    !   Um. I don't know for certain that this is unconditionally stable. 
    do l = 1, maxL - 1
      legendreP(l+1, :) = ((2*l + 1) * mus(:) * legendreP(l, :) - l * legendreP(l-1, :))/(l + 1)
    end do 
  end function computeLegendrePolynomials
  !------------------------------------------------------------------------------------------
  pure function findIndex(value, table, firstGuess)
    real,               intent( in) :: value
    real, dimension(:), intent( in) :: table
    integer, optional,  intent( in) :: firstGuess
    integer                         :: findIndex
    !
    ! Find the index i into the table such that table(i) <= value < table(i+i)
    !   This is modeled after routine "hunt" from Numerical Recipes, 2nd ed., 
    !   pg 112. Here we know that the values in the table are always increasing,
    !   that every value should be spanned by the table entries, and the firstGuess 
    !   always makes sense. 
    
    ! Local variables
    integer :: lowerBound, upperBound, midPoint
    integer :: increment

    ! Hunting; only done if a first guess is supplied
    !  Move upper and lower bounds around until the value is spanned by 
    !   table(lowerBound) and table(upperBound). Make the interval twice as
    !   big at each step
    if(present(firstGuess)) then
      lowerBound = firstGuess
      increment = 1
      huntingLoop: do
        upperBound = min(lowerBound + increment, size(table))
        if(lowerBound == size(table) .or. &
           (table(lowerBound) <= value .and. table(upperBound) > value)) exit huntingLoop
        if(table(lowerBound) > value) then
          upperBound = lowerBound
          lowerBound = max(upperBound - increment, 1)
        else 
          ! Both table(lowerBound) and table(upperBound) are <= value
          lowerBound = upperBound  
        end if
        increment = increment * 2
      end do huntingLoop
    else 
      lowerBound = 0; upperBound = size(table)
    end if 
    
    ! Bisection: figure out which half of the remaining interval holds the 
    !   desired value, discard the other half, and repeat
    bisectionLoop: do
      if(lowerBound == size(table) .or. upperBound <= lowerBound + 1) exit bisectionLoop
      midPoint = (lowerBound + upperBound)/2
      if(value >= table(midPoint)) then
        lowerBound = midPoint
      else
        upperBound = midPoint
      end if
    end do bisectionLoop
    
    findIndex = lowerBound
  end function findIndex
  !------------------------------------------------------------------------------------------
  elemental function gammln(x)
    ! Returns the natural log of Gamma(x)
    implicit none
    real, intent(in) :: x 
    real             :: gammln
    
    integer            :: j
    ! Double precision
    integer, parameter :: EightByteReal = selected_real_kind(P = 13, R = 307)
    real(kind =  EightByteReal) :: ser
    real(kind =  EightByteReal) , parameter  :: stp = 2.50662827465D0
    real(kind =  EightByteReal) , parameter, &
                    dimension(6) :: cof = (/ 76.18009173D0, -86.50532033D0, &
                                             24.01409822D0, - 1.231739516D0,&
                                               .120858003D-2,-.536382D-5 /)
                                              
    ser = 1. + sum(cof(:)/(/ (real(x + j - 1, kind = EightByteReal), j = 1, size(cof)) /))
    GAMMLN= (x - .5) * log(real(x + 4.5, kind = EightByteReal)) - &
            (x + 4.5) + log(stp * ser)
    
  end function gammln
  !------------------------------------------------------------------------------------------
end module numericUtilities
