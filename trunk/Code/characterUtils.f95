! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision: 1.3 $, $Date: 2009/03/09 19:17:22 $
! $Name: Cornish-Gilliflower $

module CharacterUtils
  implicit none
  private
  
  integer :: maxStringLength = 25
  public :: CharToInt, IntToChar, CharToReal
contains
  elemental function CharToInt(inputString)
    character(len = *), intent( in) :: inputString
    integer                         :: CharToInt
    !
    ! Reads a character string into an integer variable 
    !    No internal error checking is done - if an incorrect string
    !    is passed in the entire program will fail with a Fortran runtime error. 
    
    ! Local variables
    character(len = maxStringLength) tempString

    ! Remove trailing blanks, which might be read as 0s on some processors
    ! Use list-directed read. 
    tempString = trim(inputString)
    read(tempString, *) CharToInt
    
  end function CharToInt
  ! --------------------------------------------------------------------------------
  elemental function IntToChar(integerValue)
    integer,              intent( in) :: integerValue
    character(len = maxStringLength)  :: IntToChar
    !
    !   Creates the character representation of an integer.  
    !  

    write(IntToChar, *) integerValue
    IntToChar = AdjustL(IntToChar)
  end function IntToChar
  ! --------------------------------------------------------------------------------
  elemental function CharToReal(inputString)
    character(len = *), intent( in) :: inputString
    real                            :: CharToReal
    !
    !   Reads a character string into an real variable 
    !    No internal error checking is done - if an incorrect string
    !      is passed in the entire program will fail with a Fortran runtime error. 
    !
    ! !END
    
    ! Local variables
    character(len = maxStringLength) tempString

    ! Remove trailing blanks, which might be read as 0s on some processors
    ! Use list-directed read. 
    !
    tempString = trim(AdjustL(inputString))
    read(tempString, *) CharToReal
  end function CharToReal
  
end module CharacterUtils
