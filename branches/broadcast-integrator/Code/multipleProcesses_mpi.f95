! Copyight 2003-2009, Regents of the University of Colorado. All right reserved
! Use and duplication is permitted under the terms of the 
!   GNU public license, V2 : http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
! NASA has special rights noted in License.txt

! $Revision$, $Date$
! $URL$

module multipleProcesses
  !
  ! Module encapsulating all MPI calls needed in the I3RC model. 
  !   These can be replaced with stubs and the code compiled without MPI. 
  !
  implicit none
  include "mpif.h"
  
  logical :: MasterProc = .true. 
  
  interface sumAcrossProcesses
    module procedure sumAcrossProcesses_Real_Scalar,                         &
                     sumAcrossProcesses_Real_1D, sumAcrossProcesses_Real_2D, &
                     sumAcrossProcesses_Real_3D, sumAcrossProcesses_Real_4D
  end interface sumAcrossProcesses
  
  interface broadcastToAllProcesses
    module procedure broadcast_Logical_1D, &
                     broadcast_Int_S,      &
                     broadcast_Int_1D,   broadcast_Int_2D,  broadcast_Int_3D,  broadcast_Int_4D, & 
                     broadcast_Real_1D,  broadcast_Real_2D, broadcast_Real_3D, broadcast_Real_4D
  end interface broadcastToAllProcesses
contains
  ! -----------------------------------------------------------
  subroutine initializeProcesses(numProcs, thisProcNum)
    integer, intent(out) :: numProcs, thisProcNum
    !
    !Initial MPI calls; how many processors are being used 
    !   and which number is this one? 
    !
    integer :: ierr

    call MPI_INIT(ierr)
    call MPI_COMM_RANK( MPI_COMM_WORLD, thisProcNum, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, numProcs, ierr )
    
    MasterProc = (thisProcNum == 0)
  end subroutine initializeProcesses
  ! -----------------------------------------------------------
  subroutine synchronizeProcesses
    !
    ! Wait for all processors to get to this point
    ! 
    integer :: ierr

    CALL MPI_Barrier( MPI_COMM_WORLD, ierr )
    
  end subroutine synchronizeProcesses
  ! -----------------------------------------------------------
  subroutine finalizeProcesses
    integer :: ierr
    
    call MPI_FINALIZE(ierr)
  end subroutine finalizeProcesses
  ! -----------------------------------------------------------
  function sumAcrossProcesses_Real_Scalar(x) 
    !
    ! Add values across all processors
    !
    real, intent(in) :: x
    real             :: sumAcrossProcesses_Real_Scalar
    
    real    :: temp
    integer :: ierr
    
    call MPI_REDUCE(x, temp, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    sumAcrossProcesses_Real_Scalar = temp
    
  end function sumAcrossProcesses_Real_Scalar
  ! -----------------------------------------------------------
  function sumAcrossProcesses_Real_1D(x) 
    !
    ! Add values across all processors
    !
    real, dimension(:), intent(in) :: x
    real, dimension(size(x))       :: sumAcrossProcesses_Real_1D
    
    real, dimension(size(x)) :: temp
    integer :: ierr
    
    call MPI_REDUCE(x, temp, size(x), MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    sumAcrossProcesses_Real_1D(:) = temp(:)
    
  end function sumAcrossProcesses_Real_1D
  ! -----------------------------------------------------------
  function sumAcrossProcesses_Real_2D(x) 
    !
    ! Add values across all processors
    !
    real, dimension(:, :),       intent(in) :: x
    real, dimension(size(x, 1), size(x, 2)) :: sumAcrossProcesses_Real_2D
    
    real, dimension(size(x, 1), size(x, 2)) :: temp
    integer :: ierr
    
    call MPI_REDUCE(x, temp, size(x), MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    sumAcrossProcesses_Real_2D(:, :) = temp(:, :)
    
  end function sumAcrossProcesses_Real_2D
  ! -----------------------------------------------------------
  function sumAcrossProcesses_Real_3D(x) 
    !
    ! Add values across all processors
    !
    real, dimension(:, :, :), intent(in) :: x
    real, dimension(size(x, 1), size(x, 2), &
                    size(x, 3))          :: sumAcrossProcesses_Real_3D
    
    real, dimension(size(x, 1), size(x, 2), &
                    size(x, 3))             :: temp
    integer :: ierr
    
    call MPI_REDUCE(x, temp, size(x), MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    sumAcrossProcesses_Real_3D(:, :, :) = temp(:, :, :)
    
  end function sumAcrossProcesses_Real_3D
  ! -----------------------------------------------------------
  function sumAcrossProcesses_Real_4D(x) 
    real, dimension(:, :, :, :), intent(in) :: x
    real, dimension(size(x, 1), size(x, 2), &
                    size(x, 3), size(x, 4)) :: sumAcrossProcesses_Real_4D
    
    real, dimension(size(x, 1), size(x, 2), &
                    size(x, 3), size(x, 4)) :: temp
    integer :: ierr
    
    call MPI_REDUCE(x, temp, size(x), MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    sumAcrossProcesses_Real_4D(:, :, :, :) = temp(:, :, :, :)
    
  end function sumAcrossProcesses_Real_4D
  ! -----------------------------------------------------------
  !
  ! Broadcasting
  !
  ! -----------------------------------------------------------
  subroutine broadcast_Logical_1D(x) 
    logical, dimension(:), intent(in) :: x 
    
    integer :: ierr 
    
    call mpi_bcast(x, size(x), MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr) 
  end subroutine broadcast_Logical_1D  
  ! -----------------------------------------------------------
  ! Integers
  ! -----------------------------------------------------------
  subroutine broadcast_Int_S(x) 
    integer, intent(in) :: x 
    
    integer :: ierr 
    
    call mpi_bcast(x, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
  end subroutine broadcast_Int_S  
  ! -----------------------------------------------------------
  subroutine broadcast_Int_1D(x) 
    integer, dimension(:), intent(in) :: x 
    
    integer :: ierr 
    
    call mpi_bcast(x, size(x), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
  end subroutine broadcast_Int_1D  
  ! -----------------------------------------------------------
  subroutine broadcast_Int_2D(x) 
    integer, dimension(:, :), intent(in) :: x 
    
    integer :: ierr 
    
    call mpi_bcast(x, size(x), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
  end subroutine broadcast_Int_2D  
  ! -----------------------------------------------------------
  subroutine broadcast_Int_3D(x) 
    integer, dimension(:, :, :), intent(in) :: x 
    
    integer :: ierr 
    
    call mpi_bcast(x, size(x), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
  end subroutine broadcast_Int_3D  
  ! -----------------------------------------------------------
  subroutine broadcast_Int_4D(x) 
    integer, dimension(:, :, :, :), intent(in) :: x 
    
    integer :: ierr 
    
    call mpi_bcast(x, size(x), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
  end subroutine broadcast_Int_4D  
  ! -----------------------------------------------------------
  ! Reals
  ! -----------------------------------------------------------
  subroutine broadcast_Real_1D(x) 
    real, dimension(:), intent(in) :: x 
    
    integer :: ierr 
    
    call mpi_bcast(x, size(x), MPI_REAL, 0, MPI_COMM_WORLD, ierr) 
  end subroutine broadcast_Real_1D  
  ! -----------------------------------------------------------
  subroutine broadcast_Real_2D(x) 
    real, dimension(:, :), intent(in) :: x 
    
    integer :: ierr 
    
    call mpi_bcast(x, size(x), MPI_REAL, 0, MPI_COMM_WORLD, ierr) 
  end subroutine broadcast_Real_2D  
  ! -----------------------------------------------------------
  subroutine broadcast_Real_3D(x) 
    real, dimension(:, :, :), intent(in) :: x 
    
    integer :: ierr 
    
    call mpi_bcast(x, size(x), MPI_REAL, 0, MPI_COMM_WORLD, ierr) 
  end subroutine broadcast_Real_3D  
  ! -----------------------------------------------------------
  subroutine broadcast_Real_4D(x) 
    real, dimension(:, :, :, :), intent(in) :: x 
    
    integer :: ierr 
    
    call mpi_bcast(x, size(x), MPI_REAL, 0, MPI_COMM_WORLD, ierr) 
  end subroutine broadcast_Real_4D  
  ! -----------------------------------------------------------
end module multipleProcesses