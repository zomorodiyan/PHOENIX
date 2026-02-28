!______________________________________________________________________________
!
module mpi_info
!______________________________________________________________________________
! MPI rank/size and j-direction domain decomposition.
!
! Usage:
!   call MPI_Init(ierr)
!   call mpi_init_j_decomp(nj)
!   ...
!   call MPI_Finalize(ierr)
!
	use mpi
	implicit none

	integer :: mpi_rank = 0
	integer :: mpi_size = 1
	integer :: jlo_local = 1
	integer :: jhi_local = 1

	contains

subroutine mpi_init_j_decomp(nj)
	integer, intent(in) :: nj
	integer :: nj_per, remainder, ierr

	call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
	call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, ierr)

	nj_per    = nj / mpi_size
	remainder = mod(nj, mpi_size)

	! Spread remainder rows across the first 'remainder' ranks
	if (mpi_rank < remainder) then
		jlo_local = mpi_rank * (nj_per + 1) + 1
		jhi_local = jlo_local + nj_per
	else
		jlo_local = remainder * (nj_per + 1) + (mpi_rank - remainder) * nj_per + 1
		jhi_local = jlo_local + nj_per - 1
	end if

end subroutine mpi_init_j_decomp

end module mpi_info
