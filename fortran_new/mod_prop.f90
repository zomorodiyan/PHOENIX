!______________________________________________________________________________
!
module property
!______________________________________________________________________________
!
	use initialization
	use parameters
	use dimensions
	use mpi

	implicit none

	! Saved so j-split is computed only on the first call
	integer, save :: jlo_prop = 0
	integer, save :: jhi_prop = 0

	contains

subroutine properties
	implicit none
	integer i,j,k
	real(wp) diffs, diffl
	real(wp) visT, diffT
	integer :: ierr, mpi_rank, mpi_size
	integer :: nj_per, remainder

	! --- Get rank (needed every call for Bcast) ---
	call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)

	! --- J-decomposition (computed once and saved) ---
	if (jlo_prop == 0) then
		call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, ierr)

		nj_per    = nj / mpi_size
		remainder = mod(nj, mpi_size)

		if (mpi_rank < remainder) then
			jlo_prop = mpi_rank * (nj_per + 1) + 1
			jhi_prop = jlo_prop + nj_per
		else
			jlo_prop = remainder * (nj_per + 1) + (mpi_rank - remainder) * nj_per + 1
			jhi_prop = jlo_prop + nj_per - 1
		end if
	end if

	! --- Broadcast inputs from rank 0 to all ranks ---
	! rank0 has the up-to-date arrays; rank1/2/3 need them for their j-slice computation
	call MPI_Bcast(temp,       ni*nj*nk, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
	call MPI_Bcast(fracl,      ni*nj*nk, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
	call MPI_Bcast(solidfield, ni*nj*nk, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

	! --- Zero out arrays (each rank fills only its j-slice) ---
	vis  = 0.0_wp
	diff = 0.0_wp
	den  = 0.0_wp

	! --- Compute properties for local j-slice ---
	do k=1,nk
	do j=jlo_prop,jhi_prop
	do i=1,ni

		visT=0 !dens*dL_mix*0.3*vMag
		diffT=visT/0.9

		diffs=(thconsa*temp(i,j,k)+thconsb)/(acpa*temp(i,j,k)+acpb)
		diffl=thconl/acpl
		vis(i,j,k)=(viscos+visT)                !temperature is larger than liquidus
		diff(i,j,k)=(diffl+diffT)
		den(i,j,k)=denl

		if(temp(i,j,k).ge.tliquid) cycle
			diff(i,j,k)=diffs          !temperature is less than solidus
			vis(i,j,k)=vis_solid
			den(i,j,k)=dens

			if(z(nk)-z(k) .le. layerheight .and. solidfield(i,j,k) .le. powder_threshold)then    !powder properties
				den(i,j,k)=pden
				vis(i,j,k)=vis_solid
				diff(i,j,k)=(pthcona*temp(i,j,k)+pthconb)/(pcpa*temp(i,j,k)+pcpb)
			endif

		if(temp(i,j,k).le.tsolid) cycle   !temperature is between liquidus and solidus
			diff(i,j,k)=(fracl(i,j,k)*diffl+(1.0-fracl(i,j,k))*diffs)
			vis(i,j,k)=(viscos+visT)
			den(i,j,k)=(fracl(i,j,k)*denl+(1.0-fracl(i,j,k))*dens)

	enddo
	enddo
	enddo

	! --- Allreduce: sum each rank's j-slice into complete arrays ---
	call MPI_Allreduce(MPI_IN_PLACE, vis,  ni*nj*nk, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
	call MPI_Allreduce(MPI_IN_PLACE, diff, ni*nj*nk, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
	call MPI_Allreduce(MPI_IN_PLACE, den,  ni*nj*nk, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

end subroutine properties
end module property
