
!****************************************************************************
!
!  AM-CFD
!
!  Version 1.0
!
!  July, 2020, With English Annotation, Northwestern University
!
!****************************************************************************

program main
	use geometry
	use initialization
	use discretization
	use dimensions
	use boundary
	use source
	use residue
	use solver
	use fluxes
	use parameters
	use printing
	use constant
	use entotemp
	use convergence
	use property
	use revision
	use laserinput
	use toolpath
	use timing
	use mpi

	implicit none
	integer i,j,k
	real(wp) amaxres
	real(wp) t0, t1
	integer :: mpi_rank, ierr
	integer :: do_exit_iter   ! flag broadcast to sync iter_loop exit across ranks

	call MPI_Init(ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)

	call read_data
	call read_toolpath

	call generate_grid          ! sets ni,nj,nk and allocates geometry
	call allocate_fields(ni, nj, nk)
	call allocate_source(ni, nj, nk)
	call allocate_print(ni, nj, nk)
	call allocate_laser(ni, nj)
	if (mpi_rank == 0) call OpenFiles
	call initialize

	call StartTime

	itertot=0
	timet=small

!------time stepping loop------------------------------------
	time_loop: do while (timet.lt.timax)
		timet=timet+delt
		niter=0

!-------Move laser and calculate freesurface (rank 0 only)-----------------
		if (mpi_rank == 0) then
			call cpu_time(t0)
			call laser_beam
			call read_coordinates
			call cpu_time(t1)
			t_laser = t_laser + (t1 - t0)
		end if

!-----iteration loop within each time step----------------
		iter_loop: do while (niter.lt.maxit)
			niter=niter+1
			itertot=itertot+1

!-----solve energy equation (formerly ivar=5)-----
!-----properties: ALL ranks participate (MPI Bcast + j-split compute + Allreduce)-----
			call cpu_time(t0)
			call properties
			call cpu_time(t1)
			t_prop = t_prop + (t1 - t0)

!-----all remaining computation: rank 0 only-----
			if (mpi_rank == 0) then

				call cpu_time(t0)
				call bound_enthalpy
				call cpu_time(t1)
				t_bound = t_bound + (t1 - t0)

				call cpu_time(t0)
				call discretize_enthalpy
				call cpu_time(t1)
				t_discret = t_discret + (t1 - t0)
				t_discret_enthalpy = t_discret_enthalpy + (t1 - t0)

				call cpu_time(t0)
				call source_enthalpy
				call cpu_time(t1)
				t_sour = t_sour + (t1 - t0)
				t_sour_enthalpy = t_sour_enthalpy + (t1 - t0)

				call cpu_time(t0)
				call calc_enthalpy_residual
				call cpu_time(t1)
				t_resid = t_resid + (t1 - t0)

				call cpu_time(t0)
				call enhance_converge_speed
				call cpu_time(t1)
				t_converge = t_converge + (t1 - t0)

				call cpu_time(t0)
				call solution_enthalpy
				call cpu_time(t1)
				t_solve = t_solve + (t1 - t0)
				t_solve_enthalpy = t_solve_enthalpy + (t1 - t0)

				call cpu_time(t0)
				call enthalpy_to_temp
				call cpu_time(t1)
				t_entot = t_entot + (t1 - t0)

				call cpu_time(t0)
				call pool_size
				call cpu_time(t1)
				t_dimen = t_dimen + (t1 - t0)

				if(tpeak.gt.tsolid) then
					call cpu_time(t0)
					call cleanuvw
					call cpu_time(t1)
					t_solve = t_solve + (t1 - t0)
					t_solve_cleanuvw = t_solve_cleanuvw + (t1 - t0)

!-----solve u-momentum (formerly ivar=1)-----
					call cpu_time(t0)
					call bound_uv(1)
					call cpu_time(t1)
					t_bound = t_bound + (t1 - t0)
					call cpu_time(t0)
					call discretize_momentum(1)
					call cpu_time(t1)
					t_discret = t_discret + (t1 - t0)
					t_discret_momentum = t_discret_momentum + (t1 - t0)
					call cpu_time(t0)
					call source_momentum(1)
					call cpu_time(t1)
					t_sour = t_sour + (t1 - t0)
					t_sour_momentum = t_sour_momentum + (t1 - t0)
					call cpu_time(t0)
					call calc_momentum_residual(uVel, resoru, .true.)
					call cpu_time(t1)
					t_resid = t_resid + (t1 - t0)
					call cpu_time(t0)
					call solution_uvw(uVel)
					call cpu_time(t1)
					t_solve = t_solve + (t1 - t0)
					t_solve_uvw = t_solve_uvw + (t1 - t0)

!-----solve v-momentum (formerly ivar=2)-----
					call cpu_time(t0)
					call bound_uv(2)
					call cpu_time(t1)
					t_bound = t_bound + (t1 - t0)
					call cpu_time(t0)
					call discretize_momentum(2)
					call cpu_time(t1)
					t_discret = t_discret + (t1 - t0)
					t_discret_momentum = t_discret_momentum + (t1 - t0)
					call cpu_time(t0)
					call source_momentum(2)
					call cpu_time(t1)
					t_sour = t_sour + (t1 - t0)
					t_sour_momentum = t_sour_momentum + (t1 - t0)
					call cpu_time(t0)
					call calc_momentum_residual(vVel, resorv, .false.)
					call cpu_time(t1)
					t_resid = t_resid + (t1 - t0)
					call cpu_time(t0)
					call solution_uvw(vVel)
					call cpu_time(t1)
					t_solve = t_solve + (t1 - t0)
					t_solve_uvw = t_solve_uvw + (t1 - t0)

!-----solve w-momentum (formerly ivar=3)-----
					call cpu_time(t0)
					call bound_w
					call cpu_time(t1)
					t_bound = t_bound + (t1 - t0)
					call cpu_time(t0)
					call discretize_momentum(3)
					call cpu_time(t1)
					t_discret = t_discret + (t1 - t0)
					t_discret_momentum = t_discret_momentum + (t1 - t0)
					call cpu_time(t0)
					call source_momentum(3)
					call cpu_time(t1)
					t_sour = t_sour + (t1 - t0)
					t_sour_momentum = t_sour_momentum + (t1 - t0)
					call cpu_time(t0)
					call calc_momentum_residual(wVel, resorw, .false.)
					call cpu_time(t1)
					t_resid = t_resid + (t1 - t0)
					call cpu_time(t0)
					call solution_uvw(wVel)
					call cpu_time(t1)
					t_solve = t_solve + (t1 - t0)
					t_solve_uvw = t_solve_uvw + (t1 - t0)

!-----solve pressure correction (formerly ivar=4)-----
					call cpu_time(t0)
					call bound_pp
					call cpu_time(t1)
					t_bound = t_bound + (t1 - t0)
					call cpu_time(t0)
					call discretize_pp
					call cpu_time(t1)
					t_discret = t_discret + (t1 - t0)
					t_discret_pp = t_discret_pp + (t1 - t0)
					call cpu_time(t0)
					call source_pp
					call cpu_time(t1)
					t_sour = t_sour + (t1 - t0)
					t_sour_pp = t_sour_pp + (t1 - t0)
					call cpu_time(t0)
					call calc_pressure_residual
					call cpu_time(t1)
					t_resid = t_resid + (t1 - t0)
					call cpu_time(t0)
					call solution_uvw(pp)
					call cpu_time(t1)
					t_solve = t_solve + (t1 - t0)
					t_solve_uvw = t_solve_uvw + (t1 - t0)
					call cpu_time(t0)
					call revision_p
					call cpu_time(t1)
					t_revise = t_revise + (t1 - t0)
				endif

				call cpu_time(t0)
				call heat_fluxes
				call cpu_time(t1)
				t_flux = t_flux + (t1 - t0)
				amaxres=max(resorm, resoru,resorv,resorw)

			end if  ! mpi_rank == 0

!-----convergence criterion: rank 0 decides, broadcasts to all ranks------------
			do_exit_iter = 0
			if (mpi_rank == 0) then
				if(toolmatrix(PathNum,5) .ge. laser_on_threshold)then
					if(amaxres.lt.conv_res_heat .and. ratio.le.ratio_upper .and. ratio.ge.ratio_lower) &
						do_exit_iter = 1
				else
					if(resorh.lt.conv_res_cool) do_exit_iter = 1
				endif
			end if
			call MPI_Bcast(do_exit_iter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
			if (do_exit_iter == 1) exit iter_loop

		end do iter_loop

		call cpu_time(t0)
		if (mpi_rank == 0) call CalTime
		if (mpi_rank == 0) call outputres
		call cpu_time(t1)
		t_print = t_print + (t1 - t0)

		call cpu_time(t0)
		if (mpi_rank == 0) then
			do k=1,nk
			do j=1,nj
			do i=1,ni
				if(temp(i,j,k).le.tsolid) then
					uVel(i,j,k)=0.0
					vVel(i,j,k)=0.0
					wVel(i,j,k)=0.0
				endif
				unot(i,j,k)=uVel(i,j,k)
				vnot(i,j,k)=vVel(i,j,k)
				wnot(i,j,k)=wVel(i,j,k)
				tnot(i,j,k)=temp(i,j,k)
				hnot(i,j,k)=enthalpy(i,j,k)
				fraclnot(i,j,k)=fracl(i,j,k)
			enddo
			enddo
			enddo
			call Cust_Out
		end if
		call cpu_time(t1)
		t_other = t_other + (t1 - t0)

	end do time_loop

	call EndTime

	if (mpi_rank == 0) call write_timing_report(itertot, timet)
	call MPI_Finalize(ierr)
	stop
	end
