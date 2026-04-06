
!****************************************************************************
!
!  PHOENIX
!
!  Version 1.0
!
!  March. 2026, Arizona State University, Zhengtao Gan
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
	use omp_lib
	use defect_field
	use adaptive_mesh_mod
	use prediction
	use mechanical_solver
	use mech_io
	use species, only: allocate_species, init_species, &
		species_bc, solve_species, concentration, conc_old, &
		mix, tsolid2

	implicit none
	integer i,j,k
	integer step_idx
	integer ilo, ihi, jlo, jhi, klo, khi
	real(wp) amaxres
	real(wp) t0, t1
	real(wp) t_step_wall
	real(wp) wall_start, wall_elapsed
	integer  mech_solve_count
	real(wp) mech_res
	integer  mech_newton_iters, mech_cg_iters, n_yield
	real(wp) t_mech_cpu0, t_mech_cpu1
	logical  mech_parallel  ! .true. if mechanical runs as separate process
	character(len=32) :: run_mode  ! 'thermal', 'mechanical', or '' (default=all)

	call read_data

	! Check run mode from environment (set by run.sh for parallel)
	run_mode = ''
	block
		character(len=32) :: envval
		integer :: ios
		call get_environment_variable('PHOENIX_RUN_MODE', envval, status=ios)
		if (ios == 0) run_mode = trim(adjustl(envval))
	end block

	! Mechanical-only process: skip thermal init, go straight to mechanical loop
	if (run_mode == 'mechanical') then
		call generate_grid
		! init_mechanical and init_mech_history called inside run_mechanical_loop
		call run_mechanical_loop()
		stop
	endif

	call read_toolpath
	call generate_grid          ! sets ni,nj,nk and allocates geometry
	call allocate_fields(ni, nj, nk)
	call allocate_source(ni, nj, nk)
	call allocate_print(ni, nj, nk)
	call allocate_laser(ni, nj)
	call allocate_defect(ni, nj, nk)
	call OpenFiles
	call initialize
	if (adaptive_flag == 1) call amr_init()
	mech_solve_count = 0
	mech_parallel = (mechanical_flag == 1 .and. n_mech_threads > 0)
	if (mechanical_flag == 1 .and. .not. mech_parallel) then
		! Serial mode: init mechanical in main process
		call init_mechanical()
		call init_mech_history(Nnx, Nny, Nnz)
	endif
	call init_thermal_history
	call init_meltpool_history

	if (species_flag == 1) then
		call allocate_species
		call init_species
	endif

	call StartTime
	wall_start = omp_get_wtime()

	! Set thermal thread count (parallel mode reduces from total)
	if (mech_parallel) then
		call omp_set_num_threads(n_thermal_threads)
		write(9,'(A,I3,A,I3)') '  Parallel mechanical: thermal_threads=', &
			n_thermal_threads, '  mech_threads=', n_mech_threads
	endif

	itertot=0
	step_idx=0
	timet=small

	! Full-domain solve indices (always global)
	ilo = 2; ihi = nim1; jlo = 2; jhi = njm1; klo = 2; khi = nkm1

!------time stepping loop------------------------------------
	time_loop: do while (timet.lt.timax)
		timet=timet+delt
		step_idx = step_idx + 1
		niter=0
		call StartStepTime

!-------Move laser and calculate freesurface-----------------
		call cpu_time(t0)
		call laser_beam
		call read_coordinates
		call cpu_time(t1)
		t_laser = t_laser + (t1 - t0)

		if (adaptive_flag == 1) then
			call amr_check_remesh(step_idx)
			if (amr_needs_remesh) then
				call cpu_time(t0)
				call amr_regenerate_grid()
				call update_thermal_history_indices()
				call defect_update_map()
				if (mechanical_flag == 1 .and. .not. mech_parallel) call update_mech_grid()
				call cpu_time(t1)
				t_amr = t_amr + (t1 - t0)
				n_amr_remesh = n_amr_remesh + 1
				amr_needs_remesh = .false.
			endif
		endif

!-----pool size from previous timestep (once per timestep, before iter_loop)-----
		call cpu_time(t0)
		call pool_size(ilo, ihi, jlo, jhi, klo, khi)
		call cpu_time(t1)
		t_dimen = t_dimen + (t1 - t0)

!-----integer-cell field prediction for heating steps (no interpolation)-----
!     Extend shift region beyond melt pool in scan direction to cover where laser is moving
		if (predict_flag == 1 .and. toolmatrix(PathNum,5) .ge. laser_on_threshold .and. tpeak .gt. tsolid) then
			i = max(nint(abs(scanvelx) * delt / amr_dx_fine), 3)
			j = max(nint(abs(scanvely) * delt / amr_dx_fine), 3)
			call predict_shift_integer(scanvelx, scanvely, delt, &
				max(istat-i, 2), min(iend+i, nim1), &
				max(jstat-j, 2), min(jend+j, njm1), &
				kstat, nkm1)
			call enthalpy_to_temp( &
				max(istat-i, 2), min(iend+i, nim1), &
				max(jstat-j, 2), min(jend+j, njm1), &
				kstat, nkm1)
		endif

!-----iteration loop within each time step----------------
		iter_loop: do while (niter.lt.maxit)
			niter=niter+1
			itertot=itertot+1

!-----solve energy equation (formerly ivar=5)-----
			call cpu_time(t0)
			call properties(ilo, ihi, jlo, jhi, klo, khi)
			call cpu_time(t1)
			t_prop = t_prop + (t1 - t0)

			call cpu_time(t0)
			call bound_enthalpy(ilo, ihi, jlo, jhi, klo, khi)
			call cpu_time(t1)
			t_bound = t_bound + (t1 - t0)

			call cpu_time(t0)
			call discretize_enthalpy(ilo, ihi, jlo, jhi, klo, khi)
			call cpu_time(t1)
			t_discret = t_discret + (t1 - t0)
			t_discret_enthalpy = t_discret_enthalpy + (t1 - t0)

			call cpu_time(t0)
			call source_enthalpy(ilo, ihi, jlo, jhi, klo, khi)
			call cpu_time(t1)
			t_sour = t_sour + (t1 - t0)
			t_sour_enthalpy = t_sour_enthalpy + (t1 - t0)

			call cpu_time(t0)
			call calc_enthalpy_residual(ilo, ihi, jlo, jhi, klo, khi)
			call cpu_time(t1)
			t_resid = t_resid + (t1 - t0)

			call cpu_time(t0)
			call enhance_converge_speed(ilo, ihi, jlo, jhi, klo, khi)
			call cpu_time(t1)
			t_converge = t_converge + (t1 - t0)

			call cpu_time(t0)
			call solution_enthalpy(ilo, ihi, jlo, jhi, klo, khi)
			call cpu_time(t1)
			t_solve = t_solve + (t1 - t0)
			t_solve_enthalpy = t_solve_enthalpy + (t1 - t0)

			call cpu_time(t0)
			call enthalpy_to_temp(ilo, ihi, jlo, jhi, klo, khi)
			call cpu_time(t1)
			t_entot = t_entot + (t1 - t0)

			if(tpeak.gt.merge(min(tsolid,tsolid2), tsolid, species_flag==1)) then
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

!-----convergence criterion------------
			call cpu_time(t0)
			call heat_fluxes
			call cpu_time(t1)
			t_flux = t_flux + (t1 - t0)

			amaxres=max(resorm, resoru,resorv,resorw)

			if(toolmatrix(PathNum,5) .ge. laser_on_threshold)then
				! Laser on: transient-state criteria (heating stage)
				if(resorh.lt.conv_res_heat .and. ratio.le.ratio_upper .and. ratio.ge.ratio_lower) exit iter_loop
			else
				! Laser off: transient-state criteria (cooling stage)
				if(resorh.lt.conv_res_cool) exit iter_loop
			endif

		end do iter_loop

!-----species transport (once per timestep, after iter_loop)-----
		if (species_flag == 1) then
			call cpu_time(t0)
			call species_bc
			call solve_species
			call cpu_time(t1)
			t_species = t_species + (t1 - t0)
		endif

		call cpu_time(t0)
		call update_max_temp()
		call cpu_time(t1)
		t_defect = t_defect + (t1 - t0)

!-----mechanical solver (every mech_interval steps)-----
		if (mechanical_flag == 1 .and. mod(step_idx, mech_interval) == 0) then
			if (mech_parallel) then
				! Parallel mode: write binary input file, mechanical section picks it up
				call write_mech_input(step_idx, timet, temp, solidfield)
			else
				! Serial mode: solve in-loop
				call cpu_time(t_mech_cpu0)
				call solve_mechanical(temp, solidfield, mech_res, mech_newton_iters, mech_cg_iters)
				call get_stress_yield(sxx_out, syy_out, szz_out, vm_out, fplus_out)
				mech_solve_count = mech_solve_count + 1
				call cpu_time(t_mech_cpu1)
				t_mech = t_mech + (t_mech_cpu1 - t_mech_cpu0)
				mio_t_cpu = mio_t_cpu + (t_mech_cpu1 - t_mech_cpu0)
				mio_n_solves = mio_n_solves + 1

				n_yield = count(fplus_out > 0.0_wp)
				write(9,'(A,I6,A,es10.3,A,es10.3,A,I8)') &
					'  Mech step', mech_solve_count, &
					'  res=', mech_res, '  max_vm=', maxval(vm_out), &
					'  yield_elems=', n_yield

				if (mod(mech_solve_count, mech_output_interval) == 0) then
					call write_mech_vtk(step_idx, T_fem_last, &
						ux_mech, uy_mech, uz_mech, mech_phase, &
						sxx_out, syy_out, szz_out, vm_out, fplus_out, &
						Nnx, Nny, Nnz)
				endif

				call write_mech_history(timet, T_fem_last, &
					ux_mech, uy_mech, uz_mech, sxx_out, syy_out, &
					Nnx, Nny, Nnz)
			endif
		endif

		call cpu_time(t0)
		call CalTime
		call outputres
		call cpu_time(t1)
		t_print = t_print + (t1 - t0)

		! Accumulate heating/cooling wall-clock time
		t_step_wall = omp_get_wtime() - t_step_start
		if (toolmatrix(PathNum,5) .ge. laser_on_threshold) then
			t_heating = t_heating + t_step_wall
			n_heating = n_heating + 1
		else
			t_cooling = t_cooling + t_step_wall
			n_cooling = n_cooling + 1
		endif

		call cpu_time(t0)
		! Global step: update all cells
		!$OMP PARALLEL DO PRIVATE(i,j,k)
		do k=1,nk
		do j=1,nj
		do i=1,ni
			if(temp(i,j,k).le.merge(mix(tsolid,tsolid2,concentration(i,j,k)), tsolid, species_flag==1)) then
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
		!$OMP END PARALLEL DO
		if (species_flag == 1) conc_old = concentration
		call Cust_Out
		call write_thermal_history(timet)
		call write_meltpool_history(timet)
		call cpu_time(t1)
		t_other = t_other + (t1 - t0)

	end do time_loop

	! No mech_DONE sentinel needed — mechanical knows the last step from timax/delt

	! Post-simulation analysis (before EndTime closes output file)
	call compute_defect_determ()
	call write_defect_report()

	call EndTime
	call finalize_thermal_history
	call finalize_meltpool_history

	if (mechanical_flag == 1 .and. .not. mech_parallel) then
		! Serial mode: finalize here
		mio_t_wall = t_mech  ! approximate wall time from CPU time
		call write_mech_timing_report(file_prefix)
		call write_mech_memory_report(file_prefix, Nx, Ny, Nz, Nnx, Nny, Nnz)
		call finalize_mech_history()
		call finalize_mechanical_io()
		call cleanup_mechanical()
	endif

	wall_elapsed = omp_get_wtime() - wall_start

	call write_timing_report(itertot, timet, wall_elapsed, file_prefix)
	call write_memory_report(file_prefix)

	! In parallel mode, mechanical runs as a separate process (launched by run.sh).
	! run.sh handles waiting for both processes to complete.

	stop
	end

!****************************************************************************
subroutine run_mechanical_loop()
! Parallel mechanical solver loop: polls for binary input files,
! solves mechanical, writes VTK/history output.
	use precision
	use geometry, only: ni, nj, nk, x, y
	use parameters
	use omp_lib
	use mechanical_solver
	use mech_io
	implicit none

	real(wp), allocatable :: temp_buf(:,:,:), sf_buf(:,:,:)
	real(wp), allocatable :: x_buf(:), y_buf(:), x_prev(:), y_prev(:)
	integer  :: step_out, mech_count, n_yield, next_step, last_mech_step, nsteps_total
	real(wp) :: time_out, mech_res
	integer  :: mech_newton_iters, mech_cg_iters
	real(wp) :: t0_cpu, t1_cpu, wall_mech_start
	logical  :: found, grid_changed

	! Set thread count for mechanical section
	call omp_set_num_threads(n_mech_threads)
	wall_mech_start = omp_get_wtime()

	allocate(temp_buf(ni,nj,nk), sf_buf(ni,nj,nk))
	allocate(x_buf(ni), y_buf(nj))
	allocate(x_prev(ni), y_prev(nj))

	call init_mechanical()
	call init_mech_history(Nnx, Nny, Nnz)

	x_prev = x(1:ni)
	y_prev = y(1:nj)

	! Compute the last step that will have a mech input file
	nsteps_total = nint(timax / delt)
	last_mech_step = (nsteps_total / mech_interval) * mech_interval
	write(*,'(A,I6,A,I6)') '  [MECH] Expecting files from step', mech_interval, &
		' to', last_mech_step

	mech_count = 0
	next_step = mech_interval

	do
		! All expected files processed — exit
		if (next_step > last_mech_step) exit

		! Poll for next input file
		call read_mech_input(next_step, step_out, time_out, &
			temp_buf, sf_buf, x_buf, y_buf, found)

		if (.not. found) then
			call sleep(1)
			cycle
		endif

		! Check if grid changed (AMR)
		grid_changed = any(abs(x_buf(1:ni) - x_prev(1:ni)) > 1.0e-15_wp) .or. &
		               any(abs(y_buf(1:nj) - y_prev(1:nj)) > 1.0e-15_wp)
		if (grid_changed) then
			! Update geometry arrays so update_mech_grid can read them
			x(1:ni) = x_buf(1:ni)
			y(1:nj) = y_buf(1:nj)
			call update_mech_grid()
			x_prev = x_buf
			y_prev = y_buf
		endif

		! Solve mechanical
		call cpu_time(t0_cpu)
		call solve_mechanical(temp_buf, sf_buf, mech_res, mech_newton_iters, mech_cg_iters)
		call get_stress_yield(sxx_out, syy_out, szz_out, vm_out, fplus_out)
		mech_count = mech_count + 1
		call cpu_time(t1_cpu)
		mio_t_cpu = mio_t_cpu + (t1_cpu - t0_cpu)
		mio_n_solves = mio_n_solves + 1

		n_yield = count(fplus_out > 0.0_wp)
		write(*,'(A,I6,A,I6,A,es10.3,A,es10.3,A,I8)') &
			'  [MECH] step', step_out, '  solve', mech_count, &
			'  res=', mech_res, '  max_vm=', maxval(vm_out), &
			'  yield=', n_yield

		! VTK output
		if (mod(mech_count, mech_output_interval) == 0) then
			call write_mech_vtk(step_out, T_fem_last, &
				ux_mech, uy_mech, uz_mech, mech_phase, &
				sxx_out, syy_out, szz_out, vm_out, fplus_out, &
				Nnx, Nny, Nnz)
		endif

		! History
		call write_mech_history(time_out, T_fem_last, &
			ux_mech, uy_mech, uz_mech, sxx_out, syy_out, &
			Nnx, Nny, Nnz)

		next_step = next_step + mech_interval
	enddo

	! Finalize
	mio_t_wall = omp_get_wtime() - wall_mech_start
	call write_mech_timing_report(file_prefix)
	call write_mech_memory_report(file_prefix, Nx, Ny, Nz, Nnx, Nny, Nnz)
	call finalize_mech_history()
	call finalize_mechanical_io()
	call cleanup_mechanical()

	deallocate(temp_buf, sf_buf, x_buf, y_buf, x_prev, y_prev)

end subroutine run_mechanical_loop
