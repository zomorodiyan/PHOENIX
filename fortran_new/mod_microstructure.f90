!______________________________________________________________________________
!
module microstructure_mod
!______________________________________________________________________________
! Solidification microstructure prediction module.
! Computes thermal gradient G, cooling rate, solidification rate R,
! primary dendrite arm spacing (PDAS), and secondary dendrite arm spacing (SDAS)
! at the moment each cell solidifies (fracl transitions from >0 to 0).
!
! Input namelist (in input_param.txt):
!   &microstructure_params a1_pdas=50e-6, a2_sdas=10e-6, n1_pdas=-0.5, n2_pdas=-0.25, n3_sdas=-0.333 /
!
!   a1_pdas  (m)  - PDAS prefactor, alloy-specific (default: 50 um for IN718)
!   a2_sdas  (m)  - SDAS prefactor, alloy-specific (default: 10 um for IN718)
!   n1_pdas  (-)  - PDAS exponent for thermal gradient G (Kurz-Fisher: -0.5)
!   n2_pdas  (-)  - PDAS exponent for solidification rate R (Kurz-Fisher: -0.25)
!   n3_sdas  (-)  - SDAS exponent for cooling rate (Kattamis-Flemings: -1/3)
!
! Models:
!   PDAS: lambda_1 = a1_pdas * G^n1_pdas * R^n2_pdas
!   SDAS: lambda_2 = a2_sdas * |dT/dt|^n3_sdas
!
	use precision
	use geometry
	use initialization
	use parameters
	use sim_state
	use constant
	use defect_field, only: k_def_lo, k_def_hi, point_in_scan_region, &
		compute_scan_range, n_hull

	implicit none

	! --- Model parameters read from mod_param.f90 (a1_pdas, a2_sdas, n1_pdas, n2_pdas, n3_sdas) ---

	! --- Arrays (full X-Y, limited Z range active) ---
	real(wp), allocatable :: cool_rate_micro(:,:,:)   ! cooling rate at solidification (K/s)
	real(wp), allocatable :: therm_grad(:,:,:)        ! thermal gradient magnitude (K/m)
	real(wp), allocatable :: solid_rate(:,:,:)        ! solidification rate R (m/s)
	real(wp), allocatable :: pdas_arr(:,:,:)          ! primary dendrite arm spacing (m)
	real(wp), allocatable :: sdas_arr(:,:,:)          ! secondary dendrite arm spacing (m)
	logical, allocatable  :: micro_solidified(:,:,:)  ! solidification event mask

	contains

!********************************************************************
subroutine allocate_microstructure(nni, nnj, nnk)
	integer, intent(in) :: nni, nnj, nnk

	allocate(cool_rate_micro(nni, nnj, nnk))
	allocate(therm_grad(nni, nnj, nnk))
	allocate(solid_rate(nni, nnj, nnk))
	allocate(pdas_arr(nni, nnj, nnk))
	allocate(sdas_arr(nni, nnj, nnk))
	allocate(micro_solidified(nni, nnj, nnk))

	cool_rate_micro = 0.0_wp
	therm_grad      = 0.0_wp
	solid_rate      = 0.0_wp
	pdas_arr        = 0.0_wp
	sdas_arr        = 0.0_wp
	micro_solidified = .false.

end subroutine allocate_microstructure

!********************************************************************
subroutine update_microstructure(dt)
! Called each timestep after enthalpy solve.
! Detects solidification events and stores G, dT/dt, R, PDAS, SDAS.
	real(wp), intent(in) :: dt
	integer :: i, j, k
	real(wp) :: dTdt, dTdx, dTdy, dTdz, G, R, abs_dTdt

	!$OMP PARALLEL DO PRIVATE(i, j, k, dTdt, dTdx, dTdy, dTdz, G, R, abs_dTdt)
	do k = k_def_lo, k_def_hi
	do j = 2, njm1
	do i = 2, nim1
		! Skip cells that have already solidified
		if (micro_solidified(i,j,k)) cycle

		! Detect solidification: liquid fraction was >0 last step, now <=0
		if (fraclnot(i,j,k) > 0.0_wp .and. fracl(i,j,k) <= 0.0_wp) then

			! Cooling rate
			dTdt = (temp(i,j,k) - tnot(i,j,k)) / dt
			abs_dTdt = abs(dTdt)

			! Thermal gradient via central differences (one-sided at boundaries)
			! dT/dx
			if (i == 2) then
				dTdx = (temp(i+1,j,k) - temp(i,j,k)) / (x(i+1) - x(i))
			else if (i == nim1) then
				dTdx = (temp(i,j,k) - temp(i-1,j,k)) / (x(i) - x(i-1))
			else
				dTdx = (temp(i+1,j,k) - temp(i-1,j,k)) / (x(i+1) - x(i-1))
			endif

			! dT/dy
			if (j == 2) then
				dTdy = (temp(i,j+1,k) - temp(i,j,k)) / (y(j+1) - y(j))
			else if (j == njm1) then
				dTdy = (temp(i,j,k) - temp(i,j-1,k)) / (y(j) - y(j-1))
			else
				dTdy = (temp(i,j+1,k) - temp(i,j-1,k)) / (y(j+1) - y(j-1))
			endif

			! dT/dz
			if (k == k_def_lo) then
				dTdz = (temp(i,j,k+1) - temp(i,j,k)) / (z(k+1) - z(k))
			else if (k == k_def_hi) then
				dTdz = (temp(i,j,k) - temp(i,j,k-1)) / (z(k) - z(k-1))
			else
				dTdz = (temp(i,j,k+1) - temp(i,j,k-1)) / (z(k+1) - z(k-1))
			endif

			G = sqrt(dTdx**2 + dTdy**2 + dTdz**2)

			! Store cooling rate and thermal gradient
			cool_rate_micro(i,j,k) = abs_dTdt
			therm_grad(i,j,k) = G

			! Solidification rate R = |dT/dt| / G
			R = abs_dTdt / max(G, 1.0e-10_wp)
			solid_rate(i,j,k) = R

			! PDAS: lambda_1 = a1 * G^n1 * R^n2
			if (G > 0.0_wp .and. R > 0.0_wp) then
				pdas_arr(i,j,k) = a1_pdas * G**n1_pdas * R**n2_pdas
			else
				pdas_arr(i,j,k) = 0.0_wp
			endif

			! SDAS: lambda_2 = a2 * |dT/dt|^n3
			if (abs_dTdt > 0.0_wp) then
				sdas_arr(i,j,k) = a2_sdas * abs_dTdt**n3_sdas
			else
				sdas_arr(i,j,k) = 0.0_wp
			endif

			micro_solidified(i,j,k) = .true.
		endif
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO

end subroutine update_microstructure

!********************************************************************
subroutine report_microstructure()
! Post-simulation: compute statistics and write report + VTK files.
	integer :: i, j, k, n_solid
	integer, parameter :: lun = 94
	real(wp) :: g_min, g_max, g_sum
	real(wp) :: r_min, r_max, r_sum
	real(wp) :: cr_min, cr_max, cr_sum
	real(wp) :: pdas_min, pdas_max, pdas_sum
	real(wp) :: sdas_min, sdas_max, sdas_sum

	! Build scan region polygon if not already done
	if (n_hull == 0) call compute_scan_range()

	! Initialize statistics
	n_solid  = 0
	g_min    = great;  g_max    = 0.0_wp;  g_sum    = 0.0_wp
	r_min    = great;  r_max    = 0.0_wp;  r_sum    = 0.0_wp
	cr_min   = great;  cr_max   = 0.0_wp;  cr_sum   = 0.0_wp
	pdas_min = great;  pdas_max = 0.0_wp;  pdas_sum = 0.0_wp
	sdas_min = great;  sdas_max = 0.0_wp;  sdas_sum = 0.0_wp

	do k = k_def_lo, k_def_hi
	do j = 2, njm1
	do i = 2, nim1
		if (.not. micro_solidified(i,j,k)) cycle
		if (.not. point_in_scan_region(x(i), y(j))) cycle

		n_solid = n_solid + 1

		! Thermal gradient
		if (therm_grad(i,j,k) < g_min) g_min = therm_grad(i,j,k)
		if (therm_grad(i,j,k) > g_max) g_max = therm_grad(i,j,k)
		g_sum = g_sum + therm_grad(i,j,k)

		! Solidification rate
		if (solid_rate(i,j,k) < r_min) r_min = solid_rate(i,j,k)
		if (solid_rate(i,j,k) > r_max) r_max = solid_rate(i,j,k)
		r_sum = r_sum + solid_rate(i,j,k)

		! Cooling rate
		if (cool_rate_micro(i,j,k) < cr_min) cr_min = cool_rate_micro(i,j,k)
		if (cool_rate_micro(i,j,k) > cr_max) cr_max = cool_rate_micro(i,j,k)
		cr_sum = cr_sum + cool_rate_micro(i,j,k)

		! PDAS
		if (pdas_arr(i,j,k) < pdas_min) pdas_min = pdas_arr(i,j,k)
		if (pdas_arr(i,j,k) > pdas_max) pdas_max = pdas_arr(i,j,k)
		pdas_sum = pdas_sum + pdas_arr(i,j,k)

		! SDAS
		if (sdas_arr(i,j,k) < sdas_min) sdas_min = sdas_arr(i,j,k)
		if (sdas_arr(i,j,k) > sdas_max) sdas_max = sdas_arr(i,j,k)
		sdas_sum = sdas_sum + sdas_arr(i,j,k)
	enddo
	enddo
	enddo

	! --- Write micro_report.txt ---
	open(unit=lun, file=trim(file_prefix)//'micro_report.txt', action='write', status='replace')
	write(lun,'(a)') '============================================'
	write(lun,'(a)') '  PHOENIX Microstructure Report'
	write(lun,'(a)') '============================================'
	write(lun,'(a)')
	write(lun,'(a,i0)') '  Solidified cells in scan region: ', n_solid
	write(lun,'(a)')

	if (n_solid > 0) then
		write(lun,'(a)') '  Thermal Gradient G (K/m):'
		write(lun,'(a,es12.5)') '    Min:  ', g_min
		write(lun,'(a,es12.5)') '    Max:  ', g_max
		write(lun,'(a,es12.5)') '    Mean: ', g_sum / real(n_solid, wp)
		write(lun,'(a)')
		write(lun,'(a)') '  Solidification Rate R (m/s):'
		write(lun,'(a,es12.5)') '    Min:  ', r_min
		write(lun,'(a,es12.5)') '    Max:  ', r_max
		write(lun,'(a,es12.5)') '    Mean: ', r_sum / real(n_solid, wp)
		write(lun,'(a)')
		write(lun,'(a)') '  Cooling Rate |dT/dt| (K/s):'
		write(lun,'(a,es12.5)') '    Min:  ', cr_min
		write(lun,'(a,es12.5)') '    Max:  ', cr_max
		write(lun,'(a,es12.5)') '    Mean: ', cr_sum / real(n_solid, wp)
		write(lun,'(a)')
		write(lun,'(a)') '  Primary Dendrite Arm Spacing PDAS (m):'
		write(lun,'(a,es12.5)') '    Min:  ', pdas_min
		write(lun,'(a,es12.5)') '    Max:  ', pdas_max
		write(lun,'(a,es12.5)') '    Mean: ', pdas_sum / real(n_solid, wp)
		write(lun,'(a)')
		write(lun,'(a)') '  Secondary Dendrite Arm Spacing SDAS (m):'
		write(lun,'(a,es12.5)') '    Min:  ', sdas_min
		write(lun,'(a,es12.5)') '    Max:  ', sdas_max
		write(lun,'(a,es12.5)') '    Mean: ', sdas_sum / real(n_solid, wp)
	else
		write(lun,'(a)') '  No solidified cells found in scan region.'
	endif

	write(lun,'(a)')
	write(lun,'(a)') '  Model Parameters:'
	write(lun,'(a,es12.5)') '    a1_pdas = ', a1_pdas
	write(lun,'(a,es12.5)') '    n1_pdas = ', n1_pdas
	write(lun,'(a,es12.5)') '    n2_pdas = ', n2_pdas
	write(lun,'(a,es12.5)') '    a2_sdas = ', a2_sdas
	write(lun,'(a,es12.5)') '    n3_sdas = ', n3_sdas
	write(lun,'(a)')
	write(lun,'(a)') '  Layer Z range:'
	write(lun,'(a,i4,a,i4)') '    k indices: ', k_def_lo, ' to ', k_def_hi
	write(lun,'(a,es12.5,a,es12.5,a)') '    Z: [', z(k_def_lo), ', ', z(k_def_hi), '] m'
	write(lun,'(a)') '============================================'
	close(lun)

	! Print summary to output file
	write(9,'(a)') ''
	write(9,'(a)') '  === Microstructure Analysis ==='
	write(9,'(a,i0)') '  Solidified cells: ', n_solid
	if (n_solid > 0) then
		write(9,'(a,es12.5,a,es12.5)') '  G (K/m):   min=', g_min, '  max=', g_max
		write(9,'(a,es12.5,a,es12.5)') '  R (m/s):   min=', r_min, '  max=', r_max
		write(9,'(a,es12.5,a,es12.5)') '  PDAS (m):  min=', pdas_min, '  max=', pdas_max
		write(9,'(a,es12.5,a,es12.5)') '  SDAS (m):  min=', sdas_min, '  max=', sdas_max
	endif

	! --- Write combined VTK file ---
	call write_micro_vtk_combined()

end subroutine report_microstructure

!********************************************************************
subroutine write_micro_vtk_combined()
! Write a single VTK file with all microstructure fields as separate scalars.
	integer :: i, j, k, npts
	integer :: gridx, gridy, gridz
	real(kind=4) :: val4
	integer, parameter :: lun = 93

	gridx = nim1 - 2 + 1
	gridy = njm1 - 2 + 1
	gridz = k_def_hi - k_def_lo + 1
	npts = gridx * gridy * gridz

	! ASCII header
	open(unit=lun, file=trim(file_prefix)//'microstructure.vtk')
	write(lun,'(A)') '# vtk DataFile Version 3.0'
	write(lun,'(A)') 'PHOENIX microstructure fields'
	write(lun,'(A)') 'BINARY'
	write(lun,'(A)') 'DATASET STRUCTURED_GRID'
	write(lun,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', gridx, ' ', gridy, ' ', gridz
	write(lun,'(A,I0,A)') 'POINTS ', npts, ' float'
	close(lun)

	! Binary coordinates
	open(unit=lun, file=trim(file_prefix)//'microstructure.vtk', &
	     access='stream', form='unformatted', position='append', convert='big_endian')
	do k = k_def_lo, k_def_hi
	do j = 2, njm1
	do i = 2, nim1
		val4 = real(x(i), 4); write(lun) val4
		val4 = real(y(j), 4); write(lun) val4
		val4 = real(z(k), 4); write(lun) val4
	enddo
	enddo
	enddo
	close(lun)

	! POINT_DATA header
	open(unit=lun, file=trim(file_prefix)//'microstructure.vtk', position='append')
	write(lun,'(A,I0)') 'POINT_DATA ', npts
	close(lun)

	! Write each scalar field
	call write_micro_scalar(lun, 'cooling_rate', cool_rate_micro)
	call write_micro_scalar(lun, 'thermal_gradient', therm_grad)
	call write_micro_scalar(lun, 'solidification_rate', solid_rate)
	call write_micro_scalar(lun, 'PDAS', pdas_arr)
	call write_micro_scalar(lun, 'SDAS', sdas_arr)

end subroutine write_micro_vtk_combined

!********************************************************************
subroutine write_micro_scalar(lun, name, field)
	integer, intent(in) :: lun
	character(len=*), intent(in) :: name
	real(wp), intent(in) :: field(:,:,:)
	integer :: i, j, k
	real(kind=4) :: val4

	! ASCII scalar header
	open(unit=lun, file=trim(file_prefix)//'microstructure.vtk', position='append')
	write(lun,'(A)') 'SCALARS '//trim(name)//' float 1'
	write(lun,'(A)') 'LOOKUP_TABLE default'
	close(lun)

	! Binary data
	open(unit=lun, file=trim(file_prefix)//'microstructure.vtk', &
	     access='stream', form='unformatted', position='append', convert='big_endian')
	do k = k_def_lo, k_def_hi
	do j = 2, njm1
	do i = 2, nim1
		val4 = real(field(i,j,k), 4)
		write(lun) val4
	enddo
	enddo
	enddo
	close(lun)

end subroutine write_micro_scalar

end module microstructure_mod
