!______________________________________________________________________________
!
module crack_risk_mod
!______________________________________________________________________________
! Crack risk prediction module: computes crack susceptibility index (CSI)
! from thermal strain rate and time spent in the Brittle Temperature Range.
! CSI = accumulated thermal strain in BTR = integral of alpha * |dT/dt| dt.
!
! Input namelist (in input_param.txt):
!   &crack_params delta_t_btr=100.0 /
!
!   delta_t_btr  (K) - BTR width below T_solidus (default: 100 K for IN718)
!                      BTR = [T_solidus - delta_t_btr, T_solidus]
!                      Material has near-zero ductility in this range;
!                      thermal strain accumulated here drives cracking.
!
! Output:
!   CSI = integral of beta * |dT/dt| * dt over time spent in BTR
!   where beta = thermal expansion coefficient (from &material_properties)
!
	use precision
	use geometry
	use initialization
	use parameters
	use sim_state
	use constant
	use defect_field, only: k_def_lo, k_def_hi, point_in_scan_region, compute_scan_range, n_hull

	implicit none

	! --- Input parameter read from mod_param.f90 (delta_t_btr) ---

	! --- Arrays (full X-Y, limited Z range) ---
	real(wp), allocatable :: cool_rate_solid(:,:,:)    ! cooling rate at solidification (K/s)
	real(wp), allocatable :: strain_rate_solid(:,:,:)  ! thermal strain rate at solidification (1/s)
	real(wp), allocatable :: btr_time(:,:,:)           ! cumulative time in BTR (s)
	real(wp), allocatable :: crack_risk_arr(:,:,:)     ! CSI: accumulated thermal strain in BTR (dimensionless)
	logical, allocatable  :: crack_solidified(:,:,:)   ! has this cell solidified?

	contains

!********************************************************************
subroutine allocate_crack_risk(nni, nnj, nnk)
	integer, intent(in) :: nni, nnj, nnk

	allocate(cool_rate_solid(nni, nnj, nnk))
	allocate(strain_rate_solid(nni, nnj, nnk))
	allocate(btr_time(nni, nnj, nnk))
	allocate(crack_risk_arr(nni, nnj, nnk))
	allocate(crack_solidified(nni, nnj, nnk))

	cool_rate_solid   = 0.0_wp
	strain_rate_solid = 0.0_wp
	btr_time          = 0.0_wp
	crack_risk_arr    = 0.0_wp
	crack_solidified  = .false.

end subroutine allocate_crack_risk

!********************************************************************
subroutine update_crack_risk(dt)
! Called every timestep after enthalpy solve: track solidification and BTR.
	real(wp), intent(in) :: dt
	integer :: i, j, k
	real(wp) :: dTdt, abs_dTdt

	!$OMP PARALLEL DO PRIVATE(i, j, k, dTdt, abs_dTdt)
	do k = k_def_lo, k_def_hi
	do j = 2, njm1
	do i = 2, nim1
		dTdt = (temp(i,j,k) - tnot(i,j,k)) / dt
		abs_dTdt = abs(dTdt)

		! Solidification detection: fracl transitions from >0 to <=0
		if (fraclnot(i,j,k) > 0.0_wp .and. fracl(i,j,k) <= 0.0_wp &
		    .and. .not. crack_solidified(i,j,k)) then
			cool_rate_solid(i,j,k) = abs_dTdt
			strain_rate_solid(i,j,k) = beta * abs_dTdt
			crack_solidified(i,j,k) = .true.
		endif

		! BTR tracking: temp in [tsolid - delta_t_btr, tsolid)
		if (temp(i,j,k) < tsolid .and. temp(i,j,k) > tsolid - delta_t_btr) then
			btr_time(i,j,k) = btr_time(i,j,k) + dt
			crack_risk_arr(i,j,k) = crack_risk_arr(i,j,k) + beta * abs_dTdt * dt
		endif
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO

end subroutine update_crack_risk

!********************************************************************
subroutine compute_crack_report()
! Post-simulation: compute crack risk statistics and write report + VTK files.
	integer :: i, j, k
	integer, parameter :: lun = 91
	integer :: n_solidified, n_highrisk
	real(wp) :: csi_max, csi_sum, csi_mean
	real(wp) :: cool_max, strain_max
	real(wp) :: frac_highrisk
	real(wp), parameter :: csi_threshold = 0.01_wp   ! 1% strain

	! Ensure scan range is computed
	if (n_hull == 0) call compute_scan_range()

	! --- Compute statistics over solidified cells in scan region ---
	n_solidified = 0
	n_highrisk   = 0
	csi_max      = 0.0_wp
	csi_sum      = 0.0_wp
	cool_max     = 0.0_wp
	strain_max   = 0.0_wp

	do k = k_def_lo, k_def_hi
	do j = 2, njm1
	do i = 2, nim1
		if (.not. point_in_scan_region(x(i), y(j))) cycle
		if (.not. crack_solidified(i,j,k)) cycle

		n_solidified = n_solidified + 1
		csi_sum = csi_sum + crack_risk_arr(i,j,k)
		if (crack_risk_arr(i,j,k) > csi_max) csi_max = crack_risk_arr(i,j,k)
		if (cool_rate_solid(i,j,k) > cool_max) cool_max = cool_rate_solid(i,j,k)
		if (strain_rate_solid(i,j,k) > strain_max) strain_max = strain_rate_solid(i,j,k)
		if (crack_risk_arr(i,j,k) > csi_threshold) n_highrisk = n_highrisk + 1
	enddo
	enddo
	enddo

	if (n_solidified > 0) then
		csi_mean = csi_sum / real(n_solidified, wp)
		frac_highrisk = real(n_highrisk, wp) / real(n_solidified, wp)
	else
		csi_mean = 0.0_wp
		frac_highrisk = 0.0_wp
	endif

	! --- Write crack_report.txt ---
	open(unit=lun, file=trim(file_prefix)//'crack_report.txt', action='write', status='replace')
	write(lun,'(a)') '============================================'
	write(lun,'(a)') '  PHOENIX Crack Risk Report'
	write(lun,'(a)') '============================================'
	write(lun,'(a)')
	write(lun,'(a,i0)')      '  Solidified cells in scan region: ', n_solidified
	write(lun,'(a)')
	write(lun,'(a,es12.5)')  '  Max CSI (crack susceptibility):  ', csi_max
	write(lun,'(a,es12.5)')  '  Mean CSI:                        ', csi_mean
	write(lun,'(a)')
	write(lun,'(a,es12.5,a)') '  Max cooling rate at solidification: ', cool_max, ' K/s'
	write(lun,'(a,es12.5,a)') '  Max strain rate at solidification:  ', strain_max, ' 1/s'
	write(lun,'(a)')
	write(lun,'(a,i0)')      '  High-risk cells (CSI > 0.01):    ', n_highrisk
	write(lun,'(a,f10.6,a)') '  High-risk fraction:              ', frac_highrisk * 100.0_wp, ' %'
	write(lun,'(a)')
	write(lun,'(a)') '  Parameters:'
	write(lun,'(a,f10.2,a)') '    delta_T_BTR = ', delta_t_btr, ' K'
	write(lun,'(a,es12.5)')  '    beta        = ', beta
	write(lun,'(a,f10.2,a)') '    T_solidus   = ', tsolid, ' K'
	write(lun,'(a)')
	write(lun,'(a)') '  Layer Z range:'
	write(lun,'(a,i4,a,i4)') '    k indices: ', k_def_lo, ' to ', k_def_hi
	write(lun,'(a,es12.5,a,es12.5,a)') '    Z: [', z(k_def_lo), ', ', z(k_def_hi), '] m'
	write(lun,'(a)') '============================================'
	close(lun)

	! --- Print summary to output file ---
	write(9,'(a)') ''
	write(9,'(a)') '  === Crack Risk Analysis ==='
	write(9,'(a,es12.5)')    '  Max CSI:                         ', csi_max
	write(9,'(a,es12.5)')    '  Mean CSI:                        ', csi_mean
	write(9,'(a,es12.5,a)')  '  Max cooling rate at solidification: ', cool_max, ' K/s'
	write(9,'(a,es12.5,a)')  '  Max strain rate at solidification:  ', strain_max, ' 1/s'
	write(9,'(a,f10.6,a)')   '  High-risk fraction:              ', frac_highrisk * 100.0_wp, ' %'

	! --- Write combined VTK file ---
	call write_crack_vtk_combined()

end subroutine compute_crack_report

!********************************************************************
subroutine write_crack_vtk_combined()
! Write a single VTK file with all crack risk fields as separate scalars.
	integer :: i, j, k, npts
	integer :: gridx, gridy, gridz
	real(kind=4) :: val4
	integer, parameter :: lun = 92

	gridx = nim1 - 2 + 1
	gridy = njm1 - 2 + 1
	gridz = k_def_hi - k_def_lo + 1
	npts = gridx * gridy * gridz

	! ASCII header
	open(unit=lun, file=trim(file_prefix)//'crack_risk.vtk')
	write(lun,'(A)') '# vtk DataFile Version 3.0'
	write(lun,'(A)') 'PHOENIX crack risk fields'
	write(lun,'(A)') 'BINARY'
	write(lun,'(A)') 'DATASET STRUCTURED_GRID'
	write(lun,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', gridx, ' ', gridy, ' ', gridz
	write(lun,'(A,I0,A)') 'POINTS ', npts, ' float'
	close(lun)

	! Binary coordinates
	open(unit=lun, file=trim(file_prefix)//'crack_risk.vtk', &
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
	open(unit=lun, file=trim(file_prefix)//'crack_risk.vtk', position='append')
	write(lun,'(A,I0)') 'POINT_DATA ', npts
	close(lun)

	! Write each scalar field
	call write_crack_scalar(lun, 'crack_csi', crack_risk_arr)
	call write_crack_scalar(lun, 'cooling_rate_solid', cool_rate_solid)
	call write_crack_scalar(lun, 'strain_rate_solid', strain_rate_solid)
	call write_crack_scalar(lun, 'btr_time', btr_time)

end subroutine write_crack_vtk_combined

!********************************************************************
subroutine write_crack_scalar(lun, name, field)
	integer, intent(in) :: lun
	character(len=*), intent(in) :: name
	real(wp), intent(in) :: field(:,:,:)
	integer :: i, j, k
	real(kind=4) :: val4

	! ASCII scalar header
	open(unit=lun, file=trim(file_prefix)//'crack_risk.vtk', position='append')
	write(lun,'(A)') 'SCALARS '//trim(name)//' float 1'
	write(lun,'(A)') 'LOOKUP_TABLE default'
	close(lun)

	! Binary data
	open(unit=lun, file=trim(file_prefix)//'crack_risk.vtk', &
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

end subroutine write_crack_scalar

end module crack_risk_mod
