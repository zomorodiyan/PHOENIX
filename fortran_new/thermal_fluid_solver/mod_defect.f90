!______________________________________________________________________________
!
module defect_field
!______________________________________________________________________________
! Defect detection module: computes max_temp and defect arrays over one layer.
! When adaptive_flag==1, uses independent uniform mesh for history tracking.
!
	use precision
	use geometry
	use initialization
	use parameters
	use sim_state
	use constant

	implicit none

	! --- maxtemp_determ parameters ---
	real(wp), parameter :: k_lof = 0.9_wp    ! lack-of-fusion calibration (range 0-1)
	real(wp), parameter :: k_kep = 1.0_wp    ! keyhole porosity scaling

	! --- Arrays (on simulation mesh when adaptive_flag==0, on uniform mesh when ==1) ---
	real(wp), allocatable :: max_temp(:,:,:)
	real(wp), allocatable :: defect_arr(:,:,:)

	! --- Uniform defect mesh (only allocated when adaptive_flag==1) ---
	integer :: def_ni = 0, def_nj = 0
	integer :: def_nim1 = 0, def_njm1 = 0
	real(wp), allocatable :: def_x(:), def_y(:)   ! cell center positions
	real(wp), allocatable :: temp_def(:,:,:)       ! interpolated temp/field buffer
	real(wp), allocatable :: solidfield_def(:,:,:) ! solidID on uniform defect mesh

	! --- Index mapping: AMR mesh -> uniform defect mesh ---
	integer, allocatable  :: def_imap_lo(:), def_imap_hi(:)
	real(wp), allocatable :: def_imap_w(:)
	integer, allocatable  :: def_jmap_lo(:), def_jmap_hi(:)
	real(wp), allocatable :: def_jmap_w(:)

	! --- Z range for defect computation ---
	integer :: k_def_lo, k_def_hi

	! --- Laser scanning range (X-Y) ---
	real(wp) :: x_scan_min, x_scan_max, y_scan_min, y_scan_max

	! --- Convex hull of scanned region ---
	integer :: n_hull = 0
	real(wp) :: hull_x(200), hull_y(200)

	contains

!********************************************************************
subroutine allocate_defect(nni, nnj, nnk)
	integer, intent(in) :: nni, nnj, nnk
	integer :: k, i

	! Determine k range: from top (nkm1) down to z(nk) - layerheight
	k_def_hi = nkm1
	k_def_lo = nkm1
	do k = nkm1, 2, -1
		if (z(k) < z(nk) - layerheight) exit
		k_def_lo = k
	enddo

	if (adaptive_flag == 1) then
		! Uniform defect mesh
		def_ni = nint(dimx / amr_dx_fine) + 2
		def_nj = nint(dimy / amr_dx_fine) + 2
		def_nim1 = def_ni - 1
		def_njm1 = def_nj - 1

		allocate(def_x(def_ni), def_y(def_nj))
		do i = 1, def_ni
			def_x(i) = (real(i, wp) - 1.5_wp) * amr_dx_fine
		enddo
		do i = 1, def_nj
			def_y(i) = (real(i, wp) - 1.5_wp) * amr_dx_fine
		enddo
		! Clamp boundary nodes
		def_x(1) = 0.0_wp
		def_x(def_ni) = dimx
		def_y(1) = 0.0_wp
		def_y(def_nj) = dimy

		allocate(max_temp(def_ni, def_nj, nnk))
		allocate(defect_arr(def_ni, def_nj, nnk))
		allocate(temp_def(def_ni, def_nj, nnk))
		allocate(solidfield_def(def_ni, def_nj, nnk))
		max_temp = tempPreheat
		defect_arr = 0.0_wp
		temp_def = tempPreheat
		solidfield_def = 0.0_wp

		! Allocate mapping arrays
		allocate(def_imap_lo(def_ni), def_imap_hi(def_ni), def_imap_w(def_ni))
		allocate(def_jmap_lo(def_nj), def_jmap_hi(def_nj), def_jmap_w(def_nj))

		! Build initial mapping
		call defect_update_map()
	else
		allocate(max_temp(nni, nnj, nnk))
		allocate(defect_arr(nni, nnj, nnk))
		max_temp = tempPreheat
		defect_arr = 0.0_wp
	endif

end subroutine allocate_defect

!********************************************************************
subroutine defect_update_map()
! Rebuild index mapping from AMR simulation mesh to uniform defect mesh.
	integer :: i, j_lo
	real(wp) :: xn

	j_lo = 1
	do i = 1, def_ni
		xn = def_x(i)
		do while (j_lo < ni - 1 .and. x(j_lo + 1) < xn)
			j_lo = j_lo + 1
		enddo
		def_imap_lo(i) = j_lo
		def_imap_hi(i) = min(j_lo + 1, ni)
		if (def_imap_lo(i) == def_imap_hi(i)) then
			def_imap_w(i) = 0.0_wp
		else
			def_imap_w(i) = (xn - x(def_imap_lo(i))) / (x(def_imap_hi(i)) - x(def_imap_lo(i)))
			def_imap_w(i) = max(0.0_wp, min(1.0_wp, def_imap_w(i)))
		endif
	enddo

	j_lo = 1
	do i = 1, def_nj
		xn = def_y(i)
		do while (j_lo < nj - 1 .and. y(j_lo + 1) < xn)
			j_lo = j_lo + 1
		enddo
		def_jmap_lo(i) = j_lo
		def_jmap_hi(i) = min(j_lo + 1, nj)
		if (def_jmap_lo(i) == def_jmap_hi(i)) then
			def_jmap_w(i) = 0.0_wp
		else
			def_jmap_w(i) = (xn - y(def_jmap_lo(i))) / (y(def_jmap_hi(i)) - y(def_jmap_lo(i)))
			def_jmap_w(i) = max(0.0_wp, min(1.0_wp, def_jmap_w(i)))
		endif
	enddo
end subroutine defect_update_map

!********************************************************************
subroutine defect_interp_temp()
! Interpolate temp from AMR mesh to uniform defect mesh.
	integer :: i, j, k, i1, i2, j1, j2
	real(wp) :: wx, wy

	!$OMP PARALLEL DO PRIVATE(i, j, k, i1, i2, j1, j2, wx, wy)
	do k = k_def_lo, k_def_hi
	do j = 2, def_njm1
	do i = 2, def_nim1
		i1 = def_imap_lo(i); i2 = def_imap_hi(i); wx = def_imap_w(i)
		j1 = def_jmap_lo(j); j2 = def_jmap_hi(j); wy = def_jmap_w(j)
		temp_def(i,j,k) = (1.0_wp-wx)*(1.0_wp-wy) * temp(i1,j1,k) + &
		                   wx*(1.0_wp-wy)           * temp(i2,j1,k) + &
		                   (1.0_wp-wx)*wy           * temp(i1,j2,k) + &
		                   wx*wy                    * temp(i2,j2,k)
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO
end subroutine defect_interp_temp

!********************************************************************
subroutine update_max_temp()
! Called every time step: update max_temp and solidfield_def on defect mesh.
	use field_data, only: solidfield
	integer :: i, j, k
	integer :: i_hi, j_hi
	real(wp) :: sf_val

	if (adaptive_flag == 1) then
		i_hi = def_nim1
		j_hi = def_njm1
		! Step 1: interpolate temp → temp_def, update max_temp
		call defect_interp_temp()
		!$OMP PARALLEL DO PRIVATE(i, j, k)
		do k = k_def_lo, k_def_hi
		do j = 2, j_hi
		do i = 2, i_hi
			if (temp_def(i,j,k) > max_temp(i,j,k)) then
				max_temp(i,j,k) = temp_def(i,j,k)
			endif
		enddo
		enddo
		enddo
		!$OMP END PARALLEL DO
		! Step 2: interpolate solidfield → temp_def (reuse buffer), update solidfield_def
		call defect_interp_field(solidfield, temp_def)
		!$OMP PARALLEL DO PRIVATE(i, j, k)
		do k = k_def_lo, k_def_hi
		do j = 2, j_hi
		do i = 2, i_hi
			if (temp_def(i,j,k) > 0.5_wp) then
				solidfield_def(i,j,k) = temp_def(i,j,k)
			endif
		enddo
		enddo
		enddo
		!$OMP END PARALLEL DO
	else
		!$OMP PARALLEL DO PRIVATE(i, j, k)
		do k = k_def_lo, k_def_hi
		do j = 2, njm1
		do i = 2, nim1
			if (temp(i,j,k) > max_temp(i,j,k)) then
				max_temp(i,j,k) = temp(i,j,k)
			endif
		enddo
		enddo
		enddo
		!$OMP END PARALLEL DO
	endif

end subroutine update_max_temp

!********************************************************************
subroutine compute_defect_determ()
! Called after simulation completes: compute defect array from max_temp.
	integer :: i, j, k
	integer :: i_hi, j_hi

	if (adaptive_flag == 1) then
		i_hi = def_nim1
		j_hi = def_njm1
	else
		i_hi = nim1
		j_hi = njm1
	endif

	! Step 1: Compute defect values
	!$OMP PARALLEL DO PRIVATE(i, j, k)
	do k = k_def_lo, k_def_hi
	do j = 2, j_hi
	do i = 2, i_hi
		if (max_temp(i,j,k) < tsolid) then
			defect_arr(i,j,k) = -k_lof
		else if (max_temp(i,j,k) > tboiling) then
			defect_arr(i,j,k) = min(k_kep * (max_temp(i,j,k) - tboiling) / tboiling, 0.99_wp)
		else
			defect_arr(i,j,k) = 0.0_wp
		endif
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO

	! Step 2: Determine max laser scanning range from toolpath
	call compute_scan_range()

	! Step 3: Clean up — zero defect outside scanned region (convex hull)
	if (adaptive_flag == 1) then
		!$OMP PARALLEL DO PRIVATE(i, j, k)
		do k = k_def_lo, k_def_hi
		do j = 2, def_njm1
		do i = 2, def_nim1
			if (.not. point_in_scan_region(def_x(i), def_y(j))) then
				defect_arr(i,j,k) = 0.0_wp
			endif
		enddo
		enddo
		enddo
		!$OMP END PARALLEL DO
	else
		!$OMP PARALLEL DO PRIVATE(i, j, k)
		do k = k_def_lo, k_def_hi
		do j = 2, njm1
		do i = 2, nim1
			if (.not. point_in_scan_region(x(i), y(j))) then
				defect_arr(i,j,k) = 0.0_wp
			endif
		enddo
		enddo
		enddo
		!$OMP END PARALLEL DO
	endif

end subroutine compute_defect_determ

!********************************************************************
subroutine maxtemp_stochas()
! Placeholder for future stochastic defect method.
! Currently does nothing.
end subroutine maxtemp_stochas

!********************************************************************
subroutine compute_scan_range()
! Build scan region polygon from all track endpoints.
	integer :: n, ntracks, i
	real(wp) :: x1, x2, y1
	real(wp) :: left_x(TOOLLINES), left_y(TOOLLINES)
	real(wp) :: right_x(TOOLLINES), right_y(TOOLLINES)

	ntracks = 0
	do n = 2, TOOLLINES
		if (toolmatrix(n,1) < -0.5_wp) exit
		if (toolmatrix(n,5) >= laser_on_threshold) then
			x1 = toolmatrix(n-1,2)
			x2 = toolmatrix(n,2)
			y1 = toolmatrix(n-1,3)
			ntracks = ntracks + 1
			left_x(ntracks)  = min(x1, x2)
			left_y(ntracks)  = y1
			right_x(ntracks) = max(x1, x2)
			right_y(ntracks) = y1
		endif
	enddo

	if (ntracks == 0) return

	n_hull = 0
	do i = 1, ntracks
		n_hull = n_hull + 1
		hull_x(n_hull) = right_x(i)
		hull_y(n_hull) = right_y(i)
	enddo
	do i = ntracks, 1, -1
		n_hull = n_hull + 1
		hull_x(n_hull) = left_x(i)
		hull_y(n_hull) = left_y(i)
	enddo

	x_scan_min = minval(left_x(1:ntracks))
	x_scan_max = maxval(right_x(1:ntracks))
	y_scan_min = left_y(1)
	y_scan_max = left_y(ntracks)

end subroutine compute_scan_range

!********************************************************************
function point_in_scan_region(px, py) result(inside)
	real(wp), intent(in) :: px, py
	logical :: inside
	integer :: i, j
	real(wp) :: cross

	inside = .true.
	if (n_hull < 3) then
		inside = .false.
		return
	endif

	do i = 1, n_hull
		j = mod(i, n_hull) + 1
		cross = (hull_x(j) - hull_x(i)) * (py - hull_y(i)) - &
		        (hull_y(j) - hull_y(i)) * (px - hull_x(i))
		if (cross < -1.0e-15_wp) then
			inside = .false.
			return
		endif
	enddo

end function point_in_scan_region

!********************************************************************
subroutine write_defect_report()
	integer :: i, j, k
	integer, parameter :: lun = 89
	real(wp) :: v_cell, v_total, v_defect, v_lof, v_kep
	real(wp) :: frac_defect, frac_lof, frac_kep
	real(kind=4) :: val4
	integer :: gridx, gridy, gridk, npts
	integer :: i_hi, j_hi
	real(wp) :: dx_def, dy_def

	if (adaptive_flag == 1) then
		i_hi = def_nim1
		j_hi = def_njm1
	else
		i_hi = nim1
		j_hi = njm1
	endif

	v_total  = 0.0_wp
	v_lof    = 0.0_wp
	v_kep    = 0.0_wp

	do k = k_def_lo, k_def_hi
	do j = 2, j_hi
	do i = 2, i_hi
		if (adaptive_flag == 1) then
			if (.not. point_in_scan_region(def_x(i), def_y(j))) cycle
			v_cell = amr_dx_fine * amr_dx_fine * (zw(k+1)-zw(k))
		else
			if (.not. point_in_scan_region(x(i), y(j))) cycle
			v_cell = volume(i,j,k)
		endif
		v_total = v_total + v_cell
		if (defect_arr(i,j,k) < 0.0_wp) then
			v_lof = v_lof + abs(defect_arr(i,j,k)) * v_cell
		else if (defect_arr(i,j,k) > 0.0_wp .and. defect_arr(i,j,k) <= 1.0_wp) then
			v_kep = v_kep + defect_arr(i,j,k) * v_cell
		endif
	enddo
	enddo
	enddo

	v_defect = v_lof + v_kep
	if (v_total > 0.0_wp) then
		frac_defect = v_defect / v_total
		frac_lof    = v_lof / v_total
		frac_kep    = v_kep / v_total
	else
		frac_defect = 0.0_wp
		frac_lof    = 0.0_wp
		frac_kep    = 0.0_wp
	endif

	open(unit=lun, file=trim(file_prefix)//'defect_report.txt', action='write', status='replace')
	write(lun,'(a)') '============================================'
	write(lun,'(a)') '  PHOENIX Defect Report'
	write(lun,'(a)') '  Method: maxtemp_determ'
	write(lun,'(a)') '============================================'
	write(lun,'(a)')
	write(lun,'(a,f10.6,a)') '  Defect fraction:           ', frac_defect * 100.0_wp, ' %'
	write(lun,'(a,f10.6,a)') '  Lack-of-fusion fraction:   ', frac_lof * 100.0_wp, ' %'
	write(lun,'(a,f10.6,a)') '  Keyhole porosity fraction:  ', frac_kep * 100.0_wp, ' %'
	write(lun,'(a)')
	write(lun,'(a,es12.5,a)') '  Defect volume:             ', v_defect, ' m^3'
	write(lun,'(a,es12.5,a)') '  Lack-of-fusion volume:     ', v_lof, ' m^3'
	write(lun,'(a,es12.5,a)') '  Keyhole porosity volume:    ', v_kep, ' m^3'
	write(lun,'(a,es12.5,a)') '  Total reference volume:     ', v_total, ' m^3'
	write(lun,'(a)')
	write(lun,'(a)') '  Bounding box (X-Y plane):'
	write(lun,'(a,es12.5,a,es12.5,a)') '    X: [', x_scan_min, ', ', x_scan_max, '] m'
	write(lun,'(a,es12.5,a,es12.5,a)') '    Y: [', y_scan_min, ', ', y_scan_max, '] m'
	write(lun,'(a)')
	write(lun,'(a,i0,a)') '  Scan region convex hull (', n_hull, ' vertices):'
	do i = 1, n_hull
		write(lun,'(a,i2,a,es12.5,a,es12.5,a)') '    ', i, ':  (', hull_x(i), ', ', hull_y(i), ') m'
	enddo
	write(lun,'(a)')
	write(lun,'(a)') '  Layer Z range:'
	write(lun,'(a,i4,a,i4)') '    k indices: ', k_def_lo, ' to ', k_def_hi
	write(lun,'(a,es12.5,a,es12.5,a)') '    Z: [', z(k_def_lo), ', ', z(k_def_hi), '] m'
	write(lun,'(a,es12.5,a)') '    Layer height: ', layerheight, ' m'
	write(lun,'(a)')
	write(lun,'(a)') '  Parameters:'
	write(lun,'(a,f8.4)') '    k_lof = ', k_lof
	write(lun,'(a,f8.4)') '    k_kep = ', k_kep
	write(lun,'(a,f10.2,a)') '    T_solid   = ', tsolid, ' K'
	write(lun,'(a,f10.2,a)') '    T_boiling = ', tboiling, ' K'
	write(lun,'(a)') '============================================'
	close(lun)

	write(9,'(a)') ''
	write(9,'(a)') '  === Defect Analysis (maxtemp_determ) ==='
	write(9,'(a,f10.6,a)') '  Defect fraction:         ', frac_defect * 100.0_wp, ' %'
	write(9,'(a,f10.6,a)') '  Lack-of-fusion fraction: ', frac_lof * 100.0_wp, ' %'
	write(9,'(a,f10.6,a)') '  Keyhole porosity fraction:', frac_kep * 100.0_wp, ' %'

	call write_defect_vtk()

end subroutine write_defect_report

!********************************************************************
subroutine write_defect_vtk()
! Write single defect.vtk with multiple scalars: defect, maxtemp, solidID.
	use field_data, only: solidfield
	integer :: i, j, k, npts
	integer :: gridx, gridy, gridz, i_hi, j_hi
	real(kind=4) :: val4
	integer, parameter :: lun = 90
	character(len=256) :: vtkfile

	if (adaptive_flag == 1) then
		gridx = def_nim1 - 2 + 1
		gridy = def_njm1 - 2 + 1
		i_hi = def_nim1
		j_hi = def_njm1
	else
		gridx = nim1 - 2 + 1
		gridy = njm1 - 2 + 1
		i_hi = nim1
		j_hi = njm1
	endif
	gridz = k_def_hi - k_def_lo + 1
	npts = gridx * gridy * gridz
	vtkfile = trim(file_prefix)//'defect.vtk'

	! --- Header + coordinates ---
	open(unit=lun, file=trim(vtkfile))
	write(lun,'(A)') '# vtk DataFile Version 3.0'
	write(lun,'(A)') 'PHOENIX defect analysis'
	write(lun,'(A)') 'BINARY'
	write(lun,'(A)') 'DATASET STRUCTURED_GRID'
	write(lun,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', gridx, ' ', gridy, ' ', gridz
	write(lun,'(A,I0,A)') 'POINTS ', npts, ' float'
	close(lun)

	open(unit=lun, file=trim(vtkfile), &
	     access='stream', form='unformatted', position='append', convert='big_endian')
	do k = k_def_lo, k_def_hi
	do j = 2, j_hi
	do i = 2, i_hi
		if (adaptive_flag == 1) then
			val4 = real(def_x(i), 4); write(lun) val4
			val4 = real(def_y(j), 4); write(lun) val4
		else
			val4 = real(x(i), 4); write(lun) val4
			val4 = real(y(j), 4); write(lun) val4
		endif
		val4 = real(z(k), 4); write(lun) val4
	enddo
	enddo
	enddo
	close(lun)

	! --- POINT_DATA header ---
	open(unit=lun, file=trim(vtkfile), position='append')
	write(lun,'(A,I0)') 'POINT_DATA ', npts
	close(lun)

	! --- Scalar: defect ---
	call write_defect_scalar(lun, vtkfile, 'defect', defect_arr, i_hi, j_hi)

	! --- Scalar: maxtemp ---
	call write_defect_scalar(lun, vtkfile, 'maxtemp', max_temp, i_hi, j_hi)


end subroutine write_defect_vtk

!********************************************************************
subroutine write_defect_scalar(lun, vtkfile, name, field, i_hi, j_hi)
	integer, intent(in) :: lun, i_hi, j_hi
	character(len=*), intent(in) :: vtkfile, name
	real(wp), intent(in) :: field(:,:,:)
	integer :: i, j, k
	real(kind=4) :: val4

	open(unit=lun, file=trim(vtkfile), position='append')
	write(lun,'(A)') 'SCALARS '//trim(name)//' float 1'
	write(lun,'(A)') 'LOOKUP_TABLE default'
	close(lun)

	open(unit=lun, file=trim(vtkfile), &
	     access='stream', form='unformatted', position='append', convert='big_endian')
	do k = k_def_lo, k_def_hi
	do j = 2, j_hi
	do i = 2, i_hi
		val4 = real(field(i,j,k), 4)
		write(lun) val4
	enddo
	enddo
	enddo
	close(lun)
end subroutine write_defect_scalar

!********************************************************************
subroutine defect_interp_field(field_sim, field_def)
! Interpolate a 3D field from simulation mesh to uniform defect mesh.
	real(wp), intent(in) :: field_sim(:,:,:)
	real(wp), intent(out) :: field_def(:,:,:)
	integer :: i, j, k, i1, i2, j1, j2
	real(wp) :: wx, wy

	!$OMP PARALLEL DO PRIVATE(i, j, k, i1, i2, j1, j2, wx, wy)
	do k = k_def_lo, k_def_hi
	do j = 2, def_njm1
	do i = 2, def_nim1
		i1 = def_imap_lo(i); i2 = def_imap_hi(i); wx = def_imap_w(i)
		j1 = def_jmap_lo(j); j2 = def_jmap_hi(j); wy = def_jmap_w(j)
		field_def(i,j,k) = (1.0_wp-wx)*(1.0_wp-wy) * field_sim(i1,j1,k) + &
		                    wx*(1.0_wp-wy)           * field_sim(i2,j1,k) + &
		                    (1.0_wp-wx)*wy           * field_sim(i1,j2,k) + &
		                    wx*wy                    * field_sim(i2,j2,k)
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO
end subroutine defect_interp_field

end module defect_field
