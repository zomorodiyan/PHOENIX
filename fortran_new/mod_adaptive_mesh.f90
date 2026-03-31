!______________________________________________________________________________
!
module adaptive_mesh_mod
!______________________________________________________________________________
! Movable adaptive structured mesh in X-Y plane.
! Refined region follows laser/melt pool; coarse cells use geometric expansion.
! Z direction unchanged. Total cell count (ni, nj) conserved.
!
	use precision
	use geometry
	use parameters
	use sim_state
	use field_data
	use dimensions
	use coeff_data
	implicit none

	real(wp) :: amr_half_x, amr_half_y          ! current refined region half-sizes
	real(wp), parameter :: amr_expand_ratio = 1.3_wp  ! geometric expansion ratio
	integer  :: amr_step_counter = 0
	logical  :: amr_needs_remesh = .false.

	! Index mapping for interpolation (precomputed after each remesh)
	integer, allocatable  :: imap_lo(:), imap_hi(:)  ! bracketing old i-indices
	real(wp), allocatable :: imap_w(:)               ! interpolation weight in x
	integer, allocatable  :: jmap_lo(:), jmap_hi(:)
	real(wp), allocatable :: jmap_w(:)

	! Saved domain lengths (constant)
	real(wp) :: amr_dimx, amr_dimy
	integer  :: amr_ncvx, amr_ncvy  ! total interior cells in x, y

	contains

!********************************************************************
subroutine amr_init()
	amr_half_x = amr_local_half_x
	amr_half_y = amr_local_half_y
	amr_dimx = dimx
	amr_dimy = dimy
	amr_ncvx = nim1 - 1   ! number of interior cells = nim1-1
	amr_ncvy = njm1 - 1

	! Validate
	if (nzx /= 1 .or. nzy /= 1) then
		write(*,*) 'ERROR: adaptive_flag=1 requires nzx=1 and nzy=1'
		stop
	endif

	! Allocate mapping arrays
	allocate(imap_lo(ni), imap_hi(ni), imap_w(ni))
	allocate(jmap_lo(nj), jmap_hi(nj), jmap_w(nj))

	! Don't regenerate at init — beam position not yet set.
	! First remesh will happen at step = remesh_interval after laser_beam sets beam_pos.
end subroutine amr_init

!********************************************************************
subroutine amr_check_remesh(step_idx)
	integer, intent(in) :: step_idx
	real(wp) :: half_x_new, half_y_new
	real(wp) :: x_lo, x_hi, y_lo, y_hi
	integer  :: n_fine, n_coarse, n_min_coarse

	amr_step_counter = amr_step_counter + 1
	if (amr_step_counter < remesh_interval) return

	amr_step_counter = 0

	! Expand if melt pool exceeds initial region
	half_x_new = max(amr_local_half_x, alen)
	half_y_new = max(amr_local_half_y, width)

	! 20% check for X
	n_fine = nint(2.0_wp * half_x_new / amr_dx_fine)
	n_coarse = amr_ncvx - n_fine
	n_min_coarse = nint(0.2_wp * n_fine)
	if (n_coarse < n_min_coarse) then
		write(9,'(A,I6,A)') '  WARNING: Step ', step_idx, &
			' adaptive mesh cannot expand X further - insufficient coarse cells (melt pool too large)'
		half_x_new = amr_half_x  ! keep previous size
	endif

	! 20% check for Y
	n_fine = nint(2.0_wp * half_y_new / amr_dx_fine)
	n_coarse = amr_ncvy - n_fine
	n_min_coarse = nint(0.2_wp * n_fine)
	if (n_coarse < n_min_coarse) then
		write(9,'(A,I6,A)') '  WARNING: Step ', step_idx, &
			' adaptive mesh cannot expand Y further - insufficient coarse cells (melt pool too large)'
		half_y_new = amr_half_y
	endif

	! Clamp to domain
	x_lo = max(0.0_wp, beam_pos - half_x_new)
	x_hi = min(amr_dimx, beam_pos + half_x_new)
	y_lo = max(0.0_wp, beam_posy - half_y_new)
	y_hi = min(amr_dimy, beam_posy + half_y_new)

	! Adjust half sizes to clamped bounds
	half_x_new = min(half_x_new, beam_pos, amr_dimx - beam_pos)
	half_y_new = min(half_y_new, beam_posy, amr_dimy - beam_posy)

	amr_half_x = half_x_new
	amr_half_y = half_y_new
	amr_needs_remesh = .true.
end subroutine amr_check_remesh

!********************************************************************
subroutine amr_regenerate_grid()
	real(wp), allocatable :: x_old(:), y_old(:)
	integer :: i, j, k

	! Save old grid for interpolation
	allocate(x_old(ni), y_old(nj))
	x_old = x
	y_old = y

	! Generate new 1D grids
	call amr_generate_1d(amr_dimx, amr_ncvx, beam_pos, amr_half_x, xu, x, ni)
	call amr_generate_1d(amr_dimy, amr_ncvy, beam_posy, amr_half_y, yv, y, nj)

	! Recompute geometry metrics
	call amr_recompute_geometry()

	! Build interpolation maps
	call amr_build_map_1d(x_old, ni, x, ni, imap_lo, imap_hi, imap_w)
	call amr_build_map_1d(y_old, nj, y, nj, jmap_lo, jmap_hi, jmap_w)

	! Interpolate solution fields
	call amr_interpolate_all_fields(x_old, y_old)

	deallocate(x_old, y_old)

	! Update domain dimensions (should be unchanged)
	dimx = x(ni)
	dimy = y(nj)
end subroutine amr_regenerate_grid

!********************************************************************
subroutine amr_generate_1d(domain_len, n_total, center, half_len, vel_grid, scal_grid, n)
! Build 1D adaptive velocity grid with three regions.
! vel_grid(1)=vel_grid(2)=0, vel_grid(3..n) interior faces. n = n_total+2.
	real(wp), intent(in) :: domain_len, center, half_len
	integer, intent(in) :: n_total, n
	real(wp), intent(inout) :: vel_grid(:), scal_grid(:)

	real(wp) :: x_lo, x_hi, L_left, L_right, r_eff, dx0, gs, pos
	integer  :: n_fine, n_coarse, n_left, n_right, idx, i

	x_lo = max(0.0_wp, center - half_len)
	x_hi = min(domain_len, center + half_len)

	n_fine = nint((x_hi - x_lo) / amr_dx_fine)
	n_fine = max(1, min(n_fine, n_total))
	x_hi = x_lo + real(n_fine, wp) * amr_dx_fine
	if (x_hi > domain_len) then
		x_hi = domain_len
		n_fine = max(1, int((x_hi - x_lo) / amr_dx_fine))
	endif

	n_coarse = n_total - n_fine
	L_left = x_lo
	L_right = domain_len - x_hi

	if (L_left + L_right > 1.0e-15_wp) then
		n_left = nint(real(n_coarse, wp) * L_left / (L_left + L_right))
		n_right = n_coarse - n_left
	else
		n_left = 0; n_right = n_coarse
	endif
	if (L_left > 1.0e-15_wp .and. n_left < 1 .and. n_coarse > 1) then
		n_left = 1; n_right = n_coarse - 1
	endif
	if (L_right > 1.0e-15_wp .and. n_right < 1 .and. n_coarse > 1) then
		n_right = 1; n_left = n_coarse - 1
	endif

	vel_grid(1) = 0.0_wp
	vel_grid(2) = 0.0_wp
	idx = 2

	! --- Left coarse (largest cells near 0, smallest near x_lo) ---
	if (n_left > 0 .and. L_left > 1.0e-15_wp) then
		r_eff = amr_find_ratio(amr_dx_fine, L_left, n_left)
		! Build ascending from smallest (near x_lo) to largest (near 0), then reverse
		! Forward: face(i) = dx0 * (r^i - 1)/(r - 1), dx0=smallest cell
		if (abs(r_eff - 1.0_wp) < 1.0e-10_wp) then
			dx0 = L_left / real(n_left, wp)
			do i = 1, n_left
				idx = idx + 1
				vel_grid(idx) = real(i, wp) * dx0
			enddo
		else
			gs = (r_eff**n_left - 1.0_wp) / (r_eff - 1.0_wp)
			dx0 = L_left / gs
			! Reversed: face(i) from left side
			! face(i) = L_left - dx0 * (r^(n_left-i) - 1)/(r-1) ... nope, just reverse
			! Forward cumulative: cum(j) = dx0*(r^j-1)/(r-1)
			! Reversed face i = L_left - cum(n_left - i) for i=1..n_left-1, face(n_left)=L_left
			do i = 1, n_left
				gs = (r_eff**(n_left - i) - 1.0_wp) / (r_eff - 1.0_wp)
				pos = L_left - dx0 * gs
				if (i == n_left) pos = L_left
				idx = idx + 1
				vel_grid(idx) = pos
			enddo
		endif
	endif

	! --- Refined region (uniform) ---
	do i = 1, n_fine
		idx = idx + 1
		vel_grid(idx) = x_lo + real(i, wp) * amr_dx_fine
	enddo
	if (n_fine > 0) vel_grid(idx) = x_hi

	! --- Right coarse (smallest near x_hi, largest near domain_len) ---
	if (n_right > 0 .and. L_right > 1.0e-15_wp) then
		r_eff = amr_find_ratio(amr_dx_fine, L_right, n_right)
		if (abs(r_eff - 1.0_wp) < 1.0e-10_wp) then
			dx0 = L_right / real(n_right, wp)
			do i = 1, n_right
				idx = idx + 1
				vel_grid(idx) = x_hi + real(i, wp) * dx0
			enddo
		else
			gs = (r_eff**n_right - 1.0_wp) / (r_eff - 1.0_wp)
			dx0 = L_right / gs
			do i = 1, n_right
				gs = (r_eff**i - 1.0_wp) / (r_eff - 1.0_wp)
				pos = x_hi + dx0 * gs
				if (i == n_right) pos = domain_len
				idx = idx + 1
				vel_grid(idx) = pos
			enddo
		endif
	endif

	vel_grid(n) = domain_len

	! Scalar grid
	do i = 1, n - 1
		scal_grid(i) = (vel_grid(i) + vel_grid(i+1)) * 0.5_wp
	enddo
	scal_grid(n) = vel_grid(n)
end subroutine amr_generate_1d

!********************************************************************
function amr_find_ratio(dx_fine, length, ncells) result(r)
! Find expansion ratio r such that: dx_fine * (r^n - 1)/(r-1) = length.
! If uniform spacing suffices (n*dx_fine >= length), return 1.
! Otherwise bisection in [1, r_max].
	real(wp), intent(in) :: dx_fine, length
	integer, intent(in) :: ncells
	real(wp) :: r
	real(wp) :: r_lo, r_hi, r_mid, f_mid
	integer :: iter

	if (real(ncells, wp) * dx_fine >= length) then
		r = 1.0_wp
		return
	endif

	! Bisection: find r in [1, r_max] where dx_fine*(r^n-1)/(r-1) = length
	r_lo = 1.0_wp
	r_hi = 2.0_wp
	! Expand r_hi until it overshoots
	do iter = 1, 50
		if (dx_fine * (r_hi**ncells - 1.0_wp) / (r_hi - 1.0_wp) > length) exit
		r_hi = r_hi * 2.0_wp
	enddo

	do iter = 1, 100
		r_mid = 0.5_wp * (r_lo + r_hi)
		if (r_mid <= 1.0_wp + 1.0e-12_wp) then
			r = 1.0_wp
			return
		endif
		f_mid = dx_fine * (r_mid**ncells - 1.0_wp) / (r_mid - 1.0_wp) - length
		if (abs(f_mid) < 1.0e-12_wp * length) exit
		if (f_mid > 0.0_wp) then
			r_hi = r_mid
		else
			r_lo = r_mid
		endif
	enddo
	r = r_mid
end function amr_find_ratio

!********************************************************************
subroutine amr_recompute_geometry()
	integer :: i, j, k

	! Reciprocal node spacing
	do i = 2, ni
		dxpwinv(i) = 1.0_wp / (x(i) - x(i-1))
	enddo
	do j = 2, nj
		dypsinv(j) = 1.0_wp / (y(j) - y(j-1))
	enddo
	! Z unchanged — no need to recompute dzpbinv

	! Interpolation fractions
	do i = 1, nim1
		fracx(i) = (x(i+1) - xu(i+1)) / (x(i+1) - x(i))
	enddo
	do j = 1, njm1
		fracy(j) = (y(j+1) - yv(j+1)) / (y(j+1) - y(j))
	enddo
	! fracz unchanged

	! 3D volumes
	!$OMP PARALLEL DO PRIVATE(i, j, k)
	do k = 2, nkm1
	do j = 2, njm1
	do i = 2, nim1
		volume(i,j,k) = (xu(i+1)-xu(i)) * (yv(j+1)-yv(j)) * (zw(k+1)-zw(k))
		volume_u(i,j,k) = (x(i)-x(i-1)) * (yv(j+1)-yv(j)) * (zw(k+1)-zw(k))
		volume_v(i,j,k) = (xu(i+1)-xu(i)) * (y(j)-y(j-1)) * (zw(k+1)-zw(k))
		volume_w(i,j,k) = (xu(i+1)-xu(i)) * (yv(j+1)-yv(j)) * (z(k)-z(k-1))
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO

	! 2D face areas
	!$OMP PARALLEL DO PRIVATE(i, j)
	do j = 2, njm1
	do i = 2, nim1
		areaij(i,j) = (xu(i+1)-xu(i)) * (yv(j+1)-yv(j))
		areauij(i,j) = (x(i)-x(i-1)) * (yv(j+1)-yv(j))
		areavij(i,j) = (xu(i+1)-xu(i)) * (y(j)-y(j-1))
	enddo
	enddo
	!$OMP END PARALLEL DO

	!$OMP PARALLEL DO PRIVATE(i, k)
	do k = 2, nkm1
	do i = 2, nim1
		areaik(i,k) = (xu(i+1)-xu(i)) * (zw(k+1)-zw(k))
		areawik(i,k) = (xu(i+1)-xu(i)) * (z(k)-z(k-1))
		areauik(i,k) = (x(i)-x(i-1)) * (zw(k+1)-zw(k))
	enddo
	enddo
	!$OMP END PARALLEL DO

	!$OMP PARALLEL DO PRIVATE(j, k)
	do k = 2, nkm1
	do j = 2, njm1
		areajk(j,k) = (yv(j+1)-yv(j)) * (zw(k+1)-zw(k))
		areavjk(j,k) = (y(j)-y(j-1)) * (zw(k+1)-zw(k))
		areawjk(j,k) = (yv(j+1)-yv(j)) * (z(k)-z(k-1))
	enddo
	enddo
	!$OMP END PARALLEL DO

end subroutine amr_recompute_geometry

!********************************************************************
subroutine amr_build_map_1d(x_old, n_old, x_new, n_new, map_lo, map_hi, map_w)
! Build index mapping: for each new cell center, find bracketing old cell centers.
	real(wp), intent(in) :: x_old(:), x_new(:)
	integer, intent(in) :: n_old, n_new
	integer, intent(out) :: map_lo(:), map_hi(:)
	real(wp), intent(out) :: map_w(:)
	integer :: i, j_lo
	real(wp) :: xn

	j_lo = 1
	do i = 1, n_new
		xn = x_new(i)
		! Advance j_lo until x_old(j_lo+1) >= xn
		do while (j_lo < n_old - 1 .and. x_old(j_lo + 1) < xn)
			j_lo = j_lo + 1
		enddo
		map_lo(i) = j_lo
		map_hi(i) = min(j_lo + 1, n_old)
		if (map_lo(i) == map_hi(i)) then
			map_w(i) = 0.0_wp
		else
			map_w(i) = (xn - x_old(map_lo(i))) / (x_old(map_hi(i)) - x_old(map_lo(i)))
			map_w(i) = max(0.0_wp, min(1.0_wp, map_w(i)))
		endif
		! Reset for next search (don't reset j_lo — monotone sweep)
	enddo
end subroutine amr_build_map_1d

!********************************************************************
subroutine amr_interpolate_all_fields(x_old, y_old)
	real(wp), intent(in) :: x_old(:), y_old(:)

	! Interpolate solution fields
	call amr_interp_field(enthalpy)
	call amr_interp_field(temp)
	call amr_interp_field(fracl)
	call amr_interp_field(hnot)
	call amr_interp_field(tnot)
	call amr_interp_field(fraclnot)
	call amr_interp_field(uVel)
	call amr_interp_field(vVel)
	call amr_interp_field(wVel)
	call amr_interp_field(unot)
	call amr_interp_field(vnot)
	call amr_interp_field(wnot)
	call amr_interp_field(pressure)
	call amr_interp_field(pp)
	call amr_interp_field(solidfield)

	! Property fields (recomputed each iteration, but interpolate for smooth start)
	call amr_interp_field(den)
	call amr_interp_field(vis)
	call amr_interp_field(diff)

end subroutine amr_interpolate_all_fields

!********************************************************************
subroutine amr_interp_field(field)
! Bilinear interpolation in X-Y using precomputed maps. Z is 1:1.
	real(wp), intent(inout) :: field(:,:,:)
	real(wp), allocatable :: tmp(:,:,:)
	integer :: i, j, k
	integer :: i1, i2, j1, j2
	real(wp) :: wx, wy

	allocate(tmp(ni, nj, nk))

	!$OMP PARALLEL DO PRIVATE(i, j, k, i1, i2, j1, j2, wx, wy)
	do k = 1, nk
	do j = 1, nj
	do i = 1, ni
		i1 = imap_lo(i); i2 = imap_hi(i); wx = imap_w(i)
		j1 = jmap_lo(j); j2 = jmap_hi(j); wy = jmap_w(j)
		tmp(i,j,k) = (1.0_wp-wx)*(1.0_wp-wy) * field(i1,j1,k) + &
		             wx*(1.0_wp-wy)            * field(i2,j1,k) + &
		             (1.0_wp-wx)*wy            * field(i1,j2,k) + &
		             wx*wy                     * field(i2,j2,k)
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO

	field = tmp
	deallocate(tmp)
end subroutine amr_interp_field

!********************************************************************
subroutine amr_validate_grid()
	integer :: i, j, k
	real(wp) :: sum_x, sum_y, tol
	logical :: valid

	valid = .true.
	tol = 1.0e-10_wp

	! Check domain length
	sum_x = xu(ni) - xu(2)
	if (abs(sum_x - amr_dimx) > tol * amr_dimx) then
		write(*,'(A,es12.5,A,es12.5)') 'AMR ERROR: X domain mismatch: ', sum_x, ' vs ', amr_dimx
		valid = .false.
	endif
	sum_y = yv(nj) - yv(2)
	if (abs(sum_y - amr_dimy) > tol * amr_dimy) then
		write(*,'(A,es12.5,A,es12.5)') 'AMR ERROR: Y domain mismatch: ', sum_y, ' vs ', amr_dimy
		valid = .false.
	endif

	! Check positive volumes
	do k = 2, nkm1
	do j = 2, njm1
	do i = 2, nim1
		if (volume(i,j,k) <= 0.0_wp) then
			write(*,'(A,3I5)') 'AMR ERROR: Non-positive volume at i,j,k=', i, j, k
			valid = .false.
			goto 100
		endif
	enddo
	enddo
	enddo
100 continue

	! Check positive spacing inverses
	do i = 2, ni
		if (dxpwinv(i) <= 0.0_wp) then
			write(*,'(A,I5)') 'AMR ERROR: Non-positive dxpwinv at i=', i
			valid = .false.
			exit
		endif
	enddo
	do j = 2, nj
		if (dypsinv(j) <= 0.0_wp) then
			write(*,'(A,I5)') 'AMR ERROR: Non-positive dypsinv at j=', j
			valid = .false.
			exit
		endif
	enddo

	if (.not. valid) then
		write(*,*) 'AMR: Grid validation FAILED - aborting'
		stop
	endif
end subroutine amr_validate_grid

end module adaptive_mesh_mod
