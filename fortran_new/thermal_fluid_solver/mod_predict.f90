!______________________________________________________________________________
!
module prediction
!______________________________________________________________________________
! Field prediction module for heating steps.
! Shifts all fields (enthalpy, velocity, pressure) by integer number of cells
! along the scan direction within the melt pool region.
! Designed as a separate module for future extensibility (e.g., ML-based prediction).
!
	use precision
	use geometry
	use parameters
	use field_data
	implicit none

	contains

!********************************************************************
subroutine predict_shift_integer(vx, vy, dt, is, ie, js, je, ks, ke)
! Integer-cell shift of all fields in the melt pool region.
! No interpolation → no numerical diffusion. Sharp interfaces preserved.
	real(wp), intent(in) :: vx, vy, dt
	integer, intent(in) :: is, ie, js, je, ks, ke
	real(wp) :: dx_shift, dy_shift, dx_avg, dy_avg
	integer :: di, dj

	dx_shift = vx * dt
	dy_shift = vy * dt

	! Compute shift in cells using average local cell size in the region
	if (ie > is) then
		dx_avg = (x(ie) - x(is)) / real(ie - is, wp)
	else
		dx_avg = amr_dx_fine
	endif
	if (je > js) then
		dy_avg = (y(je) - y(js)) / real(je - js, wp)
	else
		dy_avg = amr_dx_fine
	endif

	di = nint(dx_shift / dx_avg)
	dj = nint(dy_shift / dy_avg)

	if (di == 0 .and. dj == 0) return

	call ishift_field(enthalpy, di, dj, is, ie, js, je, ks, ke)
	call ishift_field(uVel,     di, dj, is, ie, js, je, ks, ke)
	call ishift_field(vVel,     di, dj, is, ie, js, je, ks, ke)
	call ishift_field(wVel,     di, dj, is, ie, js, je, ks, ke)
	call ishift_field(pressure, di, dj, is, ie, js, je, ks, ke)

end subroutine predict_shift_integer

!********************************************************************
subroutine ishift_field(field, di, dj, is, ie, js, je, ks, ke)
! Shift field by (di, dj) cells within local region. No interpolation.
! field(i,j,k) = field(i-di, j-dj, k) where source is valid, else keep original.
! Iteration order avoids overwrite-before-read for the given shift direction.
	real(wp), intent(inout) :: field(:,:,:)
	integer, intent(in) :: di, dj, is, ie, js, je, ks, ke
	integer :: i, j, k, isrc, jsrc

	!$OMP PARALLEL DO PRIVATE(i, j, k, isrc, jsrc)
	do k = ks, ke
		if (dj >= 0) then
			do j = je, js, -1
				jsrc = j - dj
				if (jsrc < 2 .or. jsrc > njm1) cycle
				if (di >= 0) then
					do i = ie, is, -1
						isrc = i - di
						if (isrc < 2 .or. isrc > nim1) cycle
						field(i,j,k) = field(isrc, jsrc, k)
					enddo
				else
					do i = is, ie
						isrc = i - di
						if (isrc < 2 .or. isrc > nim1) cycle
						field(i,j,k) = field(isrc, jsrc, k)
					enddo
				endif
			enddo
		else
			do j = js, je
				jsrc = j - dj
				if (jsrc < 2 .or. jsrc > njm1) cycle
				if (di >= 0) then
					do i = ie, is, -1
						isrc = i - di
						if (isrc < 2 .or. isrc > nim1) cycle
						field(i,j,k) = field(isrc, jsrc, k)
					enddo
				else
					do i = is, ie
						isrc = i - di
						if (isrc < 2 .or. isrc > nim1) cycle
						field(i,j,k) = field(isrc, jsrc, k)
					enddo
				endif
			enddo
		endif
	enddo
	!$OMP END PARALLEL DO

end subroutine ishift_field

end module prediction
