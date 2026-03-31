!______________________________________________________________________________
!
module dimensions
!______________________________________________________________________________
! Melt pool dimension detection via flood-fill from tpeak location.
! Reports the connected molten region containing the hottest cell.
! Robust for AMR, multi-pool, turnaround.
!
	use initialization
	use parameters
	use laserinput
	implicit none

	real(wp) alen,depth,width,hpeak,tpeak,umax,vmax,wmax
        !   length, depth, width, peak enthalpy, peak temperature, maximum velocity at each direction

	integer	istat,jstat,kstat,iend,jend,istatp1,iendm1  !

	contains

subroutine pool_size(ilo, ihi, jlo, jhi, klo, khi)

	integer, intent(in) :: ilo, ihi, jlo, jhi, klo, khi
	integer i, j, k
	integer il, ih, jl, jh
	integer ipeak, jpeak
	integer, allocatable :: mask(:,:)  ! 2D flood-fill mask at surface
	logical :: changed

	il = max(ilo, 2)
	ih = min(ihi, nim1)
	jl = max(jlo, 2)
	jh = min(jhi, njm1)

	tpeak = maxval(temp(il:ih, jl:jh, 2:nkm1))

	alen  = 0.0_wp
	depth = 0.0_wp
	width = 0.0_wp
	imin = istart; imax = istart
	jmin = jstart; jmax = jstart
	kmin = nkm1;   kmax = nkm1

	if (tpeak .le. tsolid) then
		istat = istart; iend = istart
		jstat = jstart; jend = jstart
		kstat = nkm1
		istatp1 = istat + 1; iendm1 = iend - 1
		return
	endif

	! --- Find tpeak location on surface ---
	ipeak = il; jpeak = jl
	do j = jl, jh
	do i = il, ih
		if (temp(i,j,nk) .eq. tpeak) then
			ipeak = i; jpeak = j
			goto 10
		endif
	enddo
	enddo
	! If tpeak not at surface, search all k for surface projection
	do k = nkm1, 2, -1
	do j = jl, jh
	do i = il, ih
		if (temp(i,j,k) .eq. tpeak) then
			ipeak = i; jpeak = j
			goto 10
		endif
	enddo
	enddo
	enddo
10	continue

	! --- 2D flood fill from (ipeak, jpeak) at surface k=nk ---
	allocate(mask(ni, nj))
	mask = 0

	! Seed
	if (temp(ipeak, jpeak, nk) .gt. tsolid) then
		mask(ipeak, jpeak) = 1
	else
		! tpeak is below surface; find nearest surface cell above tsolid
		! Search outward from ipeak, jpeak
		do j = jl, jh
		do i = il, ih
			if (temp(i,j,nk) .gt. tsolid) then
				if (mask(ipeak,jpeak) == 0) then
					ipeak = i; jpeak = j
					mask(i,j) = 1
				endif
			endif
		enddo
		enddo
	endif

	! Iterative expansion: mark neighbors that are above tsolid
	changed = .true.
	do while (changed)
		changed = .false.
		do j = jl, jh
		do i = il, ih
			if (mask(i,j) == 1) cycle
			if (temp(i,j,nk) .le. tsolid) cycle
			! Check 4-connected neighbors
			if (i > il .and. mask(i-1,j) == 1) then
				mask(i,j) = 1; changed = .true.; cycle
			endif
			if (i < ih .and. mask(i+1,j) == 1) then
				mask(i,j) = 1; changed = .true.; cycle
			endif
			if (j > jl .and. mask(i,j-1) == 1) then
				mask(i,j) = 1; changed = .true.; cycle
			endif
			if (j < jh .and. mask(i,j+1) == 1) then
				mask(i,j) = 1; changed = .true.; cycle
			endif
		enddo
		enddo
	enddo

	! --- Bounding box of connected region ---
	imin = ih; imax = il
	jmin = jh; jmax = jl

	do j = jl, jh
	do i = il, ih
		if (mask(i,j) == 1) then
			imin = min(imin, i)
			imax = max(imax, i)
			jmin = min(jmin, j)
			jmax = max(jmax, j)
		endif
	enddo
	enddo

	if (imax < imin) then
		deallocate(mask)
		istat = istart; iend = istart; jstat = jstart; jend = jstart; kstat = nkm1
		istatp1 = istat + 1; iendm1 = iend - 1
		return
	endif

	! --- Sub-cell interpolation at bounding-box edges ---
	call interp_pool_length(mask, il, ih, jl, jh, alen)
	call interp_pool_width(mask, il, ih, jl, jh, width)

	! --- Depth: deepest k in connected region, with interpolation ---
	kmin = nkm1
	do k = 2, nkm1
	do j = jmin, jmax
	do i = imin, imax
		if (mask(i,j) == 1 .and. temp(i,j,k) .gt. tsolid) then
			kmin = min(kmin, k)
		endif
	enddo
	enddo
	enddo

	call interp_pool_depth(mask, depth)
	kmax = nkm1

	deallocate(mask)

	! --- Solution domain for momentum equations ---
	istat = max(imin - 3, 2)
	iend  = min(imax + 3, nim1)
	jstat = max(jmin - 3, 2)
	jend  = min(jmax + 2, njm1)
	kstat = max(kmin - 2, 3)
	istatp1 = istat + 1
	iendm1  = iend - 1

end subroutine pool_size

!********************************************************************
subroutine interp_pool_length(mask, il, ih, jl, jh, length)
! Refine x-extent by interpolating tsolid crossing at left (imin) and right (imax) edges.
	integer, intent(in) :: mask(:,:), il, ih, jl, jh
	real(wp), intent(out) :: length
	integer :: i, j
	real(wp) :: xlo, xhi, t_in, t_out, frac

	xhi = x(imin)  ! will be pushed right
	xlo = x(imax)  ! will be pushed left

	! Right edge: for each j in pool, interpolate between imax and imax+1
	do j = jmin, jmax
		if (mask(imin,j) /= 1) cycle
		if (imin > il) then
			t_in  = temp(imin, j, nk)
			t_out = temp(imin-1, j, nk)
			if (t_in > tsolid .and. t_out <= tsolid .and. abs(t_in - t_out) > 1.0_wp) then
				frac = (t_in - tsolid) / (t_in - t_out)
				xlo = min(xlo, x(imin) - frac * (x(imin) - x(imin-1)))
			else
				xlo = min(xlo, x(imin))
			endif
		else
			xlo = min(xlo, x(imin))
		endif
	enddo

	do j = jmin, jmax
		if (mask(imax,j) /= 1) cycle
		if (imax < ih) then
			t_in  = temp(imax, j, nk)
			t_out = temp(imax+1, j, nk)
			if (t_in > tsolid .and. t_out <= tsolid .and. abs(t_in - t_out) > 1.0_wp) then
				frac = (t_in - tsolid) / (t_in - t_out)
				xhi = max(xhi, x(imax) + frac * (x(imax+1) - x(imax)))
			else
				xhi = max(xhi, x(imax))
			endif
		else
			xhi = max(xhi, x(imax))
		endif
	enddo

	length = max(0.0_wp, xhi - xlo)
end subroutine interp_pool_length

!********************************************************************
subroutine interp_pool_width(mask, il, ih, jl, jh, w)
! Refine y-extent by interpolating tsolid crossing at bottom (jmin) and top (jmax) edges.
	integer, intent(in) :: mask(:,:), il, ih, jl, jh
	real(wp), intent(out) :: w
	integer :: i, j
	real(wp) :: ylo, yhi, t_in, t_out, frac

	yhi = y(jmin)
	ylo = y(jmax)

	! Bottom edge
	do i = imin, imax
		if (mask(i,jmin) /= 1) cycle
		if (jmin > jl) then
			t_in  = temp(i, jmin, nk)
			t_out = temp(i, jmin-1, nk)
			if (t_in > tsolid .and. t_out <= tsolid .and. abs(t_in - t_out) > 1.0_wp) then
				frac = (t_in - tsolid) / (t_in - t_out)
				ylo = min(ylo, y(jmin) - frac * (y(jmin) - y(jmin-1)))
			else
				ylo = min(ylo, y(jmin))
			endif
		else
			ylo = min(ylo, y(jmin))
		endif
	enddo

	! Top edge
	do i = imin, imax
		if (mask(i,jmax) /= 1) cycle
		if (jmax < jh) then
			t_in  = temp(i, jmax, nk)
			t_out = temp(i, jmax+1, nk)
			if (t_in > tsolid .and. t_out <= tsolid .and. abs(t_in - t_out) > 1.0_wp) then
				frac = (t_in - tsolid) / (t_in - t_out)
				yhi = max(yhi, y(jmax) + frac * (y(jmax+1) - y(jmax)))
			else
				yhi = max(yhi, y(jmax))
			endif
		else
			yhi = max(yhi, y(jmax))
		endif
	enddo

	w = max(0.0_wp, yhi - ylo)
end subroutine interp_pool_width

!********************************************************************
subroutine interp_pool_depth(mask, dep)
! Refine depth by interpolating tsolid crossing at kmin boundary.
	integer, intent(in) :: mask(:,:)
	real(wp), intent(out) :: dep
	integer :: i, j
	real(wp) :: t_in, t_out, frac, zbot

	dep = z(nk) - z(kmin)

	if (kmin <= 2) return   ! already at domain bottom

	! Find maximum interpolated depth across all (i,j) in connected region
	do j = jmin, jmax
	do i = imin, imax
		if (mask(i,j) /= 1) cycle
		t_in  = temp(i, j, kmin)
		t_out = temp(i, j, kmin-1)
		if (t_in > tsolid .and. t_out <= tsolid .and. abs(t_in - t_out) > 1.0_wp) then
			frac = (t_in - tsolid) / (t_in - t_out)
			zbot = z(kmin) - frac * (z(kmin) - z(kmin-1))
			dep = max(dep, z(nk) - zbot)
		endif
	enddo
	enddo
end subroutine interp_pool_depth

end module dimensions
