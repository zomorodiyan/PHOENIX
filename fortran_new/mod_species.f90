!______________________________________________________________________________
!
module species
!______________________________________________________________________________
! Species transport for dissimilar material mixing.
! When species_flag=1, solves concentration equation once per timestep
! after iter_loop exits. When species_flag=0, module is inert.
!
	use precision
	use constant
	use geometry
	use initialization
	use parameters
	use dimensions
	use cfd_utils
	use residue, only: resorc

	implicit none

	! species_flag is declared in mod_param.f90 and read from input_param.txt

	! Secondary material properties (named constants)
	real(wp), parameter :: dens2    = 8880.0_wp
	real(wp), parameter :: denl2    = 7800.0_wp
	real(wp), parameter :: viscos2  = 0.003_wp
	real(wp), parameter :: tsolid2  = 1728.0_wp
	real(wp), parameter :: tliquid2 = 1803.0_wp
	real(wp), parameter :: tboiling2= 3650.0_wp
	real(wp), parameter :: acpa2    = 0.3441_wp
	real(wp), parameter :: acpb2    = 400.0_wp
	real(wp), parameter :: acpl2    = 800.0_wp
	real(wp), parameter :: thconsa2 = 0.0205_wp
	real(wp), parameter :: thconsb2 = 10.0_wp
	real(wp), parameter :: thconl2  = 120.0_wp

	! Secondary powder properties
	real(wp), parameter :: pden2    = 7330.0_wp
	real(wp), parameter :: pcpa2    = 0.3508_wp
	real(wp), parameter :: pcpb2    = 457.7_wp
	real(wp), parameter :: pthcona2 = 0.0_wp
	real(wp), parameter :: pthconb2 = 0.795_wp

	! Species transport parameters
	real(wp), parameter :: D_m = 5.0e-9_wp         ! molecular mass diffusivity (m^2/s)
	real(wp), parameter :: urfspecies = 0.7_wp      ! under-relaxation factor
	real(wp), parameter :: dgdc_const = -0.3_wp      ! dg/dC (solutal Marangoni coefficient, N/m)
	real(wp), parameter :: diff_floor = 1.0e-30_wp  ! diffusivity floor in solid

	! Derived constants for secondary material (set in init_species)
	real(wp), save :: hsmelt2, hlcal2, hlfriz2, cpavg2, deltemp2

	! resorc (species residual) is declared in mod_resid.f90

	! Allocatable fields
	real(wp), allocatable :: concentration(:,:,:)
	real(wp), allocatable :: conc_old(:,:,:)

	contains

!********************************************************************
! Linear mixing helper: prop1*C + prop2*(1-C)
!********************************************************************
pure real(wp) function mix(prop1, prop2, C)
	real(wp), intent(in) :: prop1, prop2, C
	mix = prop1 * C + prop2 * (1.0_wp - C)
end function mix

!********************************************************************
subroutine allocate_species
	allocate(concentration(ni, nj, nk))
	allocate(conc_old(ni, nj, nk))
	concentration = 1.0_wp
	conc_old = 1.0_wp
end subroutine allocate_species

!********************************************************************
subroutine init_species
	integer :: i, j, k
	real(wp) :: z_powder_top, y_mid

	! Derived constants for secondary material
	deltemp2 = tliquid2 - tsolid2
	cpavg2   = (acpa2*tsolid2 + acpb2 + acpl2) * 0.5_wp
	hsmelt2  = acpa2*tsolid2**2/2.0_wp + acpb2*tsolid2
	hlcal2   = hsmelt2 + cpavg2*deltemp2
	hlfriz2  = hlcal2 + hlatnt

	! Initial conditions: substrate=1 (base), powder half-domain=0 (secondary)
	z_powder_top = dimz                    ! top of domain = top of powder
	y_mid = dimy * 0.5_wp                  ! half of domain in y

	do k=1,nk
	do j=1,nj
	do i=1,ni
		if (z(k) >= (dimz - layerheight) .and. y(j) < y_mid) then
			concentration(i,j,k) = 0.0_wp   ! secondary material in powder half
		else
			concentration(i,j,k) = 1.0_wp   ! base material
		endif
	enddo
	enddo
	enddo

	conc_old = concentration
end subroutine init_species

!********************************************************************
subroutine species_bc
	integer :: i, j, k

	! Zero-flux Neumann BCs on all 6 faces
	! k faces
	do j=2,njm1
	do i=2,nim1
		concentration(i,j,nk) = concentration(i,j,nk-1)
		concentration(i,j,1)  = concentration(i,j,2)
	enddo
	enddo

	! j faces
	concentration(:,1,:)  = concentration(:,2,:)
	concentration(:,nj,:) = concentration(:,nj-1,:)

	! i faces
	concentration(1,:,:)  = concentration(2,:,:)
	concentration(ni,:,:) = concentration(ni-1,:,:)
end subroutine species_bc

!********************************************************************
subroutine solve_species
! Solves species transport equation once per timestep.
! Loops over melt pool region only (istatp1:iendm1, jstat:jend, kstat:nkm1),
! same indices as momentum solver. Outside melt pool, concentration stays unchanged.
! Reuses shared coefficient arrays (an,as,ae,aw,at,ab,ap,su,sp,apnot).
	integer :: i, j, k
	real(wp) :: vn, ue, wt, fn, fe, ft
	real(wp) :: vs, uw, wb, fs, fw, fb
	real(wp) :: gamma_face, dn, de, dt, ds, dw, db, tmp1
	real(wp) :: difn, dife, dift

	! Skip if no melt pool
	if (tpeak <= tsolid) then
		resorc = 0.0_wp
		return
	endif

	! --- Discretization (power law scheme) ---
	! Main loop: north, east, top faces + transient + source
	do k=kstat,nkm1
!$OMP PARALLEL PRIVATE(vn, ue, wt, fn, fe, ft, gamma_face, difn, dife, dift, dn, de, dt, tmp1)
!$OMP DO
	do j=jstat,jend
	do i=istatp1,iendm1
		vn = vVel(i,j+1,k)
		ue = uVel(i+1,j,k)
		wt = wVel(i,j,k+1)

		! Convection coefficients (use raw velocities, no x2 factor)
		fn = den(i,j,k)*vn*areaik(i,k)
		fe = den(i,j,k)*ue*areajk(j,k)
		ft = den(i,j,k)*wt*areaij(i,j)

		! Diffusion coefficient: Gamma = den * D_m (or floor in solid)
		gamma_face = den(i,j,k) * D_m
		if (temp(i,j,k) < tsolid) gamma_face = diff_floor

		if (j == jend) then
			difn = gamma_face
		else
			difn = harmonic_mean(gamma_face, den(i,j+1,k)*D_m, fracy(j))
			if (temp(i,j+1,k) < tsolid) difn = diff_floor
		endif
		if (i == iendm1) then
			dife = gamma_face
		else
			dife = harmonic_mean(gamma_face, den(i+1,j,k)*D_m, fracx(i))
			if (temp(i+1,j,k) < tsolid) dife = diff_floor
		endif
		if (k == nkm1) then
			dift = gamma_face
		else
			dift = harmonic_mean(gamma_face, den(i,j,k+1)*D_m, fracz(k))
			if (temp(i,j,k+1) < tsolid) dift = diff_floor
		endif

		dn = difn*areaik(i,k)*dypsinv(j+1)
		de = dife*areajk(j,k)*dxpwinv(i+1)
		dt = dift*areaij(i,j)*dzpbinv(k+1)

		! Power law scheme
		tmp1 = dn*max(0.0_wp,(1.0_wp - 0.1_wp*(abs(fn)/dn))**5)
		an(i,j,k) = tmp1 + max(0.0_wp, -fn)
		as(i,j+1,k) = tmp1 + max(0.0_wp, fn)

		tmp1 = de*max(0.0_wp,(1.0_wp - 0.1_wp*(abs(fe)/de))**5)
		ae(i,j,k) = tmp1 + max(0.0_wp, -fe)
		aw(i+1,j,k) = tmp1 + max(0.0_wp, fe)

		tmp1 = dt*max(0.0_wp,(1.0_wp - 0.1_wp*(abs(ft)/dt))**5)
		at(i,j,k) = tmp1 + max(0.0_wp, -ft)
		ab(i,j,k+1) = tmp1 + max(0.0_wp, ft)

		! Transient term (use delt, not delt_eff)
		apnot(i,j,k) = den(i,j,k)/delt * volume(i,j,k)

		sp(i,j,k) = 0.0_wp
		su(i,j,k) = apnot(i,j,k) * conc_old(i,j,k)
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	! South face boundary (j=jstat)
	j = jstat
	do k=kstat,nkm1
		do i=istatp1,iendm1
			vs = vVel(i,j,k)
			fs = den(i,j,k)*vs*areaik(i,k)
			ds = den(i,j,k)*D_m*areaik(i,k)*dypsinv(j)
			if (temp(i,j,k) < tsolid) ds = diff_floor*areaik(i,k)*dypsinv(j)
			as(i,j,k) = ds*max(0.0_wp,(1.0_wp - 0.1_wp*(abs(fs)/ds))**5) + max(0.0_wp, fs)
		enddo
	enddo

	! West face boundary (i=istatp1)
	i = istatp1
	do k=kstat,nkm1
		do j=jstat,jend
			uw = uVel(i,j,k)
			fw = den(i,j,k)*uw*areajk(j,k)
			dw = den(i,j,k)*D_m*areajk(j,k)*dxpwinv(i)
			if (temp(i,j,k) < tsolid) dw = diff_floor*areajk(j,k)*dxpwinv(i)
			aw(i,j,k) = dw*max(0.0_wp,(1.0_wp - 0.1_wp*(abs(fw)/dw))**5) + max(0.0_wp, fw)
		enddo
	enddo

	! Bottom face boundary (k=kstat)
	k = kstat
	do j=jstat,jend
		do i=istatp1,iendm1
			wb = wVel(i,j,k)
			fb = den(i,j,k)*wb*areaij(i,j)
			db = den(i,j,k)*D_m*areaij(i,j)*dzpbinv(k)
			if (temp(i,j,k) < tsolid) db = diff_floor*areaij(i,j)*dzpbinv(k)
			ab(i,j,k) = db*max(0.0_wp,(1.0_wp - 0.1_wp*(abs(fb)/db))**5) + max(0.0_wp, fb)
		enddo
	enddo

	! --- Boundary coefficient transfers ---
	! k boundaries
	do j=jstat,jend
	do i=istatp1,iendm1
		su(i,j,kstat) = su(i,j,kstat) + ab(i,j,kstat)*concentration(i,j,kstat-1)
		sp(i,j,kstat) = sp(i,j,kstat) - ab(i,j,kstat)
		ab(i,j,kstat) = 0.0_wp
		su(i,j,nkm1) = su(i,j,nkm1) + at(i,j,nkm1)*concentration(i,j,nk)
		sp(i,j,nkm1) = sp(i,j,nkm1) - at(i,j,nkm1)
		at(i,j,nkm1) = 0.0_wp
	enddo
	enddo

	! j boundaries
	do k=kstat,nkm1
	do i=istatp1,iendm1
		su(i,jstat,k) = su(i,jstat,k) + as(i,jstat,k)*concentration(i,jstat-1,k)
		sp(i,jstat,k) = sp(i,jstat,k) - as(i,jstat,k)
		as(i,jstat,k) = 0.0_wp
		su(i,jend,k) = su(i,jend,k) + an(i,jend,k)*concentration(i,jend+1,k)
		sp(i,jend,k) = sp(i,jend,k) - an(i,jend,k)
		an(i,jend,k) = 0.0_wp
	enddo
	enddo

	! i boundaries
	do k=kstat,nkm1
	do j=jstat,jend
		su(istatp1,j,k) = su(istatp1,j,k) + aw(istatp1,j,k)*concentration(istat,j,k)
		sp(istatp1,j,k) = sp(istatp1,j,k) - aw(istatp1,j,k)
		aw(istatp1,j,k) = 0.0_wp
		su(iendm1,j,k) = su(iendm1,j,k) + ae(iendm1,j,k)*concentration(iendm1+1,j,k)
		sp(iendm1,j,k) = sp(iendm1,j,k) - ae(iendm1,j,k)
		ae(iendm1,j,k) = 0.0_wp
	enddo
	enddo

	! --- Assembly and under-relaxation ---
	do k=kstat,nkm1
!$OMP PARALLEL
!$OMP DO
	do j=jstat,jend
	do i=istatp1,iendm1
		ap(i,j,k) = an(i,j,k)+as(i,j,k)+ae(i,j,k)+aw(i,j,k)+at(i,j,k)+ab(i,j,k)+apnot(i,j,k)-sp(i,j,k)
		ap(i,j,k) = ap(i,j,k)/urfspecies
		su(i,j,k) = su(i,j,k) + (1.0_wp - urfspecies)*ap(i,j,k)*concentration(i,j,k)
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	! --- TDMA solve (same bounds as momentum) ---
	call solution_species_tdma

	! --- Block correction ---
	call enhance_species_speed

	! --- Concentration clipping ---
	do k=kstat,nkm1
!$OMP PARALLEL
!$OMP DO
	do j=jstat,jend
	do i=istatp1,iendm1
		concentration(i,j,k) = max(0.0_wp, min(1.0_wp, concentration(i,j,k)))
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	! --- Residual ---
	call calc_species_residual

end subroutine solve_species

!********************************************************************
subroutine solution_species_tdma
! Line-by-line TDMA solver for concentration field.
! Uses melt pool bounds (same as solution_uvw).
	use omp_lib
	integer :: i, j, k, ksweep, nthreads, tid
	real(wp) :: d, denom
	real(wp), allocatable :: prbuf(:,:), qrbuf(:,:)
	integer :: nxi

	nxi = ni

	!$OMP PARALLEL
	!$OMP MASTER
	nthreads = omp_get_num_threads()
	!$OMP END MASTER
	!$OMP END PARALLEL
	if (nthreads < 1) nthreads = 1

	allocate(prbuf(nxi, nthreads))
	allocate(qrbuf(nxi, nthreads))

	!$OMP PARALLEL DEFAULT(NONE) &
	!$OMP SHARED(concentration, prbuf, qrbuf, &
	!$OMP        kstat, nkm1, jstat, jend, istat, istatp1, iendm1, &
	!$OMP        at, ab, an, as, ae, aw, ap, su) &
	!$OMP PRIVATE(tid, k, j, i, ksweep, d, denom)
	tid = omp_get_thread_num() + 1
	do ksweep = 1, 2
		do k = nkm1, kstat, -1
			!$OMP DO SCHEDULE(STATIC)
			do j = jstat, jend
				prbuf(istat, tid) = 0.0_wp
				qrbuf(istat, tid) = concentration(istat, j, k)

				do i = istatp1, iendm1
					d = at(i,j,k)*concentration(i,j,k+1) + ab(i,j,k)*concentration(i,j,k-1) &
						+ an(i,j,k)*concentration(i,j+1,k) + as(i,j,k)*concentration(i,j-1,k) + su(i,j,k)
					denom = ap(i,j,k) - aw(i,j,k)*prbuf(i-1, tid)
					if (abs(denom) < 1.0e-12_wp) denom = sign(1.0e-12_wp, denom)
					prbuf(i, tid) = ae(i,j,k) / denom
					qrbuf(i, tid) = (d + aw(i,j,k)*qrbuf(i-1, tid)) / denom
				enddo

				do i = iendm1, istatp1, -1
					concentration(i,j,k) = prbuf(i, tid)*concentration(i+1,j,k) + qrbuf(i, tid)
				enddo
			enddo
			!$OMP END DO
			!$OMP BARRIER
		enddo
	enddo
	!$OMP END PARALLEL

	deallocate(prbuf, qrbuf)
end subroutine solution_species_tdma

!********************************************************************
subroutine enhance_species_speed
! Block correction for species: accumulate j-k planes into 1D system
! along x (melt pool region), solve with TDMA, broadcast correction.
	integer :: i, j, k
	real(wp), allocatable :: bl(:), blp(:), blm(:), blc(:), delh(:), pib(:), qib(:)
	real(wp) :: denom

	allocate(bl(ni), blp(ni), blm(ni), blc(ni), delh(ni), pib(ni), qib(ni))

	bl  = 0.0_wp
	blp = 0.0_wp
	blm = 0.0_wp
	blc = 0.0_wp

	do k=kstat,nkm1
!$OMP PARALLEL
!$OMP DO REDUCTION(+: bl, blp, blm, blc)
	do j=jstat,jend
	do i=istatp1,iendm1
		bl(i)  = bl(i) + ap(i,j,k) - an(i,j,k) - as(i,j,k) - at(i,j,k) - ab(i,j,k)
		blp(i) = blp(i) + ae(i,j,k)
		blm(i) = blm(i) + aw(i,j,k)
		blc(i) = blc(i) + ae(i,j,k)*concentration(i+1,j,k) + aw(i,j,k)*concentration(i-1,j,k) &
			+ an(i,j,k)*concentration(i,j+1,k) + as(i,j,k)*concentration(i,j-1,k) &
			+ at(i,j,k)*concentration(i,j,k+1) + ab(i,j,k)*concentration(i,j,k-1) &
			+ su(i,j,k) - ap(i,j,k)*concentration(i,j,k)
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	! Forward sweep
	pib(istatp1) = blp(istatp1)/bl(istatp1)
	qib(istatp1) = blc(istatp1)/bl(istatp1)

	do i=istatp1+1,iendm1
		denom = bl(i) - blm(i)*pib(i-1)
		pib(i) = blp(i)/denom
		qib(i) = (blc(i) + blm(i)*qib(i-1))/denom
	enddo

	! Back substitution
	delh(iendm1) = qib(iendm1)
	do i=iendm1-1,istatp1,-1
		delh(i) = pib(i)*delh(i+1) + qib(i)
	enddo

	! Apply correction
	do k=kstat,nkm1
!$OMP PARALLEL
!$OMP DO
	do j=jstat,jend
	do i=istatp1,iendm1
		concentration(i,j,k) = concentration(i,j,k) + delh(i)
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	deallocate(bl, blp, blm, blc, delh, pib, qib)
end subroutine enhance_species_speed

!********************************************************************
subroutine calc_species_residual
	integer :: i, j, k
	real(wp) :: sumc, sumd, resor

	sumc = 0.0_wp
	sumd = 0.0_wp

	do k=kstat,nkm1
!$OMP PARALLEL PRIVATE(resor)
!$OMP DO REDUCTION(+: sumd, sumc)
	do j=jstat,jend
	do i=istatp1,iendm1
		resor = (an(i,j,k)*concentration(i,j+1,k) + as(i,j,k)*concentration(i,j-1,k) &
			+ ae(i,j,k)*concentration(i+1,j,k) + aw(i,j,k)*concentration(i-1,j,k) &
			+ at(i,j,k)*concentration(i,j,k+1) + ab(i,j,k)*concentration(i,j,k-1) &
			+ su(i,j,k))/ap(i,j,k) - concentration(i,j,k)
		sumd = sumd + abs(resor)
		sumc = sumc + abs(concentration(i,j,k))
	enddo
	enddo
!$OMP END PARALLEL
	enddo

	resorc = sumd / (sumc + small)
end subroutine calc_species_residual

!********************************************************************
subroutine write_species_vtk(outputnum)
! Write standalone species VTK file: {file_prefix}species{N}.vtk
	integer, intent(in) :: outputnum
	integer :: i, j, k, gridx, gridy, gridz, npts
	character(len=3) :: cTemp
	character(len=512) :: fname
	real(kind=4) :: val4

	write(cTemp, '(i3)') outputnum
	fname = trim(file_prefix) // 'species' // trim(adjustl(cTemp)) // '.vtk'

	gridx = nim1 - 2 + 1  ! 2:nim1
	gridy = njm1 - 2 + 1
	gridz = nkm1 - 2 + 1
	npts = gridx * gridy * gridz

	! ASCII header
	open(unit=42, file=trim(fname))
	write(42,'(A)') '# vtk DataFile Version 3.0'
	write(42,'(A)') 'AMCFD Species Concentration'
	write(42,'(A)') 'BINARY'
	write(42,'(A)') 'DATASET STRUCTURED_GRID'
	write(42,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', gridx, ' ', gridy, ' ', gridz
	write(42,'(A,I0,A)') 'POINTS ', npts, ' float'
	close(42)

	! Binary coordinates
	open(unit=42, file=trim(fname), access='stream', form='unformatted', &
	     position='append', convert='big_endian')
	do k=2,nkm1
	do j=2,njm1
	do i=2,nim1
		val4 = real(x(i), 4); write(42) val4
		val4 = real(y(j), 4); write(42) val4
		val4 = real(z(k), 4); write(42) val4
	enddo
	enddo
	enddo
	close(42)

	! POINT_DATA header (ASCII)
	open(unit=42, file=trim(fname), position='append')
	write(42,'(A,I0)') 'POINT_DATA ', npts
	write(42,'(A)') 'SCALARS concentration float 1'
	write(42,'(A)') 'LOOKUP_TABLE default'
	close(42)

	! Binary concentration data
	open(unit=42, file=trim(fname), access='stream', form='unformatted', &
	     position='append', convert='big_endian')
	do k=2,nkm1
	do j=2,njm1
	do i=2,nim1
		val4 = real(concentration(i,j,k), 4)
		write(42) val4
	enddo
	enddo
	enddo
	close(42)
end subroutine write_species_vtk

end module species
