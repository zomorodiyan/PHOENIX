!______________________________________________________________________________
!
module mechanical_solver
!______________________________________________________________________________
! EBE (Element-By-Element) FEM mechanical solver for residual stress.
! Ported from fortran-ebe-lpbf with adaptations for PHOENIX grid system.
! Supports non-uniform grid spacing (AMR compatible).
! One-way coupling: T -> thermal strain -> displacement -> stress.
!
	use precision
	use geometry, only: x, y, z, ni, nj, nk, nim1, njm1, nkm1
	use parameters, only: wp, tempPreheat
	use mech_material
	implicit none

	private
	public :: init_mechanical, cleanup_mechanical, solve_mechanical, get_stress_yield
	public :: Nnx, Nny, Nnz, Nx, Ny, Nz

	! Grid dimensions for FEM (set in init_mechanical)
	! FEM nodes = PHOENIX interior cell centers: x(2:nim1), y(2:njm1), z(2:nkm1)
	! Nnx = nim1-1, Nny = njm1-1, Nnz = nkm1-1 (number of FEM nodes)
	! Nx = Nnx-1, Ny = Nny-1, Nz = Nnz-1 (number of FEM elements)
	integer :: Nnx, Nny, Nnz  ! FEM node counts
	integer :: Nx, Ny, Nz      ! FEM element counts
	integer :: inode           ! temp variable for shape function loop

	! Precomputed B matrices at each GP for each unique dz layer: B_layer(6,24,8,n_dz_layers)
	! For uniform x-y (AMR refined region may have varying dx/dy — use per-element B if needed)
	! Simplified: precompute per-element Jacobian using actual node spacing
	real(wp) :: N_gp(8,8)           ! Shape functions at 8 GPs

	! Precomputed element stiffness per Z-layer (dx/dy uniform, dz varies per layer)
	real(wp), allocatable :: Ke_solid_z(:,:,:)  ! (24,24,Nz)
	real(wp), allocatable :: Ke_soft_z(:,:,:)   ! (24,24,Nz)
	real(wp) :: Ke_solid(24,24), Ke_soft(24,24) ! reference (kept for compatibility)

	! GP-level state arrays: (6 Voigt components, 8 GPs, Nx, Ny, Nz)
	real(wp), allocatable :: sig_gp(:,:,:,:,:)
	real(wp), allocatable :: eps_gp(:,:,:,:,:)

	! Temperature reference for incremental thermal strain
	real(wp), allocatable :: T_old_mech(:,:,:)

	! Yield function output (node-centered)
	real(wp), allocatable :: f_plus(:,:,:)

	! Element phase cache
	integer, allocatable :: elem_phase(:,:,:)

	! Displacement fields (persistent across solves)
	real(wp), allocatable, public :: ux_mech(:,:,:), uy_mech(:,:,:), uz_mech(:,:,:)

	! Stress output arrays (node-centered, smoothed)
	real(wp), allocatable, public :: sxx_out(:,:,:), syy_out(:,:,:), szz_out(:,:,:)
	real(wp), allocatable, public :: vm_out(:,:,:), fplus_out(:,:,:)

	! Mechanical phase (public for VTK output)
	integer, allocatable, public :: mech_phase(:,:,:)

	! Last extracted FEM temperature (public for VTK/history output)
	real(wp), allocatable, public :: T_fem_last(:,:,:)

	! FEM node coordinates (coarsened from PHOENIX grid by mech_mesh_ratio)
	real(wp), allocatable, public :: fem_x(:), fem_y(:), fem_z(:)
	integer :: mratio  ! stored copy of mech_mesh_ratio

	contains

!********************************************************************
subroutine init_mechanical()
	use parameters, only: mech_mesh_ratio
	integer :: i, j, k, g
	real(wp) :: gp_c(2), xi, eta_q, zeta
	real(wp) :: xi_n(8), eta_n(8), zeta_n(8)
	real(wp) :: dx_ref, dy_ref, dz_ref
	integer :: n_thermal_x, n_thermal_y, n_thermal_z

	mratio = mech_mesh_ratio

	! Thermal interior cell count
	n_thermal_x = nim1 - 1   ! = ncvx(1)
	n_thermal_y = njm1 - 1
	n_thermal_z = nkm1 - 1

	! FEM grid: coarsen by mratio (take every mratio-th thermal cell center)
	Nnx = n_thermal_x / mratio + 1
	Nny = n_thermal_y / mratio + 1
	Nnz = n_thermal_z / mratio + 1
	Nx = Nnx - 1
	Ny = Nny - 1
	Nz = Nnz - 1

	! Ensure FEM nodes stay within PHOENIX interior (2 to nim1, etc.)
	! Last FEM node index in PHOENIX: 2 + (Nnx-1)*mratio must be <= nim1
	if (2 + (Nnx-1)*mratio > nim1) Nnx = (nim1 - 2) / mratio + 1
	if (2 + (Nny-1)*mratio > njm1) Nny = (njm1 - 2) / mratio + 1
	if (2 + (Nnz-1)*mratio > nkm1) Nnz = (nkm1 - 2) / mratio + 1
	Nx = Nnx - 1; Ny = Nny - 1; Nz = Nnz - 1

	! Build FEM node coordinates from PHOENIX grid
	! Force last node to be at the interior boundary (nim1, njm1, nkm1)
	! to capture the surface layer
	allocate(fem_x(Nnx), fem_y(Nny), fem_z(Nnz))
	do i = 1, Nnx
		fem_x(i) = x(min(2 + (i-1)*mratio, nim1))
	enddo
	do j = 1, Nny
		fem_y(j) = y(min(2 + (j-1)*mratio, njm1))
	enddo
	do k = 1, Nnz
		fem_z(k) = z(min(2 + (k-1)*mratio, nkm1))
	enddo
	! Ensure last node IS at the surface
	fem_x(Nnx) = x(nim1)
	fem_y(Nny) = y(njm1)
	fem_z(Nnz) = z(nkm1)

	! Allocate state arrays
	allocate(sig_gp(6, 8, Nx, Ny, Nz))
	allocate(eps_gp(6, 8, Nx, Ny, Nz))
	allocate(T_old_mech(Nnx, Nny, Nnz))
	allocate(f_plus(Nnx, Nny, Nnz))
	allocate(elem_phase(Nx, Ny, Nz))
	allocate(ux_mech(Nnx, Nny, Nnz))
	allocate(uy_mech(Nnx, Nny, Nnz))
	allocate(uz_mech(Nnx, Nny, Nnz))
	allocate(sxx_out(Nnx, Nny, Nnz))
	allocate(syy_out(Nnx, Nny, Nnz))
	allocate(szz_out(Nnx, Nny, Nnz))
	allocate(vm_out(Nnx, Nny, Nnz))
	allocate(fplus_out(Nnx, Nny, Nnz))
	allocate(mech_phase(Nnx, Nny, Nnz))
	allocate(T_fem_last(Nnx, Nny, Nnz))

	sig_gp = 0.0_wp
	eps_gp = 0.0_wp
	T_old_mech = tempPreheat
	f_plus = 0.0_wp
	elem_phase = MECH_POWDER
	ux_mech = 0.0_wp; uy_mech = 0.0_wp; uz_mech = 0.0_wp
	sxx_out = 0.0_wp; syy_out = 0.0_wp; szz_out = 0.0_wp
	vm_out = 0.0_wp; fplus_out = 0.0_wp
	mech_phase = MECH_POWDER
	T_fem_last = tempPreheat

	! Precompute shape functions at 8 GPs
	xi_n   = (/ -1.0_wp, 1.0_wp,-1.0_wp, 1.0_wp,-1.0_wp, 1.0_wp,-1.0_wp, 1.0_wp /)
	eta_n  = (/ -1.0_wp,-1.0_wp, 1.0_wp, 1.0_wp,-1.0_wp,-1.0_wp, 1.0_wp, 1.0_wp /)
	zeta_n = (/ -1.0_wp,-1.0_wp,-1.0_wp,-1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp /)

	gp_c(1) = -1.0_wp / sqrt(3.0_wp)
	gp_c(2) =  1.0_wp / sqrt(3.0_wp)

	N_gp = 0.0_wp
	g = 0
	do k = 1, 2
	do j = 1, 2
	do i = 1, 2
		g = g + 1
		xi = gp_c(i); eta_q = gp_c(j); zeta = gp_c(k)
		do inode = 1, 8
			N_gp(inode, g) = (1.0_wp + xi_n(inode)*xi) * (1.0_wp + eta_n(inode)*eta_q) &
			                * (1.0_wp + zeta_n(inode)*zeta) / 8.0_wp
		enddo
	enddo
	enddo
	enddo

	! Precompute Ke for each Z-layer (dx/dy uniform within layer, dz varies)
	dx_ref = fem_x(2) - fem_x(1)
	dy_ref = fem_y(2) - fem_y(1)
	dz_ref = fem_z(2) - fem_z(1)
	call compute_Ke_uniform(E_solid, nu_mech, dx_ref, dy_ref, dz_ref, Ke_solid)
	call compute_Ke_uniform(E_soft,  nu_mech, dx_ref, dy_ref, dz_ref, Ke_soft)

	allocate(Ke_solid_z(24,24,Nz), Ke_soft_z(24,24,Nz))
	do k = 1, Nz
		dz_ref = fem_z(k+1) - fem_z(k)
		call compute_Ke_uniform(E_solid, nu_mech, dx_ref, dy_ref, dz_ref, Ke_solid_z(:,:,k))
		call compute_Ke_uniform(E_soft,  nu_mech, dx_ref, dy_ref, dz_ref, Ke_soft_z(:,:,k))
	enddo

end subroutine init_mechanical

!********************************************************************
subroutine cleanup_mechanical()
	if (allocated(sig_gp))      deallocate(sig_gp)
	if (allocated(eps_gp))      deallocate(eps_gp)
	if (allocated(T_old_mech))  deallocate(T_old_mech)
	if (allocated(f_plus))      deallocate(f_plus)
	if (allocated(elem_phase))  deallocate(elem_phase)
	if (allocated(ux_mech))     deallocate(ux_mech, uy_mech, uz_mech)
	if (allocated(sxx_out))     deallocate(sxx_out, syy_out, szz_out, vm_out, fplus_out)
	if (allocated(mech_phase))  deallocate(mech_phase)
	if (allocated(fem_x))       deallocate(fem_x, fem_y, fem_z)
	if (allocated(Ke_solid_z))  deallocate(Ke_solid_z, Ke_soft_z)
	if (allocated(T_fem_last))  deallocate(T_fem_last)
end subroutine cleanup_mechanical

!********************************************************************
subroutine compute_Ke_uniform(E_val, nu_val, dxe, dye, dze, Ke)
! Compute 24x24 element stiffness matrix for a uniform hex element.
	real(wp), intent(in)  :: E_val, nu_val, dxe, dye, dze
	real(wp), intent(out) :: Ke(24,24)
	real(wp) :: C(6,6), gp_c(2)
	real(wp) :: xi, eta_q, zeta, w, det_J
	real(wp) :: xi_n(8), eta_n(8), zeta_n(8)
	real(wp) :: dN_dx(8), dN_dy(8), dN_dz(8)
	real(wp) :: dN_dxi(8), dN_deta(8), dN_dzeta(8)
	real(wp) :: B_a(6,3), CB(6,3)
	integer  :: i, j, k, a, b, p, q

	call build_C_matrix(E_val, nu_val, C)

	xi_n   = (/ -1.0_wp, 1.0_wp,-1.0_wp, 1.0_wp,-1.0_wp, 1.0_wp,-1.0_wp, 1.0_wp /)
	eta_n  = (/ -1.0_wp,-1.0_wp, 1.0_wp, 1.0_wp,-1.0_wp,-1.0_wp, 1.0_wp, 1.0_wp /)
	zeta_n = (/ -1.0_wp,-1.0_wp,-1.0_wp,-1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp /)

	gp_c(1) = -1.0_wp / sqrt(3.0_wp)
	gp_c(2) =  1.0_wp / sqrt(3.0_wp)

	det_J = (dxe/2.0_wp) * (dye/2.0_wp) * (dze/2.0_wp)

	Ke = 0.0_wp

	do k = 1, 2
	do j = 1, 2
	do i = 1, 2
		xi = gp_c(i); eta_q = gp_c(j); zeta = gp_c(k)
		w = 1.0_wp  ! weight = 1 for 2-point Gauss

		do a = 1, 8
			dN_dxi(a)   = xi_n(a)   * (1.0_wp + eta_n(a)*eta_q) * (1.0_wp + zeta_n(a)*zeta) / 8.0_wp
			dN_deta(a)  = (1.0_wp + xi_n(a)*xi) * eta_n(a)  * (1.0_wp + zeta_n(a)*zeta) / 8.0_wp
			dN_dzeta(a) = (1.0_wp + xi_n(a)*xi) * (1.0_wp + eta_n(a)*eta_q) * zeta_n(a) / 8.0_wp
		enddo

		dN_dx = dN_dxi   * (2.0_wp / dxe)
		dN_dy = dN_deta  * (2.0_wp / dye)
		dN_dz = dN_dzeta * (2.0_wp / dze)

		do b = 1, 8
			B_a = 0.0_wp
			B_a(1,1) = dN_dx(b); B_a(2,2) = dN_dy(b); B_a(3,3) = dN_dz(b)
			B_a(4,1) = dN_dy(b); B_a(4,2) = dN_dx(b)
			B_a(5,1) = dN_dz(b); B_a(5,3) = dN_dx(b)
			B_a(6,2) = dN_dz(b); B_a(6,3) = dN_dy(b)

			CB = 0.0_wp
			do q = 1, 3
			do p = 1, 6
				CB(p,q) = C(p,1)*B_a(1,q) + C(p,2)*B_a(2,q) + C(p,3)*B_a(3,q) &
				         + C(p,4)*B_a(4,q) + C(p,5)*B_a(5,q) + C(p,6)*B_a(6,q)
			enddo
			enddo

			do a = 1, 8
				B_a = 0.0_wp
				B_a(1,1) = dN_dx(a); B_a(2,2) = dN_dy(a); B_a(3,3) = dN_dz(a)
				B_a(4,1) = dN_dy(a); B_a(4,2) = dN_dx(a)
				B_a(5,1) = dN_dz(a); B_a(5,3) = dN_dx(a)
				B_a(6,2) = dN_dz(a); B_a(6,3) = dN_dy(a)

				do q = 1, 3
				do p = 1, 3
					Ke(3*(a-1)+p, 3*(b-1)+q) = Ke(3*(a-1)+p, 3*(b-1)+q) &
						+ (B_a(1,p)*CB(1,q) + B_a(2,p)*CB(2,q) + B_a(3,p)*CB(3,q) &
						 + B_a(4,p)*CB(4,q) + B_a(5,p)*CB(5,q) + B_a(6,p)*CB(6,q)) &
						* det_J * w
				enddo
				enddo
			enddo
		enddo
	enddo
	enddo
	enddo
end subroutine compute_Ke_uniform

!********************************************************************
subroutine get_B_at_gp(dxe, dye, dze, g, B_gp)
! Compute B matrix (6x24) at Gauss point g for element with spacing dxe,dye,dze.
	real(wp), intent(in) :: dxe, dye, dze
	integer, intent(in)  :: g
	real(wp), intent(out) :: B_gp(6,24)
	real(wp) :: gp_c(2), xi, eta_q, zeta
	real(wp) :: xi_n(8), eta_n(8), zeta_n(8)
	real(wp) :: dN_dx(8), dN_dy(8), dN_dz(8)
	real(wp) :: dN_dxi(8), dN_deta(8), dN_dzeta(8)
	integer :: i, j, k, a, gi, gj, gk

	xi_n   = (/ -1.0_wp, 1.0_wp,-1.0_wp, 1.0_wp,-1.0_wp, 1.0_wp,-1.0_wp, 1.0_wp /)
	eta_n  = (/ -1.0_wp,-1.0_wp, 1.0_wp, 1.0_wp,-1.0_wp,-1.0_wp, 1.0_wp, 1.0_wp /)
	zeta_n = (/ -1.0_wp,-1.0_wp,-1.0_wp,-1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp /)
	gp_c(1) = -1.0_wp / sqrt(3.0_wp)
	gp_c(2) =  1.0_wp / sqrt(3.0_wp)

	gi = mod(g-1,2) + 1; gj = mod((g-1)/2,2) + 1; gk = (g-1)/4 + 1
	xi = gp_c(gi); eta_q = gp_c(gj); zeta = gp_c(gk)

	do a = 1, 8
		dN_dxi(a)   = xi_n(a)   * (1.0_wp + eta_n(a)*eta_q) * (1.0_wp + zeta_n(a)*zeta) / 8.0_wp
		dN_deta(a)  = (1.0_wp + xi_n(a)*xi) * eta_n(a)  * (1.0_wp + zeta_n(a)*zeta) / 8.0_wp
		dN_dzeta(a) = (1.0_wp + xi_n(a)*xi) * (1.0_wp + eta_n(a)*eta_q) * zeta_n(a) / 8.0_wp
	enddo

	dN_dx = dN_dxi   * (2.0_wp / dxe)
	dN_dy = dN_deta  * (2.0_wp / dye)
	dN_dz = dN_dzeta * (2.0_wp / dze)

	B_gp = 0.0_wp
	do a = 1, 8
		B_gp(1, 3*(a-1)+1) = dN_dx(a)
		B_gp(2, 3*(a-1)+2) = dN_dy(a)
		B_gp(3, 3*(a-1)+3) = dN_dz(a)
		B_gp(4, 3*(a-1)+1) = dN_dy(a); B_gp(4, 3*(a-1)+2) = dN_dx(a)
		B_gp(5, 3*(a-1)+1) = dN_dz(a); B_gp(5, 3*(a-1)+3) = dN_dx(a)
		B_gp(6, 3*(a-1)+2) = dN_dz(a); B_gp(6, 3*(a-1)+3) = dN_dy(a)
	enddo
end subroutine get_B_at_gp

!********************************************************************
subroutine compute_elem_phases(phase)
	integer, intent(in) :: phase(Nnx, Nny, Nnz)
	integer :: ie, je, ke, di, dj, dk, n_solid

	!$OMP PARALLEL DO PRIVATE(ie,je,ke,di,dj,dk,n_solid) COLLAPSE(3)
	do ke = 1, Nz
	do je = 1, Ny
	do ie = 1, Nx
		n_solid = 0
		do dk = 0, 1
		do dj = 0, 1
		do di = 0, 1
			if (phase(ie+di, je+dj, ke+dk) == MECH_SOLID) n_solid = n_solid + 1
		enddo
		enddo
		enddo
		if (n_solid > 0) then
			elem_phase(ie,je,ke) = MECH_SOLID
		else
			elem_phase(ie,je,ke) = MECH_POWDER
		endif
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO
end subroutine compute_elem_phases

!********************************************************************
subroutine fem_node_x(ix, xval)
! Get x-coordinate of FEM node ix (1-based). Maps to PHOENIX x(ix+1).
	integer, intent(in)   :: ix
	real(wp), intent(out) :: xval
	xval = x(ix + 1)
end subroutine fem_node_x

!********************************************************************
subroutine get_elem_dx(ie, je, ke, dxe, dye, dze)
! Get element spacing from FEM node coordinates.
	integer, intent(in) :: ie, je, ke
	real(wp), intent(out) :: dxe, dye, dze
	dxe = fem_x(ie+1) - fem_x(ie)
	dye = fem_y(je+1) - fem_y(je)
	dze = fem_z(ke+1) - fem_z(ke)
end subroutine get_elem_dx

!********************************************************************
subroutine extract_temp_to_fem(temp_phoenix, T_fem)
! Extract PHOENIX cell-center temperature to coarsened FEM node array.
! Last node forced to nim1/njm1/nkm1 to capture surface.
	real(wp), intent(in)  :: temp_phoenix(:,:,:)
	real(wp), intent(out) :: T_fem(Nnx, Nny, Nnz)
	integer :: i, j, k, ip, jp, kp

	!$OMP PARALLEL DO PRIVATE(i,j,k,ip,jp,kp)
	do k = 1, Nnz
	do j = 1, Nny
	do i = 1, Nnx
		ip = min(2 + (i-1)*mratio, nim1)
		jp = min(2 + (j-1)*mratio, njm1)
		kp = min(2 + (k-1)*mratio, nkm1)
		T_fem(i,j,k) = temp_phoenix(ip, jp, kp)
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO
end subroutine extract_temp_to_fem

!********************************************************************
subroutine extract_solidfield_to_fem(sf_phoenix, sf_fem)
	real(wp), intent(in)  :: sf_phoenix(:,:,:)
	real(wp), intent(out) :: sf_fem(Nnx, Nny, Nnz)
	integer :: i, j, k, ip, jp, kp

	!$OMP PARALLEL DO PRIVATE(i,j,k,ip,jp,kp)
	do k = 1, Nnz
	do j = 1, Nny
	do i = 1, Nnx
		ip = min(2 + (i-1)*mratio, nim1)
		jp = min(2 + (j-1)*mratio, njm1)
		kp = min(2 + (k-1)*mratio, nkm1)
		sf_fem(i,j,k) = sf_phoenix(ip, jp, kp)
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO
end subroutine extract_solidfield_to_fem

!********************************************************************
subroutine compute_dT_gp(T_new, dT_gp_arr)
	real(wp), intent(in)  :: T_new(Nnx, Nny, Nnz)
	real(wp), intent(out) :: dT_gp_arr(8, Nx, Ny, Nz)
	integer :: ie, je, ke, di, dj, dk, ln, g, a
	real(wp) :: dT_nodes(8)

	!$OMP PARALLEL DO PRIVATE(ie,je,ke,di,dj,dk,ln,g,a,dT_nodes) COLLAPSE(3)
	do ke = 1, Nz
	do je = 1, Ny
	do ie = 1, Nx
		do dk = 0, 1
		do dj = 0, 1
		do di = 0, 1
			ln = 1 + di + 2*dj + 4*dk
			dT_nodes(ln) = T_new(ie+di, je+dj, ke+dk) - T_old_mech(ie+di, je+dj, ke+dk)
		enddo
		enddo
		enddo
		do g = 1, 8
			dT_gp_arr(g, ie, je, ke) = 0.0_wp
			do a = 1, 8
				dT_gp_arr(g, ie, je, ke) = dT_gp_arr(g, ie, je, ke) + N_gp(a, g) * dT_nodes(a)
			enddo
		enddo
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO
end subroutine compute_dT_gp

!********************************************************************
subroutine compute_residual(ux, uy, uz, dT_gp_arr, T_fem, phase, Rx, Ry, Rz)
	real(wp), intent(in)  :: ux(Nnx,Nny,Nnz), uy(Nnx,Nny,Nnz), uz(Nnx,Nny,Nnz)
	real(wp), intent(in)  :: dT_gp_arr(8, Nx, Ny, Nz)
	real(wp), intent(in)  :: T_fem(Nnx, Nny, Nnz)
	integer,  intent(in)  :: phase(Nnx,Nny,Nnz)
	real(wp), intent(out) :: Rx(Nnx,Nny,Nnz), Ry(Nnx,Nny,Nnz), Rz(Nnx,Nny,Nnz)

	real(wp) :: u_e(24), R_e(24), eps_curr(6), eps_inc(6), eps_th(6)
	real(wp) :: s_trial(6), sigma(6), BtSig(24), C_loc(6,6), B_gp(6,24)
	real(wp) :: av, dT_g, f_dum, dxe, dye, dze, det_J, T_gp_val
	integer  :: ie, je, ke, di, dj, dk, ln, dof, g, p, q, a
	integer  :: color, ic, jc, kc, ph_g

	Rx = 0.0_wp; Ry = 0.0_wp; Rz = 0.0_wp

	do color = 0, 7
		ic = mod(color, 2); jc = mod(color/2, 2); kc = mod(color/4, 2)

		!$OMP PARALLEL DO PRIVATE(ie,je,ke,u_e,R_e,eps_curr,eps_inc,eps_th,s_trial,sigma,BtSig) &
		!$OMP   PRIVATE(C_loc,B_gp,av,dT_g,f_dum,dxe,dye,dze,det_J,di,dj,dk,ln,dof,g,p,q,ph_g) COLLAPSE(3)
		do ke = 1 + kc, Nz, 2
		do je = 1 + jc, Ny, 2
		do ie = 1 + ic, Nx, 2

			call get_elem_dx(ie, je, ke, dxe, dye, dze)
			det_J = (dxe/2.0_wp) * (dye/2.0_wp) * (dze/2.0_wp)

			! Gather 24 DOFs
			do dk = 0, 1
			do dj = 0, 1
			do di = 0, 1
				ln = 1 + di + 2*dj + 4*dk
				dof = 3*(ln-1)
				u_e(dof+1) = ux(ie+di, je+dj, ke+dk)
				u_e(dof+2) = uy(ie+di, je+dj, ke+dk)
				u_e(dof+3) = uz(ie+di, je+dj, ke+dk)
			enddo
			enddo
			enddo

			R_e = 0.0_wp

			do g = 1, 8
				call get_B_at_gp(dxe, dye, dze, g, B_gp)

				! Strain at GP: eps = B * u_e
				do p = 1, 6
					eps_curr(p) = 0.0_wp
					do q = 1, 24
						eps_curr(p) = eps_curr(p) + B_gp(p, q) * u_e(q)
					enddo
				enddo

				eps_inc = eps_curr - eps_gp(:, g, ie, je, ke)

				! Phase at GP
				di = mod(g-1, 2); dj = mod((g-1)/2, 2); dk = (g-1) / 4
				ph_g = phase(ie+di, je+dj, ke+dk)

				if (ph_g == MECH_SOLID) then
					call build_C_matrix(E_solid, nu_mech, C_loc)
					av = alpha_V
				else
					call build_C_matrix(E_soft, nu_mech, C_loc)
					av = 0.0_wp
				endif

				dT_g = dT_gp_arr(g, ie, je, ke)
				eps_th = 0.0_wp
				eps_th(1) = av * dT_g; eps_th(2) = av * dT_g; eps_th(3) = av * dT_g

				! Trial stress
				do p = 1, 6
					s_trial(p) = sig_gp(p, g, ie, je, ke)
					do q = 1, 6
						s_trial(p) = s_trial(p) + C_loc(p,q) * (eps_inc(q) - eps_th(q))
					enddo
				enddo

				! Temperature at GP for yield strength
				T_gp_val = 0.0_wp
				do a = 1, 8
					di = mod(a-1,2); dj = mod((a-1)/2,2); dk = (a-1)/4
					T_gp_val = T_gp_val + N_gp(a, g) * T_fem(ie+di, je+dj, ke+dk)
				enddo
				call j2_return_map(s_trial, sigma, f_dum, T_gp_val)

				! R_e += B^T * sigma * detJ
				do p = 1, 24
					BtSig(p) = 0.0_wp
					do q = 1, 6
						BtSig(p) = BtSig(p) + B_gp(q, p) * sigma(q)
					enddo
				enddo
				R_e = R_e + BtSig * det_J
			enddo

			! Scatter
			do dk = 0, 1
			do dj = 0, 1
			do di = 0, 1
				ln = 1 + di + 2*dj + 4*dk
				dof = 3*(ln-1)
				Rx(ie+di, je+dj, ke+dk) = Rx(ie+di, je+dj, ke+dk) + R_e(dof+1)
				Ry(ie+di, je+dj, ke+dk) = Ry(ie+di, je+dj, ke+dk) + R_e(dof+2)
				Rz(ie+di, je+dj, ke+dk) = Rz(ie+di, je+dj, ke+dk) + R_e(dof+3)
			enddo
			enddo
			enddo
		enddo
		enddo
		enddo
		!$OMP END PARALLEL DO
	enddo

	! Dirichlet BC: bottom face (k=1) clamped
	Rx(:,:,1) = 0.0_wp; Ry(:,:,1) = 0.0_wp; Rz(:,:,1) = 0.0_wp
end subroutine compute_residual

!********************************************************************
subroutine ebe_matvec_mech(ux, uy, uz, Aux, Auy, Auz, phase)
	real(wp), intent(in)  :: ux(Nnx,Nny,Nnz), uy(Nnx,Nny,Nnz), uz(Nnx,Nny,Nnz)
	real(wp), intent(out) :: Aux(Nnx,Nny,Nnz), Auy(Nnx,Nny,Nnz), Auz(Nnx,Nny,Nnz)
	integer,  intent(in)  :: phase(Nnx,Nny,Nnz)

	real(wp) :: xe(24), Axe(24), Ke_loc(24,24)
	real(wp) :: dxe, dye, dze
	integer  :: ie, je, ke, a, b, di, dj, dk, ln, dof
	integer  :: color, ic, jc, kc

	Aux = 0.0_wp; Auy = 0.0_wp; Auz = 0.0_wp

	do color = 0, 7
		ic = mod(color, 2); jc = mod(color/2, 2); kc = mod(color/4, 2)

		!$OMP PARALLEL DO PRIVATE(ie,je,ke,xe,Axe,Ke_loc,dxe,dye,dze,a,b,di,dj,dk,ln,dof) COLLAPSE(3)
		do ke = 1 + kc, Nz, 2
		do je = 1 + jc, Ny, 2
		do ie = 1 + ic, Nx, 2
			! Use precomputed per-layer Ke (dx/dy uniform, dz varies per layer)
			if (elem_phase(ie,je,ke) /= MECH_POWDER) then
				Ke_loc = Ke_solid_z(:,:,ke)
			else
				Ke_loc = Ke_soft_z(:,:,ke)
			endif

			! Gather
			do dk = 0, 1
			do dj = 0, 1
			do di = 0, 1
				ln = 1 + di + 2*dj + 4*dk
				dof = 3*(ln-1)
				xe(dof+1) = ux(ie+di, je+dj, ke+dk)
				xe(dof+2) = uy(ie+di, je+dj, ke+dk)
				xe(dof+3) = uz(ie+di, je+dj, ke+dk)
			enddo
			enddo
			enddo

			! Ke * xe
			Axe = 0.0_wp
			do b = 1, 24
			do a = 1, 24
				Axe(a) = Axe(a) + Ke_loc(a,b) * xe(b)
			enddo
			enddo

			! Scatter
			do dk = 0, 1
			do dj = 0, 1
			do di = 0, 1
				ln = 1 + di + 2*dj + 4*dk
				dof = 3*(ln-1)
				Aux(ie+di, je+dj, ke+dk) = Aux(ie+di, je+dj, ke+dk) + Axe(dof+1)
				Auy(ie+di, je+dj, ke+dk) = Auy(ie+di, je+dj, ke+dk) + Axe(dof+2)
				Auz(ie+di, je+dj, ke+dk) = Auz(ie+di, je+dj, ke+dk) + Axe(dof+3)
			enddo
			enddo
			enddo
		enddo
		enddo
		enddo
		!$OMP END PARALLEL DO
	enddo

	! Dirichlet BC
	Aux(:,:,1) = ux(:,:,1); Auy(:,:,1) = uy(:,:,1); Auz(:,:,1) = uz(:,:,1)
end subroutine ebe_matvec_mech

!********************************************************************
subroutine solve_mech_cg(ux, uy, uz, fx, fy, fz, phase, cg_iters_out)
	real(wp), intent(inout) :: ux(Nnx,Nny,Nnz), uy(Nnx,Nny,Nnz), uz(Nnx,Nny,Nnz)
	real(wp), intent(in)    :: fx(Nnx,Nny,Nnz), fy(Nnx,Nny,Nnz), fz(Nnx,Nny,Nnz)
	integer,  intent(in)    :: phase(Nnx,Nny,Nnz)
	integer,  intent(out)   :: cg_iters_out

	real(wp), allocatable :: rx(:,:,:), ry(:,:,:), rz_arr(:,:,:)
	real(wp), allocatable :: zx(:,:,:), zy(:,:,:), zz_arr(:,:,:)
	real(wp), allocatable :: px(:,:,:), py(:,:,:), pz(:,:,:)
	real(wp), allocatable :: Apx(:,:,:), Apy(:,:,:), Apz(:,:,:)
	real(wp), allocatable :: diag_inv(:,:,:)
	real(wp) :: rz_old, rz_new, pAp, alpha_cg, beta_cg, rnorm, bnorm
	real(wp) :: dxe, dye, dze, Ke_loc(24,24)
	integer  :: iter, ie, je, ke, di, dj, dk, ln, dof

	allocate(rx(Nnx,Nny,Nnz), ry(Nnx,Nny,Nnz), rz_arr(Nnx,Nny,Nnz))
	allocate(zx(Nnx,Nny,Nnz), zy(Nnx,Nny,Nnz), zz_arr(Nnx,Nny,Nnz))
	allocate(px(Nnx,Nny,Nnz), py(Nnx,Nny,Nnz), pz(Nnx,Nny,Nnz))
	allocate(Apx(Nnx,Nny,Nnz), Apy(Nnx,Nny,Nnz), Apz(Nnx,Nny,Nnz))
	allocate(diag_inv(Nnx,Nny,Nnz))

	! Jacobi preconditioner: assemble diagonal using precomputed per-layer Ke
	diag_inv = 0.0_wp
	do ke = 1, Nz
	do je = 1, Ny
	do ie = 1, Nx
		if (elem_phase(ie,je,ke) /= MECH_POWDER) then
			Ke_loc = Ke_solid_z(:,:,ke)
		else
			Ke_loc = Ke_soft_z(:,:,ke)
		endif
		do dk = 0, 1; do dj = 0, 1; do di = 0, 1
			ln = 1 + di + 2*dj + 4*dk
			dof = 3*(ln-1) + 1
			diag_inv(ie+di,je+dj,ke+dk) = diag_inv(ie+di,je+dj,ke+dk) + Ke_loc(dof,dof)
		enddo; enddo; enddo
	enddo
	enddo
	enddo
	where (diag_inv > 1.0e-30_wp)
		diag_inv = 1.0_wp / diag_inv
	elsewhere
		diag_inv = 1.0_wp
	end where
	diag_inv(:,:,1) = 1.0_wp  ! Dirichlet

	! Enforce Dirichlet
	ux(:,:,1) = 0.0_wp; uy(:,:,1) = 0.0_wp; uz(:,:,1) = 0.0_wp

	! r = f - A*u
	call ebe_matvec_mech(ux, uy, uz, Apx, Apy, Apz, phase)
	rx = fx - Apx; ry = fy - Apy; rz_arr = fz - Apz

	! z = M^{-1} r
	zx = diag_inv * rx; zy = diag_inv * ry; zz_arr = diag_inv * rz_arr
	px = zx; py = zy; pz = zz_arr

	rz_old = sum(rx*zx) + sum(ry*zy) + sum(rz_arr*zz_arr)
	bnorm  = sqrt(sum(fx*fx) + sum(fy*fy) + sum(fz*fz))
	if (bnorm < 1.0e-30_wp) then
		cg_iters_out = 0
		deallocate(rx,ry,rz_arr,zx,zy,zz_arr,px,py,pz,Apx,Apy,Apz,diag_inv)
		return
	endif

	do iter = 1, cg_maxiter_mech
		call ebe_matvec_mech(px, py, pz, Apx, Apy, Apz, phase)
		pAp = sum(px*Apx) + sum(py*Apy) + sum(pz*Apz)
		if (abs(pAp) < 1.0e-30_wp) exit
		alpha_cg = rz_old / pAp

		ux = ux + alpha_cg * px; uy = uy + alpha_cg * py; uz = uz + alpha_cg * pz
		rx = rx - alpha_cg * Apx; ry = ry - alpha_cg * Apy; rz_arr = rz_arr - alpha_cg * Apz

		zx = diag_inv * rx; zy = diag_inv * ry; zz_arr = diag_inv * rz_arr

		rnorm = sqrt(sum(rx*rx) + sum(ry*ry) + sum(rz_arr*rz_arr))
		if (rnorm / bnorm < cg_tol_mech) exit

		rz_new = sum(rx*zx) + sum(ry*zy) + sum(rz_arr*zz_arr)
		beta_cg = rz_new / rz_old
		px = zx + beta_cg * px; py = zy + beta_cg * py; pz = zz_arr + beta_cg * pz
		rz_old = rz_new
	enddo

	cg_iters_out = iter
	ux(:,:,1) = 0.0_wp; uy(:,:,1) = 0.0_wp; uz(:,:,1) = 0.0_wp

	deallocate(rx,ry,rz_arr,zx,zy,zz_arr,px,py,pz,Apx,Apy,Apz,diag_inv)
end subroutine solve_mech_cg

!********************************************************************
subroutine update_gp_state(ux, uy, uz, dT_gp_arr, T_fem, phase)
	real(wp), intent(in) :: ux(Nnx,Nny,Nnz), uy(Nnx,Nny,Nnz), uz(Nnx,Nny,Nnz)
	real(wp), intent(in) :: dT_gp_arr(8, Nx, Ny, Nz)
	real(wp), intent(in) :: T_fem(Nnx, Nny, Nnz)
	integer,  intent(in) :: phase(Nnx,Nny,Nnz)

	real(wp) :: u_e(24), eps_curr(6), eps_inc(6), eps_th(6)
	real(wp) :: s_trial(6), sigma(6), C_loc(6,6), B_gp(6,24)
	real(wp) :: av, dT_g, f_yield_g, dxe, dye, dze, T_gp_val
	integer  :: ie, je, ke, di, dj, dk, ln, dof, g, p, q, ph_g, a

	f_plus = 0.0_wp

	!$OMP PARALLEL DO PRIVATE(ie,je,ke,u_e,eps_curr,eps_inc,eps_th,s_trial,sigma,C_loc,B_gp) &
	!$OMP   PRIVATE(av,dT_g,f_yield_g,dxe,dye,dze,di,dj,dk,ln,dof,g,p,q,ph_g) COLLAPSE(3)
	do ke = 1, Nz
	do je = 1, Ny
	do ie = 1, Nx
		call get_elem_dx(ie, je, ke, dxe, dye, dze)

		do dk = 0, 1
		do dj = 0, 1
		do di = 0, 1
			ln = 1 + di + 2*dj + 4*dk
			dof = 3*(ln-1)
			u_e(dof+1) = ux(ie+di, je+dj, ke+dk)
			u_e(dof+2) = uy(ie+di, je+dj, ke+dk)
			u_e(dof+3) = uz(ie+di, je+dj, ke+dk)
		enddo
		enddo
		enddo

		do g = 1, 8
			call get_B_at_gp(dxe, dye, dze, g, B_gp)

			do p = 1, 6
				eps_curr(p) = 0.0_wp
				do q = 1, 24
					eps_curr(p) = eps_curr(p) + B_gp(p, q) * u_e(q)
				enddo
			enddo

			eps_inc = eps_curr - eps_gp(:, g, ie, je, ke)

			di = mod(g-1, 2); dj = mod((g-1)/2, 2); dk = (g-1) / 4
			ph_g = phase(ie+di, je+dj, ke+dk)

			if (ph_g == MECH_SOLID) then
				call build_C_matrix(E_solid, nu_mech, C_loc)
				av = alpha_V
			else
				call build_C_matrix(E_soft, nu_mech, C_loc)
				av = 0.0_wp
			endif

			dT_g = dT_gp_arr(g, ie, je, ke)
			eps_th = 0.0_wp
			eps_th(1) = av * dT_g; eps_th(2) = av * dT_g; eps_th(3) = av * dT_g

			do p = 1, 6
				s_trial(p) = sig_gp(p, g, ie, je, ke)
				do q = 1, 6
					s_trial(p) = s_trial(p) + C_loc(p,q) * (eps_inc(q) - eps_th(q))
				enddo
			enddo

			! Temperature at GP
		T_gp_val = 0.0_wp
		do a = 1, 8
			di = mod(a-1,2); dj = mod((a-1)/2,2); dk = (a-1)/4
			T_gp_val = T_gp_val + N_gp(a, g) * T_fem(ie+di, je+dj, ke+dk)
		enddo
		call j2_return_map(s_trial, sigma, f_yield_g, T_gp_val)
			sig_gp(:, g, ie, je, ke) = sigma
			eps_gp(:, g, ie, je, ke) = eps_curr
		enddo
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO

	! Scatter f_plus using 8-color
	call scatter_fplus_to_nodes(T_fem)
end subroutine update_gp_state

!********************************************************************
subroutine scatter_fplus_to_nodes(T_fem)
	real(wp), intent(in) :: T_fem(Nnx, Nny, Nnz)
	integer :: ie, je, ke, g, di, dj, dk, color, ic, jc, kc
	real(wp) :: sm, sd(6), sn, fy, sy

	f_plus = 0.0_wp
	do color = 0, 7
		ic = mod(color, 2); jc = mod(color/2, 2); kc = mod(color/4, 2)
		!$OMP PARALLEL DO PRIVATE(ie,je,ke,g,di,dj,dk,sm,sd,sn,fy,sy) COLLAPSE(3)
		do ke = 1+kc, Nz, 2
		do je = 1+jc, Ny, 2
		do ie = 1+ic, Nx, 2
			do g = 1, 8
				di = mod(g-1,2); dj = mod((g-1)/2,2); dk = (g-1)/4
				sm = (sig_gp(1,g,ie,je,ke)+sig_gp(2,g,ie,je,ke)+sig_gp(3,g,ie,je,ke))/3.0_wp
				sd(1)=sig_gp(1,g,ie,je,ke)-sm; sd(2)=sig_gp(2,g,ie,je,ke)-sm; sd(3)=sig_gp(3,g,ie,je,ke)-sm
				sd(4)=sig_gp(4,g,ie,je,ke); sd(5)=sig_gp(5,g,ie,je,ke); sd(6)=sig_gp(6,g,ie,je,ke)
				sn = sqrt(1.5_wp*(sd(1)**2+sd(2)**2+sd(3)**2+2.0_wp*(sd(4)**2+sd(5)**2+sd(6)**2)))
				sy = get_sig_yield(T_fem(ie+di, je+dj, ke+dk))
				fy = max(sn - sy, 0.0_wp)
				f_plus(ie+di,je+dj,ke+dk) = max(f_plus(ie+di,je+dj,ke+dk), fy)
			enddo
		enddo
		enddo
		enddo
		!$OMP END PARALLEL DO
	enddo
end subroutine scatter_fplus_to_nodes

!********************************************************************
subroutine solve_mechanical(T_phoenix, sf_phoenix, res_out, newton_iters_out, cg_iters_out)
! Main public interface. Extracts T from PHOENIX, solves mechanical, updates stress.
	real(wp), intent(in)  :: T_phoenix(:,:,:)    ! PHOENIX temp(ni,nj,nk)
	real(wp), intent(in)  :: sf_phoenix(:,:,:)   ! PHOENIX solidfield(ni,nj,nk)
	real(wp), intent(out) :: res_out
	integer,  intent(out) :: newton_iters_out, cg_iters_out

	real(wp), allocatable :: T_fem(:,:,:), sf_fem(:,:,:)
	real(wp), allocatable :: Rx(:,:,:), Ry(:,:,:), Rz(:,:,:)
	real(wp), allocatable :: dux(:,:,:), duy(:,:,:), duz(:,:,:)
	real(wp), allocatable :: dT_gp_arr(:,:,:,:)
	real(wp) :: R_norm, R_norm0
	integer  :: newton_iter, cg_iter_this

	allocate(T_fem(Nnx,Nny,Nnz), sf_fem(Nnx,Nny,Nnz))
	allocate(Rx(Nnx,Nny,Nnz), Ry(Nnx,Nny,Nnz), Rz(Nnx,Nny,Nnz))
	allocate(dux(Nnx,Nny,Nnz), duy(Nnx,Nny,Nnz), duz(Nnx,Nny,Nnz))
	allocate(dT_gp_arr(8, Nx, Ny, Nz))

	! Extract fields from PHOENIX grid to FEM grid
	call extract_temp_to_fem(T_phoenix, T_fem)
	call extract_solidfield_to_fem(sf_phoenix, sf_fem)
	call update_mech_phase(mech_phase, T_fem, sf_fem, Nnx, Nny, Nnz, fem_z)
	call compute_elem_phases(mech_phase)
	call compute_dT_gp(T_fem, dT_gp_arr)

	! Newton iteration
	cg_iters_out = 0
	R_norm0 = 0.0_wp
	do newton_iter = 1, newton_maxiter
		call compute_residual(ux_mech, uy_mech, uz_mech, dT_gp_arr, T_fem, mech_phase, Rx, Ry, Rz)

		R_norm = sqrt(sum(Rx**2) + sum(Ry**2) + sum(Rz**2))
		if (newton_iter == 1) R_norm0 = R_norm
		if (R_norm / max(R_norm0, 1.0e-30_wp) < newton_tol) exit

		Rx = -Rx; Ry = -Ry; Rz = -Rz
		dux = 0.0_wp; duy = 0.0_wp; duz = 0.0_wp
		call solve_mech_cg(dux, duy, duz, Rx, Ry, Rz, mech_phase, cg_iter_this)
		cg_iters_out = cg_iters_out + cg_iter_this
		ux_mech = ux_mech + dux; uy_mech = uy_mech + duy; uz_mech = uz_mech + duz
	enddo

	newton_iters_out = newton_iter
	res_out = R_norm / max(R_norm0, 1.0e-30_wp)

	! Update GP state
	call update_gp_state(ux_mech, uy_mech, uz_mech, dT_gp_arr, T_fem, mech_phase)

	! Save FEM temperature for VTK/history output
	T_fem_last = T_fem

	! Update temperature reference
	T_old_mech = T_fem

	deallocate(T_fem, sf_fem, Rx, Ry, Rz, dux, duy, duz, dT_gp_arr)
end subroutine solve_mechanical

!********************************************************************
subroutine get_stress_yield(sxx, syy, szz, vonmises, fplus_arr)
! Smooth GP stresses to nodes and compute von Mises.
	real(wp), intent(out) :: sxx(Nnx,Nny,Nnz), syy(Nnx,Nny,Nnz), szz(Nnx,Nny,Nnz)
	real(wp), intent(out) :: vonmises(Nnx,Nny,Nnz), fplus_arr(Nnx,Nny,Nnz)
	integer :: ie, je, ke, g, di, dj, dk
	integer :: cnt(Nnx,Nny,Nnz)
	real(wp) :: sxx_e, syy_e, szz_e

	sxx = 0.0_wp; syy = 0.0_wp; szz = 0.0_wp
	cnt = 0

	do ke = 1, Nz
	do je = 1, Ny
	do ie = 1, Nx
		sxx_e = 0.0_wp; syy_e = 0.0_wp; szz_e = 0.0_wp
		do g = 1, 8
			sxx_e = sxx_e + sig_gp(1, g, ie, je, ke)
			syy_e = syy_e + sig_gp(2, g, ie, je, ke)
			szz_e = szz_e + sig_gp(3, g, ie, je, ke)
		enddo
		sxx_e = sxx_e / 8.0_wp; syy_e = syy_e / 8.0_wp; szz_e = szz_e / 8.0_wp

		do dk = 0, 1
		do dj = 0, 1
		do di = 0, 1
			sxx(ie+di, je+dj, ke+dk) = sxx(ie+di, je+dj, ke+dk) + sxx_e
			syy(ie+di, je+dj, ke+dk) = syy(ie+di, je+dj, ke+dk) + syy_e
			szz(ie+di, je+dj, ke+dk) = szz(ie+di, je+dj, ke+dk) + szz_e
			cnt(ie+di, je+dj, ke+dk) = cnt(ie+di, je+dj, ke+dk) + 1
		enddo
		enddo
		enddo
	enddo
	enddo
	enddo

	where (cnt > 0)
		sxx = sxx / real(cnt, wp)
		syy = syy / real(cnt, wp)
		szz = szz / real(cnt, wp)
	end where

	! Von Mises stress
	vonmises = sqrt(0.5_wp * ((sxx-syy)**2 + (syy-szz)**2 + (szz-sxx)**2))

	fplus_arr = f_plus
end subroutine get_stress_yield

end module mechanical_solver
