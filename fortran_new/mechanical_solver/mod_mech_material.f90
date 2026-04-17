!______________________________________________________________________________
!
module mech_material
!______________________________________________________________________________
! Phase-dependent mechanical material properties for EBE FEM solver.
! Parameters identical to fortran-ebe-lpbf.
!
	use precision
	use parameters, only: tsolid, wp, layerheight
	use geometry, only: z, nk
	use constant, only: powder_threshold
	implicit none

	! Phase constants
	integer, parameter :: MECH_POWDER = 0
	integer, parameter :: MECH_LIQUID = 1
	integer, parameter :: MECH_SOLID  = 2

	! Mechanical material parameters (read from input_param_mechanical.txt)
	real(wp) :: E_solid
	real(wp) :: E_soft
	real(wp) :: nu_mech
	real(wp) :: sig_yield_0
	real(wp) :: T_ref_yield
	real(wp) :: alpha_V

	! Solver parameters (read from input_param_mechanical.txt)
	integer  :: cg_maxiter_mech
	real(wp) :: cg_tol_mech
	integer  :: newton_maxiter
	real(wp) :: newton_tol

	namelist / mechanical_material / E_solid, E_soft, nu_mech, sig_yield_0, T_ref_yield, alpha_V
	namelist / mechanical_numerics / cg_maxiter_mech, cg_tol_mech, newton_maxiter, newton_tol

	contains

!********************************************************************
subroutine update_mech_phase(mech_phase, temp_arr, solidfield_arr, nnx, nny, nnz, fem_z_arr)
! Update mechanical phase from PHOENIX temperature and solidfield.
! MECH_LIQUID: currently above solidus (in melt pool)
! MECH_SOLID:  substrate OR solidified (was once melted)
! MECH_POWDER: in powder layer (z > z_top - layerheight) and never melted
	integer,  intent(out) :: mech_phase(nnx, nny, nnz)
	real(wp), intent(in)  :: temp_arr(nnx, nny, nnz)
	real(wp), intent(in)  :: solidfield_arr(nnx, nny, nnz)
	real(wp), intent(in)  :: fem_z_arr(nnz)
	integer, intent(in)   :: nnx, nny, nnz
	integer :: i, j, k
	real(wp) :: z_top, z_powder_min

	z_top = z(nk)              ! top of domain
	z_powder_min = z_top - layerheight  ! powder layer starts here

	! DEBUG: check solidfield at domain center surface
	i = nnx/2; j = nny/2; k = nnz
	write(9,'(A,3I5,A,F8.1,A,F6.2,A,F10.6)') &
		'  PHASE DBG: ijk=', i,j,k, &
		' T=', temp_arr(i,j,k), ' sf=', solidfield_arr(i,j,k), &
		' z=', fem_z_arr(k)*1e3

	!$OMP PARALLEL DO PRIVATE(i,j,k)
	do k = 1, nnz
	do j = 1, nny
	do i = 1, nnx
		if (temp_arr(i,j,k) >= tsolid) then
			! Currently molten
			mech_phase(i,j,k) = MECH_LIQUID
		else if (solidfield_arr(i,j,k) > powder_threshold) then
			! Was melted and now solidified
			mech_phase(i,j,k) = MECH_SOLID
		else if (fem_z_arr(k) >= z_powder_min) then
			! In powder layer height range and never melted → powder
			mech_phase(i,j,k) = MECH_POWDER
		else
			! Below powder layer, never melted → substrate (treated as solid)
			mech_phase(i,j,k) = MECH_SOLID
		endif
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO
end subroutine update_mech_phase

!********************************************************************
function get_sig_yield(T) result(sy)
! Temperature-dependent yield strength.
! Linear decrease from sig_yield_0 at T_ref to 0 at tsolid.
	real(wp), intent(in) :: T
	real(wp) :: sy

	if (T <= T_ref_yield) then
		sy = sig_yield_0
	else if (T >= tsolid) then
		sy = 0.0_wp
	else
		sy = sig_yield_0 * (tsolid - T) / (tsolid - T_ref_yield)
	endif
end function get_sig_yield

!********************************************************************
subroutine build_C_matrix(E_val, nu_val, C)
! Build 6x6 isotropic elasticity matrix (Voigt notation).
	real(wp), intent(in)  :: E_val, nu_val
	real(wp), intent(out) :: C(6,6)
	real(wp) :: lam, mu

	lam = E_val * nu_val / ((1.0_wp + nu_val) * (1.0_wp - 2.0_wp * nu_val))
	mu  = E_val / (2.0_wp * (1.0_wp + nu_val))

	C = 0.0_wp
	C(1,1) = lam + 2.0_wp*mu; C(1,2) = lam;             C(1,3) = lam
	C(2,1) = lam;             C(2,2) = lam + 2.0_wp*mu; C(2,3) = lam
	C(3,1) = lam;             C(3,2) = lam;             C(3,3) = lam + 2.0_wp*mu
	C(4,4) = mu
	C(5,5) = mu
	C(6,6) = mu
end subroutine build_C_matrix

!********************************************************************
subroutine j2_return_map(s_trial, s_mapped, f_yield_out, T_gp)
! J2 von Mises radial return mapping with temperature-dependent yield.
	real(wp), intent(in)  :: s_trial(6)
	real(wp), intent(out) :: s_mapped(6)
	real(wp), intent(out) :: f_yield_out
	real(wp), intent(in)  :: T_gp    ! temperature at this Gauss point
	real(wp) :: s_mean, s_dev(6), s_norm, f_yield, sy

	sy = get_sig_yield(T_gp)

	s_mean = (s_trial(1) + s_trial(2) + s_trial(3)) / 3.0_wp
	s_dev(1) = s_trial(1) - s_mean
	s_dev(2) = s_trial(2) - s_mean
	s_dev(3) = s_trial(3) - s_mean
	s_dev(4) = s_trial(4)
	s_dev(5) = s_trial(5)
	s_dev(6) = s_trial(6)

	s_norm = sqrt(1.5_wp * (s_dev(1)**2 + s_dev(2)**2 + s_dev(3)**2 &
	         + 2.0_wp * (s_dev(4)**2 + s_dev(5)**2 + s_dev(6)**2)))

	f_yield = s_norm - sy
	f_yield_out = max(f_yield, 0.0_wp)

	if (f_yield > 0.0_wp .and. s_norm > 1.0e-30_wp) then
		s_mapped(1) = s_mean + s_dev(1) * (1.0_wp - f_yield / s_norm)
		s_mapped(2) = s_mean + s_dev(2) * (1.0_wp - f_yield / s_norm)
		s_mapped(3) = s_mean + s_dev(3) * (1.0_wp - f_yield / s_norm)
		s_mapped(4) = s_dev(4) * (1.0_wp - f_yield / s_norm)
		s_mapped(5) = s_dev(5) * (1.0_wp - f_yield / s_norm)
		s_mapped(6) = s_dev(6) * (1.0_wp - f_yield / s_norm)
	else
		s_mapped = s_trial
	endif
end subroutine j2_return_map

!********************************************************************
subroutine read_mech_params()
	integer :: iu, ios
	logical :: fexist

	inquire(file='./mechanical_solver/inputfile/input_param_mechanical.txt', exist=fexist)
	if (.not. fexist) then
		write(*,'(A)') 'ERROR: mechanical_solver/inputfile/input_param_mechanical.txt not found'
		stop 1
	endif

	iu = 88
	open(unit=iu, file='./mechanical_solver/inputfile/input_param_mechanical.txt', &
	     form='formatted', status='old', iostat=ios)
	if (ios /= 0) then
		write(*,'(A)') 'ERROR: cannot open input_param_mechanical.txt'
		stop 1
	endif

	read(iu, NML=mechanical_material, iostat=ios)
	if (ios /= 0) then
		write(*,'(A)') 'ERROR: failed reading &mechanical_material from input_param_mechanical.txt'
		stop 1
	endif

	read(iu, NML=mechanical_numerics, iostat=ios)
	if (ios /= 0) then
		write(*,'(A)') 'ERROR: failed reading &mechanical_numerics from input_param_mechanical.txt'
		stop 1
	endif

	close(iu)

	write(*,'(A)')        '  [mech] Parameters from input_param_mechanical.txt:'
	write(*,'(A,ES10.3)') '    E_solid      = ', E_solid
	write(*,'(A,ES10.3)') '    E_soft       = ', E_soft
	write(*,'(A,F6.3)')   '    nu_mech      = ', nu_mech
	write(*,'(A,ES10.3)') '    sig_yield_0  = ', sig_yield_0
	write(*,'(A,F8.1)')   '    T_ref_yield  = ', T_ref_yield
	write(*,'(A,ES10.3)') '    alpha_V      = ', alpha_V
	write(*,'(A,I8)')     '    cg_maxiter   = ', cg_maxiter_mech
	write(*,'(A,ES10.3)') '    cg_tol       = ', cg_tol_mech
	write(*,'(A,I8)')     '    newton_max   = ', newton_maxiter
	write(*,'(A,ES10.3)') '    newton_tol   = ', newton_tol
end subroutine read_mech_params

end module mech_material
