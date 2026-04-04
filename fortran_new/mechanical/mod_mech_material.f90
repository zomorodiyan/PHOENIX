!______________________________________________________________________________
!
module mech_material
!______________________________________________________________________________
! Phase-dependent mechanical material properties for EBE FEM solver.
! Parameters identical to fortran-ebe-lpbf.
!
	use precision
	use parameters, only: tsolid, wp
	use constant, only: powder_threshold
	implicit none

	! Phase constants
	integer, parameter :: MECH_POWDER = 0
	integer, parameter :: MECH_LIQUID = 1
	integer, parameter :: MECH_SOLID  = 2

	! Mechanical material parameters (identical to fortran-ebe-lpbf)
	real(wp), parameter :: E_solid   = 70.0e9_wp      ! Young's modulus, solid (Pa)
	real(wp), parameter :: E_soft    = 0.7e9_wp        ! Young's modulus, powder/liquid (Pa)
	real(wp), parameter :: nu_mech   = 0.3_wp          ! Poisson's ratio
	real(wp), parameter :: sig_yield = 250.0e6_wp      ! J2 yield stress (Pa)
	real(wp), parameter :: alpha_V   = 1.0e-5_wp       ! Volumetric thermal expansion (1/K)

	! Solver parameters
	integer, parameter  :: cg_maxiter_mech = 20000
	real(wp), parameter :: cg_tol_mech     = 1.0e-4_wp
	integer, parameter  :: newton_maxiter  = 10
	real(wp), parameter :: newton_tol      = 1.0e-4_wp

	contains

!********************************************************************
subroutine update_mech_phase(mech_phase, temp_arr, solidfield_arr, nnx, nny, nnz)
! Update mechanical phase from PHOENIX temperature and solidfield.
! MECH_POWDER: never melted and below solidus
! MECH_LIQUID: currently above solidus (in melt pool)
! MECH_SOLID:  was once melted and now below solidus
	integer,  intent(out) :: mech_phase(nnx, nny, nnz)
	real(wp), intent(in)  :: temp_arr(nnx, nny, nnz)
	real(wp), intent(in)  :: solidfield_arr(nnx, nny, nnz)
	integer, intent(in)   :: nnx, nny, nnz
	integer :: i, j, k

	!$OMP PARALLEL DO PRIVATE(i,j,k)
	do k = 1, nnz
	do j = 1, nny
	do i = 1, nnx
		if (temp_arr(i,j,k) >= tsolid) then
			mech_phase(i,j,k) = MECH_LIQUID
		else if (solidfield_arr(i,j,k) > powder_threshold) then
			mech_phase(i,j,k) = MECH_SOLID
		else
			mech_phase(i,j,k) = MECH_POWDER
		endif
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO
end subroutine update_mech_phase

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
subroutine j2_return_map(s_trial, s_mapped, f_yield_out)
! J2 von Mises radial return mapping.
	real(wp), intent(in)  :: s_trial(6)
	real(wp), intent(out) :: s_mapped(6)
	real(wp), intent(out) :: f_yield_out
	real(wp) :: s_mean, s_dev(6), s_norm, f_yield

	s_mean = (s_trial(1) + s_trial(2) + s_trial(3)) / 3.0_wp
	s_dev(1) = s_trial(1) - s_mean
	s_dev(2) = s_trial(2) - s_mean
	s_dev(3) = s_trial(3) - s_mean
	s_dev(4) = s_trial(4)
	s_dev(5) = s_trial(5)
	s_dev(6) = s_trial(6)

	s_norm = sqrt(1.5_wp * (s_dev(1)**2 + s_dev(2)**2 + s_dev(3)**2 &
	         + 2.0_wp * (s_dev(4)**2 + s_dev(5)**2 + s_dev(6)**2)))

	f_yield = s_norm - sig_yield
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

end module mech_material
