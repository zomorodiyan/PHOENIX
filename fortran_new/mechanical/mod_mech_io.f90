!______________________________________________________________________________
!
module mech_io
!______________________________________________________________________________
! VTK output, timing report, and history for mechanical solver.
!
	use precision
	use geometry, only: x, y, z, nim1, njm1, nkm1
	use parameters, only: wp, file_prefix, case_name, result_dir
	use mech_material, only: sig_yield
	use mechanical_solver, only: fem_x, fem_y, fem_z
	implicit none

	! Timing accumulators
	real(wp), save :: mio_t_cpu  = 0.0_wp
	real(wp), save :: mio_t_wall   = 0.0_wp
	integer,  save :: mio_n_solves = 0
	integer,  save :: n_mech_vtk    = 0

	! History file initialized flag
	logical, save :: mech_hist_init = .false.

	contains

!********************************************************************
subroutine write_mech_vtk(step_idx, T_fem, ux, uy, uz, phase, sxx, syy, szz, vonmises, fplus, &
                           Nnx, Nny, Nnz)
! Write mechanical VTK file with all fields.
	integer, intent(in)   :: step_idx, Nnx, Nny, Nnz
	real(wp), intent(in)  :: T_fem(Nnx,Nny,Nnz)
	real(wp), intent(in)  :: ux(Nnx,Nny,Nnz), uy(Nnx,Nny,Nnz), uz(Nnx,Nny,Nnz)
	integer,  intent(in)  :: phase(Nnx,Nny,Nnz)
	real(wp), intent(in)  :: sxx(Nnx,Nny,Nnz), syy(Nnx,Nny,Nnz), szz(Nnx,Nny,Nnz)
	real(wp), intent(in)  :: vonmises(Nnx,Nny,Nnz), fplus(Nnx,Nny,Nnz)

	integer, parameter :: lun = 91
	integer :: i, j, k, npts
	real(kind=4) :: val4
	character(len=10) :: cstep
	character(len=512) :: vtkfile

	write(cstep, '(I5.5)') step_idx
	vtkfile = trim(file_prefix)//'mech_'//trim(adjustl(cstep))//'.vtk'
	npts = Nnx * Nny * Nnz

	! ASCII header
	open(unit=lun, file=trim(vtkfile))
	write(lun,'(A)') '# vtk DataFile Version 3.0'
	write(lun,'(A)') 'PHOENIX Mechanical Output'
	write(lun,'(A)') 'BINARY'
	write(lun,'(A)') 'DATASET STRUCTURED_GRID'
	write(lun,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', Nnx, ' ', Nny, ' ', Nnz
	write(lun,'(A,I0,A)') 'POINTS ', npts, ' float'
	close(lun)

	! Binary coordinates (from FEM node arrays)
	open(unit=lun, file=trim(vtkfile), access='stream', form='unformatted', &
	     position='append', convert='big_endian')
	do k = 1, Nnz
	do j = 1, Nny
	do i = 1, Nnx
		val4 = real(fem_x(i), 4); write(lun) val4
		val4 = real(fem_y(j), 4); write(lun) val4
		val4 = real(fem_z(k), 4); write(lun) val4
	enddo
	enddo
	enddo
	close(lun)

	! POINT_DATA header
	open(unit=lun, file=trim(vtkfile), position='append')
	write(lun,'(A,I0)') 'POINT_DATA ', npts
	close(lun)

	! Scalar fields
	call write_mech_scalar(lun, vtkfile, 'Temperature', T_fem, Nnx, Nny, Nnz)
	call write_mech_scalar(lun, vtkfile, 'ux', ux, Nnx, Nny, Nnz)
	call write_mech_scalar(lun, vtkfile, 'uy', uy, Nnx, Nny, Nnz)
	call write_mech_scalar(lun, vtkfile, 'uz', uz, Nnx, Nny, Nnz)
	call write_mech_scalar_int(lun, vtkfile, 'phase', phase, Nnx, Nny, Nnz)
	call write_mech_scalar(lun, vtkfile, 'sxx', sxx, Nnx, Nny, Nnz)
	call write_mech_scalar(lun, vtkfile, 'syy', syy, Nnx, Nny, Nnz)
	call write_mech_scalar(lun, vtkfile, 'szz', szz, Nnx, Nny, Nnz)
	call write_mech_scalar(lun, vtkfile, 'von_mises', vonmises, Nnx, Nny, Nnz)
	call write_mech_scalar(lun, vtkfile, 'fplus', fplus, Nnx, Nny, Nnz)

	n_mech_vtk = n_mech_vtk + 1
end subroutine write_mech_vtk

!********************************************************************
subroutine write_mech_scalar(lun, vtkfile, name, field, Nnx, Nny, Nnz)
	integer, intent(in) :: lun, Nnx, Nny, Nnz
	character(len=*), intent(in) :: vtkfile, name
	real(wp), intent(in) :: field(Nnx, Nny, Nnz)
	integer :: i, j, k
	real(kind=4) :: val4

	open(unit=lun, file=trim(vtkfile), position='append')
	write(lun,'(A)') 'SCALARS '//trim(name)//' float 1'
	write(lun,'(A)') 'LOOKUP_TABLE default'
	close(lun)

	open(unit=lun, file=trim(vtkfile), access='stream', form='unformatted', &
	     position='append', convert='big_endian')
	do k = 1, Nnz
	do j = 1, Nny
	do i = 1, Nnx
		val4 = real(field(i,j,k), 4); write(lun) val4
	enddo
	enddo
	enddo
	close(lun)
end subroutine write_mech_scalar

!********************************************************************
subroutine write_mech_scalar_int(lun, vtkfile, name, field, Nnx, Nny, Nnz)
	integer, intent(in) :: lun, Nnx, Nny, Nnz
	character(len=*), intent(in) :: vtkfile, name
	integer, intent(in) :: field(Nnx, Nny, Nnz)
	integer :: i, j, k
	integer(kind=4) :: val4

	open(unit=lun, file=trim(vtkfile), position='append')
	write(lun,'(A)') 'SCALARS '//trim(name)//' int 1'
	write(lun,'(A)') 'LOOKUP_TABLE default'
	close(lun)

	open(unit=lun, file=trim(vtkfile), access='stream', form='unformatted', &
	     position='append', convert='big_endian')
	do k = 1, Nnz
	do j = 1, Nny
	do i = 1, Nnx
		val4 = int(field(i,j,k), 4); write(lun) val4
	enddo
	enddo
	enddo
	close(lun)
end subroutine write_mech_scalar_int

!********************************************************************
subroutine write_mech_timing_report(file_prefix_in)
	character(len=*), intent(in) :: file_prefix_in
	integer, parameter :: lun = 92

	open(unit=lun, file=trim(file_prefix_in)//'mech_timing_report.txt', action='write', status='replace')
	write(lun,'(a)') '============================================'
	write(lun,'(a)') '  PHOENIX Mechanical Solver Timing Report'
	write(lun,'(a)') '============================================'
	write(lun,'(a,i0)')     '  Solves performed:    ', mio_n_solves
	write(lun,'(a,f12.3,a)') '  Total CPU time:      ', mio_t_cpu, ' s'
	write(lun,'(a,f12.3,a)') '  Total wall time:     ', mio_t_wall, ' s'
	if (mio_n_solves > 0) then
		write(lun,'(a,f12.3,a)') '  Avg wall per solve:  ', mio_t_wall / real(mio_n_solves, wp), ' s'
	endif
	write(lun,'(a,i0)')     '  VTK files written:   ', n_mech_vtk
	write(lun,'(a)') '============================================'
	close(lun)
end subroutine write_mech_timing_report

!********************************************************************
subroutine write_mech_memory_report(file_prefix_in, Nx, Ny, Nz, Nnx, Nny, Nnz)
	character(len=*), intent(in) :: file_prefix_in
	integer, intent(in) :: Nx, Ny, Nz, Nnx, Nny, Nnz
	integer, parameter :: lun = 92
	real(wp) :: gp_mem, disp_mem, stress_mem, total_mem

	gp_mem = 2.0_wp * 6.0_wp * 8.0_wp * real(Nx,wp) * real(Ny,wp) * real(Nz,wp) * 8.0_wp / (1024.0_wp**2)
	disp_mem = 3.0_wp * real(Nnx,wp) * real(Nny,wp) * real(Nnz,wp) * 8.0_wp / (1024.0_wp**2)
	stress_mem = 5.0_wp * real(Nnx,wp) * real(Nny,wp) * real(Nnz,wp) * 8.0_wp / (1024.0_wp**2)
	total_mem = gp_mem + disp_mem + stress_mem

	open(unit=lun, file=trim(file_prefix_in)//'mech_memory_report.txt', action='write', status='replace')
	write(lun,'(a)') '============================================'
	write(lun,'(a)') '  PHOENIX Mechanical Memory Report'
	write(lun,'(a)') '============================================'
	write(lun,'(a,i0,a,i0,a,i0)') '  FEM elements: ', Nx, ' x ', Ny, ' x ', Nz
	write(lun,'(a,i0,a,i0,a,i0)') '  FEM nodes:    ', Nnx, ' x ', Nny, ' x ', Nnz
	write(lun,'(a)')
	write(lun,'(a,f10.1,a)') '  GP state (sig_gp + eps_gp):  ', gp_mem, ' MB'
	write(lun,'(a,f10.1,a)') '  Displacement (ux,uy,uz):     ', disp_mem, ' MB'
	write(lun,'(a,f10.1,a)') '  Stress output (sxx..fplus):  ', stress_mem, ' MB'
	write(lun,'(a,f10.1,a)') '  Total mechanical:            ', total_mem, ' MB'
	write(lun,'(a)') '============================================'
	close(lun)
end subroutine write_mech_memory_report

!********************************************************************
subroutine init_mech_history(Nnx, Nny, Nnz)
	use printing, only: thist_px, thist_py, n_thist
	integer, intent(in) :: Nnx, Nny, Nnz
	integer :: p
	integer, parameter :: lun = 93

	open(unit=lun, file=trim(file_prefix)//'mech_history.txt', status='replace')
	write(lun,'(a)') '# Mechanical History - 10 Monitoring Points'
	write(lun,'(a)') '# Columns: time(s)  T1..T10(K)  ux1..ux10(m)  uy1..uy10(m)  uz1..uz10(m)  sxx1..sxx10(Pa)  syy1..syy10(Pa)'
	write(lun,'(a)') '# Point physical coordinates:'
	do p = 1, n_thist
		write(lun,'(a,i2,a,2(es10.3,a))') '#  P', p, ':  x=', thist_px(p), 'm  y=', thist_py(p), 'm'
	enddo
	close(lun)
	mech_hist_init = .true.
end subroutine init_mech_history

!********************************************************************
subroutine write_mech_history(t, T_fem, ux, uy, uz, sxx, syy, Nnx, Nny, Nnz)
	use printing, only: thist_px, thist_py, n_thist
	real(wp), intent(in) :: t
	real(wp), intent(in) :: T_fem(Nnx,Nny,Nnz), ux(Nnx,Nny,Nnz), uy(Nnx,Nny,Nnz), uz(Nnx,Nny,Nnz)
	real(wp), intent(in) :: sxx(Nnx,Nny,Nnz), syy(Nnx,Nny,Nnz)
	integer, intent(in)  :: Nnx, Nny, Nnz
	integer, parameter :: lun = 93
	integer :: p, i1, i2, j1, j2
	real(wp) :: wx, wy, xp, yp

	if (.not. mech_hist_init) return

	open(unit=lun, file=trim(file_prefix)//'mech_history.txt', status='old', position='append')
	write(lun,'(es12.5)', advance='no') t

	do p = 1, n_thist
		xp = thist_px(p); yp = thist_py(p)
		! Find bracketing FEM node indices (FEM node i corresponds to x(i+1))
		call find_fem_bracket(xp, yp, Nnx, Nny, i1, i2, j1, j2, wx, wy)
		! Interpolate and write: T, ux, uy, uz, sxx, syy
		write(lun,'(6(2x,es12.5))', advance='no') &
			interp2d(T_fem, i1,i2,j1,j2,wx,wy,Nnx,Nny,Nnz), &
			interp2d(ux,    i1,i2,j1,j2,wx,wy,Nnx,Nny,Nnz), &
			interp2d(uy,    i1,i2,j1,j2,wx,wy,Nnx,Nny,Nnz), &
			interp2d(uz,    i1,i2,j1,j2,wx,wy,Nnx,Nny,Nnz), &
			interp2d(sxx,   i1,i2,j1,j2,wx,wy,Nnx,Nny,Nnz), &
			interp2d(syy,   i1,i2,j1,j2,wx,wy,Nnx,Nny,Nnz)
	enddo
	write(lun,*)
	close(lun)
end subroutine write_mech_history

!********************************************************************
function interp2d(field, i1, i2, j1, j2, wx, wy, Nnx, Nny, Nnz) result(val)
! Bilinear interpolation at surface (k=Nnz) of FEM field.
	real(wp), intent(in) :: field(Nnx, Nny, Nnz), wx, wy
	integer, intent(in)  :: i1, i2, j1, j2, Nnx, Nny, Nnz
	real(wp) :: val
	integer :: kk
	kk = Nnz  ! surface
	val = (1.0_wp-wx)*(1.0_wp-wy) * field(i1,j1,kk) + &
	      wx*(1.0_wp-wy)           * field(i2,j1,kk) + &
	      (1.0_wp-wx)*wy           * field(i1,j2,kk) + &
	      wx*wy                    * field(i2,j2,kk)
end function interp2d

!********************************************************************
subroutine find_fem_bracket(xp, yp, Nnx, Nny, i1, i2, j1, j2, wx, wy)
! Find bracketing FEM node indices for physical coordinate (xp, yp).
	real(wp), intent(in)  :: xp, yp
	integer,  intent(in)  :: Nnx, Nny
	integer,  intent(out) :: i1, i2, j1, j2
	real(wp), intent(out) :: wx, wy
	integer :: i

	i1 = 1
	do i = 1, Nnx - 1
		if (fem_x(i+1) > xp) exit
		i1 = i + 1
	enddo
	if (i1 >= Nnx) i1 = Nnx - 1
	i2 = i1 + 1
	if (abs(fem_x(i2) - fem_x(i1)) > 1.0e-15_wp) then
		wx = (xp - fem_x(i1)) / (fem_x(i2) - fem_x(i1))
	else
		wx = 0.0_wp
	endif
	wx = max(0.0_wp, min(1.0_wp, wx))

	j1 = 1
	do i = 1, Nny - 1
		if (fem_y(i+1) > yp) exit
		j1 = i + 1
	enddo
	if (j1 >= Nny) j1 = Nny - 1
	j2 = j1 + 1
	if (abs(fem_y(j2) - fem_y(j1)) > 1.0e-15_wp) then
		wy = (yp - fem_y(j1)) / (fem_y(j2) - fem_y(j1))
	else
		wy = 0.0_wp
	endif
	wy = max(0.0_wp, min(1.0_wp, wy))
end subroutine find_fem_bracket

!********************************************************************
subroutine finalize_mech_history()
! Write Python plot script for mechanical history and execute it.
	integer, parameter :: lun = 93

	if (.not. mech_hist_init) return

	open(unit=lun, file=trim(file_prefix)//'plot_mech_history.py', status='replace')
	write(lun,'(a)') 'import numpy as np'
	write(lun,'(a)') 'import matplotlib'
	write(lun,'(a)') 'matplotlib.use("Agg")'
	write(lun,'(a)') 'import matplotlib.pyplot as plt'
	write(lun,'(a)') ''
	write(lun,'(a)') 'data = np.loadtxt("'//trim(adjustl(case_name))//'_mech_history.txt", comments="#")'
	write(lun,'(a)') 'if data.ndim == 1: data = data[np.newaxis, :]'
	write(lun,'(a)') 't = data[:, 0] * 1e3  # ms'
	write(lun,'(a)') 'n = 10  # monitoring points'
	write(lun,'(a)') 'T   = data[:, 1:1+n]'
	write(lun,'(a)') 'ux  = data[:, 1+n:1+2*n]'
	write(lun,'(a)') 'uy  = data[:, 1+2*n:1+3*n]'
	write(lun,'(a)') 'uz  = data[:, 1+3*n:1+4*n]'
	write(lun,'(a)') 'sxx = data[:, 1+4*n:1+5*n]'
	write(lun,'(a)') 'syy = data[:, 1+5*n:1+6*n]'
	write(lun,'(a)') 'disp_mag = np.sqrt(ux**2 + uy**2 + uz**2) * 1e6  # um'
	write(lun,'(a)') ''
	write(lun,'(a)') 'fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)'
	write(lun,'(a)') 'for i in range(n):'
	write(lun,'(a)') '    axes[0].plot(t, T[:, i], label=f"P{i+1}")'
	write(lun,'(a)') 'axes[0].set_ylabel("Temperature (K)")'
	write(lun,'(a)') 'axes[0].legend(fontsize=7, ncol=5)'
	write(lun,'(a)') 'axes[0].grid(True, alpha=0.3)'
	write(lun,'(a)') ''
	write(lun,'(a)') 'for i in range(n):'
	write(lun,'(a)') '    axes[1].plot(t, disp_mag[:, i], label=f"P{i+1}")'
	write(lun,'(a)') 'axes[1].set_ylabel("Displacement (um)")'
	write(lun,'(a)') 'axes[1].legend(fontsize=7, ncol=5)'
	write(lun,'(a)') 'axes[1].grid(True, alpha=0.3)'
	write(lun,'(a)') ''
	write(lun,'(a)') 'for i in range(n):'
	write(lun,'(a)') '    axes[2].plot(t, sxx[:, i] * 1e-6, label=f"P{i+1}")'
	write(lun,'(a)') 'axes[2].axhline(250, color="red", linestyle="--", linewidth=0.8, label="sig_yield")'
	write(lun,'(a)') 'axes[2].axhline(-250, color="red", linestyle="--", linewidth=0.8)'
	write(lun,'(a)') 'axes[2].set_ylabel("Stress sxx (MPa)")'
	write(lun,'(a)') 'axes[2].set_xlabel("Time (ms)")'
	write(lun,'(a)') 'axes[2].legend(fontsize=7, ncol=5)'
	write(lun,'(a)') 'axes[2].grid(True, alpha=0.3)'
	write(lun,'(a)') ''
	write(lun,'(a)') 'plt.tight_layout()'
	write(lun,'(a)') 'plt.savefig("'//trim(adjustl(case_name))//'_mech_history.png", dpi=150)'
	write(lun,'(a)') 'print("Saved '//trim(adjustl(case_name))//'_mech_history.png")'
	close(lun)

	call execute_command_line('cd '//trim(result_dir)// &
		' && python3 '//trim(adjustl(case_name))//'_plot_mech_history.py', wait=.true.)
end subroutine finalize_mech_history

end module mech_io
