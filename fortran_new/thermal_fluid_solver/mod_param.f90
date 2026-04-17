!______________________________________________________________________________
!
module parameters
!______________________________________________________________________________
!
	use constant

	implicit none

	!--------------------------------------------------------------
	! Derived type: material_t (encapsulates material properties)
	!--------------------------------------------------------------
	type :: material_t
		real(wp) :: dens, denl, viscos, tsolid, tliquid, tboiling
		real(wp) :: hsmelt, hlfriz, acpa, acpb, acpl
		real(wp) :: thconsa, thconsb, thconl, beta, emiss, dgdtp
	end type material_t

	!--------------------------------------------------------------
	! Module-level variables (backward compatibility)
	!--------------------------------------------------------------
	real(wp) alaspow, alaseta, alasrb, alasfact, scanvelx, scanvely
	real(wp) alaspowvol, alasetavol, sourcerad, sourcedepth
	real(wp) dens, denl, viscos, tsolid, tliquid, tboiling, hsmelt, hlfriz, hlatnt, acpa, acpb, acpl, &
		thconsa, thconsb, thconl, beta, emiss, dgdtp
	real(wp) layerheight, pden, pcpa, pcpb, pthcona, pthconb
	real(wp) delt, timax, urfu, urfv, urfw, urfp, urfh
	real(wp) xzone(nx1),yzone(ny1),zzone(nz1),powrx(nx1),powry(ny1),powrz(nz1)
	real(wp) htci, htcj, htck1, htckn, tempWest, tempEast, tempNorth, tempBottom, tempPreheat, tempAmb	
	integer nzx, nzy, nzz, maxit, outputintervel
	integer ncvx(nx1),ncvy(ny1),ncvz(nz1)
	character(len=256) :: toolpath_file = './ToolFiles/B26.crs'
	character(len=128) :: case_name = 'default'
	character(len=256) :: result_dir = './result/default/'
	character(len=256) :: file_prefix = './result/default/default_'

	namelist / process_parameters /alaspow, alaseta, alasrb, alasfact
	namelist / volumetric_parameters/ alaspowvol, alasetavol, sourcerad, sourcedepth
	namelist / material_properties /dens, denl, viscos, tsolid, tliquid, tboiling, hlatnt, acpa,acpb, acpl,  &
		thconsa, thconsb, thconl,  beta, emiss, dgdtp
	namelist / powder_properties / layerheight, pden, pcpa, pcpb, pthcona, pthconb
	namelist / numerical_relax / maxit, delt, timax, urfu, urfv, urfw, urfp, urfh
	namelist / boundary_conditions / htci, htcj, htck1, htckn, tempWest, tempEast, tempNorth, &
		tempBottom, tempPreheat, tempAmb
	integer species_flag

	! Adaptive mesh parameters
	integer  :: adaptive_flag = 0
	real(wp) :: amr_local_half_x = 1.0e-3_wp
	real(wp) :: amr_local_half_y = 2.0e-4_wp
	real(wp) :: amr_dx_fine = 10.0e-6_wp
	integer  :: remesh_interval = 20

	! Enthalpy prediction flag
	integer :: predict_flag = 0   ! 0=off, 1=integer-cell shift prediction during heating

	! Mechanical solver parameters
	integer :: mechanical_flag = 0        ! 0=off, 1=on
	integer :: mech_interval = 10         ! solve mechanical every N thermal steps
	integer :: mech_output_interval = 5   ! output mech VTK every N mechanical solves
	integer :: mech_mesh_ratio = 2        ! mechanical grid coarsening ratio vs thermal (1=same, 2=half cells, etc.)
	integer :: n_thermal_threads = 0      ! 0=use all OMP threads (set from env PHOENIX_MECH_THREADS)
	integer :: n_mech_threads = 0         ! 0=serial mechanical (in-loop), >0=parallel

	namelist / output_control / outputintervel, case_name, toolpath_file, species_flag, predict_flag
	namelist / adaptive_mesh / adaptive_flag, amr_local_half_x, amr_local_half_y, amr_dx_fine, remesh_interval
	namelist / mechanical_params / mechanical_flag, mech_interval, mech_output_interval, mech_mesh_ratio

	contains

subroutine read_data

	integer i,j,k

	species_flag = 0  ! default: no species transport

	open(unit=10,file='./inputfile/input_param.txt',form='formatted')
	
!-----geometrical parameters-------------------------------
	
	read(10,*)			  !
	read(10,*) nzx	                  !
	read(10,*) (xzone(i),i=1,nzx)	  !
	read(10,*) (ncvx(i),i=1,nzx)	  !
	read(10,*) (powrx(i),i=1,nzx)	  !
	read(10,*) nzy			  !
	read(10,*) (yzone(i),i=1,nzy)	  !
	read(10,*) (ncvy(i),i=1,nzy)      !
	read(10,*) (powry(i),i=1,nzy)	  !
	read(10,*) nzz	
	read(10,*) (zzone(i),i=1,nzz)	
	read(10,*) (ncvz(i),i=1,nzz)	
	read(10,*) (powrz(i),i=1,nzz)	
 
	READ (10, NML=process_parameters)  
        !read process parameters
	READ (10, NML=volumetric_parameters) 

	READ (10, NML=material_properties) 
	!read material properties

	READ (10, NML=powder_properties) 
	!read powder_properties

 	READ (10, NML=numerical_relax)
	!read numerical parameters

	READ (10, NML=boundary_conditions)
	!read boudary parameters

	READ (10, NML=output_control)
	! read output control parameters

	READ (10, NML=adaptive_mesh)
	! read adaptive mesh parameters

	READ (10, NML=mechanical_params)
	! read mechanical solver parameters

	close(10)    ! close file 10

	! Build result directory and file prefix from case_name
	result_dir = './result/' // trim(adjustl(case_name)) // '/'
	file_prefix = trim(result_dir) // trim(adjustl(case_name)) // '_'
	call execute_command_line('mkdir -p ' // trim(result_dir))

	! Read thread allocation from environment (set by run.sh)
	block
		character(len=32) :: envval
		integer :: ios
		call get_environment_variable('PHOENIX_THERMAL_THREADS', envval, status=ios)
		if (ios == 0) read(envval, *, iostat=ios) n_thermal_threads
		call get_environment_variable('PHOENIX_MECH_THREADS', envval, status=ios)
		if (ios == 0) read(envval, *, iostat=ios) n_mech_threads
	end block

	return
end subroutine read_data
end module parameters
