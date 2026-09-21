program eBDIMS2_mod
	use omp_lib
	implicit none

	interface
		subroutine get_CA(pdb_id, chain, x_coord, y_coord, z_coord, bead_mass, bead_name, &
                    bead_chain, bead_number)
			character (len = 100), intent(in) :: pdb_id, chain
			double precision, intent(out), allocatable :: x_coord(:), y_coord(:), z_coord(:), bead_mass(:)
			character (len = 3), intent(out), allocatable :: bead_name(:)
			character (len = 1), intent(out), allocatable :: bead_chain(:)
			integer, intent(out), allocatable :: bead_number(:)
		end subroutine get_CA

		subroutine get_force_res_list(step,n_beads,x,y,z,chain_id_ref,&
				num_force_list,list,num_step_print,res_numb_ref)
			integer, intent(in) :: n_beads, step
			double precision, intent(in) :: x(n_beads), y(n_beads), z(n_beads)
			character (len = 1), intent(in) :: chain_id_ref(:)
			integer, intent(out) :: num_force_list
			integer, allocatable, intent(out) :: list(:,:)
			integer, intent(in) :: num_step_print
			integer, intent(in) :: res_numb_ref(n_beads)
		end subroutine get_force_res_list
			
		subroutine get_DIMS_res_list(n_beads_ref,r_numb_ref,r_chain_ref,&
		n_beads_tar,r_numb_tar,r_chain_tar,num_interaction_list, DIMS_list,&
		num_correspondent_residues,list,chains_ref,chains_tar)
			integer, intent(in) :: n_beads_ref, n_beads_tar
			integer, intent(in) :: r_numb_ref(n_beads_ref), r_numb_tar(n_beads_tar)
			integer, intent(out) :: num_interaction_list
			integer, allocatable, intent(out) :: DIMS_list(:,:)
			character (len = 1), intent(in) :: r_chain_ref(n_beads_ref), r_chain_tar(n_beads_tar)
			integer, intent(out) :: num_correspondent_residues
			integer, allocatable, intent(out) :: list(:,:)
			character (len = 100), intent(in) :: chains_ref, chains_tar
		end subroutine get_DIMS_res_list

	end interface	

	character (len = 100) :: pdb_ref, pdb_tar
	character (len = 10) :: save_frames_frequency_str
	character (len = 100) :: chain_ID_ref, chain_ID_tar
	integer :: count, n_beads_ref, n_beads_tar, i, j, seq_dist, index
	integer :: num_interactions, num_force_interactions
	double precision, allocatable :: x_ref(:), y_ref(:), z_ref(:), mass_ref(:)
	double precision, allocatable :: x_tar(:), y_tar(:), z_tar(:), mass_tar(:)
	double precision :: rcmx, rcmy, rcmz
	double precision :: DIMS_parameter, DIMS_parameter_curr
	double precision :: DIMS_parameter0, DIMS_progress
	double precision :: dist_ref, dist_tar, dist_curr
	integer, allocatable :: DIMS_list(:,:), force_list(:,:)
	integer :: time(3), iseed, time_fin(3), time_elaps, time_curr(3), time_prev(3), time_cum
	double precision, allocatable :: v(:), vx(:), vy(:), vz(:)
	double precision :: theta, phi, dran_u
	double precision :: exp_const
	double precision :: spring_cost, dist, dist0, delta_dist
	integer :: n_frames, runs, n_improvements, save_frames_frequency
	double precision, allocatable :: randx(:), randy(:), randz(:)
	double precision, allocatable :: rx(:), ry(:), rz(:)
	double precision, allocatable :: fx(:), fy(:), fz(:)
	double precision, allocatable :: gamma(:), c_v(:), c_r(:)
	double precision, allocatable :: rx_ckp(:), ry_ckp(:), rz_ckp(:)
	double precision, allocatable :: vx_ckp(:), vy_ckp(:), vz_ckp(:)
	character (len = 3), allocatable :: res_name_ref(:), res_name_tar(:)
	character (len = 1), allocatable :: res_chain_ref(:), res_chain_tar(:)
	integer, allocatable :: res_numb_ref(:), res_numb_tar(:)
	character (len = 18) :: output_file
	character (len = 7) :: auxstr
	logical :: exist
	double precision :: DIMS_convergence
	character (len = 10) :: DIMS_convergence_str, itempsmax_str
	integer :: flag_DIMS_conv, itempsmax
	character (len = 7) :: couple
	integer :: num_correspondent_residues
	integer, allocatable :: corr_res_list(:,:)


	double precision, parameter :: pi = 3.141593d0					!PI constant
	double precision, parameter :: k_B = 1.380649d-23				!Boltzmann constant
	double precision, parameter :: uma_kg_conv = 1.66054d-27			!Conversion factor between uma and kg
	double precision, parameter :: conv_kcal_N_stiff = 0.69477d0			!Conversion factor between kcal/mol/Å^2 and N/m
		
	double precision, parameter :: dt = 2d-15 					!BD time step: 2 fs (in seconds)
	double precision, parameter :: tau = 2.5d-12					!BD solvent time constant : 2.5 ps (in seconds)
	double precision, parameter :: temp = 300.0d0 					!Temperature - K

	double precision, parameter :: Cseq = 600, Ccart = 6   				!edENM spring constants short- (x10) and long-range interactions
	integer, parameter :: exp_seq = 2, exp_cart = 6					!edENM pwer decays for short- and long-range interactions
	integer, parameter :: Slim = 3, cutoff = 10 					!edENM cutoffs for short- and long-range interactions

	integer, parameter :: num_step_print_force_list = 100000 			!Number of steps force list results are printed
	integer, parameter :: num_step_force_list_update = 1000 			!Number of steps the force list is updated
	integer, parameter :: k_steps = 10  						!eBDIMS biasing frequency	

	integer, parameter :: num_threads = 16 						!Number of threads for OpenMP parallelization


	!$ call omp_set_num_threads(num_threads)

	!!!! READ INPUTS !!!!

	call getarg(1, pdb_ref)
	call getarg(2, chain_ID_ref)
	call getarg(3, pdb_tar)
	call getarg(4, chain_ID_tar)
	call getarg(5, save_frames_frequency_str)
	call getarg(6, DIMS_convergence_str)
	call getarg(7, itempsmax_str)

	read (save_frames_frequency_str, *) save_frames_frequency
	read (DIMS_convergence_str, *) DIMS_convergence
	read (itempsmax_str, *) itempsmax

	print *, "Computing eBDIMS2 transition"
	print *, "Reference: "//trim(pdb_ref)//" - chain(s) "//trim(chain_ID_ref)
	print *, "Target: "//trim(pdb_tar)//" - chain(s) "//trim(chain_ID_tar)
	print *

	call get_CA(pdb_ref, chain_ID_ref, x_ref, y_ref, z_ref, mass_ref, res_name_ref, res_chain_ref, res_numb_ref)
	call get_CA(pdb_tar, chain_ID_tar, x_tar, y_tar, z_tar, mass_tar, res_name_tar, res_chain_tar, res_numb_tar)

	mass_ref = mass_ref * uma_kg_conv
	mass_tar = mass_tar * uma_kg_conv
	
	n_beads_ref = size(mass_ref)
	n_beads_tar = size(mass_tar)

	rcmx = 0.0d0
	rcmy = 0.0d0
	rcmz = 0.0d0
	!$omp parallel do private(i) reduction(+:rcmx,rcmy,rcmz)
	do i = 1, n_beads_ref
		rcmx = rcmx + x_ref(i)*mass_ref(i)
		rcmy = rcmy + y_ref(i)*mass_ref(i)
		rcmz = rcmz + z_ref(i)*mass_ref(i)
	end do
	!$omp end parallel do
	rcmx = rcmx/sum(mass_ref)
	rcmy = rcmy/sum(mass_ref)
	rcmz = rcmz/sum(mass_ref)

	!$omp parallel do private(i)
	do i = 1, n_beads_ref
		x_ref(i) = x_ref(i) - rcmx
		y_ref(i) = y_ref(i) - rcmy
		z_ref(i) = z_ref(i) - rcmz
	end do
	!$omp end parallel do
	
	!$omp parallel do private(i)
	do i = 1, n_beads_tar
		x_tar(i) = x_tar(i) - rcmx
		y_tar(i) = y_tar(i) - rcmy
		z_tar(i) = z_tar(i) - rcmz
	end do
	!$omp end parallel do

	call get_DIMS_res_list(n_beads_ref, res_numb_ref, res_chain_ref, n_beads_tar, res_numb_tar, res_chain_tar, &
	num_interactions, DIMS_list, num_correspondent_residues, corr_res_list, chain_ID_ref, chain_ID_tar)

	DIMS_parameter = 0.0d0

	!$omp parallel do private(i, dist_ref, dist_tar) reduction(+:DIMS_parameter)
	do i = 1,num_interactions
		dist_ref = sqrt((x_ref(DIMS_list(i,1)) - x_ref(DIMS_list(i,2)))**2 &
		+ (y_ref(DIMS_list(i,1)) - y_ref(DIMS_list(i,2)))**2 &
		+ (z_ref(DIMS_list(i,1)) - z_ref(DIMS_list(i,2)))**2)

		dist_tar = sqrt((x_tar(DIMS_list(i,3)) - x_tar(DIMS_list(i,4)))**2 &
		+ (y_tar(DIMS_list(i,3)) - y_tar(DIMS_list(i,4)))**2 &
	   	+ (z_tar(DIMS_list(i,3)) - z_tar(DIMS_list(i,4)))**2)

		DIMS_parameter = DIMS_parameter + (dist_ref - dist_tar)**2
	end do
	!$omp end parallel do

	DIMS_parameter0 = DIMS_parameter

    	call get_force_res_list(0,n_beads_ref,x_ref,y_ref,z_ref, res_chain_ref,num_force_interactions,&
    				force_list,num_step_print_force_list,res_numb_ref)

	!!! INITIALIZATION OF VELOCITIES !!!

	call itime(time)
	time_prev = time
	time_cum = 0

	open(20,file='log_time.txt',status='new')
	write(20,'(a29)') "#Time (s)    DIMS conv (%)   "

	open(53,file='log.txt',status='new')
900	format(a28)
	write(53,900) "Welcome: eBDIMS is starting!"

910 format(a14,1x,i2,1x,i2,1x,i2)
    write(53,910) "Starting time:", time(1),time(2),time(3)

	iseed = time(2) + time(3)
	call amrset(iseed)

	allocate(v(n_beads_ref))
	allocate(vx(n_beads_ref))
	allocate(vy(n_beads_ref))
	allocate(vz(n_beads_ref))

	call ranvel(n_beads_ref,v,mass_ref,temp) ! Velocity outputs in m/s
    	iseed = iseed + 20
    	call dran_ini(iseed)

	do i = 1, n_beads_ref
		theta = -pi/2.d0 + pi*dran_u()
	    	phi = 2.0*pi*dran_u()
	    	vx(i) = abs(v(i))*cos(theta)*cos(phi)
	    	vy(i) = abs(v(i))*cos(theta)*sin(phi)
	    	vz(i) = abs(v(i))*sin(theta)
	end do

	!!! START BROWNIAN SIMULATION !!!
	exp_const = exp(-dt/tau)
	
	allocate(gamma(n_beads_ref))
	allocate(c_v(n_beads_ref))
	allocate(c_r(n_beads_ref))
	
	!$omp parallel do private(i)
	do i = 1, n_beads_ref
		gamma(i) = mass_ref(i)/tau						! Gamma in kg/s
		c_v(i) = exp_const*(sqrt(2*k_B*temp*gamma(i)*dt)/mass_ref(i))		! Velocity coefficient in m/s
		c_r(i) = (1 - exp_const)*sqrt(2*k_B*temp*dt/gamma(i)) * (1.d10)		! Displacement coefficient in Å (E10 is applied to convert m to Å)
	end do
	!$omp end parallel do

	allocate(randx(n_beads_ref))
	allocate(randy(n_beads_ref))
	allocate(randz(n_beads_ref))

	allocate(rx(n_beads_ref))
	allocate(ry(n_beads_ref))
	allocate(rz(n_beads_ref))

	!$omp parallel do private(i)
	do i = 1, n_beads_ref
		rx(i) = x_ref(i)
		ry(i) = y_ref(i)
		rz(i) = z_ref(i)
	end do
	!$omp end parallel do

	allocate(fx(n_beads_ref))
	allocate(fy(n_beads_ref))
	allocate(fz(n_beads_ref))

	allocate(rx_ckp(n_beads_ref))
	allocate(ry_ckp(n_beads_ref))
	allocate(rz_ckp(n_beads_ref))

	allocate(vx_ckp(n_beads_ref))
	allocate(vy_ckp(n_beads_ref))
	allocate(vz_ckp(n_beads_ref))

	n_improvements = 0
	n_frames = 0

	rx_ckp = rx
	ry_ckp = ry
	rz_ckp = rz

	vx_ckp = vx
	vy_ckp = vy
	vz_ckp = vz

	inquire(file='eBDIMS_confrms_list.txt',exist=exist)
	if (exist) then
		open(54,file = 'eBDIMS_confrms_list.txt',status = 'old')
		close(54,status='delete')
	end if
	open(55,file='eBDIMS_confrms_list.txt',status='new')

	DIMS_progress = 0.0d0
	flag_DIMS_conv = 0

	do runs = 1,itempsmax

		call dran_gv(randx,n_beads_ref)
		call dran_gv(randy,n_beads_ref)
		call dran_gv(randz,n_beads_ref)

		fx = 0.0d0
		fy = 0.0d0
		fz = 0.0d0

		if (mod(runs,num_step_force_list_update) == 0) then
			deallocate(force_list)

			call get_force_res_list(runs,n_beads_ref,rx,ry,rz, &
				res_chain_ref,num_force_interactions,force_list,num_step_print_force_list,res_numb_ref)
		end if

		!$omp parallel do private(index,i,j,dist,dist0,seq_dist,spring_cost,delta_dist,couple) reduction(+:fx,fy,fz)
		do index = 1, num_force_interactions
			i = force_list(index,1)
			j = force_list(index,2)
			dist = sqrt((rx(i) - rx(j))**2 + (ry(i) - ry(j))**2 + (rz(i) - rz(j))**2)
			dist0 = sqrt((x_ref(i) - x_ref(j))**2 + (y_ref(i) - y_ref(j))**2 &
			+ (z_ref(i) - z_ref(j))**2)
			
		        if (res_chain_ref(j).eq.res_chain_ref(i)) then
		        	seq_dist = abs(res_numb_ref(j) - res_numb_ref(i))
		        	if (seq_dist.le.Slim) then
					spring_cost = Cseq/((dble(seq_dist))**(exp_seq))
				else
		        	        if (dist.le.cutoff) then
		        			spring_cost = (Ccart/dist)**exp_cart
		                	else
		                    		spring_cost = 0.0d0
		                	end if
				end if
		        else
		        	if (dist.le.cutoff) then
		        	        spring_cost = (Ccart/dist)**exp_cart
		        	else
		                	spring_cost = 0.0d0
		            	end if
		        end if

            		spring_cost = spring_cost*conv_kcal_N_stiff 			! Spring const. in N/m (conversion from kcal/molÅ² to N/m)
			delta_dist = (dist - dist0)*1.d-10				! Displacement in m (conversion from Å to m)
			
			fx(i) = fx(i) + spring_cost*delta_dist*(-rx(i) + rx(j))/dist	! All forces in N
			fy(i) = fy(i) + spring_cost*delta_dist*(-ry(i) + ry(j))/dist
			fz(i) = fz(i) + spring_cost*delta_dist*(-rz(i) + rz(j))/dist
			fx(j) = fx(j) + spring_cost*delta_dist*(rx(i) - rx(j))/dist
			fy(j) = fy(j) + spring_cost*delta_dist*(ry(i) - ry(j))/dist
			fz(j) = fz(j) + spring_cost*delta_dist*(rz(i) - rz(j))/dist
		end do
		!$omp end parallel do

		!$omp parallel do private(i)
		do i = 1, n_beads_ref
			vx(i) = vx(i)*exp_const + fx(i)*(1-exp_const)/gamma(i) + c_v(i)*randx(i)	! All velocities in m/s
			vy(i) = vy(i)*exp_const + fy(i)*(1-exp_const)/gamma(i) + c_v(i)*randy(i)
			vz(i) = vz(i)*exp_const + fz(i)*(1-exp_const)/gamma(i) + c_v(i)*randz(i)

			rx(i) = rx(i) + (1-exp_const)*tau*vx(i)*1.d10 &
			+ (1-(tau/dt)*(1-exp_const))*fx(i)*(dt/gamma(i))*1.d10 + c_r(i)*randx(i)	! Coordinates in Å (E10 to convert from m to Å)
			ry(i) = ry(i) + (1-exp_const)*tau*vy(i)*1.d10 &
			+ (1-(tau/dt)*(1-exp_const))*fy(i)*(dt/gamma(i))*1.d10 + c_r(i)*randy(i)
			rz(i) = rz(i) + (1-exp_const)*tau*vz(i)*1.d10 &
			+ (1-(tau/dt)*(1-exp_const))*fz(i)*(dt/gamma(i))*1.d10 + c_r(i)*randz(i)
		end do
		!$omp end parallel do

		if (mod(runs,k_steps) == 0) then
			DIMS_parameter_curr = 0.0d0
			!$omp parallel do private(i, dist_curr, dist_tar) reduction(+:DIMS_parameter_curr)
			do i = 1,num_interactions
				dist_curr = sqrt((rx(DIMS_list(i,1)) - rx(DIMS_list(i,2)))**2 &
				+ (ry(DIMS_list(i,1)) - ry(DIMS_list(i,2)))**2 &
				+ (rz(DIMS_list(i,1)) - rz(DIMS_list(i,2)))**2)

				dist_tar = sqrt((x_tar(DIMS_list(i,3)) - x_tar(DIMS_list(i,4)))**2 &
				+ (y_tar(DIMS_list(i,3)) - y_tar(DIMS_list(i,4)))**2 &
				+ (z_tar(DIMS_list(i,3)) - z_tar(DIMS_list(i,4)))**2)

				DIMS_parameter_curr = DIMS_parameter_curr + (dist_curr - dist_tar)**2
			end do
			!$omp end parallel do

			if (DIMS_parameter_curr.le.DIMS_parameter) then
				DIMS_parameter = DIMS_parameter_curr
				DIMS_progress = (1 - DIMS_parameter/DIMS_parameter0)*1.d2
				rx_ckp = rx
				ry_ckp = ry
				rz_ckp = rz
				vx_ckp = vx
				vy_ckp = vy
				vz_ckp = vz
				n_improvements = n_improvements + 1
				if (mod(n_improvements,save_frames_frequency).eq.0) then
					n_frames = n_frames + 1

210					format(i7.7)
					write(auxstr,210) runs/10
					output_file = "DIMS_MD"//auxstr//".pdb"

230					format(a18)
					write(55,230) output_file

220					format(a4,i7,2x,a3,1x,a3,1x,a1,i4,4x,3f8.3)
					open(60+n_frames,file=output_file,status='new')
					do i = 1,n_beads_ref
						write(60+n_frames,220) "ATOM",i,"CA ",res_name_ref(i),res_chain_ref(i),&
						res_numb_ref(i),rx(i),ry(i),rz(i)
					end do
					close(60+n_frames)

240					format(a18,2x,a6,f18.3,2x,a11,f6.2,a1)
					write(53,240) output_file, "DIMS: ", DIMS_parameter, "DIMS conv: ", DIMS_progress,"%"
					print '(a18,2x,a6,f18.3,2x,a11,f6.2,a1)', output_file, "DIMS: ", DIMS_parameter,&
					"DIMS conv: ", DIMS_progress, "%"

					call itime(time_curr)
					time_elaps = 3600*(time_curr(1) - time_prev(1)) + 60*(time_curr(2) - time_prev(2)) + &
					(time_curr(3) - time_prev(3))
					if (time_elaps < 0) then
						time_elaps = time_elaps + 86400
					end if
					time_cum = time_cum + time_elaps
					time_prev = time_curr

250					format(i9,4x,f13.2)
					write(20,250) time_cum, DIMS_progress
				end if
				
				if (DIMS_progress.ge.DIMS_convergence) then
					flag_DIMS_conv = 1

					n_frames = n_frames + 1
					write(auxstr,210) runs/10
					output_file = "DIMS_MD"//auxstr//".pdb"
					write(55,230) output_file

					open(60+n_frames,file=output_file,status='new')
					do i = 1,n_beads_ref
						write(60+n_frames,220) "ATOM",i,"CA ",res_name_ref(i),res_chain_ref(i),&
						res_numb_ref(i),rx(i),ry(i),rz(i)
					end do
					close(60+n_frames)

					write(53,240) output_file, "DIMS: ", DIMS_parameter, "DIMS conv: ", DIMS_progress,"%"
					print '(a18,2x,a6,f18.3,2x,a11,f6.2,a1)', output_file, "DIMS: ", DIMS_parameter,&
					"DIMS conv: ", DIMS_progress, "%"

					call itime(time_curr)
					time_elaps = 3600*(time_curr(1) - time_prev(1)) + 60*(time_curr(2) - time_prev(2)) + &
					(time_curr(3) - time_prev(3))
					if (time_elaps < 0) then
						time_elaps = time_elaps + 86400
					end if
					time_cum = time_cum + time_elaps
					time_prev = time_curr
					write(20,250) time_cum, DIMS_progress

					exit
				end if
			else
				rx = rx_ckp
				ry = ry_ckp
				rz = rz_ckp

				iseed = iseed + 1
				call amrset(iseed)
				call ranvel(n_beads_ref,v,mass_ref,temp)
    				iseed = iseed + 20
    				call dran_ini(iseed)
    			
				do i = 1,n_beads_ref
					theta = -pi/2.d0 + pi*dran_u()
    					phi = 2.0*pi*dran_u()
    					vx(i) = abs(v(i))*cos(theta)*cos(phi)
    					vy(i) = abs(v(i))*cos(theta)*sin(phi)
    					vz(i) = abs(v(i))*sin(theta)
				end do
			end if
		end if
	end do
	
	close(55)

	if (flag_DIMS_conv.eq.1) then
920 	format(a39,f5.2,a13)
		write(53,920) "Nice job: eBDIMS transition completed (", DIMS_convergence, "% DIMS conv)!"
	else
970 	format(a85)
		write(53,970) "Convergence not reached in the number of steps used, but eBDIMS transition completed!"
	end if

930 format(a27,i8)
	write(53,930) "Total number of used steps:", runs

	call itime(time_fin)
940	format(a15,1x,i2,1x,i2,1x,i2)
	write(53,940) "Finishing time:", time_fin(1),time_fin(2),time_fin(3)

950	format(a13,i8,a4)
	write(53,950) "Elapsed time:", time_cum, " sec"

    close(53)
    close(20)
end program eBDIMS2_mod

subroutine get_CA(pdb_id, chain, x_coord, y_coord, z_coord, bead_mass, bead_name, bead_chain, bead_number)
	implicit none
	
	character (len = 100), intent(in) :: pdb_id, chain
	character (len = 1) :: alt_loc, chain_id
	character (len = 3) :: atom_type
	character (len = 4) :: label
	character (len = 3) :: res_type
    	character (len = 80) :: string
	integer :: num_beads, n, atom_num, res_num, pdb_iostat, num_chains, num_atoms, i
	double precision :: x, y, z
	double precision, intent(out), allocatable :: x_coord(:), y_coord(:), z_coord(:), bead_mass(:)
	character (len = 3), intent(out), allocatable :: bead_name(:)
	character (len = 1), intent(out), allocatable :: bead_chain(:)
	integer, intent(out), allocatable :: bead_number(:)
	logical :: exist

	open(11,file = trim(pdb_id),status = 'old',iostat = pdb_iostat)
	if (pdb_iostat .ne. 0) then
		print *, "I couldn't open the file "//trim(pdb_id)//" or the file does not exist!"
		stop
	end if
	
	inquire(file = "ATOM.pdb", exist=exist)
	if (exist) then
		open(12,file = "ATOM.pdb",status = 'old')
		close(12,status='delete')
	end if
	
	open(13,file = "ATOM.pdb",status = 'new')
	num_atoms = 0
	do
	    read(11,'(a80)',end = 10) string
	    if (trim(string(1:4)) == 'ATOM') then
	    	num_atoms = num_atoms + 1
	    	write (13,'(a80)') string
	    end if
	end do
10  close(11)
	close(13)
	
	! format :ATOM      1  CA  LYS A   2      -0.921  45.314  37.381 
20  format(a4,i7,2x,a3,a1,a3,1x,a1,i4,4x,3f8.3)

	num_chains = len(trim(chain))
	num_beads = 0
	do i = 1, num_chains
		open(16,file = "ATOM.pdb",status = 'old')
		do n = 1, num_atoms
			read(16,20) label, atom_num, atom_type, alt_loc, res_type, chain_id, res_num, x, y, z
			if (atom_type(1:2) == 'CA') then
				if (alt_loc == ' ' .or. alt_loc == 'A') then
					if (chain_id == chain(i:i)) then
						num_beads = num_beads + 1
					end if
				end if
			end if
		end do
		close(16)
	end do
	close(15)
	
	allocate (x_coord(num_beads))
	allocate (y_coord(num_beads))
	allocate (z_coord(num_beads))
	allocate (bead_mass(num_beads))
	allocate (bead_name(num_beads))
	allocate (bead_chain(num_beads))
	allocate (bead_number(num_beads))

	num_beads = 0
	do i = 1, num_chains
		open(12,file = "ATOM.pdb",status = 'old')
		do n = 1, num_atoms
			read(12,20) label, atom_num, atom_type, alt_loc, res_type, chain_id, res_num, x, y, z
			if (atom_type(1:2) == 'CA') then
				if (alt_loc == ' ' .or. alt_loc == 'A') then
					if (chain_id == chain(i:i)) then
						num_beads = num_beads + 1
						x_coord(num_beads) = x
						y_coord(num_beads) = y
						z_coord(num_beads) = z
						bead_chain(num_beads) = chain_id
						bead_number(num_beads) = res_num
						bead_name(num_beads) = res_type
						if (res_type == 'ALA') then
							bead_mass(num_beads) = 71
						else if (res_type == 'ARG') then
							bead_mass(num_beads) = 156
						else if (res_type == 'ASN') then
							bead_mass(num_beads) = 114
						else if (res_type == 'ASP') then
							bead_mass(num_beads) = 115
						else if (res_type == 'CYS') then
							bead_mass(num_beads) = 103
						else if (res_type == 'GLU') then
							bead_mass(num_beads) = 129
						else if (res_type == 'GLH') then
							bead_mass(num_beads) = 129
						else if (res_type == 'GLN') then
							bead_mass(num_beads) = 128
						else if (res_type == 'GLY') then
							bead_mass(num_beads) = 57
						else if (res_type == 'HIS') then
							bead_mass(num_beads) = 137
						else if (res_type == 'HIE') then
							bead_mass(num_beads) = 137
						else if (res_type == 'HID') then
							bead_mass(num_beads) = 137
						else if (res_type == 'HIP') then
							bead_mass(num_beads) = 137
						else if (res_type == 'ILE') then
							bead_mass(num_beads) = 113
						else if (res_type == 'LEU') then
							bead_mass(num_beads) = 113
						else if (res_type == 'LYS') then
							bead_mass(num_beads) = 128
						else if (res_type == 'LYN') then
							bead_mass(num_beads) = 128
						else if (res_type == 'MET') then
							bead_mass(num_beads) = 131
						else if (res_type == 'PHE') then
							bead_mass(num_beads) = 147
						else if (res_type == 'PRO') then
							bead_mass(num_beads) = 97
						else if (res_type == 'SER') then
							bead_mass(num_beads) = 87
						else if (res_type == 'THR') then
							bead_mass(num_beads) = 101
						else if (res_type == 'TRP') then
							bead_mass(num_beads) = 186
						else if (res_type == 'TYR') then
							bead_mass(num_beads) = 163
						else if (res_type == 'VAL') then
							bead_mass(num_beads) = 99
						else
							bead_mass(num_beads) = 110
						end if
					end if
				end if
			end if
		end do
		close(12)
	end do
	
	open(17,file = "ATOM.pdb",status = 'old')
	close(17,status='delete')
	
end subroutine get_CA

subroutine amrset(iseed)
	implicit none

	integer :: iseed
	
	double precision :: u(97), c, cd, cm
	integer :: i97, j97
	logical :: set
	common/raset1/u,c,cd,cm,i97,j97,set
	
	integer :: is1, is2, is1max, is2max
	integer :: i, j, k, l, m
	integer :: ii, jj
	double precision :: s, t

	data is1max, is2max /31328, 30081/

	is1 = max((iseed/is2max)+1,1)
    	is1 = min(is1,is1max)

    	is2 = max(1,mod(iseed,is2max)+1)
    	is2 = min(is2,is2max)

    	i = mod(is1/177,177) + 2
    	j = mod(is1,177) + 2
    	k = mod(is2/169,178) + 1
    	l = mod(is2,169)

    	do ii = 1,97
    		s = 0.0d0
    		t = 0.5d0
    		do jj = 1,24
            		m = mod(mod(i*j,179)*k,179)
            		i = j
            		j = k
            		k = m
            		l = mod(53*l+1,169)
            		if(mod(l*m,64).ge.32) then
            			s = s+t
            		end if
            		t = 0.5d0*t
    		end do
    		u(ii) = s
   	end do

    	c = 362436.0d0/16777216.0d0
    	cd = 7654321.0d0/16777216.0d0
    	cm = 16777213.0d0/16777216.0d0

    	i97 = 97
    	j97 = 33

    	set = .true.
    	return
end subroutine amrset

subroutine gauss(am,sd,v)
	implicit none

	double precision :: u(97), c, cd, cm
	integer :: i97, j97
	logical :: set
	common/raset1/u,c,cd,cm,i97,j97,set
	double precision :: a, uni, am, sd, v, zero, six
	integer :: i

	data zero, six /0.0d0,6.0d0/
	
	if (.not. set) then
		print *, "amrset not called!"
		stop
	end if

	a = zero
	do i = 1,12
		uni = u(i97) - u(j97)
		if (uni.lt.0.0d0) then
			uni = uni+1.0d0
		end if
		u(i97) = uni

		i97 = i97 - 1
        	if (i97.eq.0) then
        		i97 = 97
        	end if
        	j97 = j97 - 1
        	if (j97.eq.0) then
        		j97 = 97
        	end if

        	c = c - cd
        	if (c.lt.0.0d0) then
        		c = c + cm
        	end if
        	uni = uni - c
        	if (uni.lt.0.0d0) then
        		uni=uni+1.0d0
        	end if
        	a=a+uni      
	end do

	v = (a - six)*sd + am
	return
end subroutine gauss

subroutine ranvel(nrp,v,mass,temp)
	implicit none

	integer :: nrp, i, j
	double precision :: v(nrp), mass(nrp)
	double precision :: temp, y, sd
	double precision, parameter :: k_B = 1.380649d-23

	do i = 1,nrp
		sd = sqrt(k_B*temp/mass(i))
		call gauss(0.0d0,sd,v(i))	! Velocity in m/s
	end do
	return
end subroutine ranvel


function dran_u()
	implicit none

    	integer, parameter :: ip = 1279
    	integer, parameter :: iq = 418
    	integer, parameter :: is = ip - iq
    	double precision, parameter :: rmax = 2147483647.0d0
    	integer :: ix(ip), ic
    	double precision :: dran_u

    	common /ixx/ ix
    	common /icc/ ic
    
    	ic = ic + 1
    	if(ic.gt.ip) then
    		ic = 1
   		end if
    	if(ic.gt.iq) then
        	ix(ic) = ieor(ix(ic),ix(ic-iq))
    	else
        	ix(ic) = ieor(ix(ic),ix(ic+is))
    	endif
    	dran_u = dble(ix(ic))/rmax
    	return	
end function dran_u

subroutine dran_gv(u,n)
	implicit none

	integer, parameter :: ip = 1279
    	integer, parameter :: iq = 418
    	integer, parameter :: is = ip - iq
    	integer, parameter :: np = 14
    	integer, parameter :: nbit = 31
    	integer, parameter :: m = 2**np
    	integer, parameter :: np1 = nbit - np
    	integer, parameter :: nn = 2**np1 - 1
    	integer, parameter :: nn1 = nn + 1
    	integer :: ix(ip), ic, k, n, i, i2
    	double precision :: g(0:m), u(n)

    	common /ixx/ ix
    	common /icc/ ic
    	common /gg/ g

    	do k = 1,n
        	ic = ic+1
        	if (ic.gt.ip) then
        		ic = 1
       		end if
        	if (ic.gt.iq) then
           		ix(ic) = ieor(ix(ic),ix(ic-iq))
        	else
            		ix(ic) = ieor(ix(ic),ix(ic+is))
        	endif
        	i = ishft(ix(ic),-np1)
        	i2 = iand(ix(ic),nn)
        	u(k) = i2*g(i+1) + (nn1-i2)*g(i)
    	end do
    	return
end subroutine dran_gv

subroutine dran_ini(iseed0)
	implicit none

	integer, parameter :: ip = 1279
	integer, parameter :: np = 14
    	integer, parameter :: nbit = 31
	integer, parameter :: m = 2**np
    	integer, parameter :: np1 = nbit - np
    	integer, parameter :: nn = 2**np1 - 1
    	integer, parameter :: nn1 = nn + 1
    	integer :: ix(ip), ic, i, j, iseed0
    	double precision :: c0, c1, c2, d1, d2, d3
    	double precision :: dseed, pi, p, t, x, u2th
  	double precision :: g(0:m), rand_xx

    	data c0,c1,c2/2.515517d0,0.802853d0,0.010328d0/
    	data d1,d2,d3/1.432788d0,0.189269d0,0.001308d0/

	common /ixx/ ix
    	common /icc/ ic
    	common /gg/ g

    	dseed = iseed0
    	do i = 1,ip
    		ix(i) = 0
    		do j = 0, nbit - 1
    			if (rand_xx(dseed).lt.0.5) then
    				ix(i) = ibset(ix(i),j)
    			end if
    		end do
    	end do
    	ic = 0

    	pi = 4.0d0*datan(1.0d0)
    	do i = m/2,m
		p = 1.0d0 - dble(i+1)/(m+2)
	    	t = sqrt(-2.0d0*log(p))
	    	x = t - (c0 + t*(c1 + c2*t))/(1.0 + t*(d1+t*(d2 + t*d3)))
	    	g(i) = x
	    	g(m - i) = -x
    	end do

    	u2th = 1.0d0 - dble(m + 2)/m*sqrt(2.0d0/pi)*g(m)*exp(-g(m)*g(m)/2)
    	u2th = nn1 * sqrt(u2th)
    	do i = 0,m
    		g(i) = g(i)/u2th
    	end do
	return
end subroutine dran_ini

function rand_xx(dseed)
    	implicit none

    	double precision :: dseed, rand_xx
    	double precision, parameter :: xmm = 2.d0**32
    	double precision, parameter :: rm = 1.d0/xmm
    	double precision, parameter :: a = 69069.d0
    	double precision, parameter :: c = 1.d0

    	dseed = mod(dseed*a+c,xmm)
    	rand_xx = dseed*rm
    	return
end function rand_xx

subroutine get_force_res_list(step,n_beads,x,y,z,chain_id_ref,num_force_list,list,num_step_print,res_numb_ref)
	implicit none
	
	integer, intent(in) :: n_beads, step
	double precision, intent(in) :: x(n_beads), y(n_beads), z(n_beads)
	character (len = 1), intent(in) :: chain_id_ref(:)
	integer :: i, j, seq_dist
	integer, intent(out) :: num_force_list
	integer, allocatable, intent(out) :: list(:,:)
	double precision :: dist
	integer, intent(in) :: num_step_print
	integer, intent(in) :: res_numb_ref(n_beads)

	integer, parameter :: Slim = 3, cutoff = 10
	
	num_force_list = 0
	do i = 1, n_beads - 1
		do j = i + 1, n_beads
			dist = sqrt((x(i) - x(j))**2 + (y(i) - y(j))**2 + (z(i) - z(j))**2)
			if (chain_id_ref(j).eq.chain_id_ref(i)) then
				seq_dist = abs(res_numb_ref(j) - res_numb_ref(i))
				if (seq_dist.le.Slim) then
					num_force_list = num_force_list + 1
				else
					if (dist.le.cutoff) then
				        	num_force_list = num_force_list + 1
				        end if
				end if
			else
				if (dist.le.cutoff) then
					num_force_list = num_force_list + 1
				end if
			end if
		end do
	end do

	allocate(list(num_force_list,2))

	num_force_list = 0
	do i = 1, n_beads - 1
		do j = i + 1, n_beads
			dist = sqrt((x(i) - x(j))**2 + (y(i) - y(j))**2 + (z(i) - z(j))**2)
			if (chain_id_ref(j).eq.chain_id_ref(i)) then
				seq_dist = abs(res_numb_ref(j) - res_numb_ref(i))
				if (seq_dist.le.Slim) then
					num_force_list = num_force_list + 1
					list(num_force_list,1) = i
					list(num_force_list,2) = j
				else
					if (dist.le.cutoff) then
						num_force_list = num_force_list + 1
						list(num_force_list,1) = i
						list(num_force_list,2) = j
					end if
				end if
			else
				if (dist.le.cutoff) then
					num_force_list = num_force_list + 1
					list(num_force_list,1) = i
					list(num_force_list,2) = j
				end if
			end if
		end do
	end do

	if (mod(step,num_step_print) == 0) then
		print *, "Step: ", step, " Now we have only ", num_force_list, &
		" force interactions, instead of", (n_beads**2 - n_beads)/2
	end if
end subroutine get_force_res_list

subroutine get_DIMS_res_list(n_beads_ref, r_numb_ref, r_chain_ref,&
n_beads_tar, r_numb_tar, r_chain_tar, num_interaction_list, DIMS_list,&
num_correspondent_residues, list, chains_ref, chains_tar)
	implicit none
	
	integer, intent(in) :: n_beads_ref, n_beads_tar
	integer, intent(in) :: r_numb_ref(n_beads_ref), r_numb_tar(n_beads_tar)
	integer :: i, j, k
	integer, intent(out) :: num_interaction_list
	integer, allocatable, intent(out) :: DIMS_list(:,:)
	character (len = 1), intent(in) :: r_chain_ref(n_beads_ref), r_chain_tar(n_beads_tar)
	integer, intent(out) :: num_correspondent_residues
	integer, allocatable, intent(out) :: list(:,:)
	character (len = 100), intent(in) :: chains_ref, chains_tar

	num_correspondent_residues = 0
	do i = 1, n_beads_ref
		do j = 1, n_beads_tar
			if ((r_numb_ref(i).eq.r_numb_tar(j))) then
				do k = 1, len(chains_ref)
					if ((r_chain_ref(i).eq.chains_ref(k:k)).and.(r_chain_tar(j).eq.chains_tar(k:k))) then
						num_correspondent_residues = num_correspondent_residues + 1
					end if
				end do
			end if
		end do
	end do

	allocate(list(num_correspondent_residues,2))
	num_correspondent_residues = 0
	do i = 1, n_beads_ref
		do j = 1, n_beads_tar
			if ((r_numb_ref(i).eq.r_numb_tar(j))) then
				do k = 1, len(chains_ref)
					if ((r_chain_ref(i).eq.chains_ref(k:k)).and.(r_chain_tar(j).eq.chains_tar(k:k))) then
						num_correspondent_residues = num_correspondent_residues + 1
						list(num_correspondent_residues,1) = i
						list(num_correspondent_residues,2) = j
					end if
				end do
			end if
		end do
	end do
	
	num_interaction_list = (num_correspondent_residues**2-num_correspondent_residues)/2
	allocate(DIMS_list(num_interaction_list,4))
	num_interaction_list = 0
	do i = 1, num_correspondent_residues - 1
		do j = i + 1, num_correspondent_residues
			num_interaction_list = num_interaction_list + 1
			DIMS_list(num_interaction_list,1) = list(i,1)
			DIMS_list(num_interaction_list,2) = list(j,1)
			DIMS_list(num_interaction_list,3) = list(i,2)
			DIMS_list(num_interaction_list,4) = list(j,2)
		end do
	end do

end subroutine get_DIMS_res_list
