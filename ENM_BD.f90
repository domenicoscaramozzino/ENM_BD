program ENM_BD
    	use omp_lib
    	implicit none

    	interface
        	subroutine get_CA(pdb_id, chain, x_coord, y_coord, z_coord, res_mass, res_name, res_chain, res_number)
            		character (len = 100), intent(in) :: pdb_id, chain
            		double precision, intent(out), allocatable :: x_coord(:), y_coord(:), z_coord(:), res_mass(:)
            		character (len = 3), intent(out), allocatable :: res_name(:)
            		character (len = 1), intent(out), allocatable :: res_chain(:)
            		integer, intent(out), allocatable :: res_number(:)
		end subroutine get_CA

        	subroutine get_force_res_list(step,n_CA,x,y,z,Slim,cutoff,chain_id_ref,&
        			num_force_list,list,num_step_print,res_numb_ref)
            		integer, intent(in) :: n_CA, step
            		double precision, intent(in) :: x(n_CA), y(n_CA), z(n_CA)
            		integer, intent(in) :: cutoff
            		integer, intent(in) :: Slim
            		character (len = 1), intent(in) :: chain_id_ref(:)
            		integer, intent(out) :: num_force_list
            		integer, allocatable, intent(out) :: list(:,:)
            		integer, intent(in) :: num_step_print
            		integer, intent(in) :: res_numb_ref(n_CA)
        	end subroutine get_force_res_list
	end interface	

    	character (len = 100) :: pdb_ref, chain_ID_ref
    	character (len = 10) :: itempsmax_str, save_frames_frequency_str
    	integer :: count, n_CA, i, j, seq_dist, index
    	integer :: num_force_interactions
    	double precision, allocatable :: x_ref(:), y_ref(:), z_ref(:), mass_ref(:)
    	integer, allocatable :: force_list(:,:)
    	integer :: time(3), iseed, time_fin(3), time_elaps, time_curr(3), time_prev(3), time_cum
    	double precision, allocatable :: v(:), vx(:), vy(:), vz(:)
    	double precision :: theta, phi, dran_u, exp_const
    	double precision :: spring_cost, dist, dist0, delta_dist
    	integer :: itempsmax, runs, save_frames_frequency
    	double precision, allocatable :: randx(:), randy(:), randz(:)
    	double precision, allocatable :: rx(:), ry(:), rz(:)
    	double precision, allocatable :: fx(:), fy(:), fz(:)
    	double precision, allocatable :: gamma(:), c_v(:), c_r(:)
    	character (len = 3), allocatable :: res_ref(:)
    	character (len = 1), allocatable :: res_chain_ref(:)
    	integer, allocatable :: res_numb_ref(:)
    	character (len = 16) :: output_file
    	character (len = 7) :: auxstr
    	logical :: exist
    	double precision :: RMSD_curr

    	double precision, parameter :: pi = 3.141593d0				!PI constant
    	double precision, parameter :: k_B = 1.380649d-23			!Boltzmann constant J/K
    	double precision, parameter :: uma_kg_conv = 1.66054d-27		!Conversion factor between uma and kg
    	double precision, parameter :: conv_kcal_N_stiff = 0.69477d0		!Conversion factor between kcal/mol/Å^2 and N/m
    
    	double precision, parameter :: dt = 2d-15				!BD time-step = 2 fs (in seconds)
    	double precision, parameter :: tau = 2.5d-12				!BD solvent time costant = 2.5 ps (in seconds)
    	double precision, parameter :: temp = 300.0d0 				!Temperature - K
    	
    	double precision, parameter :: Ccart = 6, Cseq = 600			!Spring constant parameters (in kcal/molÅ²)
    	integer, parameter :: Slim = 3, cutoff = 10				!Sequential cutoff (imensionless)
    	integer, parameter :: exp_cart = 6, exp_seq = 2				!Power-decay parameters (dimensionless)
    	
	integer, parameter :: num_step_print_force_list = 100000 		!Number of steps the results from the force list are printed
	integer, parameter :: num_step_force_list_update = 100 			!Number of steps the force list is updated    
    	
    	integer, parameter :: num_threads = 16 					!Number of threads for OpenMP parallelization    	
    	
    	!$ call omp_set_num_threads(num_threads)

    	!!!! READ INPUTS !!!!

    	call getarg(1, pdb_ref)
    	call getarg(2, chain_ID_ref)
    	call getarg(3, itempsmax_str)
    	call getarg(4, save_frames_frequency_str)

    
    	read (itempsmax_str, *) itempsmax
    	read (save_frames_frequency_str, *) save_frames_frequency


    	!!!! READ C-ALPHAS OF REFERENCE STRUCTURE
    	!!!! READ INTERACTION LIST FOR FORCE COMPUTATION

    	call get_CA(pdb_ref, chain_ID_ref, x_ref, y_ref, z_ref, mass_ref, res_ref, res_chain_ref, res_numb_ref)
    	n_CA = size(mass_ref)
    	mass_ref = mass_ref * uma_kg_conv
    	
    	print *
    	print *, "Welcome to BD! Computing BD of "//trim(pdb_ref)//" - chain(s) "//trim(chain_ID_ref)
    	print *
    	write(*,fmt='(a32)') " ... using edENM_mod version ..."
    	write(*,fmt='(a9,f3.1)') " Ccart = ", Ccart
    	write(*,fmt='(a8,f5.1)') " Cseq = ", Cseq
    	write(*,fmt='(a10,i2,a2)') " Cutoff = ", cutoff, " A"
    	write(*,fmt='(a8,i1)') " Slim = ", Slim
    	write(*,fmt='(a12,i1)') " Exp-cart = ", exp_cart
    	write(*,fmt='(a11,i1)') " Exp-seq = ", exp_seq
    	print *
    	write(*,fmt='(a16,i10,a22,f10.3,a3)') " Simulating for ", itempsmax, " steps. Equivalent to ", itempsmax*dt*1E9, " ns"
    	write(*,fmt='(a21,f10.3,a3)') " Saving frames every ", save_frames_frequency*dt*1E9, " ns"
    	print *
    	print *, "Note: RMSD values also include rotations and translations! Perform alignment after the simulation!"
    	print *
    	
    	call get_force_res_list(0,n_CA,x_ref,y_ref,z_ref,Slim,cutoff, res_chain_ref, num_force_interactions, &
    			force_list,num_step_print_force_list,res_numb_ref)

    	!!! INITIALIZATION OF VELOCITIES !!!

    	call itime(time)
    	time_prev = time
    	time_cum = 0

    	open(20,file='log_time.txt',status='new')
    	write(20,'(a20)') "#Time (s)   RMSD (A)"

    	open(53,file='log.txt',status='new')
900	format(a24)
    	write(53,900) "Welcome: BD is starting!"

910 	format(a14,1x,i2,1x,i2,1x,i2)
    	write(53,910) "Starting time:", time(1),time(2),time(3)

    	iseed = time(2) + time(3)
    	call amrset(iseed)

    	allocate(v(n_CA))
    	allocate(vx(n_CA))
    	allocate(vy(n_CA))
    	allocate(vz(n_CA))

    	call ranvel(n_CA,v,mass_ref,temp) ! velocity output v in m/s
    	iseed = iseed + 20
    	call dran_ini(iseed)

    	do i = 1,n_CA
        	theta = -pi/2.d0 + pi*dran_u()
        	phi = 2.0*pi*dran_u()
        	vx(i) = abs(v(i))*cos(theta)*cos(phi)
        	vy(i) = abs(v(i))*cos(theta)*sin(phi)
        	vz(i) = abs(v(i))*sin(theta)
    	end do
	
    	!!! START BROWNIAN SIMULATION !!!

    	exp_const = exp(-dt/tau)
    
    	allocate(gamma(n_CA))
    	allocate(c_v(n_CA))
    	allocate(c_r(n_CA))
    
    	!$omp parallel do private(i)
    	do i = 1,n_CA
        	gamma(i) = mass_ref(i)/tau    !gamma in kg/s!
        	c_v(i) = exp_const*(sqrt(2*k_B*temp*gamma(i)*dt)/(mass_ref(i))) ! Velocity coefficient in m/s
        	c_r(i) = (1 - exp_const)*sqrt(2*k_B*temp*dt/gamma(i)) * (1.d10) ! Displacement coefficient in Å (*1E10 is applied to convert m to Å)
    	end do
    	!$omp end parallel do
	
    	allocate(randx(n_CA))
    	allocate(randy(n_CA))
    	allocate(randz(n_CA))

    	allocate(rx(n_CA))
    	allocate(ry(n_CA))
    	allocate(rz(n_CA))

    	!$omp parallel do private(i)
    	do i = 1, n_CA
        	rx(i) = x_ref(i)
        	ry(i) = y_ref(i)
        	rz(i) = z_ref(i)
    	end do
    	!$omp end parallel do

    	allocate(fx(n_CA))
    	allocate(fy(n_CA))
    	allocate(fz(n_CA))

    	inquire(file='BD_confrms_list.txt',exist=exist)
    	if (exist) then
        	open(54,file = 'BD_confrms_list.txt',status = 'old')
        	close(54,status='delete')
    	end if
    	open(55,file='BD_confrms_list.txt',status='new')

210 format(i7.7)
230 format(a16)
220 format(a4,i7,2x,a3,1x,a3,1x,a1,i4,4x,3f8.3,2f6.2,11x,a1)
240 format(a16,2x,a6,f6.2,a1)

	runs = 0
        write(auxstr,210) runs/10
        output_file = "BD_MD"//auxstr//".pdb"
        write(55,230) output_file
        open(60,file=output_file,status='new')
        do i = 1,n_CA
                write(60,220) "ATOM",i,"CA ",res_ref(i),res_chain_ref(i),res_numb_ref(i),rx(i),ry(i),rz(i),1.0d0,0.0d0,"C"
        end do
        close(60)
            		
    	do runs = 1,itempsmax

        	call dran_gv(randx,n_CA)
        	call dran_gv(randy,n_CA)
        	call dran_gv(randz,n_CA)

        	fx = 0.0d0
        	fy = 0.0d0
        	fz = 0.0d0

        	if (mod(runs,num_step_force_list_update) == 0) then
        		deallocate(force_list)
			call get_force_res_list(runs,n_CA,rx,ry,rz,Slim,cutoff, &
            			res_chain_ref,num_force_interactions,force_list,num_step_print_force_list,res_numb_ref)
        	end if

        	!$omp parallel do private(index,i,j,dist,dist0,seq_dist,spring_cost,delta_dist) reduction(+:fx,fy,fz)
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
            		spring_cost = spring_cost * conv_kcal_N_stiff ! Spring const. in N/m (conversion from kcal/molÅ² to N/m for spring constants)
            		delta_dist = (dist - dist0)*1.d-10            ! Displacement in m (conversion from Å to m)
            		
            		fx(i) = fx(i) + spring_cost*delta_dist*(-rx(i) + rx(j))/dist  ! All forces in N
            		fy(i) = fy(i) + spring_cost*delta_dist*(-ry(i) + ry(j))/dist
            		fz(i) = fz(i) + spring_cost*delta_dist*(-rz(i) + rz(j))/dist
            		fx(j) = fx(j) + spring_cost*delta_dist*(rx(i) - rx(j))/dist
            		fy(j) = fy(j) + spring_cost*delta_dist*(ry(i) - ry(j))/dist
            		fz(j) = fz(j) + spring_cost*delta_dist*(rz(i) - rz(j))/dist
		end do
        	!$omp end parallel do

        	!$omp parallel do private(i)
        	do i = 1, n_CA
        	
            		vx(i) = vx(i)*exp_const + fx(i)*(1-exp_const)/gamma(i) + c_v(i)*randx(i) ! All velocities in m/s
            		vy(i) = vy(i)*exp_const + fy(i)*(1-exp_const)/gamma(i) + c_v(i)*randy(i)
           		vz(i) = vz(i)*exp_const + fz(i)*(1-exp_const)/gamma(i) + c_v(i)*randz(i)
		
            		rx(i) = rx(i) + (1-exp_const)*tau*vx(i)*1.d10 &
            			+ (1-(tau/dt)*(1-exp_const))*fx(i)*(dt/gamma(i))*1.d10 + c_r(i)*randx(i) ! All coordinates in Å (*1E10 to convert m to Å)
            		ry(i) = ry(i) + (1-exp_const)*tau*vy(i)*1.d10 &
            			+ (1-(tau/dt)*(1-exp_const))*fy(i)*(dt/gamma(i))*1.d10 + c_r(i)*randy(i)
            		rz(i) = rz(i) + (1-exp_const)*tau*vz(i)*1.d10 &
            			+ (1-(tau/dt)*(1-exp_const))*fz(i)*(dt/gamma(i))*1.d10 + c_r(i)*randz(i)
        	end do
        	!$omp end parallel do

        	RMSD_curr = 0.0d0
        	!$omp parallel do private(i) reduction(+:RMSD_curr)
        	do i = 1, n_CA
            		RMSD_curr = RMSD_curr + (x_ref(i) - rx(i))**2 + (y_ref(i) - ry(i))**2 + (z_ref(i) - rz(i))**2
        	end do
        	!$omp end parallel do
        	RMSD_curr = sqrt(RMSD_curr/dble(n_CA))

        	if (mod(runs,save_frames_frequency).eq.0) then
            		write(auxstr,210) runs/10
            		output_file = "BD_MD"//auxstr//".pdb"
            		write(55,230) output_file
            		open(60,file=output_file,status='new')
            		do i = 1,n_CA
                		write(60,220) "ATOM",i,"CA ",res_ref(i),res_chain_ref(i),&
                		res_numb_ref(i),rx(i),ry(i),rz(i),1.0d0,0.0d0,"C"
            		end do
            		close(60)
            		write(53,240) output_file, "RMSD: ", RMSD_curr, "A"
            		print '(a22,a16,2x,a6,f6.2,a1)', " Saving conformation: ", output_file, "RMSD: ", RMSD_curr, "A"
            		call itime(time_curr)
            		time_elaps = 3600*(time_curr(1) - time_prev(1)) + 60*(time_curr(2) - time_prev(2)) + &
            			(time_curr(3) - time_prev(3))
            		if (time_elaps < 0) then
                		time_elaps = time_elaps + 86400
            		end if
            		time_cum = time_cum + time_elaps
            		time_prev = time_curr

250 format(i9,3x,f8.2)
            		write(20,250) time_cum, RMSD_curr
        	end if
	end do
	
	close(55)

970 format(a54)
    	write(53,970) "BD simulation completed in the chosen number of steps!"
    
930 format(a27,i8)
    	write(53,930) "Total number of used steps:", runs

    	call itime(time_fin)
940 format(a15,1x,i2,1x,i2,1x,i2)
    	write(53,940) "Finishing time:", time_fin(1),time_fin(2),time_fin(3)

950 format(a13,i8,a4)
    	write(53,950) "Elapsed time:", time_cum, " sec"

    	close(53)
    	close(20)

end program ENM_BD

subroutine get_CA(pdb_id, chain, x_coord, y_coord, z_coord, res_mass, res_name, res_chain, res_number)
	implicit none
    
    	character (len = 100), intent(in) :: pdb_id, chain
    	character (len = 1) :: alt_loc, chain_id
    	character (len = 3) :: atom_type
    	character (len = 4) :: label
    	character (len = 3) :: res_type
    	character (len = 80) :: string
    	integer :: num_res, i, n, atom_num, res_num, pdb_iostat, num_atoms, num_chains
    	double precision :: x, y, z
    	double precision, intent(out), allocatable :: x_coord(:), y_coord(:), z_coord(:), res_mass(:)
    	character (len = 3), intent(out), allocatable :: res_name(:)
    	character (len = 1), intent(out), allocatable :: res_chain(:)
    	integer, intent(out), allocatable :: res_number(:)
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
    
    	! format :ATOM     47  CA  VAL A   7      24.703 -13.450  -0.420  1.00 13.94           C
20  format(a4,i7,2x,a3,a1,a3,1x,a1,i4,4x,3f8.3)
	
	num_chains = len(trim(chain))
	num_res = 0
	do i = 1, num_chains
		open(16,file = "ATOM.pdb",status = 'old')
		do n = 1, num_atoms
			read(16,20) label, atom_num, atom_type, alt_loc, res_type, chain_id, res_num, x, y, z
			if (atom_type(1:2) == 'CA') then
				if (alt_loc == ' ' .or. alt_loc == 'A') then
					if (chain_id == chain(i:i)) then
						num_res = num_res + 1
					end if
				end if
			end if
		end do
		close(16)
	end do

    	allocate (x_coord(num_res))
    	allocate (y_coord(num_res))
    	allocate (z_coord(num_res))
    	allocate (res_mass(num_res))
    	allocate (res_name(num_res))
    	allocate (res_chain(num_res))
    	allocate (res_number(num_res))

	num_res = 0
	do i = 1, num_chains
		open(16,file = "ATOM.pdb",status = 'old')
		do n = 1, num_atoms
			read(16,20) label, atom_num, atom_type, alt_loc, res_type, chain_id, res_num, x, y, z
			if (atom_type(1:2) == 'CA') then
				if (alt_loc == ' ' .or. alt_loc == 'A') then
					if (chain_id == chain(i:i)) then
						num_res = num_res + 1
						x_coord(num_res) = x
        					y_coord(num_res) = y
        					z_coord(num_res) = z
        					res_name(num_res) = res_type
						res_chain(num_res) = chain_id
						res_number(num_res) = res_num
						if (res_type == 'ALA') then
					    		res_mass(num_res) = 71
						else if (res_type == 'ARG') then
					    		res_mass(num_res) = 156
						else if (res_type == 'ASN') then
					    		res_mass(num_res) = 114
						else if (res_type == 'ASP') then
					    		res_mass(num_res) = 115
						else if (res_type == 'CYS') then
					    		res_mass(num_res) = 103
						else if (res_type == 'GLU') then
					    		res_mass(num_res) = 129
					    	else if (res_type == 'GLH') then
					    		res_mass(num_res) = 129
						else if (res_type == 'GLN') then
					    		res_mass(num_res) = 128
						else if (res_type == 'GLY') then
					    		res_mass(num_res) = 57
						else if (res_type == 'HIS') then
					    		res_mass(num_res) = 137
					    	else if (res_type == 'HIE') then
					    		res_mass(num_res) = 137
					    	else if (res_type == 'HID') then
					    		res_mass(num_res) = 137
					    	else if (res_type == 'HIP') then
					    		res_mass(num_res) = 137
						else if (res_type == 'ILE') then
					    		res_mass(num_res) = 113
						else if (res_type == 'LEU') then
					    		res_mass(num_res) = 113
						else if (res_type == 'LYS') then
					    		res_mass(num_res) = 128
					    	else if (res_type == 'LYN') then
					    		res_mass(num_res) = 128
						else if (res_type == 'MET') then
					    		res_mass(num_res) = 131
						else if (res_type == 'PHE') then
					    		res_mass(num_res) = 147
						else if (res_type == 'PRO') then
					    		res_mass(num_res) = 97
						else if (res_type == 'SER') then
					    		res_mass(num_res) = 87
						else if (res_type == 'THR') then
					    		res_mass(num_res) = 101
						else if (res_type == 'TRP') then
					    		res_mass(num_res) = 186
						else if (res_type == 'TYR') then
					    		res_mass(num_res) = 163
						else if (res_type == 'VAL') then
					    		res_mass(num_res) = 99
						else
					    		res_mass(num_res) = 110
						end if
					end if
				end if
			end if
		end do
		close(16)
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
        	call gauss(0.0d0,sd,v(i)) ! v(i) m/s
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

subroutine get_force_res_list(step,n_CA,x,y,z,Slim,cutoff,chain_id_ref,num_force_list,list,num_step_print,res_numb_ref)
    	implicit none
    
    	integer, intent(in) :: n_CA, step
    	double precision, intent(in) :: x(n_CA), y(n_CA), z(n_CA)
    	integer, intent(in) :: cutoff
    	character (len = 1), intent(in) :: chain_id_ref(:)
    	integer :: i, j, seq_dist
    	integer, intent(out) :: num_force_list
    	integer, intent(in) :: Slim 
    	integer, allocatable, intent(out) :: list(:,:)
    	double precision :: dist
    	integer, intent(in) :: num_step_print
    	integer, intent(in) :: res_numb_ref(n_CA)

    
   	num_force_list = 0
    	!$omp parallel do private(i,j,dist,seq_dist) reduction(+:num_force_list)
    	do i = 1, n_CA - 1
        	do j = i + 1, n_CA
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
    	!$omp end parallel do

    	allocate(list(num_force_list,2))

    	num_force_list = 0
    	do i = 1, n_CA - 1
        	do j = i + 1, n_CA
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
    
end subroutine get_force_res_list
