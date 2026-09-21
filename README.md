# ENM_BD
Brownian Dynamics (BD) with Elastic Network Models (ENMs) for sampling of protein conformational dynamics

You can use the ENM_BD.f90 (source code in Fortran) to run BD and sample new protein conformations from an input PDB structure.
First compile the code with: gfortran ENM_BD.f90 -o ENM.f90 -fopenmp (by default, we use 16 OpenMP threads. You can change this number by modifying the "num_threads" parameter in the source code).
Then, simply run providing 4 inputs from the terminal with: ./ENM_BD $pdb_file_name $chain_IDs $tot_number_of_steps $save_frequency . Input 1 = name of the PDB file (with .pdb extension and formatting); input 2 = chain identifiers (e.g., "A", "AB", "ABCD", etc.); input 3 = total length of the BD simulations as number of steps (note that dt = 2 fs, so if you e.g. want to run 1 ns, you have to type 500000 steps); input 4 = number of steps after which you will write an output conformation (if you want to save every 10 ps, type 5000).
