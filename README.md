# ENM_BD
Brownian Dynamics (BD) with Elastic Network Models (ENMs) for sampling of protein conformational dynamics

You can use the ENM_BD.f90 (source code in Fortran) to run BD and sample new protein conformations from an input PDB structure.

First compile the code with: gfortran ENM_BD.f90 -o ENM.f90 -fopenmp (by default, we use 16 OpenMP threads. You can change this by modifying the "num_threads" parameter in the source code).

Then, simply run the code providing 4 inputs from the terminal with: ./ENM_BD $pdb_file $chain_IDs $tot_number_steps $save_frequency

Input 1 = name of the PDB file (with .pdb extension and standard PDB formatting. Note that only CA atoms with empty or "A" altLoc will be read); input 2 = chain identifiers (e.g., "A", "AB", "ABCD", etc.); input 3 = total number of steps of the simulation (note that dt = 2 fs, so if you want to run 1 ns, you have to type 500,000 steps); input 4 = number of steps after which the code saves conformation coordinates (if you want to save e.g. every ps, type 500). The folder RBP_example/BD_1ns_1000conformers reports BD conformations for RBP obtained from 1ba2.pdb using: ./ENM_BD 1ba2.pdb A 500000 500


