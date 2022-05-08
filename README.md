# LGS-PET

This code documents the image pre-processing steps relevant to the following paper:

_Frontoparietal 18FDG-PET hypo-metabolism in Lennox- Gastaut syndrome: further evidence for a shared epileptic network._

By: Tom Balfroid, Aaron E.L Warren, Linda J. Dalic, Alec Aeby, Salvatore U. Berlangieri, John S. Archer

Analyses were performed using the High Performance Computing (HPC) system ([Spartan](https://dashboard.hpc.unimelb.edu.au)) operated by Research Computing Services at The University of Melbourne:

Lev Lafayette, Greg Sauter, Linh Vu, Bernard Meade, "Spartan Performance and Flexibility: An HPC-Cloud Chimera", OpenStack Summit, Barcelona, October 27, 2016. doi.org/10.4225/49/58ead90dceaaa

This bash script creates SLURM submission scripts for each subject. After per-subject scripts are created, they are ready to submit via sbatch.

_Usage: bash LGS_PET_PP.sh_
