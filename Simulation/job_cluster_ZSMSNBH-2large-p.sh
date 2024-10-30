#!/bin/bash
#SBATCH --array=40-80
#SBATCH --account=def-aschmidt  # replace this with your supervisors account
#SBATCH --ntasks=1              # number of processes
#SBATCH --mem-per-cpu=8500M      # memory; default unit is megabytes
#SBATCH --time=08:00:00         # time (HH:MM:SS)
#SBATCH --output=/home/mingchi/logs/rps_result/%x-%j.out

module load gcc/9.3.0 r/4.0.2

# Export the nodes names.
# If all processes are allocated on the same node, NODESLIST contains : node1 node1 node1 node1
# Cut the domain name and keep only the node name
export NODESLIST=$(echo $(srun hostname | cut -f 1 -d '.'))
R -f rps_sim_ZSMSNBH_2large_p.R --args $SLURM_ARRAY_TASK_ID
