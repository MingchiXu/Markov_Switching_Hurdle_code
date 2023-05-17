#!/bin/bash
#SBATCH --array=280-380
#SBATCH --account=def-aschmidt  # replace this with your supervisors account
#SBATCH --ntasks=1              # number of processes
#SBATCH --mem-per-cpu=20000M      # memory; default unit is megabytes
#SBATCH --time=110:00:00         # time (HH:MM:SS)
#SBATCH --output=/home/mingchi/logs/rps_result/%x-%j.out

module load gcc/9.3.0 r/4.0.2

# Export the nodes names.
# If all processes are allocated on the same node, NODESLIST contains : node1 node1 node1 node1
# Cut the domain name and keep only the node name
export NODESLIST=$(echo $(srun hostname | cut -f 1 -d '.'))
R -f rps-NB-cluster.R --args $SLURM_ARRAY_TASK_ID
