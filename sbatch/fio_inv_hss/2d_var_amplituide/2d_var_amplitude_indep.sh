#!/bin/bash

#SBATCH --job-name=FIO2d_indep
#SBATCH --output=FIO2d_indep_%j.out
#SBATCH --error=FIO2d_indep_%j.err
#SBATCH --nodelist=bigMem0
#SBATCH --exclusive
#SBATCH --time=36:00:00
#SBATCH --qos=double

module unload MATLAB
module load MATLAB/R2023b

echo "=========================================="
echo "        SLURM JOB INFORMATION"
echo "=========================================="
echo "Job ID:           $SLURM_JOB_ID"
echo "Job Name:         $SLURM_JOB_NAME"
echo "Submit Directory: $SLURM_SUBMIT_DIR"
echo "Submit Host:      $SLURM_SUBMIT_HOST"
echo "Partition:        $SLURM_JOB_PARTITION"
echo "Node List:        $SLURM_JOB_NODELIST"
echo "Nodes:            $SLURM_NNODES"
echo "Tasks:            $SLURM_NTASKS"
echo "CPUs per Task:    $SLURM_CPUS_PER_TASK"
echo "Total CPUs:       $SLURM_NPROCS"
echo "Working Directory: $(pwd)"
echo "=========================================="
echo ""
echo "        RESOURCE ALLOCATION"
echo "=========================================="
echo "Memory per Node:  $SLURM_MEM_PER_NODE"
echo "Time Limit:       $SLURM_TIMELIMIT"
echo "=========================================="
echo ""
echo "        ENVIRONMENT VARIABLES"
echo "=========================================="
echo "PATH:             $PATH"
echo "LD_LIBRARY_PATH:  $LD_LIBRARY_PATH"
echo "LIBRARY_PATH:     $LIBRARY_PATH"
echo "CPATH:            $CPATH"
echo "=========================================="
echo ""
echo "        MODULE INFORMATION"
echo "=========================================="
module list 2>&1
echo "=========================================="
echo ""
echo "Job started at: $(date)"
echo "=========================================="
echo ""

matlab -r 'cd /home/jyliu/FIO; fio_startup; cd experiments/fio_inv_hss/2d_var_amplitude; exp_2d_var_amplitude_indep;'

echo ""
echo "MATLAB finished at: $(date)"
echo "=========================================="
echo "Job completed"
echo "=========================================="

