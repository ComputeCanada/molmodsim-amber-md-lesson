---
title: "Hands-on 6: A Comparative Performance Ranking of the Molecular Dynamics Software"
teaching: 30
exercises: 5
questions:
- "How to evaluate CPU efficiency of a simulation?"
- "How fast and efficient my simulation will run with different programs and computing resources?"
objectives:
- "Learn how to request the right computational resources"
keypoints:
- "To assess CPU efficiency, you need to know how fast a serial simulation runs"
---

Requesting the right computational resources is essential for fast and efficient simulations. Submitting a simulation with more CPUs does not necessarily mean that it will run faster. In some cases, a simulation will run slower with more CPUs. There is also a choice between using CPU or GPU versions. When deciding on the number of CPUs, it is crucial to consider both simulation speed and CPU efficiency. If CPU efficiency is low, you will be wasting resources. This will negatively impact your priority, and as a result, you will not be able to run as many jobs as you would if you used CPUs more efficiently. To assess CPU efficiency, you need to know how fast a serial simulation runs and then compare the expected 100% efficient speedup (speed on 1CPU x N) with the actual speedup on N CPUs.

Here is the [chart of the maximum simulation speed](https://mdbench.ace-net.ca/mdbench/) of all MD engines tested on the Alliance systems. These results may give you valuable insight into how fast and efficient you can expect your simulation to run with different packages/resources.

## Submission scripts for running the benchmarks.
### GROMACS
Extend simulation for 10000 steps
~~~
gmx convert-tpr -s topol.tpr -nsteps 10000 -o next.tpr
~~~
{: .language-bash}

#### Submission script for a CPU simulation
~~~
#SBATCH --mem-per-cpu=4000 --time=10:0:0 -c4 --ntasks=2

module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 gromacs/2023.2
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

srun gmx mdrun -s next.tpr -cpi state.cpt
~~~
{: .file-content}

[Benchmark](https://mdbench.ace-net.ca/mdbench/bform/?software_contains=GROMACS.mpi&software_id=&module_contains=&module_version=&site_contains=Narval&gpu_model=&cpu_model=&arch=&dataset=6n4o)

#### Submission script for a single GPU simulation 
~~~
#!/bin/bash
#SBATCH  --mem-per-cpu 2000 --time 1:0:0   
#SBATCH -c 12 --gpus-per-node 1  

module load StdEnv/2020 gcc/9.3.0 cuda/11.4 openmpi/4.0.3 gromacs/2023.2

gmx mdrun -ntomp ${SLURM_CPUS_PER_TASK:-1} \
-nb gpu -pme gpu -update gpu -bonded cpu -s topol.tpr
~~~
{: .file-content}

#### Submission script for a multiple GPU simulation 
~~~
#!/bin/bash
#SBATCH  --mem-per-cpu 2000 --time 1:0:0   
#SBATCH  --nodes 1 --ntasks 2 -c 12 --gpus-per-task 1  

module load StdEnv/2020 gcc/9.3.0 cuda/11.4 openmpi/4.0.3 gromacs/2023.2

srun gmx mdrun -ntomp ${SLURM_CPUS_PER_TASK:-1} \
-nb gpu -pme gpu -update gpu -bonded cpu -s topol.tpr
~~~
{: .file-content}

[Benchmark](https://mdbench.ace-net.ca/mdbench/bform/?software_contains=GROMACS.cuda.mpi&software_id=&module_contains=&module_version=&site_contains=Narval&gpu_model=&cpu_model=&arch=&dataset=6n4o)

### PMEMD
#### Submission script for a single GPU simulation
~~~
#!/bin/bash
#SBATCH -c1 --gpus 1 --mem-per-cpu=2000  --time=1:0:0

module --force purge
ml StdEnv/2020  gcc/9.3.0 cuda/11.4 openmpi/4.0.3 amber/20.12-20.15
pmemd.cuda -O -i pmemd.in -o production.log -p prmtop.parm7 -c restart.rst7

~~~
{: .file-content}

#### Submission script for a multiple GPU simulation
Multiple GPU pmemd version is meant to be used only for AMBER methods running multiple simulations, such as replica exchange. A single simulation does not scale beyond 1 GPU.

~~~
#!/bin/bash
#SBATCH --nodes=1 --ntasks=2 --gpus-per-node=2
#SBATCH --mem-per-cpu=2000 --time=1:0:0

module --force purge
ml StdEnv/2020  gcc/9.3.0 cuda/11.4 openmpi/4.0.3 amber/20.12-20.15

srun pmemd.cuda.MPI -O -i pmemd_prod.in -o production.log \
-p prmtop.parm7 -c restart.rst7
~~~
{: .file-content}

[Benchmark](https://mdbench.ace-net.ca/mdbench/bform/?software_contains=PMEMD.cuda.MPI&software_id=&module_contains=&module_version=&site_contains=Narval&gpu_model=&cpu_model=&arch=&dataset=6n4o)

### NAMD 3
#### Submission script for a GPU simulation
~~~
#!/bin/bash
#SBATCH -c2 --gpus-per-node=a100:2  
#SBATCH --mem-per-cpu=2000 --time=1:0:0
NAMDHOME=$HOME/NAMD_3.0b3_Linux-x86_64-multicore-CUDA

$NAMDHOME/namd3 +p${SLURM_CPUS_PER_TASK} +idlepoll namd3_input.in  
~~~
[Benchmark](https://mdbench.ace-net.ca/mdbench/bform/?software_contains=&software_id=36&module_contains=&module_version=&site_contains=Narval&gpu_model=&cpu_model=&arch=&dataset=6n4o)

### How to make your simulation run faster?
It is possible to increase time step to 4 fs with hydrogen mass re-partitioning. The idea is that hydrogen masses are increased and at the same time masses of the atoms to which these hydrogens are bonded are decreased to keep the total mass constant. Hydrogen masses can be automatically re-partitioned with the *parmed* program.
~~~
module --force purge
module load StdEnv/2020 gcc ambertools python scipy-stack
source $EBROOTAMBERTOOLS/amber.sh
parmed prmtop.parm7
~~~
{: .language-bash}

~~~
ParmEd: a Parameter file Editor

Loaded Amber topology file prmtop.parm7
Reading input from STDIN...
~~~
{: .output}

~~~
> hmassrepartition
> outparm prmtop_hmass.parm7
> quit
~~~
{: .language-bash}

References:  
1.[Lessons learned from comparing molecular dynamics engines on the SAMPL5 dataset](https://link.springer.com/article/10.1007%2Fs10822-016-9977-1)   
2.[Delivering up to 9X the Throughput with NAMD v3 and NVIDIA A100 GPU](https://developer.nvidia.com/blog/  delivering-up-to-9x-throughput-with-namd-v3-and-a100-gpu/)  
3.[AMBER GPU Docs](https://ambermd.org/GPUHowTo.php)  
4.[Long-Time-Step Molecular Dynamics through Hydrogen Mass Repartitioning](https://pubs.acs.org/doi/abs/10.1021/ct5010406)
