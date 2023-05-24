---
title: "Running Molecular Dynamics Simulations with AMBER"
teaching: 30
exercises: 5
questions:
- "What simulation programs are available in the AMBER package?"
- "How to minimize energy?"
- "How to heat up a simulation system?"
- "How to equilibrate a simulation system?"
objectives:
- "Learn how to minimize energy of a system."
- "Learn how to heat up and equilibrate a simulation."
keypoints:
- " "
---

### AMBER MD engines.
Amber package includes two MD engines: SANDER and PMEMD. Both programs are available in serial and parallel versions. 


#### SANDER

{: .instructor_notes} 
SANDER is a free simulation engine distributed with the AmberTools package. For parallel distributed simulations, it uses the MPI (message passing interface). The parallel version of Sander implements replicated data structure. Each CPU computes a portion of the potential energy and corresponding gradients for a set of atoms assigned to it.  A global part of the code then sums the force vector and sends the result to each CPU. The processors then perform a molecular dynamics update step for the assigned atoms and communicate the updated positions to all CPUs in preparation for the subsequent molecular dynamics step.
{: .instructor_notes} 

{: .instructor_notes} 
This model provides a convenient programming environment, but the main problem is that the communication required at each step grows with the number of processors limiting parallel scaling. 
{: .instructor_notes} 

{: .self_study_text} 
- SANDER is a free simulation engine distributed with the AmberTools package.
{: .self_study_text} 


#### PMEMD

{: .instructor_notes} 
PMEMD is an extensively revised version of SANDER available only in the commercial AMBER package. Developers made many optimizations to improve both single-processor performance and parallel scaling. To avoid data transfer bottleneck, PMEMD communicates to each processor only the coordinate information necessary for computing the pieces of the potential energy assigned to it. This code, however, does not support all of the options found in the SANDER.  
{: .instructor_notes} 

{: .self_study_text} 
- PMEMD is an extensively revised version of SANDER available only in the commercial AMBER package.
{: .self_study_text} 


#### GPU-Accelerated PMEMD
GPU - accelerated PMEMD version of PMEMD (pmemd.cuda) leverages NVIDIA GPUs to perform MD simulations.  It is significantly faster than the CPU version achieving high simulation speed by executing all calculations on a single GPU within its memory. This approach eliminates the bottleneck of moving data between CPU and GPU and allows very efficient GPU usage.  

**<font color="red">Modern GPUs are so fast that communication overhead between GPUs does not allow for efficient parallel scaling of an MD simulation to two or more GPUs.</font>**

{: .instructor_notes} 
While you can run a single simulation on several GPUs using the parallel PMEMD GPU version (pmemd.cuda.MPI) it will run not run much faster than on a single GPU. Parallel GPU version is useful only for specific simulations such as thermodynamic integration and replica-exchange MD. These types of jobs involve several completely independent simulations that can be executed concurrently on different GPUs. PMEMD allows running multiple copies of simulations within a single parallel run via the multi-pmemd mechanism described below. 
{: .instructor_notes} 


[PMEMD parallel scaling, A100](https://mdbench.ace-net.ca/mdbench/bform/?software_contains=PMEMD.cuda.MPI&software_id=&module_contains=&module_version=&site_contains=Narval&gpu_model=&cpu_model=&arch=&dataset=6n4o)  
[PMEMD parallel scaling, P100](https://mdbench.ace-net.ca/mdbench/bform/?software_contains=PMEMD.cuda.MPI&software_id=&module_contains=&module_version=&site_contains=Cedar&gpu_model=P100-PCIE&cpu_model=&arch=&dataset=6n4o)

#### Multi-sander and multi-pmemd simulations
Multi-sander and multi-pmemd are wrappers around parallel versions of these programs. These wrappers are invoked by preparing special input files. The wrappers allow running multiple copies of simulations within a single parallel job. The multi-sander and multi-pmemd mechanisms are also utilized for methods requiring multiple simulations to communicate with one another, such as thermodynamic integration and replica exchange molecular dynamics.


<br>
#### Summary of available AMBER MD executables: 

|  Verion     | SANDER |  PMEMD |  
|:---|:---|:---|
|Serial | sander  | pmemd|      |
|Parallel | sander.MPI  | pmemd.MPI|   
|Single GPU, default link| -  | **pmemd.cuda** -> pmemd.cuda_SPFP|
|Single GPU, single precision| -  | pmemd.cuda_SPFP|   
|Single GPU, double precision| -  | pmemd.cuda_DPFP|
|Multi GPU, default link| -  | **pmemd.cuda.MPI** -> pmemd.cuda_SPFP.MPI|
|Multi GPU, single precision| -  | pmemd.cuda_SPFP.MPI|   
|Multi GPU, double precision| -  | pmemd.cuda_DPFP.MPI|


### Energy minimization.
Before simulating a system we need to relax it. Any atomic clashes must be resolved, and potential energy minimized to avoid unphysically large forces that can crash a simulation. 

{: .instructor_notes} 
The general minimization strategy is first to restrict all solute atoms with the experimental coordinates and relax all atoms that were added. (solvent, ions and missing fragments). This will help to stabilize the native conformation. There are no strict rules defining how many minimization steps are necessary. The choice will depend on the composition of a simulation system. For a big systems with a large amount of missing residues it is safer to carry out several minimization steps gradually releasing restraints. For example, you can first relax only solvent and ions, then lipid bilayer (if this is a membrane protein), then added fragments, then the original protein side-chains. Having more steps may be unnecessary, but it will not cause any problems. 
{: .instructor_notes} 

{: .self_study_text} 
- it is safer to carry out several minimization steps gradually releasing restraints.
{: .self_study_text} 

For example, we could do a two stage minimization. In the first stage we restrain all original atoms. In the second stage we restrain only the original backbone atoms. Our example protein is very small and we have limited time, so we skip the first step and restrain only protein backbone.

~~~
cd ~/workshop/pdb/1RGG/AMBER/1_minimization
~~~
{: .language-bash}

MD programs can do a lot of different things, and every type of calculation has a number of parameters that allow us to control what will be done. To run a minimization we need to make an input file describing exactly what we want to do and how we want to do it:

- instruct a simulation program to minimize energy  
- choose a method of minimization  
- specify the maximum number of cycles of minimization  
- apply restraints to a subset of atoms (optionally)

AMBER MD programs read simulation parameters from an input file. Simulation parameters in AMBER are called "FLAGS". The Table below lists some important minimization FLAGS.   

#### AMBER minimization parameters

| Flag        | Value     | Description
|-------------|-----------|-----------
|imin         |     1     | Turn on minimization
|ntmin        | 0,1,2,3   | Flag for the method of minimization    
|maxcyc       | integer   | The maximum number of cycles of minimization  
|ncyc         | integer   | If NTMIN=1 switch from SD to CG after NCYC cycles  
|ntpr         | integer n | Print energies every n steps    
|ntr          |    1      | Use cartesian harmonic restraints   
|restraint_wt | float     | Restraint force kcal/mol   
|restraintmask| ambermask | Specifies restrained atoms   

#### Methods of minimization

|--|
|0 |Steepest descent+conjugate gradient. The first 4 cycles are steepest descent at the start of the run and after every non-bonded pair-list update.
|1 | For NCYC cycles the steepest descent method is used then conjugate gradient is switched on.
|2 | Steepest descent only
|3 | XMIN family methods. The default is LBFGS (Limited-memory Broyden-Fletcher-Goldfarb-Shanno). It is a popular algorithm in machine learning. The method incrementally learns from previous steps, so that it can make the next step more accurate. It converges considerably faster than CG, but requires more memory.

Minimization input file *min.in*:
~~~
Energy minimization 
&cntrl 
imin=1, ntmin=0, maxcyc=1000,   ! Minimization, method, number of cycles 
ntpr=5,                         ! Print energies every ntpr steps
ntr=1,                          ! Use harmonic cartesian restraints   
restraint_wt=10.0,              ! Restraint force kcal/mol
restraintmask="(:1-96)&(@CA,N,O)",
&end
END
~~~
{: .file-content}

- There are several molecular dynamics programs in AMBER package: *sander, sander.MPI, pmemd, pmemd.MPI, pmemd.cuda*, and *pmemd.cuda.MPI*.  
- *Sander* is distributed free, *pmemd* is high performance binary which requires a license.

Submission script *submit.sh*:
~~~
#!/bin/bash
#SBATCH --mem-per-cpu=1000 --time=3:0:0 --ntasks=4
ml --force purge
ml StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 amber/20.9-20.15
srun pmemd.MPI -O -i min.in -p prmtop -c inpcrd -ref inpcrd -r minimized.nc -o mdout
~~~
{: .file-content}

This job runs about 30 sec on 4 CPUs.

- The option -O means: overwrite the output files if present.  
- The output from the minimization goes into the file *mdout*. The total energy of the system is printed in the lines beginning with "EAMBER =". 

If minimization is successful we expect to see large negative energies.

### Heating
~~~
cd ~/workshop/pdb/1RGG/AMBER/2_heating
~~~
{: .language-bash}

#### Molecular dynamics parameters

Flag        | Value  | Description
-------------|--------|-----------
dt           | 0.001  | Time step, ps. Default 0.001 
ntf          | 2      | Force evaluation. Omit interactions involving bonds with H-atoms. Default 1 (complete interaction)
ntc          | 2      | Flag for SHAKE. Bonds involving hydrogen are constrained.
ntt          | 1      | Constant temperature, using the Berendsen weak-coupling algorithm.
tempi        | 150    | Initial temperature. The velocities are assigned from a Maxwellian distribution at TEMPI 
temp0        | 300    | Reference temperature at which the system is to be kept
tautp        | 1      | Time constant, in ps, for heat bath coupling, default is 1 ps. 
ntp          | 1      | Flag for constant pressure dynamics. 1 - MD with isotropic position scaling
barostat     | 1      | Berendsen (default)
pres0        | 1      | Reference pressure, default 1 bar
taup         | 4      | Pressure relaxation time (in ps), default 1 
ntwx         | 1000   | Every ntwx steps, the coordinates will be written to the mdcrd file
ntpr         | 100    | Print energies in the log every 100 steps, default 50 

In the examples bonds with hydrogens are not constrained and the default timestep of 1 fs in used. To turn on SHAKE use ntf=2 and ntc=2.

#### Heating input file:
~~~
Heating 
&cntrl 
ntt=1,                               ! Use Berendsen thermostat
tempi=150,temp0=300,tautp=1,         ! Initial and reference temperature, time constant
ntp=0,                               ! No barostat
ntpr=100,                            ! Print energies every ntpr steps
ntwx=1000,                           ! Write coordinates every ntws steps
nstlim=10000,                        ! Simulate nstlim steps
ntr=1,                               ! Use harmonic cartesian restraints 
restraint_wt=10,                     ! Restraint force kcal/mol
restraintmask="(:1-96)&(@CA,N,O)",
&end
END
~~~
{: .file-content}

#### Heating submission script:
~~~
#!/bin/bash
#SBATCH --mem-per-cpu=4000 --time=3:0:0 --ntasks=4
ml --force purge
ml StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 amber/20.9-20.15
srun pmemd.MPI -O -i heat.in -p prmtop -c minimized.nc -ref inpcrd -r heated.nc -o mdout
~~~
{: .language-bash}

This job runs about 2 min on 4 CPUs.

### Equilibration
#### Constrained equilibration
~~~
cd ~/workshop/pdb/1RGG/AMBER/3_equilibration
~~~
{: .language-bash}

- Turn on restart flag.

Flag         | Value  | Description
-------------|--------|-----------
ntx          | 5      | Coordinates and velocities will be read from a restart file
irest        | 1      | Restart simulations


Input file *equilibrate_1.in*:
~~~
&cntrl 
irest=1,ntx=5,                      ! Restart using coordinates and velocities
ntt=1,temp0=300,tautp=1,            ! Use Berendsen thermostat
ntp=1,barostat=1,pres0=1,taup=1,    ! Use Berendsen barostat
ntpr=100,                           ! Print energies every ntpr steps   
ntwx=1000,                          ! Save coordinates every ntwx steps
nstlim=5000,                        ! Simulate nstlim steps
ntr=1,                              ! Turn on restraints
restraint_wt=10,                    ! Restraint force, kcal/mol
restraintmask="(:1-96)&(@CA,N,O)",
&end
END
~~~
{: .file-content}

Submission script for CPU-only job *submit_1.sh*:
~~~
#!/bin/bash
#SBATCH --mem-per-cpu=4000 --time=3:0:0 --ntasks=4
module --force purge
ml StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 amber/20.9-20.15
srun pmemd.MPI -O -i equilibrate_1.in -p prmtop -c heated.nc -ref inpcrd -r equilibrated_1.nc -o equilibration_1.log
~~~
{: .language-bash}

Submission script for GPU job *submit_1.sh*:
~~~
#!/bin/bash
#SBATCH --mem-per-cpu=4000 --time=3:0:0 --ntasks=1 --gres=gpu:v100:1 --partition=all_gpus
ml --force purge
ml StdEnv/2020 gcc/9.3.0 cuda/11.0 openmpi/4.0.3 amber/20.9-20.15
pmemd.cuda -O -i equilibrate_1.in -p prmtop -c heated.nc -ref inpcrd -r equilibrated_1.nc -o equilibration_1.log
~~~
{: .language-bash}

This job runs about 3 min on 4 CPUs.

#### Unconstrained equilibration

- Use Langevin thermostat + Berendsen barostat.
- Simulate on GPU

Input file *equilibrate_2.in*:
~~~
Equilibration, Langevin thermostat + Berendsen barostat
&cntrl 
irest=1,ntx=5,                     ! Restart using coordinates and velocities
ntt=3,temp0=300,gamma_ln=1.0,      ! Langevin thermostat, collision frequency
ntp=1,barostat=1,pres0=1,taup=1.0, ! Berendsen barostat
ntpr=100,                          ! Print energies every ntpr steps   
ntwx=1000,                         ! Save coordinates every ntwx steps
nstlim=10000000,                   ! Simulate nstlim steps
&end
END
~~~
{: .file-content}

Submission script *submit_2.sh*
~~~
#!/bin/bash
#SBATCH --mem-per-cpu=4000 --time=3:0:0 --ntasks=1 --gres=gpu:v100:1 --partition=all_gpus
ml --force purge
ml StdEnv/2020 gcc/9.3.0 cuda/11.0 openmpi/4.0.3 amber/20.9-20.15

pmemd.cuda -O  -i equilibrate_2.in -p prmtop -c equilibrated_1.nc  -r equilibrated_2.nc -o equilibration_2.log
~~~
{: .language-bash}

### Simulation with NAMD
Because NAMD natively supports AMBER topology files, simulating a system prepared with AMBER tools requires only NAMD simulation input files and NAMD - compatible coordinate files such as pdb or NAMD binary coordinates.

In the worksop data, you will find example simulation input files for minimization, heating and equilibration:
~~~
ls ~/workshop/pdb/6N4O/simulation/sim_namd
1-minimization  2-heating  3-equilibration  4-production
~~~
{: .language-bash}

### Analyzing simulation logs
#### Extract selected energy components from MD log and save in a table using *cpptraj*.

Use the script `~/bin/extract_energies.sh`:
~~~
mkdir ~/bin
nano ~/bin/extract_energies.sh
~~~
{: .language-bash}

The contents of the script *~/bin/extract_energies.sh*:
~~~
#!/bin/bash
echo "Usage: extract_energies simulation_log_file" 

log=$1
cpptraj << EOF
readdata $log
writedata energy.dat $log[Etot] $log[TEMP] $log[PRESS] $log[VOLUME] time 0.1
EOF
~~~
{: .file-content}

Extract selected energy components. 
~~~
cd ~/workshop/pdb/1RGG/AMBER/3_equilibration/
ml purge
ml StdEnv/2020 gcc/9.3.0 cuda/11.4 ambertools/22
bash extract_energies.sh equilibration_2.log
~~~
{: .language-bash}

- Modify the script to add or remove energy components.
- Plot data table with any plotting program.

#### Plot energy components with python

Get interactive allocation and start vncserver. Connect VNC viewer to the node running vncserver.

On the compute node load python and scipy-stack modules:
~~~
cd ~/workshop/pdb/1RGG/AMBER/3_equilibration
ml StdEnv/2020 python scipy-stack
python
~~~
{: .language-bash}

Read table from the file `energies.dat` into pandas dataframe and plot it:

~~~
python ~/bin/plot_energies.py
~~~
{: .language-bash}

File `plot_energies.py`:

~~~
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_table('energy.dat', delim_whitespace=True)
df.columns=["Time", "Etot", "Temp", "Press", "Volume"]

df.plot(subplots=True, x="Time", figsize=(6, 8))
plt.legend(loc='best')
plt.savefig('energies.png')
# plt.show()
~~~
{: .language-python}

### Managing trajectories
You can remove from a trajectory all components that are not essential for analysis, for example water and ions. The following command will remove everything except residues from 1 to 96 and save every second frame.
~~~
cpptraj<<EOF
parm prmtop
trajin mdcrd.nc 1 last 2
strip !(:1-96)
trajout mdcrd_nowat.nc 
run
EOF
cpptraj<<EOF
parm prmtop
parmstrip !(:1-96)
parmwrite out prmtop_nowat.parm7
run
EOF
~~~
{: .language-bash}

## Transferring equilibrated system between simulation packages.
Simulation packages have different methods and performance. It is useful to be able to transfer a running simulation from one software to another. 

{: .instructor_notes} 
Imagine that you started your project with GROMACS, but later realized that you need to run a constant pH simulation. You need to switch to AMBER. Want to study conformational transitions? Gaussian accelerated MD is not available in GROMACS. Another reason to move to AMBER/NAMD. Want to apply custom forces? It is easy to with Tcl scripting in NAMD.
{: .instructor_notes} 

### Moving simulation from AMBER to GROMACS.
To transfer simulation to GROMACS we need to convert topology and restart files.

~~~
cd ~/workshop/pdb/1RGG/AMBER_to_GROMACS
~~~
{: .language-bash}

#### Convert AMBER topology to GROMACS
~~~
module --force purge
module load StdEnv/2020 gcc/9.3.0 cuda/11.0 openmpi/4.0.3 
module load amber/20.9-20.15 gromacs/2021.2 scipy-stack python
~~~
{: .language-bash}
~~~
import parmed as pmd
amber = pmd.load_file("prmtop.parm7", "inpcrd.rst7")
amber.save('topol.top')
amber.save('inpcrd.gro')
~~~
{: .language-python}

#### Make index file
The default index groups are OK:
~~~
gmx make_ndx -f inpcrd.gro <<EOF
q
EOF
~~~
{: .language-bash}

#### Generate positional restraints file for the protein backbone.
~~~
gmx genrestr -f inpcrd.gro -fc 500.0 -n index.ndx -o backbone.itp<<EOF
Backbone
EOF
~~~
{: .language-bash}

Add definitions of the position restraints to the topology "gromacs.top". Use a text editor of your choice to insert the following lines at the end of the syste1 molecule block:
~~~
#ifdef POSRES
#include "backbone.itp"
#endif

[ moleculetype ]
; Name            nrexcl
Na+          3
~~~
{: .file-content}

Define position restraints in the input file min.mdp:

~~~
; Turn on position restraints
define = -D_POSRES
~~~
{: .file-content}

#### Convert AMBER restart to GROMACS restart.

~~~
import parmed as pmd

amber = pmd.load_file("prmtop.parm7", "rest.rst7")
ncrest=pmd.amber.Rst7("rest.rst7")
amber.velocities=ncrest.vels
amber.save("restart.gro")
~~~
{: .language-python}

#### Convert GROMACS restart to portable trajectory and make binary topology
~~~
gmx trjconv -f restart.gro -o restart.trr
gmx grompp -p topol.top  -c restart.gro -t restart.trr -f gromacs_production.mdp
~~~
- Tested with gromacs/2021 and gromacs/2022.- This procedure does not work with gromacs/2023. 
{: .language-bash}
