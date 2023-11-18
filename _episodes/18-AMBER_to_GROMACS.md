---
title: "Hands-on 2: Transferring AMBER simulation to GROMACS"
teaching: 20
exercises: 5
questions:
- "How to use GROMACS to continue a simulation equilibrated with AMBER?"
objectives:
- "Restart a running AMBER simulation with GROMACS"
keypoints:
- " "
---

## Moving simulations from AMBER to GROMACS
A specific MD software may provide you with a method that you cannot find in other software. What to do if you already have prepared and equilibrated a simulation? MD packages are reasonably compatible with each other, which makes transferring simulations possible. During this hands-on activity, you will transfer the equilibrated bilayer simulation you prepared in hands-on 1 to GROMACS.

Typically GROMACS simulations are restarted from checkpoint files in cpt format. By default, GROMACS writes a checkpoint file of your system every 15 minutes.AMBER restart files cannot be used to create a GROMACS checkpoint file, so we will prepare a GROMACS run input file (tpr).

- Change directory to 4-equilibration.
- Create a directory for GROMACS simulation
- Load ambertools and gromacs modules and start python

~~~
mkdir 5-amb2gro && cd 5-amb2gro
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 ambertools/23 gromacs/2023.2
python
~~~
{: .language-bash}


### Converting AMBER restart and topology to GROMACS.
Using parmed, convert AMBER topology and restart files to GROMACS.
~~~
import parmed as pmd
amber = pmd.load_file("../bilayer.parm7", "../4-equilibration/eq-2.rst7")
ncrest=pmd.amber.Rst7("../4-equilibration/eq-2.rst7")
amber.velocities=ncrest.vels
amber.save("restart.gro")
amber.save("bilayer.top")
quit()
~~~
{: .language-python}

These parmed commands will create two GROMACS files: restart.gro (coordinates and velocities) and bilayer.top.

### Creating and editing index file

The default index file can be created using the command:

~~~
gmx make_ndx -f restart.gro
~~~
{: .language-bash}

It has 13 groups:

~~~
...
  0 System              : 113855 atoms
  1 Other               : 47155 atoms
  2 PA                  : 16882 atoms
  3 SPM                 :  3441 atoms
  4 SA                  :  4092 atoms
  5 PC                  :  3762 atoms
  6 OL                  : 13700 atoms
  7 PE                  :  3422 atoms
  8 PS                  :  1767 atoms
  9 Na+                 :    73 atoms
 10 Cl-                 :    16 atoms
 11 Water               : 66700 atoms
 12 SOL                 : 66700 atoms
 13 non-Water           : 47155 atoms
...
~~~
{: .output}

Groups SOL and non-Water can be deleted. There will be two thermostats: one thermostat will be used for lipids, and another one for everything else. For these thermostats we have to create two new groups: Lipids and Water+Ions.

Here are the commands to do it:
~~~
gmx make_ndx -f restart.gro << EOF
del 12 
del 12
1&!rNa+&!rCl- 
name 12 Lipids
11|rNa+|rCl- 
name 13 WaterIons   
q
EOF
~~~
{: .language-bash}

Check created index file:
~~~
gmx make_ndx -f restart.gro -n index.ndx
~~~
{: .language-bash}

You should have two new groups:  
~~~
12 Lipids              : 53850 atoms  
13 WaterIons           : 66792 atoms 
~~~
{: .output}

As packing bilayer and solvating it is non-deterministic, the number of atoms may be slightly different in your case. 

### Creating a portable binary run input file. 

GROMACS run input files (.tpr) contain everything needed to run a simulation: topology, coordinates, velocities and MD parameters. 

~~~
gmx grompp -p bilayer.top -c restart.gro -f production.mdp -n index.ndx -maxwarn 1
~~~
{: .language-bash}


GROMACS simulation input file production.mdp
~~~
Title                   = bilayer
; Run parameters
integrator              = md
nsteps                  = 400000
dt                      = 0.001
; Output control
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstenergy               = 100
nstlog                  = 10000
nstxout-compressed      = 5000
compressed-x-grps       = System
; Bond parameters
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
; Neighborsearching
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = Ener ; anaytic VDW correction
; Electrostatics
coulombtype             = PME
pme_order               = 4
; Temperature coupling is on
tcoupl                  = V-rescale
tc-grps                 = WaterIons Lipids
tau_t                   = 0.1 0.1
ref_t                   = 300 300
; Pressure coupling is on
pcoupl                  = Parrinello-Rahman
pcoupltype              = semiisotropic
tau_p                   = 5.0
ref_p                   = 1.0 1.0
compressibility         = 4.5e-5 4.5e-5
; Periodic boundary conditions
pbc                     = xyz
; Velocity generation
gen_vel                 = no
~~~

### Running GROMACS simulation
~~~
#!/bin/bash
#SBATCH  --mem-per-cpu 4000 --time 1:0:0   
#SBATCH  -c10 --gpus-per-node=v100:1
  
module load StdEnv/2020 gcc/9.3.0 cuda/11.4 openmpi/4.0.3 gromacs/2023.2

srun gmx mdrun -ntomp ${SLURM_CPUS_PER_TASK:-1} \
-nb gpu -pme gpu -update gpu -bonded cpu -s topol.tpr
~~~
{: .language-bash}


