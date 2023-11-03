---
title: "How to pack a complex mixture of different lipid species into a bilayer and get it moving?"
teaching: 30
exercises: 5
questions:
- "How to pack a complex mixture of different lipid species into a bilayer?"
- "How to minimize energy?"
- "How to heat up a simulation system?"
- "How to equilibrate a simulation system?"
objectives:
- "Prepare a lipid bilayer"
- "Energy minimize, heat up and equilibrate the system"
keypoints:
- " "
---

## Introduction

This in-depth tutorial guides you through the process of building a bilayer, preparing simulation input files, minimizing energy, equilibrating the system, and running an equilibrium molecular dynamics simulation. You will need to follow a number of steps to complete the tutorial. If you are already familiar with some of these topics, you can skip them and focus on the ones you don't know about.

### Step 1: Packing lipids into a bilayer.
~~~
#!/bin/bash
#SBATCH -c1 --mem-per-cpu=4000 --time=6:0:0

module purge
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 ambertools/23

rm -f bilayer* *.log
packmol-memgen \
    --lipids CHL1:POPC:POPE:PSM//CHL1:POPC:POPE:PSM \
    --salt --salt_c Na+ --salt_a Cl- --saltcon 0.15 \
    --dist_wat 15 \
    --ratio 4:3:1:2//4:3:1:2 \
    --distxy_fix 100 \
    --parametrize
~~~
{: .language-bash}


This command generates several output files: AMBER topology (bilayer_only_lipids.top), coordinates (bilayer_only_lipids.crd), and PDB formatted system (bilayer_only_lipids.pdb). 

### Step 2: Changing force fields.
The file bilayer_only_lipids.pdb can be used to create topology files with new AMBER force fields not yet available in packmol-memgen. There is only one additional thing that we need to know before we can proceed: the periodic box size. It can be extracted from the file bilayer_only_lipid.top:

~~~
grep -A2 BOX bilayer_only_lipid.top
~~~
{: .language-bash}

~~~
...
%FORMAT(5E16.8)                                                                 
  9.00000000E+01  1.04281300E+02  1.03995000E+02  1.00909000E+02
~~~
{: .output}

The last three numbers refer to the size of the box in the X, Y, and Z axes. We can now parameterize the system using the Lipid21 force field and OPC3-Pol polarizable water as follows:

~~~
tleap -f <(cat << EOF
source leaprc.water.opc3pol
source leaprc.lipid21 
sys=loadpdb bilayer_only_lipid.pdb
set sys box {104.2 103.9 100.9}
saveamberparm sys bilayer.parm7 bilayer.rst7
quit
EOF)
~~~
{: .language-bash}

Don't forget to replace {104.2 103.9 100.9} with your own numbers.

### Step 3: Energy minimization
Minimization input file (minimization.in).
~~~
Minimization
 &cntrl
  imin=1,       ! Minimization 
  ntmin=3,      ! Steepest Descent 
  maxcyc=1000,  ! Maximum number of cycles for minimization
  ntpr=10,      ! Print to mdout every ntpr steps        
  ntwr=1000,    ! Write a restart file every ntwr steps
  cut=10.0,     ! Nonbonded cutoff
 /
 ~~~
 {: .file-content}

 Submission script
 ~~~
 #!/bin/bash
#SBATCH --mem-per-cpu=4000 --ntasks=10 --time=6:0:0 

module purge
module load StdEnv/2020  gcc/9.3.0  cuda/11.4  openmpi/4.0.3 amber/20.12-20.15

srun sander.MPI -O -i min.in \
-o minimization.out \
-p ../bilayer.parm7 \
-c ../bilayer.rst7 \
-r minimized.rst7 
 ~~~
{: .language-bash}

### Step 4: Heating
~~~
Heating 
 &cntrl
  imin=0,         	! Molecular dynamics
  ntx=1,          	! Read coordinates with no initial velocities
  ntc=2,          	! SHAKE on for bonds with hydrogen
  ntf=2,         	! No force evaluation for bonds with hydrogen
  tol=0.0000001,  	! SHAKE tolerance
  nstlim=20000,   	! Number of MD steps
  ntt=1,          	! Berendsen thermostat
  tempi=100,       	! Initial temperature
  temp0=300,        ! Target temperature
  tautp=5.0,        ! Time constant for heat bath coupling
  ntpr=100,         ! Write energy info every ntwr steps 
  ntwr=10000,       ! Write restart file every ntwr steps
  ntwx=2000,      	! Write to trajectory file every ntwx steps
  dt=0.001,       	! Timestep 
  ntb=2,            ! Constant pressure periodic bc.
  barostat=1,       ! Berendsen barostat
  taup=5,           ! Pressure relaxation time
  ntp=3,            ! Semiisotropic pressure scaling
  csurften=3,       ! Constant surface tension with interfaces in the xy plane
  cut=10.0,         ! Nonbonded cutoff
 /
~~~
{: .file-content}

Submission script
~~~
#!/bin/bash
#SBATCH --mem-per-cpu=4000 --gpus-per-node=v100:1 --time=6:0:0 

module purge
module load StdEnv/2020  gcc/9.3.0  cuda/11.4  openmpi/4.0.3 amber/20.12-20.15

pmemd.cuda -O -i heating.in \
-o heating.out \
-p ../bilayer.parm7 \
-c minimized.rst7 \
-r heated.rst7 \
~~~
{: .language-bash}

### Step 5: The first phase of the equilibration process.
~~~
Equilibration 1 
 &cntrl
  imin=0,         	! Molecular dynamics
  irest=1,          ! Restart
  ntx=5,          	! Read positions and velocities
  ntc=2,          	! SHAKE on for bonds with hydrogen
  ntf=2,         	! No force evaluation for bonds with hydrogen
  tol=0.0000001,  	! SHAKE tolerance
  nstlim=100000,   	! Number of MD steps
  ntt=1,          	! Berendsen thermostat
  temp0=300,        ! Set temperature
  tautp=5.0,        ! Time constant for heat bath coupling
  ntpr=100,         ! Write energy info every ntwr steps 
  ntwr=10000,       ! Write restart file every ntwr steps
  ntwx=2000,      	! Write to trajectory file every ntwx steps
  dt=0.001,       	! Timestep (ps)
  ntb=2,            ! Constant pressure periodic bc.
  barostat=1,       ! Berendsen barostat
  taup=5,           ! Pressure relaxation time
  ntp=3,            ! Semiisotropic pressure scaling
  csurften=3,       ! Constant surface tension with interfaces in the xy plane
  cut=10.0,         ! Nonbonded cutoff 
 /
~~~
{: .file-content}

Submission script
~~~
#!/bin/bash
#SBATCH --mem-per-cpu=4000 --gpus-per-node=v100:1 --time=6:0:0 

module purge
module load StdEnv/2020  gcc/9.3.0  cuda/11.4  openmpi/4.0.3 amber/20.12-20.15

pmemd.cuda -O -i equilibration.in \
-o eq-1.out \
-p ../bilayer.parm7 \
-c heated.rst7 \
-r eq-1.rst7 \
~~~
{: .language-bash}

### Step 6: The second phase of the equilibration process.
~~~
Equilibration 2 
 &cntrl
  imin=0,         	! Molecular dynamics
  irest=1           ! Restart
  ntx=5,            ! Read  positions and  velocities
  ntc=2,          	! SHAKE on for bonds with hydrogen
  ntf=2,         	! No force evaluation for bonds with hydrogen
  tol=0.0000001,  	! SHAKE tolerance
  nstlim=2500000,   ! Number of MD steps
  ntt=3,            ! Langevin thermostat 
  temp0=300,        ! Target temperature
  gamma_ln=2.0,     ! Collision frquency 
  ntpr=100,         ! Write energy info every ntwr steps 
  ntwr=10000,       ! Write restart file every ntwr steps
  ntwx=10000,      	! Write to trajectory file every ntwx steps
  dt=0.002,       	! Timestep (ps)
  ntb=2,            ! Constant pressure periodic bc.
  barostat=1,       ! Berendsen barostat
  taup=1.0,         ! Pressure relaxation time
  ntp=3,            ! Semiisotropic pressure scaling
  csurften=3,       ! Constant surface tension with interfaces in the xy plane
  cut=10.0,         ! Nonbonded cutoff 
 /
~~~
{: .file-content}

Submission script
~~~
#!/bin/bash
#SBATCH --mem-per-cpu=4000 --gpus-per-node=v100:1 --time=6:0:0 

module purge
module load StdEnv/2020  gcc/9.3.0  cuda/11.4  openmpi/4.0.3 amber/20.12-20.15

pmemd.cuda -O -i equilibration-2.in \
-o eq-2.out \
-p ../bilayer.parm7 \
-c eq-1.rst7 \
-r eq-2.rst7 \
~~~
{: .language-bash}
