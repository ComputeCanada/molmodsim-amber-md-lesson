---
title: "Solvating a System, Adding Ions and Generating Input Files"
teaching: 30
exercises: 5
questions:
- "Why the simulation system should be neutralized?"
- "How to add water and ions to a simulation system?"
- "How to choose size and shape of a solvent box?"
objectives:
- "Understand why it is necessary to neutralize the simulation system."
- "Neutralize a system."
- "Solvate a macromolecule."
- "Add ions to a simulation system to mimic a salt solution."
- "Generate molecular topology for simulation with GROMACS and NAMD."
keypoints:
- "Simulation system must be neutralized by adding counter-ions to obtain the correct electrostatic energy."
- "Ions are added to a simulations system to mimic composition of a local macromolecule environment."
- "Solvent box should be large enough to allow for unrestricted conformational dynamics of a macromolecule."
---
### Why ions are added to simulations?

When periodic boundary conditions are applied and grid-based methods are used to compute the Coulomb energy, the simulation box interacts with the infinite number of its periodic images. As a result, if a simulation system is charged, the electrostatic energy will essentially add to infinity. To solve this issue, we need to add counter-ions to neutralize the system so that the electrostatic energy can be correctly calculated during simulation.

Another reason is that the conformations, dynamics, and functions of biological macromolecules are influenced by ion concentration and composition of the local environment.

## Neutralizing a system

Fist we will add enough counter-ions to neutralize the system. Ions can be added using two approaches:
1. Solvate the system and then replace random solvent molecules with ions.
2. Place ions according to the electrostatic potential of the macromolecule before solvation.

Adding more ions to a neutralized system will be necessary to represent physiological salt concentrations.

### Caveats and limitations of the random ion placement

{: .instructor_notes} 
Random placement of ions will generate a system in the completely dissociated, energetically unfavorable state. The random placement of ions is problematic if the electric charge of the macromolecule is big (for example DNA) because ions tend to form screening clouds around charged molecules rather than being distributed randomly. Random placement of ions will negatively impact the time required for the system equilibration and may affect structural stability of a macromolecule. A better approach is to place ions according to the electrostatic potential of the macromolecule. Such method is implemented in the *leap* module of the *AmberTools*. The *addions* command adds ions to simulation cells near the minima of the solute's electrostatic potential field.
{: .instructor_notes} 

{: .self_study_text} 
- Random placement of ions will require longer equilibration and may affect structural stability of a macromolecule.
- A better approach is to place ions according to the electrostatic potential of the macromolecule. 
{: .self_study_text} 

Let's neutralize 1RGG protein using the *leap* module. We will add ions prior to solvation so that the potential from un-equilibrated water will not interfere with ion placement:

~~~
cd ~/workshop/pdb/1RGG/AMBER
module load gcc/9.3.0 cuda/11.4 ambertools/22
tleap
~~~
{: .language-bash}

~~~
source leaprc.water.opc
source leaprc.protein.ff19SB
s = loadpdb ../1RGG_chain_A_prot.pdb
charge s
addions s Na+ 0
~~~
{: .leap}

## Adding Ions to Mimic the Macroscopic Salt Concentration
To mimic the macroscopic salt concentration in the environment being studied we will need to add more cation/anion pairs to the simulation system. The number of ion pairs can be estimated using the formula:

$N_{Ions}=0.0187\cdot[Molarity]\cdot{N_{WaterMol}}$


{: .instructor_notes} 
The drawback of this calculation is that it does not take into account the charge of a macromolecule. As charged solute perturbs the local solvent environment by depleting ions from the bulk this method generates a system with the salt concentration that is too high. For more accurate salt concentration you can calculate the number of ions corrected for screening effects using the [*SLTCAP*](https://www.phys.ksu.edu/personal/schmit/SLTCAP/SLTCAP.html) server.

As you can see from the equation above, to calculate the number of ions we need to know the number of water molecules in the simulation system. So we continue our *leap* session and solvate the simulation system. In this lesson we will create a simple cubic solvent box. As we discussed in the episode "Periodic Boundary Conditions", a properly solvated simulation system should have at least 10 <span>&#8491;</span> distance between the solute and the edge of the periodic box after equilibration. Standard practice is to tile a pre-equilibrated solvent box across the system and eliminate solvent molecules which clash with the solute.
{: .instructor_notes} 

{: .self_study_text} 
- This calculation does not take into account the charge of a macromolecule.
- For more accurate salt concentration you can calculate the number of ions corrected for screening effects using the [*SLTCAP*](https://www.phys.ksu.edu/personal/schmit/SLTCAP/SLTCAP.html) server.
- To calculate the number of ions we need to know the number of water molecules in the simulation system.
{: .self_study_text} 

{: .instructor_notes} 
When water is added to the system in this way, some VDW voids at the macromolecule and the box interfaces are inevitable because packing is not perfect. In the process of equilibration the water molecules will move to fill the voids and minimize the interaction energy. The box will shrink and the distance between the solute and the edge of the periodic box will become smaller. To compensate for this box shrinkage we need to start with a larger box size than the desired. The rule of thumb is that you need to add at least 2.5 <span>&#8491;</span> to the box size.
{: .instructor_notes} 

We will use the *solvateBox* command to create the periodic solvent box around the macromolecule. The *solvateBox* command has many options. Let's create a cuboid water box around the 1RGG protein. We will use the pre-equilibrated box of SPCE water (SPCBOX), set the minimum distance between the solute and the edge of the box to 15 <span>&#8491;</span>, and request an isometric (cubic) box:
~~~
solvatebox s SPCBOX 15 iso
~~~
{: .leap}

~~~
  Solute vdw bounding box:              40.514 32.235 37.352
  Total bounding box for atom centers:  70.514 70.514 70.514
      (box expansion for 'iso' is  18.6%)
  Solvent unit box:                     18.774 18.774 18.774
  Volume: 399256.044 A^3
  Total mass 194369.824 amu,  Density 0.808 g/cc
  Added 10202 residues.
 ~~~
 {: .output}

Now that we know the number of water molecules in the simulation system, we can add salt to the desired concentration.

> ## Preparing an Aqueous Salt Solution
> How many Na+ and Cl- ions do we need to add to the simulation box with 1RGG protein and 10202 water molecules to prepare 0.15 M salt solution?
> Calculate the number of ions using two methods: the formula above and the [*SLTCAP*](https://www.phys.ksu.edu/personal/schmit/SLTCAP/SLTCAP.html) server. For this calculation you need to know molecular weight of the protein. You can calculate it [here](https://www.bioinformatics.org/sms/prot_mw.html). FASTA sequence of the protein is available [here](https://www.rcsb.org/fasta/entry/1RGG/display).
>
>>## Solution
>> 1. N_ions = 0.0187 x 0.15 x 10202 = 29. We need to add 35 Na+ and 29 Cl- ions
>> 2. *SLTCAP* calculation with the following input: (MW 11 KDa, 10202 water molecules, charge 6, 150 mM salt) yields 30 Na+ and 24 Cl- ions.
>>
> {: .solution}
{: .challenge}

We already have the neutralized and solvated simulation system, and in the exercise above we determined that we need to add 24 ion pairs to prepare 150 mM salt solution. Let's replace 48 randomly selected water molecules with 24 Na+ and 24 Cl- ions:
~~~
addionsrand s Na+ 24 Cl- 24
~~~
{: .leap}


## Generating Molecular Topology for Simulation with *AMBER* or *NAMD*
Setup of our simulation is almost complete. Our protein has cross-linked cysteine residues, so the last thing to do is to make disulfide bond between Cys7 and Cys96:

~~~
bond s.7.SG s.96.SG
~~~
{: .leap}

We can now save the molecular topology (*parm7*) file and *AMBER* coordinates (*rst7*). To build *GROMACS* topology later we will also save the solvated system in PDB format:

~~~
saveamberparm s 1RGG_chain_A.parm7 1RGG_chain_A.rst7
savepdb s 1RGG_chain_A_solvated.pdb
quit
~~~
{: .leap}


### A complete script for preparing coordinates and the topology in LEAP
Save the following commands in a file, e.g. solvate_1RRG.leap
~~~
source leaprc.water.opc
source leaprc.protein.ff19SB
s = loadpdb ../1RGG_chain_A_prot.pdb
addions s Na+ 0
solvatebox s SPCBOX 15 iso
addionsrand s Na+ 24 Cl- 24
bond s.7.SG s.96.SG
saveamberparm s 1RGG_chain_A.parm7 1RGG_chain_A.rst7
savepdb s 1RGG_chain_A_solvated.pdb
quit
~~~
{: .leap}
Execute the script: 
~~~
tleap -f solvate_1RRG.leap
~~~
{: .language-bash}


### Automation and Reproducibility of a Simulation Setup

{: .instructor_notes} 
The process of molecular dynamics system setup can be automated by saving the whole sequence of commands into a text file. This file (shell script) can be easily modified to prepare simulation systems from other PDB structure files or used to reproduce your simulation setup. All system preparation steps can be completed in seconds with the script.
{: .instructor_notes} 

{: .self_study_text} 
- The process of molecular dynamics system setup can be automated by saving the whole sequence of commands into a text file.
{: .self_study_text} 


You can download an example shell script that performs all preparation steps [here]({{ page.root }}/code/run_setup.sh). The script downloads the molecular structure file from PDB and generates input files for simulation with AMBER, NAMD, and GROMACS.

### How to create ligand topology

[Automated Topology Builder](https://atb.uq.edu.au/index.py)



