---
title: "Preparing a Protein-Membrane Simulation System"
teaching: 30
exercises: 5
questions:
- "How to prepare a membrane simulation system?"
- "How to pack a protein in a lipid bilayer?"
objectives:
- "Learn how to prepare a membrane simulation system?"
- "Learns to pack a protein in a lipid bilayer?"
keypoints:
- " "
---

### Lipid force field
Amber currently includes Lipid21 as its main membrane force field.  In this modular force field, lipids are modeled as polymers composed of a headgroup and acyl tails. Essentially, this means that each headgroup and tail are independent modules, analogous to protein residues. Each lipid molecule is composed of a tail analogous to an "N-terminal", a central headgroup and another tail analogous to a "C-terminal". You can combine any headgroup with any pair of tails in this force field. Ok, now you should have an idea how the lipids are represented in the Lipid21 force field, and what lipids you can include in a simulation. 

- There is a great deal of information provided in [Dickson, 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9007451/) regarding best practices when using AMBER for lipid simulations.
- Full list of supported lipids is in Section 3.4 of [Amber22 manual](https://ambermd.org/doc12/Amber22.pdf)

Known Issues: 
Using MC barostat with hard LJ cutoff is known to cause bilayer deformation. It is recommended to use an LJ force switch when running simulations with the MC barostat.[Gomez, 2021](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.26798)

### Creating simulation systems with packmol-memgen
- What lipids are available?

~~~
ml StdEnv/2020 gcc cuda/11.4 ambertools/22
packmol-memgen --available_lipids 
~~~
{: .language-bash}

- To see all available lipids use  `--available_lipids_all`, but the list will have thousands of items!

#### Creating a membrane-only simulation system
Submission script for creating a membrane-only simulation system:
~~~
#!/bin/bash
#SBATCH -c1 --mem-per-cpu=2000 --time=3:0:0

module purge
module load StdEnv/2020 gcc cuda/11.4 ambertools/22
packmol-memgen \
    --lipids DOPE:DOPG \
    --ratio 3:1 \
    --distxy_fix 100 \
    --parametrize
~~~
{: .language-bash}

- if the option --parametrize is given the solvated system is bilayer_only_lipid.pdb
- without --parametrize the solvated system is  bilayer_only.pdb

Another example:
~~~
#!/bin/bash
#SBATCH -c1 --mem-per-cpu=4000 --time=3:0:0

module purge
module load StdEnv/2020 gcc/9.3.0 cuda/11.4 openmpi/4.0.3 ambertools/22
rm -f bilayer* *.log slurm*
packmol-memgen \
    --lipids CHL1:POPC:POPE:PSM \
    --salt --salt_c Na+ --saltcon 0.15 \
    --dist_wat 15 \
    --ratio 4:2:1:1 \
    --distxy_fix 100 \
    --parametrize
~~~
{: .language-bash}

This generates bilayer_only_lipids.pdb which can be used to create topology files with AMBER force fields not available in packmol-memgen. Example using OPC3-Pol polarizable water:

~~~
leap_input=$(cat << EOF
source leaprc.water.opc3pol
source leaprc.lipid21 
sys=loadpdb bilayer_only_lipid.pdb
set sys box {111.3 111.3 81.5}
saveamberparm sys bilayer.parm7 bilayer.rst7
quit
EOF)
tleap -f <(echo "$leap_input")
~~~
{: .language-bash}

Box size can be extracted from the file bilayer_only_lipid.top:

~~~
grep -A2 BOX bilayer_only_lipid.top
~~~
{: .language-bash}

### Embedding a protein into a bilayer
We will use PDB file 6U9P (wild-type MthK pore in ~150 mM K+) for this exercise.

1. Use [PPM server](https://opm.phar.umich.edu/ppm_server) to orient a protein. PPM server will also take care of assembling the complete tetrameric pore.
2. Use vmd to remove ligands and conformers B.

~~~
module load StdEnv/2020 vmd
vmd
~~~
{: .language-bash}

~~~
mol new 6u9pout.pdb
set sel [atomselect top "protein and not altloc B"]
$sel writepdb 6U9P-clean.pdb
quit
~~~
{: .vmd}

3. Submission script to embed a protein into a lipid bilayer

 ~~~
 #!/bin/bash
 #SBATCH -c1 --mem-per-cpu=2000 --time=1:0:0
 
 module purge
 module load StdEnv/2020 gcc cuda/11.4 ambertools/22
 packmol-memgen \
	--pdb 6U9P-clean.pdb \ 
	--lipids DOPE:DOPG \
	--ratio 3:1 \
    --preoriented 
 ~~~
 {: .language-bash}

- Use `--parametrize` to make simulation input files. The default protein force field is FF14SB, water model TIP3P. 
- To use a different force field:
`--fprot ff19SB --ffwat opc --gaff2`
- To add salt (default K+, 0.15M): 
`--salt`


### Preparing ligands 

- Download hexanediol HEZ.cif: wget https://www.rcsb.org/ligand/HEZ
- Convert .cif to .pdb: https://mmcif.pdbj.org/converter/index.php?l=en

1. Make mol2 file with antechamber:

~~~
 module load StdEnv/2020 gcc cuda/11.4 ambertools/22
 antechamber -i HEZ.pdb -fi pdb -o HEZ.mol2 -fo mol2 -c bcc -s 2
~~~
{: .language-bash}

The file HEZ.mol2 contains the definition of our HEZ residue including connectivity, all of the charges and atom types. 

2. Run parmchk2 to find out if there are any missing force field parameters:

~~~
 parmchk2 -i HEZ.mol2 -f mol2 -o HEZ.frcmod
~~~
{: .language-bash}

If it can antechamber will fill in these missing parameters by analogy to a similar parameter.

3. Create the library file for HEZ using tleap:

~~~
source leaprc.gaff2
HEZ=loadmol2 HEZ.mol2
check HEZ
saveoff HEZ hez.lib 
~~~
{: .leap}

Of course point charges are not very accurate because they are derived using semi-empirical method, but antechamber can also use results of gaussian QM calculations.  

### Links to advanced AMBER tutorials
- [Placing waters and ions using 3D-RISM and MOFT](http://ambermd.org/tutorials/advanced/tutorial34/index.html)
- [Building and Equilibrating a Membrane System with PACKMOL-Memgen](https://ambermd.org/tutorials/advanced/tutorial38/index.php#Lipid_System)
- [Setup and simulation of a membrane protein with AMBER Lipid21 and PACKMOL-Memgen](https://github.com/callumjd/AMBER-Membrane_protein_tutorial). 