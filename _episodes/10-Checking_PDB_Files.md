---
title: "Checking PDB Files"
teaching: 20
exercises: 10
questions:
- "What problems are commonly found in PDB files?"
- "Why fixing errors in PDB files is essential for a simulation?"
objectives:
- "Understand why it is necessary to check PDB files before simulation."
- "Learn how to look for problems in PDB files."
- "Learn how to fix the common errors in PDB files."
keypoints:
- "Small errors in the input structure may cause MD simulations to became unstable or give unrealistic results."
---
For this lesson, we will go over the steps for setting up a fully solvated protein system for simulation with AMBER/NAMD.

Many commercial programs and interactive graphical interfaces are available to help prepare the system. These tools are easy to use and do not require as much learning effort as command-line tools, however, the functionality is limited, and results obtained with WEB/GUI tools are not reproducible and prone to human error. Therefore, we will focus on preparing the system using only scriptable command-line-driven tools. This lesson is intended to expose you to various methods that can be used to create a reproducible molecular modeling workflow by automating preparation and simulation steps. One benefit of this approach is that once a workflow script has been developed, it can be easily adapted to other systems or conditions (for example, if a new pdb file is released, you can prepare a new simulation system with one click).

Before we can successfully import a PDB file into LEAP and produce the system topology file, we need to ensure that the original files are free from errors and the molecules we want to simulate are chemically correct.


## Important Things to Check in a PDB File
Small errors in the input structure may cause MD simulations to become unstable or give unrealistic results. The most common problems in PDB files include:

- non-protein molecules (crystallographic waters, ligands, modified amino acids, etc.)
- alternate conformations
- missing side-chain atoms
- missing fragments
- clashes between atoms
- multiple copies of the same protein chains
- di-sulfide bonds
- wrong assignment of the N and O atoms in the amide groups of ASN and GLN, and the N and C atoms in the imidazole ring of HIS

Several problems can be identified and corrected automatically (such as missing atoms and some steric clashes), while others may have more than one solution  and require your decision. In this section, you will learn how to recognize and correct problems associated with multiple chains, alternate conformations, non-protein molecules, and disulphide bonds.

#### Downloading a Protein Structure File from PDB
Let's start with downloading a protein structure file:
~~~
cd ~/scratch/workshop/pdb/1ERT
wget http://files.rcsb.org/view/1ERT.pdb
~~~
{: .language-bash}

#### Checking PDB Files for Presence of Non-Protein Molecules
PDB files are essentially text files containing important structural information, as well as information about crystallographic experiments, secondary structures, and residues that are missing. To set up an MD simulation system, we will only need the coordinate section, which includes ATOM, TER, and HETATM records. 

The lines beginning with "ATOM" represent the atomic coordinates for standard amino acids and nucleotides. "TER" records indicate which atoms are at the terminal of a protein chain. For chemical compounds other than proteins or nucleic acids, the "HETATM" record type is used. Record types of both types use a simple fixed-column format explained [here](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM). 

A PDB file containing any molecules other than proteins or nucleic acids needs special treatment. It is common for PDB files to contain solvents, ions, lipid molecules, protein co-factors, e.t.c. In some cases, these extra components are necessary for the protein function and should be included in the simulation. It is common to add compounds to facilitate crystallization. These compounds are usually not necessary for simulation. In this introductory lesson, we won't consider them.

Let's select only protein atoms and save them in a new file called "protein.pdb". We will be using [VMD](https://www.ks.uiuc.edu/Research/vmd/) for this.


Load vmd module and start the program:
~~~
ml vmd
vmd
~~~
{: .language-bash}

~~~
mol new 1ERT.pdb
set prot [atomselect top "protein"]
$prot writepdb protein.pdb
quit
~~~
{: .vmd}
The first line of code loads a new molecule from 1ERT.pdb. Using the **atomselect** method, we then select all protein atoms from the top molecule. 

The Atom Selection Language has many capabilities.You can learn more about it by visiting the following [webpage](https://www.ks.uiuc.edu/Research/vmd/vmd-1.3/ug/node132.html). 

Finally, we store the selection in a file named "protein.pdb".


>## Selecting ATOM Records Using Shell Commands
>Standard Linux text searching utility `grep` can find and print all "ATOM" records from a PDB file. This is a good example of using Unix command line, and `grep` is very useful for many other purposes such as finding important things in log files. `Grep` searches for a given pattern in files and prints out each line that matches a pattern. 
>1. Check if a PDB file has "HETATM" records using `grep` command.
>2. Select only protein atoms from the file `1ERT.pdb` and save them in the new file `protein.pdb` using `grep` command to select protein atoms (the "ATOM" and the "TER" records). 
>
>Hint: the `OR` operator in *grep* is `\|`. The output from a command can be redirected into a file using the output redirection operator `>`.
>
>>## Solution
>> 1.
>>
>>~~~
>>grep "^HETATM " 1ERT.pdb | wc -l
>>~~~
>>{: .language-bash}
>>~~~
>>      46
>>~~~
>>{: .output}
>>The `^` expression matches beginning of line. We used the `grep` command to find all lines beginning with the word "HETATM" and then we sent these lines to the character counting command `wc`. The output tells us that the downloaded PDB file contains 46 non-protein atoms. In this case, they are just oxygen atoms of the crystal water molecules.
>> 
>> 2.
>>~~~
>> grep "^ATOM\|^TER " 1ERT.pdb > protein.pdb
>>~~~
>>{: .language-bash}
>{: .solution}
{: .challenge}

#### Checking PDB Files for Alternate Conformations.
Some PDB files may contain alternate positions of residues. Only one conformation is acceptable for molecular dynamics simulation. Standard simulation preparation programs such as `pdb2gmx` or `pdb4amber` will automatically select the first conformation labeled "A" in the "altLoc" column (column 17). 

Sometimes you may want to compare simulations starting from different initial conformations. If you want to select a particular conformation, all conformations except the desired one must be removed from a PDB file.

>## Selecting Alternate Conformations with VMD
>1. Check if the file 1ERT.pdb has any alternate conformations. 
>2. Select conformation A for residues 43, 90. Select conformation B for residue 20. Save the selection in the file "protein_20B_43A_90A.pdb". 
>
>> ## Solution
>>1.
>>
>>~~~
>>mol new 1ERT.pdb
>>set s [atomselect top "altloc A"]
>>$s get resid
>>set s [atomselect top "altloc B"]
>>$s get resid
>>$s get {resid resname name} 
>>set s [atomselect top "altloc C"]
>>$s get resid
>>quit
>>~~~
>>{: .vmd}
>>The output of the commands tells us that residues 20, 43 and 90 have alternate conformations A and B.    
>>  
>>2.
>>  
>>~~~
>>mol new 1ERT.pdb
>>set s [atomselect top "(protein and altloc '') or (altloc B and resid 20) or (altloc A and resid 43 90)"]
>>$s writepdb protein_20B_43A_90A.pdb
>>quit
>>~~~
>>{: .vmd}
>>
>{:.solution}
{:.challenge}

#### Checking PDB Files for cross-linked cysteines.
Disulfide bonds are covalent bonds between the sulfur atoms of two cysteine residues. They are very important for the stabilization of protein structure.
Disulfide bonds are easy to spot in PDB files with any visualization program. For example, [MDWeb](http://mmb.irbbarcelona.org/MDWeb2) server can identify disulfide bonds and many other problems in PDB files. GROMACS `pdb2gmx` utility can automatically add S-S bonds to the topology based on the distance between sulfur atoms (option *-ss*).  
For simulation preparation with the AMBER `tleap` program, cross-linked cysteines must be renamed from "CYS" to "CYX" to distinguish them from normal cysteines. Then the bond between them must be made manually using the `bond` command.

#### Check_structure utility from BioExcel building blocks project
[Check_structure](https://pypi.org/project/biobb-structure-checking/) is a command-line utility from [BioBB project](https://github.com/bioexcel/biobb) for exhaustive structure quality checking (residue chirality, amide orientation, vdw clashes, etc.).  Using this utility, you can perform manipulations with structures, such as selecting chains or conformations, removing components, mutating residues, adding missing atoms, adding hydrogens, etc.  

Installation
~~~
~/scratch/workshop/scripts/install_check_structure.sh 
~~~
{: .language-bash}

Usage
~~~
module load StdEnv/2020 python
source ~/env-biobb/bin/activate

check_structure commands
cd ~/scratch/workshop/pdb/1ERT
check_structure -i 1ERT.pdb checkall
check_structure -i 1ERT.pdb -o output.pdb altloc --select A20:A,A43:B,A90:B 
~~~
{: .language-bash}

#### Useful Links
[MDWeb](http://mmb.irbbarcelona.org/MDWeb2) server can help to identify problems with PDB files and visually inspect them. It can also perform complete simulation setup, but options are limited and waiting time in the queue may be quite long.

[CHARMM-GUI](http://www.charmm-gui.org) can be used to generate input files for simulation with CHARMM force fields. CHARMM-GUI offers useful features, for example the "Membrane Builder" and the "Multicomponent Assembler".
