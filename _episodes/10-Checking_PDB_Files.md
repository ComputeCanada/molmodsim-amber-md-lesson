---
title: "Checking and Preparing PDB Files"
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
{% if site.show_comments %}
Many commercial programs and interactive graphical interfaces such as CHARMM-GUI are available to help prepare a simulation system. These tools are easy to use and do not require as much learning effort as command-line tools, however, the functionality is limited, and results obtained with WEB/GUI tools are not reproducible and prone to human error.  Therefore, we will focus on preparing the system using only scriptable command-line-driven tools. This lesson is intended to expose you to various methods that can be used to create a reproducible molecular modeling workflow by automating preparation and simulation steps. One benefit of this approach is that once a workflow script has been developed, it can be easily adapted to other systems or conditions (for example, if a new pdb file is released, you can prepare a new simulation system with one click).
{% endif %}

#### What data is needed to setup a simulation?
{% if site.show_comments %}
Let's have a closer look at what data is needed to setup a simulation. Molecular simulation systems are typically built using PDB files. PDB files are essentially plain text files. The PDB format defines many types of records, describing structural information, crystallographic experiments, secondary structures, missing residues and other information. For setting up a simulation system, preparation programs only need the coordinate section, which consists of ATOM, HETATM, and TER records.
{% else %}
 - Molecular simulation systems are typically prepared from PDB files. 
 - For a simulation to be setup, only the coordinate section consisting of ATOM, HETATM, and TER records is required.
{% endif %}

```
         atomName  chain            coordinates                temperatureFactor (beta)
             |       |            x       y       z             |
ATOM      1  N   MET A   1      39.754  15.227  24.484  1.00 46.61
                 |       |                               |
           residueName  residueID                       occupancy
 
ATOM    147  CG AASP A  20      53.919   7.536  24.768  0.50 31.95          
ATOM    148  CG BASP A  20      55.391   5.808  23.334  0.50 32.16    
                |                                     
               conformation                         

TER  - indicates the end of a chain

HETATM  832  O   HOH A 106      32.125   6.262  24.443  1.00 21.18    
```
{% if site.show_comments %}
The lines beginning with "ATOM" represent the atomic coordinates for standard amino acids and nucleotides. "TER" records indicate the end of a chain. For chemical compounds other than proteins or nucleic acids, the "HETATM" record type is used. Records of both types use a simple fixed-column format explained [here](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM). 
{% else %}
- The lines beginning with "ATOM" represent the atomic coordinates for standard amino acids and nucleotides.
- For chemical compounds other than proteins or nucleic acids, the "HETATM" record type is used.
- Records of both types use a simple fixed-column format explained [here](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM). 
- "TER" records indicate which atoms are at the terminal of a protein chain.
{% endif %} 
{% if site.show_comments %}
Before we can successfully import a PDB file into LEAP and produce the system topology file, we need to ensure that the original PDB files are free from errors and the molecules we want to simulate are chemically correct.
{% endif %} 

## Important Things to Check in a PDB File
{% if site.show_comments %}
To simulate molecules correctly we need to ensure that:
- the original input PDB files are error-free 
- the molecules we want to simulate are chemically correct
{% else %}
A correct simulation of molecules requires error-free input PDB files.
{% endif %}

{% if site.show_comments %}
Small errors in the input structure may cause MD simulations to become unstable or give unrealistic results.
{% endif %}

There are several common problems with PDB files, including:
- presence of non-protein molecules (crystallographic waters, ligands, modified amino acids, etc.)
- alternate conformations
- missing side-chain atoms
- missing fragments
- clashes between atoms
- multiple copies of the same protein chains
- di-sulfide bonds
- wrong assignment of the N and O atoms in the amide groups of ASN and GLN, and the N and C atoms in the imidazole ring of HIS
{% if site.show_comments %}
Some problems can be identified and corrected automatically (such as missing atoms and some steric clashes), while others may have more than one solution and require your decision. In this section, you will learn how to recognize and correct problems associated with multiple chains, alternate conformations, non-protein molecules, and disulphide bonds.
{% endif %}

#### Connect to the training cluster
Let's consider some example protein PDB files. The first step is to connect to the training cluster. Sign in to Jupyter Hub **jupyter.moledyn.ace-net.training** and start a server with the following arguments:  
- **4 CPUs** 
- **4 hours** 
- **default RAM**
- **1 GPUs**
<br>

#### Workshop data: 

On the training cluster copy archive in your home directory and unpack:
~~~
cd
cp /project/def-sponsor00/workshop_amber_2024.tar.gz .
tar xf workshop_amber_2024.tar.gz
~~~
{: .language-bash}

[Download link](https://github.com/ComputeCanada/molmodsim-amber-md-lesson/releases/download/workshop-2021-04/workshop_amber_2024.tar.gz)


### Checking a molecular structure 
{% if site.show_comments %}
[Check_structure](https://pypi.org/project/biobb-structure-checking/) is a command-line utility from [BioBB project](https://github.com/bioexcel/biobb) for exhaustive structure quality checking (residue chirality, amide orientation, vdw clashes, etc.).  Using this utility, you can perform manipulations with structures, such as selecting chains or conformations, removing components, mutating residues, adding missing atoms, adding hydrogens, etc. 
{% else %}
- [Check_structure](https://pypi.org/project/biobb-structure-checking/) is a command-line utility from [BioBB project](https://github.com/bioexcel/biobb) for exhaustive structure quality checking.
{% endif %}
{% if site.show_comments %}
Check_structure is a python module. Python modules are installed in user accounts in a python virtual environment.
{% else %}
Installing check_structure. 
{% endif %}
~~~
module load StdEnv/2023 python scipy-stack
virtualenv ~/env-biobb
source ~/env-biobb/bin/activate
pip install biobb-structure-checking
~~~
{: .language-bash}

Using check_structure. 
~~~
cd ~/workshp_amber/example_01
check_structure commands # print help on commands
check_structure -i 2qwo.pdb checkall 
~~~
{: .language-bash}

~~~
...
Detected no residues with alternative location labels
...
Found 154 Residues requiring selection on adding H atoms
...
Detected 348 Water molecules
...
Detected 8 Ligands
...
Detected 1 Possible SS Bonds
...
No severe clashes detected
~~~
{: .output}

#### Removing Non-Protein Molecules
{% if site.show_comments %}
A PDB file containing any molecules other than proteins or nucleic acids needs special treatment. It is common for PDB files to contain solvents, ions, lipid molecules, protein co-factors, e.t.c. In some cases, these extra components are necessary for the protein function and should be included in the simulation. It is common to add compounds to facilitate crystallization. These compounds are usually not necessary for simulation. In this introductory lesson, we won't consider them.
{% endif %}

Let's remove ligands and save the output in a new file called "protein.pdb".
~~~
check_structure -i 2qwo.pdb -o protein.pdb ligands --remove all  
~~~
{: .language-bash}

>## Selecting protein atoms using VMD
>Select only protein atoms from the file `2qwo.pdb` and save them in the new file `protein.pdb` using VMD.
>
>>## Solution
>>Load vmd module and start the program:
>>~~~
>>module load StdEnv/2023 vmd
>>vmd
>>~~~
>>{: .language-bash}
>>
>>~~~
>>mol new 2qwo.pdb
>>set prot [atomselect top "protein"]
>>$prot writepdb protein.pdb
>>quit
>>~~~
>>{: .vmd}
>>The first line of code loads a new molecule from 2qwo.pdb. Using the **atomselect** method, we then select all protein atoms from the top molecule. Finally, we save the selection in the file "protein.pdb".  
>>The Atom Selection Language has many capabilities. You can learn more about it by visiting the following [webpage](https://www.ks.uiuc.edu/Research/vmd/vmd-1.3/ug/node132.html). 
>>
>{: .solution}
{: .challenge}

>## Selecting protein atoms using shell commands
>Standard Linux text searching utility `grep` can find and print all "ATOM" records from a PDB file. This is a good example of using Unix command line, and `grep` is very useful for many other purposes such as finding important things in log files. `Grep` searches for a given pattern in files and prints out each line that matches a pattern. 
>1. Check if a PDB file has "HETATM" records using `grep` command.
>2. Select only protein atoms from the file `2qwo.pdb` and save them in the new file `protein.pdb` using `grep` command to select protein atoms (the "ATOM" and the "TER" records). 
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


#### Checking PDB Files for alternate conformations.
{% if site.show_comments %}
Some PDB files may contain alternate positions of residue side chains. Only one conformation is acceptable for molecular dynamics simulation. Standard simulation preparation programs such as `pdb2gmx` or `pdb4amber` will automatically select the first conformation labeled "A" in the "altLoc" column (column 17).  
<br>
Sometimes you may want to compare simulations starting from different initial conformations. If you want to select a particular conformation, all conformations except the desired one must be removed from a PDB file.
{% else %}
Check conformations
{% endif %}

~~~
cd ~/workshop_amber/example_02
check_structure -i 1ert.pdb checkall
~~~
{: .language-bash}

~~~
ASP A20
  CG   A (0.50) B (0.50)
  OD1  A (0.50) B (0.50)
  OD2  A (0.50) B (0.50)
HIS A43
  CG   A (0.50) B (0.50)
  ND1  A (0.50) B (0.50)
  CD2  A (0.50) B (0.50)
  CE1  A (0.50) B (0.50)
  NE2  A (0.50) B (0.50)
SER A90
  OG   A (0.50) B (0.50)
~~~
{: .output}

Select conformers A ASP20 and B HIS43. 
~~~
check_structure -i 1ert.pdb -o output.pdb altloc --select A20:A,A43:B,A90:B 
~~~
{: .language-bash}


>## Selecting Alternate Conformations with VMD
>1. Check if the file 1ERT.pdb has any alternate conformations. 
>2. Select conformation A for residues 43, 90. Select conformation B for residue 20. Save the selection in the file "protein_20B_43A_90A.pdb". 
>
>> ## Solution
>>1.
>>
>>~~~
>>mol new 1ert.pdb
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
{% if site.show_comments %}
Disulfide bonds are covalent bonds between the sulfur atoms of two cysteine residues. They are very important for the stabilization of protein structure.
Disulfide bonds are easy to spot in PDB files with any visualization program.  <br>

For simulation preparation with the AMBER `tleap` program, cross-linked cysteines must be renamed from "CYS" to "CYX" to distinguish them from normal cysteines. Check_structure can detect and mark disulphide bonds.
{% else %}
- For simulation preparation with the AMBER, cross-linked cysteines must be renamed from "CYS" to "CYX" 
{% endif %} 

~~~
cd ~/workshop_amber/example_01
check_structure -i 2qwo.pdb -o output.pdb getss --mark all
grep CYX output.pdb
~~~
{: .language-bash}

- GROMACS `pdb2gmx` utility can automatically add S-S bonds to the topology based on the distance between sulfur atoms (option *-ss*).

#### Useful Links
[MDWeb](http://mmb.irbbarcelona.org/MDWeb2) server can help to identify problems with PDB files and visually inspect them. It can also perform complete simulation setup, but options are limited and waiting time in the queue may be quite long.

[CHARMM-GUI](http://www.charmm-gui.org) can be used to generate input files for simulation with CHARMM force fields. CHARMM-GUI offers useful features, for example the "Membrane Builder" and the "Multicomponent Assembler".
