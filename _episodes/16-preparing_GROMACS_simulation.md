---
title: "Preparing a system for simulation with GROMACS"
teaching: 30
exercises: 5
questions:
- "How to Prepare Input Files for Simulation with GROMACS?"
objectives:
- "Learn to Generate Input Files for Simulation with GROMACS"

keypoints:
- " "
---

## Generating Input Files for Simulation with GROMACS.
>## What force fields are available in the loaded GROMACS module?
>When the *GROMACS* module is loaded the environment variable *EBROOTGROMACS* will be set. This variable is pointing to the GROMACS installation directory. Knowing where the *GROMACS* installation is we can find out what force fields are available:
>~~~
>module load gromacs
>ls -d $EBROOTGROMACS/share/gromacs/top/*.ff | xargs -n1 basename
>~~~
>{: .language-bash}
>~~~
>amber03.ff		amber99sb-ildn.ff	gromos45a3.ff
>amber94.ff		amberGS.ff		gromos53a5.ff
>amber96.ff		charmm27.ff		gromos53a6.ff
>amber99.ff		gromos43a1.ff		gromos54a7.ff
>amber99sb.ff		gromos43a2.ff		oplsaa.ff
>~~~
>{: .output}
{: .callout}

### Generate GROMACS Topology and Coordinate Files from the Solvated System.

~~~
cd ~/scratch/workshop/pdb/1RGG/GROMACS
~~~
{: .language-bash}

We can generate GROMACS topology from the complete simulation system prepared previously and saved in the file 1RGG_chain_A_solvated.pdb. For *pdb2gmx* to work correctly we need to rename ions from (Na+, Cl-) to (NA, CL), and rename CYX to CYS:

~~~
ATOM   1444 Na+  Na+    97      -5.058 -11.031  -0.206  1.00  0.00
ATOM   1450 Cl-  Cl-   103      19.451  -3.022   8.361  1.00  0.00
~~~
{: .file-content}

We will also assign ions to chain B.  

~~~
mol new ../1RGG_chain_A_solvated.pdb 
set s [atomselect top "resname 'Na+'"]
$s set resname "NA"
$s set name "NA"
$s set chain "B"
set s [atomselect top "resname 'Cl-'"]
$s set resname "CL"
$s set name "CL"
$s set chain "B"
set s [atomselect top "resname CYX"]
$s set resname "CYS"
set s [atomselect top all]
$s writepdb 1RGG_chain_A_solvated_gro.pdb
quit
~~~
{:.vmd}

- Do it using the global substitution function of the stream editor (sed).

~~~
cat ../1RGG_chain_A_solvated.pdb |\
sed s/"Cl-  Cl-  "/" CL  CL  B"/g |\
sed s/"Na+  Na+  "/" NA  NA  B"/g |\
sed s/CYX/CYS/g > 1RGG_chain_A_solvated_gro.pdb
~~~
{: .language-bash}


Let's make the topology using the *AMBER ff99SBildn* force field and the *tip3* water model:
~~~
gmx pdb2gmx -f 1RGG_chain_A_solvated_gro.pdb -ff amber99sb-ildn -water tip3 -ignh -chainsep id -ss << EOF 
y
EOF
~~~
{: .language-bash}

Option       |Value| Description                         
--------------|:---|:-----------------------
ignh          | -  | Ignore hydrogens in file|  
chainsep      |id  | Separate chains by chain ID. Since we assigned ions to chain B pbb2gmx will ignore TER records and put them in a separate chain  
ss            | -  | Interactive S-S bridge selection. Detect potential S-S bonds, and ask for confirmation.

The construct
~~~
<< EOF
y
EOF
~~~
{: .language-bash}

at the end of the command is to automatically confirm 'y' S-S bond. 

By default *pdb2gmx* program saved three output files:

Type                | Filename 
--------------------|---------- 
topology            | topol.top
coordinates         | conf.gro
position restraints | posre.itp
 
The names of the output files can be changed by using output options *-p*, *-o* and *-i*.

### Prepare the System Using *GROMACS* Module *pdb2gmx*.
To demonstrate how to solvate protein and add ions using *pdb2gmx* we can go back to the protein structure file 1RGG_chain_A_prot.pdb saved before solvation and repeat all system preparation steps with this GROMACS utility. Note that in this case the neutralizing ions will be added in randomly selected positions.

First we generate the topology and the coordinate file using the *AMBER ff99SBildn* force field and the *spc/e* water model:
~~~
gmx pdb2gmx -f ../1RGG_chain_A_prot.pdb -ff amber99sb-ildn -water tip3 -ignh -chainsep id -ss << EOF 
y
EOF
~~~
{: .language-bash}

Once the gromacs coordinate file *conf.gro* is created we add a periodic box to it:
~~~
gmx editconf -f conf.gro -o boxed.gro -c -d 1.5 -bt cubic
~~~
{: .language-bash}
The option '-c' positions solute in the middle of the box, the option -d specifies the distance (in nm) between the solute and the box.
Add water
~~~
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
~~~
{: .language-bash}
 Next, we need to create a binary topology file. For this we need a MD run parameters file. An empy file will be sufficient for now. Using the empty file *'mdrun.mdp'* we can generate the binary topology *'solvated.tpr'* file:
~~~
touch mdrun.mdp
gmx grompp -f mdrun.mdp -p topol.top -c solvated.gro -o solvated.tpr >& grompp.log
~~~
{: .language-bash}
When the grompp program runs it generates a lot of diagnostic messages and prints out the net charge. We saved the output in the file *'grompp.log'* so that we can find out what is the total charge of the system:
~~~
grep "total charge" grompp.log
~~~
{: .language-bash}
~~~
System has non-zero total charge: -5.000000
~~~
{: .output}

Finally we can use the *GROMACS genion* command to replace random solvent molecules with ions. We will first add cation/anion pairs to mimic a desired salt concentration and then neutralize the system by adding sodium ions (the options *-conc* [Mol/L] and *-neutral*). By default genion uses Na+ and Cl- ions. Other ions can be chosen by selecting options *-pname* [positive ion] and *-nname* [negative ion]. We also need to select a target group of solvent molecules to be replaced with ions. We will chose the 'SOL' group which is the default name of the solvent group in *GROMACS*:
~~~
$ echo "SOL" | gmx genion -s solvated.tpr -p topol.top -neutral -conc 0.15 -o neutralized.pdb
~~~
{: .language-bash}

Let's inspect the last section of the updated topology file:
~~~
tail -n 4 topol.top
~~~
{: .language-bash}
~~~
Protein_chain_A     1
SOL         11482
NA               38
CL               33
~~~
{: .output}

Create a binary topology file with the new topology and run parameters. 
~~~
gmx grompp -f minimization.mdp -p topol.top -c neutralized.pdb -o solvated_neutral.tpr
~~~
{: .language-bash}


You can see that 38 sodium and 33 chloride ions were added to the system.

>## Why PDB to GROMACS conversion fails?
>Download the structure file 2F4K from PDB and try to generate the molecular topology with *pdb2gmx*:
>~~~
>wget http://files.rcsb.org/view/2F4K.pdb
>gmx pdb2gmx -f 2F4K.pdb -ff amber99sb-ildn -water spce
>~~~
>{: .language-bash}
>Why this command fails and how to fix it?
>> ## Solution
>> The file contains two noncanonical norleucine aminoacids (NLE65 and 70). GROMACS does not have parameters for this aminoacid. You have two options: replace norleucine with leucine or add norleucine parameters. 
>>
> {: .solution}
{: .challenge}
