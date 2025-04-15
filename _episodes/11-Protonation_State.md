---
title: "Assigning Protonation States to Residues in a Protein"
teaching: 30
exercises: 10
questions:
- "Why titratable aminoacid residues can have different protonation states?"
- "How to determine protonation state of a residue in a protein?"
- "What are the weaknesses of fixed protonation state simulations?"
objectives:
- "Understand why it is necessary to assign the correct protonation state"
- "Learn how to determine protonation state of a protein"
- "Learn how to assign protonation state to a residue"
keypoints:
- "Assigning correct protonation states of aminoacids in proteins is crucial for realistic MD simulations"
- "Conformational changes in proteins may be accompanied by changes in protonation pattern."
---


### It is important to consider amino acid protonation states

{% if site.show_comments %}
Setting up a simulation system requires assigning protonation states and, possibly, tautomers to the HIS residues. The protonation states of titratable amino acids (Arg, Lys, Tyr, Cys, His, Glu, Asp) depend on the local micro-environment and pH. A highly polar microenvironment will stabilize the charged form, while a less polar microenvironment will favor the neutral form. At physiological pH, TYR, LYS, CYS, and ARG are almost always in their standard protonation states, while GLU, ASP, and HIS can be in non-standard forms.  

The protonation pattern of proteins is crucial for their catalytic function and structural stability. For a classic MD simulation, protonation states must be determined and assigned before the simulation can begin. For realistic MD simulations, it is essential to assign the right protonation states. The wrong protonation states can have dramatic effects and invalidate the results. Numerous MD simulation studies have demonstrated the importance of protein protonation states [[1](http://doi.org/10.1529/biophysj.105.059329), [2](http://doi.org/10.7554/eLife.16616), [3](https://doi.org/10.1016/j.cplett.2018.12.039), [4](https://doi.org/10.1021/jacs.9b06064), [5](https://doi.org/10.1016/j.dib.2016.07.040)].

{% else %}
- The protonation pattern of proteins is crucial for their catalytic function and structural stability. 
-  Numerous MD simulation studies have demonstrated the importance of protein protonation states [[1](http://doi.org/10.1529/biophysj.105.059329), [2](http://doi.org/10.7554/eLife.16616), [3](https://doi.org/10.1016/j.cplett.2018.12.039), [4](https://doi.org/10.1021/jacs.9b06064), [5](https://doi.org/10.1016/j.dib.2016.07.040)]. 
{% endif %}

### How to Determine Protonation States of Residues in a Protein?

{% if site.show_comments %}
For predicting the pKa values of protein residues, several web servers and standalone programs are available. The underlying methods used by all of them fall into two broad categories: empirical (such as propka) and Poisson-Boltzmann methods. The disadvantage of the empirical method is that it works well only for proteins similar to those in the training set. It is not guaranteed to work well for all proteins. Programs implementing Poisson-Boltzmann continuum electrostatics have a solid physical basis. Most of them, however, only sample a single protein conformation, which is the primary source of inaccuracy. One exception is MCCE, which samples sidechain rotamers, giving a more precise picture of coupled ionization and position changes. MCCE also has limitations because it does not sample backbone conformations. The most rigorous method is constant pH simulations in which the ionization states are dynamically adjusted in MD simulations. It is a computationally expensive method, but a very efficient GPU implementation of constant pH molecular dynamics has recently been developed and [implemented in AMBER](https://pubs.acs.org/doi/10.1021/acs.jcim.9b00754).

- [*H++*](http://newbiophysics.cs.vt.edu/H++/) server. *H++* calculations are based on the classical continuum electrostatic model (a continuous medium where solvent molecules are replaced with the average properties of the solvent). The *H++* program uses AmberTools modules to preprocess a PDB file, and it is able to generate basic topology and coordinate files for MD simulations in AMBER format. *H++* is a single-configuration method, and it automatically selects the "A" conformation without any user intervention. Details of the methodology can be found in the [reference](https://doi.org/10.1093/nar/gks375).
- [*PROPKA3.0*](https://github.com/jensengroup/propka-3.0) is the empirical pKa prediction software. 
- [PlayMolecule-ProteinPrepare](https://www.playmolecule.org/proteinPrepare/) provides complete system preparation. Among the preparation steps are determination of protonation states with PROPKA 3.1, the addition of missing atoms, and optimization of the H-bond network with PDB2PQR 2.1.
- [*PDB2PQR*](https://server.poissonboltzmann.org/pdb2pqr) server. *PDB2PQR* solves Poisson-Boltzmann equation using the APBS solver. 
- [*MCCE*](https://sites.google.com/site/mccewiki/). For more rigorous calculations, try *MCCE* program. *MCCE* uses the same classical continuum electrostatic approach as H++. Besides, *MCCE* calculations consider protein conformational degrees of freedom giving a more accurate picture of coupled ionization and position changes. Taking into account conformational flexibility improves *pKa* prediction significantly. For details of the methodology, see [[ref]](https://doi.org/10.1002/jcc.21222).
- [*PKAD*](http://compbio.clemson.edu/pkad) database. *PKAD* is a database of experimentally determined pKa values. This database currently contains *pKa* values for residues in 157 wild-type and 45 mutant proteins (https://doi.org/10.1093/database/baz024).
 
{% else %}
For predicting the pKa values of protein residues, several web servers and standalone programs are available. 

- [*H++*](http://newbiophysics.cs.vt.edu/H++/). Continuum electrostatics model
- [*PROPKA3.0*](https://github.com/jensengroup/propka-3.0). Empirical pKa prediction. 
- [PlayMolecule-ProteinPrepare](https://www.playmolecule.org/proteinPrepare/). Based on PROPKA3. 
- [*PDB2PQR*](https://server.poissonboltzmann.org/pdb2pqr). Solves Poisson-Boltzmann equation. 
- [*MCCE*](https://sites.google.com/site/mccewiki/). Takes into account conformational flexibility. More accurate and more challenging to use. 
- [*PKAD*](http://compbio.clemson.edu/pkad). A database of experimentally determined pKa values. 
{% endif %}

There is no perfect *pKa* prediction method. Deviations from experimental values can sometimes be significant. The best practice is to compare the results obtained from different techniques and, if possible, to use experimentally measured values.

> ## Calculating *pKa*'s
>1. Calculate *pKa*'s of residues in the PDB entry 1RGG using *H++* server.
>2. What protonation states of Asp79 and His53 are appropriate for simulation at *pH* 6?
>3. Repeat calculations using *PDB2PQR* server and compare the results.
>4. Compare calculated *pKa*'s with the experimental. How accurate are the predicted *pKa* values?
>
What protonation states are appropriate for simulating Asp79 and His53 at pH 6?
>> ## Solution
>>If pKa > pH the probability that the residue is protonated is > 50%, and we use the protonated form.  
>>If pKa < pH the probability that the residue is protonated is < 50% and we use the deprotonated form.
>>
>>ASP79 has pKa 7.2 (experimental 7.37), it is protonated at pH 6 and we rename it to ASH  
>>HIS53 has pKa 8.3 (experimental 8.27), it is also protonated at pH 6 and we rename it to HIP
> {: .solution}
{: .challenge}

## How to select protonation state of a residue?

### Assigning protonation states in structure files
{% if site.show_comments %}
A more consistent and convenient way to select the desired form of amino acid is to change its name in the structure file before generating a topology. The neutral states of LYS, ASP, and GLU can be chosen by renaming them  LYN, ASH, and GLH, respectively.  You can select the appropriate form of HIS by renaming HIS to HIE (proton on NE2), HID (proton on ND1), or HIP (both protons).
{% endif %}

-  To change the form of an amino acid, change its name in the structure file

```
LYS (+1) - LYN  (0)  
ASP (-1) - ASH  (0)  
GLU (-1) - GLH  (0)  
HIS  (0) - HIE  (0)  
HIS  (0) - HID  (0)  
HIS  (0) - HIP (+1)    
```

### Selecting protonation states with check_structure.
Let's change ASP20 and ASP26 in the file 1ert.pdb to the neutral form ASH. 
~~~
cd ~/workshop_amber/example_02
check_structure -i 1ert.pdb -o 1ert_protonated.pdb \
command_list --list \
"add_hydrogen --add_mode list --list A:asp20ash,A:asp26ash; \
water --remove yes"
~~~
{: .language-bash}


You can verify that residues are changed by grepping ASH.

{% if site.show_comments %}
You can read commands from a file:
~~~
check_structure -i 1ert.pdb -o 1ert_protonated.pdb command_list --list check.in
~~~
{: .language-bash}

~~~
# File: check.in
add_hydrogen --add_mode list --list A:asp20ash,A:asp26ash
water --remove yes
~~~
{: .file-content}
{% endif %}

>## Selecting protonation states with *VMD*.
>
>Change ASP20 and ASP26 in the file 1ert.pdb to the neutral form ASH and remove water.
>>## Solution
>>~~~
>>ml StdEnv/2023 vmd
>>vmd
>>~~~
>>{: .language-bash}
>>
>>~~~
>>mol new 1ert.pdb
>>set s [atomselect top "resid 20 26"]
>>$s set resname ASH
>>set s [atomselect top "protein"]
>>$s writepdb 1ert_protonated.pdb
>>quit
>>~~~
>>{: .vmd}
>{: .solution}
{: .challenge}

### Assigning protonation states with the *GROMACS pdb2gmx* module
{% if site.show_comments %}
Protonation states can be assigned using the *GROMACS pdb2gmx* program. By default, *pdb2gmx* will select charged forms of LYS, ASP, or GLU. For HIS, it will try to place the proton on either ND1 or NE2 based on an optimal hydrogen bonding conformation. You can override the default behavior and select protonation manually using options  -lys, -asp, -glu, -his.
{% endif %}

- By default, *pdb2gmx* will select charged forms of LYS, ASP, or GLU
- For HIS, it will try to place the proton optimally

{% if site.show_comments %}
The downside of this method is that it can not be scripted. The manual selection is cumbersome because *pdb2gmx* will be prompting to select the protonation state for each of the residues in a PDB file. Besides, *pdb2gmx* changes residue names only in the output topology file. Residue names are left unchanged in the output .pdb and .gro files.  
{% else %}
Downsides:
- Interactive 
- Residue names are changed only in the topology
{% endif %}

### Limitations of Fixed Protonation State Simulations
{% if site.show_comments %}
The use of constant protonation states in molecular dynamics simulations has its disadvantages. When conformational changes occur in natural systems, protonation patterns often change as well. In molecular dynamics simulations with fixed states, these processes are not coupled, making it difficult to fully understand proton-coupled conformational dynamics. Consider using constant pH simulations if proton-coupled dynamics are essential to your research. MD with constant pH is currently implemented in both AMBER and NAMD. At the moment, it is not officially implemented in *GROMACS*, but a modified version is available [[6](https://pubs.acs.org/doi/10.1021/ct200061r)].
{% else %}
- Difficult to understand proton-coupled conformational dynamics.
- Consider using constant pH simulations if proton-coupled dynamics are essential to your research.
{% endif %}

### Combining all structure preparation steps in one check_structure script
~~~
cd ~/workshop_amber/example_03
check_structure -i 1rgg.pdb -o 1RGG_chain_A_prot.pdb \
    command_list --list "\
        chains --select A;\
        add_hydrogen --list A:his53hip,A:asp79ash --add_mode list;\
        altloc --select A5:B,A54:B,A6:A,A13:A,A42:A,A85:A,A91:A;\
        getss --mark all;\
        ligands --remove all;\
        water --remove yes"
~~~
{: .language-bash}

You can also save commands in a file and pass it as an argument to "command_list --list":

~~~
check_structure -i 1rgg.pdb -o 1RGG_chain_A_prot.pdb command_list --list prep_1rgg.chk
~~~
{: .language-bash}

#### References

1. [Constant-pH Molecular Dynamics Simulations for Large Biomolecular Systems](https://pubs.acs.org/doi/10.1021/acs.jctc.7b00875)

2. [ GPU-Accelerated Implementation of Continuous Constant pH Molecular Dynamics in Amber: pKa Predictions with Single-pH Simulations](https://pubs.acs.org/doi/10.1021/acs.jcim.9b00754)

> ## Combining all structure preparation steps in one *VMD* script
> Combine all previous steps together and create *VMD* script to prepare MD simulation system for the hydrolaze PDB structure 1RGG. The script should perform the following steps:
>
> 1. Select molecule A
> 2. Remove non-protein molecules
> 3. Select location 'B' for residues 5, 54 and location 'A' for all other residues with alternative locations
> 4. Protonate Asp79 and His53
> 5. Rename CYS 7 and 96 into CYX (cross-linked cystein)
> 6. Save the resulting structure as 1RGG_chain_A_prot.pdb
>
>>## Solution
>>~~~
>> cd ~/workshop_amber/example_03
>>~~~
>>{: .language-bash}
>> Save the following commands in a file,  e.g. prep_1RGG.vmd
>> ~~~
>># Load 1rgg.pdb into a new (top) molecule
>>mol pdbload 1rgg
>># Select and save all chain A protein atoms
>>set s [atomselect top "protein and chain A"]
>>$s writepdb 1RGG_chain_A.pdb
>># Delete the top molecule
>>mol delete top
>># Load chain A into a new molecule 
>># Loading only one chain will simplify selections commands 
>>mol new 1RGG_chain_A.pdb
>># Protonate ASP79
>>set s [atomselect top "resid 79"]
>>$s set resname ASH
>># Protonate HIS53
>>set s [atomselect top "resid 53"]
>>$s set resname HIP
>># Rename cross-linked cysteins
>>set s [atomselect top "resid 7 96"]
>>$s set resname CYX
>># Select the base and the alternate locations
>>set s [atomselect top "(altloc '') or (altloc A and resid 6 13 42 85 91) or (altloc B and resid 5 54)"]
>># Save the selection
>>$s writepdb 1RGG_chain_A_prot.pdb
>>quit
>>~~~
>>{: .vmd}
>> Execute the script
>>~~~
>> vmd -e prep_1RGG.vmd
>>~~~
>>{: .language-bash}
> {: .solution}
{: .challenge}


