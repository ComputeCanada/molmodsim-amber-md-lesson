---
title: "Hands-on 5: Generating topologies and parameters for small molecules."
teaching: 30
exercises: 5
questions:
- "How to parameterize small molecules?"
objectives:
- "Parameterize hexanediol"
keypoints:
- " "
---
## Generating topologies and parameters for small molecules.

### ANTECHAMBER

#### Download hexanediol 3D structure

- Find hexanediol (compound 147023) in [PubChem](https://pubchem.ncbi.nlm.nih.gov) and download 3D SDF file. 
- Convert SDF to PDB using Open Babel and change residue name: 

~~~
cd ~/workshop_amber/example_05
module purge
module load StdEnv/2023 openbabel
obabel Conformer3D_COMPOUND_CID_147023.sdf -O hexanediol.pdb
sed  "s/UNL/HEZ/g" hexanediol.pdb > HEZ.pdb
~~~
{: .language-bash}

#### Make mol2 file with antechamber

~~~
module purge
module load StdEnv/2023 ambertools/23
antechamber -i HEZ.pdb -fi pdb -o HEZ.mol2 -fo mol2 -c bcc -s 2
~~~
{: .language-bash}

The file HEZ.mol2 contains the definition of our HEZ residue including connectivity, all of the charges and atom types. 

#### Run parmchk2 to find out if there are any missing force field parameters

~~~
parmchk2 -i HEZ.mol2 -f mol2 -o HEZ.frcmod
~~~
{: .language-bash}

If it can antechamber will fill in these missing parameters by analogy to a similar parameter.

#### Create the library file for HEZ using tleap:

~~~
source leaprc.gaff2
HEZ=loadmol2 HEZ.mol2
check HEZ
saveoff HEZ hez.lib 
~~~
{: .leap}

Of course point charges are not very accurate because they are derived using semi-empirical method, but antechamber can also use results of gaussian QM calculations.  

#### Deriving accurate point charges
User can derive charges using RESP and supply them in mol2 file. 

[Electrostatic Parameterization with py_resp.py](https://ambermd.org/tutorials/basic/tutorial19/index.php)

1. QM Geometry Optimization (gaussian)
2. Electrostatic Potential Calculation (gaussian)
3. Convert the Gaussian ESP data format for PyRESP (ambertools:espgen)
4. Generate input for py_resp.py (ambertools:pyresp_gen.py)
4. RESP Parameterization (ambertools:py_resp.py)


### Charge derivation methods

Activity 1: Derive RESP, CM5 and AM1-BCC2 (sqm) charges and compare them   
Activity 2: Compare binding  free energy calculated using different charge sets.

[Comparison of Charge Derivation Methods Applied to Amino Acid Parameterization](https://pubs.acs.org/doi/10.1021/acsomega.8b00438). - Derivation does not matter much for aminoacids?

[Molecular Insights into the Covalent Binding of Zoxamide to the β-Tubulin of Botrytis cinerea](https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00911) - Some ligands for exercise on parameterization: carbendazim (CBZ), diethofencarb (DEF), zoxamide (ZOX)   
- Retrieve from Pubmed, Optimize using QUICK (B3LYP/6-31G*)


### AnteChamber PYthon Parser interfacE (ACPYPE)

[ANTECHAMBER](https://www.sciencedirect.com/science/article/abs/pii/S1093326305001737?via%3Dihub), a module of the AmberTools package, is the main tool for creating topological parameters in AMBER force fields. It can be used to generate topologies for most organic molecules. 

[ACPYPE - AnteChamber PYthon Parser interfacE](https://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-5-367). 

- Simplifies and automates generation of topology and parameters in different formats for different molecular mechanics programs.
- Translates force field to GROMACS 


[Installing and using ACPYPE on the Aliiance clusters](https://docs.alliancecan.ca/wiki/ACPYPE)

~~~
module load  python  openbabel
virtualenv ~/venv_acpype
source ~/venv_acpype/bin/activate
pip install --no-index acpype
~~~
{: .language-bash}


Create force field files:  
~~~
acpype -i HEZ.pdb -n 0
~~~

## Free energy calculations
### MMPBSA
Prepare tolopogies:
~~~
ante-MMPBSA.py -p ../../start.prmtop -s '!(:214-456,669-1029)' -c complex.prmtop
ante-MMPBSA.py -p complex.prmtop -n ':1-243' -l ligand.prmtop -r receptor.prmtop 
~~~

~~~
Input file for running PB and GB in serial
&general
   startframe=0, endframe=10, interval=1,
   keep_files=2, verbose=1, use_sander=1,
   strip_mask=!(:214-456,669-1029),
   ligand_mask=:1-243, receptor_mask=:244-604,
/
&gb
  igb=2, saltcon=0.150,
/
~~~

~~~
#!/bin/bash
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=3:00:00
module purge
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 ambertools/23
MMPBSA.py.MPI -O -i mmpbsa.in -o FINAL_RESULTS_MMPBSA.dat -sp ../../start.prmtop -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop -y ../../../min/min.dcd
~~~

#### Memory requirements

Generalized Born: 3N (atom positions) + 2N (atom parameters) + data structures for evaluating the full energy. GB memory requirements should scale more or less linearly with the number of atoms in the system

Normal modes: slightly more than (3N * 3N)/2 to store Hessian matrix + data structures for evaluating the full energy. The major expense here is the N^2 scaling of the Hessian storage. 

Poisson Boltzmann: the memory is dominated primarily by the grid, it depends strongly on the grid spacing. 

3D-RISM: also requires a grid. The 3D-RISM grid needs to be denser than the corresponding grid for PB, so RAM requirements for 3D-RISM are typically a bit higher than PB.

