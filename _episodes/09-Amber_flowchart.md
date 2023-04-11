---
title: "An Overview of Information Flow in AMBER"
teaching: 30
exercises: 5
questions:
- "What files are needed to run an MD simulation with AMBER?"
objectives:
- "Learn what types of files are needed to run an MD simulation with AMBER"
keypoints:
- "To run an MD simulation with AMBER 3 files are needed: an input file, a parameter file, and a file describing coordinates/velocities . "
---
### Introduction
In this lesson we will go through the steps of setting up a fully solvated protein system for simulation with AMBER/NAMD/OPENMM. There are many commercial programs and interactive graphical interfaces designed to assist with system preparation. While these tools are easy to use and don't require as much learning efforts as command line tools, they are offer only a limited functionality, and most importantly results obtained with WEB/GUI tools are not reproducible and prone to human error. When a user needs to generate multiple systems with different proteins or membrane compositions or requires different starting configurations, relying on GUI/WEB interface becomes a challenge, as the process becomes time-consuming. Therefore, we will focus on system preparation using only scriptable command line driven tools. The emphasis of this lesson is to expose you to the various methods that can be used to create a reproducible molecular modeling workflow by automating preparation and simulation steps. One of the advantages of such approach is that once a workflow script have been developed it can be easily modified for other systems or conditions (for example if an updated version of pdb file is released, you can prepare a new simulation system with a single click).

### Information flow in AMBER
Simulation workflow in MD usually involves three steps: system preparation, simulation, and analysis. Let's take a closer look at these steps.
<div class="mermaid" style="height: 30%">
flowchart TB
style Q stroke-dasharray: 5 5
style T stroke-dasharray: 5 5
subgraph "Prepare"
    A(["PDB files"]) ==> |Load <br/>coordinates| C{TLEAP, <br/>XLEAP}
    B([FF files]) --> |Load <br/>parameters|C
    H([LEaP commands])-->|Build <br/>simulation system|C
end
subgraph "Run MD"
    C-->|Save <br/>topology|D(["prmtop"])
    C==>|Save <br/>coordinates|E(["inpcrd"])
    E==>|Load <br/>coordinates|F
    Q([disang])-.->|Load <br/>NMR restraints|F
    T([groupfile])-.->|Setup <br/> multiple  <br/> simulations|F
    G(["mdin"])-->|Load <br/>run options|F
    D-->|Load <br/>topology|F{SANDER, <br/>PMEMD}
    F==>|Save <br/>restart|P([restrt])
    P==>|Load <br/>restart|F
end   
subgraph Analyze 
    N(["CPPTRAJ <br/>commands"])-->|Run <br/>analysis|J 
    F===>|Save <br/>trajectory|I([mdcrd, mdvel])
    F--->|Print <br/>energies|L([mdout,<br/> mdinfo])
    L-->|Load <br/>energies|J
    I==>|Load <br/>frames|J{CPPTRAJ, <br/> PYTRAJ}
    I==>|Load <br/>frames|K{MMPBSA}
    R([MM-PBSA <br/>commands])-.->|Compute <br/> PB energy|K
    J-->S[Results]
    K--->S((Results))
end
    D---->|Load <br/>topology|J
    linkStyle 1,3,9,21 fill:none,stroke:blue,stroke-width:3px;
    linkStyle 2,8,12 fill:none,stroke:red,stroke-width:3px;   
</div>
### System Preparation
AMBER package includes two utilities for simulation preparation: tLEaP and xLEaP. The command line version (tLeap) is very efficient for executing scripts automatically. The graphical version (xLEaP) is ideal for building simulation systems interactively and visualizing them. The first step in creating a simulation is to load coordinates in PDB format and a Force Field of your choice. Adding ions, setting up periodic boundary conditions, adding solvent are the following steps. As needed, you can modify the system, for example, by adding extra bonds, changing charges, or modifying some amino acids' ionization states. Finish by creating a system topology and coordinate files. You can perform each of these actions by executing LEaP commands interactively or importing them from a file.

### Simulation
A package of AMBER includes two MD engines: SANDER and PMEMD. Serial and parallel versions are available for both programs. AMBER also provides a GPU-accelerated version of PMEMD. Along with the topology file and the coordinate or resume file, you will also need an input file describing the parameters of your simulation, such as integrator, thermostat, barostat, time step, cut-off distance, etc. It's also possible to use bond and distance constraints derived from NMR experiments. PMEMD and SANDER can run multiple simulations, including replica exchange and constant pH. Parameters for such simulations are contained in a special groupfile. A simulation program will save energy components and MD trajectories. 

### Analysis
Analysis of simulation output files can be done using CPPTRAJ, PYTRAJ, MMPBSA and other tools included in the AMBER package.

