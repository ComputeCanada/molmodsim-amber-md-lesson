---
title: "Information flow and files needed to run an MD simulation with AMBER"
teaching: 30
exercises: 5
questions:
- "What files are needed to run an MD simulation with AMBER?"
objectives:
- "Learn what types of files are needed to run an MD simulation with AMBER"
keypoints:
- "To run an MD simulation with AMBER 3 files are needed: an input file, a parameter file, and a file describing coordinates/velocities . "
---

## AMBER MD simulation engines.
Amber package includes two MD endines: sander and pmemd.

### SANDER
Sander is a free simulation engine distributed with AmerTools package. It uses the MPI (message passing interface) to communicate between CPUs (by CPU we mean CPU core, not a physical chip). Sander uses a replicated data structure, in which each CPU computes interactions between certain atoms, but all CPUs know the coordinates of all atoms. At each step, processors compute a portion of the potential energy and corresponding gradients for a set of atoms assigned to each of them. A global portion of code then sums the force vector, and sends the result to each of the CPUs so that each processor gets the full force vector components. The processors then perform a molecular dynamics update step for the assigned to them atoms, and communicate the updated positions to all processors, in preparation for the next molecular dynamics step.

This model provides a convenient programming environment, but the main problem is that the communication required at each step grows with the number of processors limiting parallel scaling. 

### PMEMD
PMEMD is an extensively revised version of SANDER. Many optimizations were made to improve both single-processor performance and parallel scaling. To improve performance, PMEMD communicates to each processor only the coordinate information necessary for computing the pieces of the potential energy assigned to it. This code, however, does not support all of the options found in sander.

[Improving the Efficiency of Free Energy Calculations in the Amber Molecular Dynamics Package](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3811123/)

#### GPU - accelerated PMEMD
GPU - accelerated PMEMD version of PMEMD (pmemd.cuda) leverages NVIDIA GPUs to perform MD simulations. Pmemd.cuda executable is available only in the commercial AMBER package. It is significantly faster than CPU version. Pmemd.cuda achieves high simulation speed by executing all calculations on a single GPU within its memory. This approach eliminates the bottleneck of moving data between CPU and GPU and allows very efficient GPU usage for systems of practical sizes. Modern GPUs are so fast that communication overhead between GPUs does not allow for parallel scaling to multiple GPUs. Thus the parallel version (pmemd.cuda.MPI) is useful only for some specific types of simulations such as thermodynamic integration and replica-exchange where simulations in different lambda windows (different replicas) can be executed independently on different GPUs. In these simulations the number of lambda windows must be a multiple of the available GPU because individual lambda windows cannot span multiple GPUs but one GPU can run multiple windows.   

### Multi-sander and multi-pmemd simulations
Multi-sander and multi-pmemd are wrappers around parallel versions of these programs. The wrappers allow to run multiple copies of simulations within a single parallel run. 

### Summary of AMBER MD engine executables: 

|  Verion     | sander |  pmemd|  
|:---|:---|:---|
|Serial | sander  | pmemd|      |
|Parallel | sander.MPI  | pmemd.MPI|   
|Single GPU, default| -  | pmemd.cuda -> pmemd.cuda_SPFP|
|Single GPU, single precision| -  | pmemd.cuda_SPFP|   
|Single GPU, double precision| -  | pmemd.cuda_DPFP|
|Multi GPU, default| -  | pmemd.cuda.MPI -> pmemd.cuda_SPFP.MPI|
|Multi GPU, single precision| -  | pmemd.cuda_SPFP.MPI|   
|Multi GPU, double precision| -  | pmemd.cuda_DPFP.MPI|


### Information flow in AMBER
<div class="mermaid" style="height: 30%">
flowchart TD
subgraph "Prepare"
    A(["PDB files"]) ==> |Load coordinates| C{TLEAP, XLEAP}
    B([FF files]) --> |Load parameters|C
    H([LEaP commands])-..->|Build simulation system|C
end
subgraph "Run"
    C-->|Save topology|D(["prmtop"])
    C==>|Save coordinates|E(["inpcrd"])
    D-->|Load topology|F{SANDER, PMEMD}
    E==>|Load coordinates|F
    G(["mdin"])-.->|Load simulation parameters|F
end   
subgraph "Analyze"
    N(["CPPTRAJ commands"])-.->|Analyze data|J 
    F===>|Save trajectory|I([mdcrd])
    F--->|Print energies|L([mdout,mdinfo])
    I==>|Load frames|J{CPPTRAJ, PYTRAJ}
    I===>|Load frames|K{MM-PBSA}
    L-->|Load energies|J
    R([MM-PBSA commands])-.->|Compute PB energy|K
end
    F==>|Save restart|P([restrt])
    P==>|Load restart|F
    D---->|Load|J
</div>


## AMBER input files
&cntrl
&ewald
&qmmm




CPhMD
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6934042/
https://pubs.acs.org/doi/10.1021/acs.jcim.8b00227
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8062021/ charmm
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7795291/