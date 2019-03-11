# PolyOnTheFly

PolyOnTheFly is a Molecular Dynamics code (MD) written in FORTRAN 90/95.
PolyOnTheFly is a project by Matteo Bonfanti (@matbonfanti) @ [CDGT Group](http://users.unimi.it/cdtg)

### Capabilities
* Ring Polymer Molecular Dynamics
* Equilibration with Langevin Dynamics 
* Symplectic propagators: [Velocity-Verlet](http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet), [Vanden-Eijnden Ciccotti Integrator](http://dx.doi.org/10.1016/j.cplett.2006.07.086), [Path Integral Langevin Equation (PILE)](http://dx.doi.org/10.1063/1.3489925)
* Integrated with [SIESTA](http://www.icmab.es/siesta/) and [DFTB+](http://www.dftb-plus.info/) for on-the-fly computation of the forces
* MPI parallelization over the trajectories

### Acknowledgments
We thanks CINECA for the support in developing the MPI parallelization (project Si-RPMD within the LISA initiative)

