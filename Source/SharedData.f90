!***************************************************************************************
!*                              MODULE SharedData
!***************************************************************************************
!
!>  \brief     Common data
!>  \details   This module include the data to be shared by the
!>             all the other modules of the code
!
!***************************************************************************************
MODULE SharedData
   USE ErrorTrap
   ! Use the following modules to define the necessary derived data types
   USE RandomNumberGenerator
   USE ClassicalEqMotion
   USE FFTWrapper
   
   IMPLICIT NONE

   PUBLIC          ! This module contain data that is supposed to be shared

!=============================================================================================================

   ! PARAMETERS 
   

!=============================================================================================================

   ! VARIABLES SET FROM INPUT

   ! variables of the ring polymer dynamics 
   INTEGER :: NBeads                    !< Nr of replicas of the system + bath 
   REAL    :: BeadsForceConst           !< Force constant of the harmonic spring between the beads
   REAL    :: BeadsFrequency            !< Harmonic frequency of the spring between the beads

   ! Initial conditions of the bath
   REAL    :: Temperature               !< Temperature of the simulation

   ! Nr of trajectories
   INTEGER :: NrTrajs                   !< Nr of trajectories
   
   ! Info about output
   REAL    :: PrintTimeStep             !< time between writing the output
   
   ! Variables of the equilbration
   REAL    :: LangevinGamma             !< Langevin friction coefficient for the PILE thermalization
   REAL    :: EquilTotalTime            !< Total time of the equilibration
   REAL    :: EquilTimeStep             !< Time step of the equilibration
   INTEGER :: EquilNrSteps              !< Total nr of equilibration time step per trajectory
   INTEGER :: EquilPrintStepInterval    !< Step interval between each printing step during the equilibration
   
   ! Variables of the propagation
   REAL    :: TimeStep                  !< Time step for the integration of the classical EOM
   REAL    :: DynamicsTotalTime         !< Total time of the dynamics
   INTEGER :: NrSteps                   !< Total nr of time step per trajectory
   INTEGER :: PrintStepInterval         !< Step interval between each printing step during the dynamics

   ! Time evolution dataset
   TYPE(Evolution), SAVE :: MolecularDynamics   !< Propagate in micro/canonical ensamble to extract results
   TYPE(Evolution), SAVE :: InitialConditions   !< Propagate in microcanonical ensamble to generate init cond

   ! Transform to ring normal modes 
   TYPE(FFTHalfComplexType), SAVE :: RingNormalModes   !< Transform ring coords to normal modes (and viceversa)

   ! Internal state of the random number generator
   TYPE(RNGInternalState), SAVE :: RandomNr       !< Internal state of the random number generator

   ! Control additional output 
   LOGICAL, SAVE :: Out_VelDistrib        !< Print kinetic energy distribution for each trajectory
   INTEGER, SAVE :: Out_VelDistrib_nV     !< Number of intervals of the kinetic energy binning
   REAL, SAVE    :: Out_VelDistrib_DV     !< Energy interval of the kinetic energy binning

   
!=============================================================================================================

   ! INFORMATION ON THE SYSTEM

   INTEGER :: NAtoms                                     !< Nr of atoms of the system
   INTEGER :: NDim                                       !< Nr of dimension of the system
   REAL, DIMENSION(:), ALLOCATABLE :: MassVector         !< Vector with the masses of the system

!=============================================================================================================

   ! INFORMATION ON CURRENT TRAJECTORY STEP

   INTEGER :: iTraj                                      !< Number of current trajectory
   REAL    :: Time                                       !< Time of the current trajectory step
   INTEGER :: kStep                                      !< Printing steps counter

!=============================================================================================================

   ! POSITION, VELOCITY, ACCELERATION 

   REAL, DIMENSION(:), ALLOCATABLE, TARGET :: X    !< Position at given timestep
   REAL, DIMENSION(:), ALLOCATABLE, TARGET :: V    !< Velocity at given timestep
   REAL, DIMENSION(:), ALLOCATABLE, TARGET :: A    !< Acceleration at ginve timestep

!=============================================================================================================

   ! ISTANTANEOUS AVERAGES 

   REAL ::  PotEnergy                              !< Potential energy
   REAL ::  KinEnergy                              !< Total kinetic energy (in AIMD computed from the virial average)
   REAL, DIMENSION(:), ALLOCATABLE :: KinPerCoord  !< Kinetic energy per coordinate ( still virial )
   REAL ::  TotEnergy                              !< Total energy
   REAL ::  RPKinEnergy                            !< Kinetic energy of the ring polymer
   REAL ::  RPPotEnergy                            !< Potential energy of the ring polymer
   REAL ::  RPTotEnergy                            !< Total mechanical energy of the ring polymer

   REAL, DIMENSION(:), ALLOCATABLE :: CentroidPos  !< Position of the centroid of the RP
   REAL, DIMENSION(:), ALLOCATABLE :: CentroidVel  !< Velocity of the centroid of the RP

!=============================================================================================================

END MODULE SharedData
