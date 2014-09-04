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
   
!    !> Gamma of the relaxation during dynamics (its meaning depends on the bath representation)
!    REAL    :: DynamicsGamma             

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

!    ! Averages computed during propagation
!    REAL, DIMENSION(:), ALLOCATABLE      :: PositionCorrelation       !< Position correlation function
!    REAL, DIMENSION(:,:), ALLOCATABLE    :: AverageCoord              !< Average i-th coordinate vs time
!    INTEGER, PARAMETER                   :: NrOfEnergyAverages = 5    !< Nr of average values computed by the function EnergyAverages
!    REAL, DIMENSION(:,:), ALLOCATABLE    :: AverageE                  !< Traj averages of the energy in time

   
!=============================================================================================================

   ! INFORMATION ON THE SYSTEM

   INTEGER :: NAtoms                                     !< Nr of atoms of the system
   INTEGER :: NDim                                       !< Nr of dimension of the system
   REAL, DIMENSION(:), ALLOCATABLE :: MassVector         !< Vector with the masses of the system

   
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

   REAL, DIMENSION(:), ALLOCATABLE :: CentroidPos  !< Position of the centroid of the RP
   REAL, DIMENSION(:), ALLOCATABLE :: CentroidVel  !< Velocity of the centroid of the RP

!=============================================================================================================


! CONTAINS

!    !> Variable to define which kind of calculation is required 
!    INTEGER :: RunType
!    INTEGER, PARAMETER :: EQUILIBRIUM     = 1,  & ! Equilibrium calculation with H already adsorbed
!                          HARMONICMODEL   = 2,  & ! Test the parameters with a 1D harmonic model 
!                          RELAXATION      = 3,  & ! Relaxation dynamics of a CH bound state, with the bath at 0K
!                          SCATTERING      = 4,  & ! Scattering calculation with H coming from gas-phase
!                          RPMD_RELAXATION = 5,  & ! Relaxation dynamics with RING POLYMER DYNAMICS
!                          RPMD_EQUILIBRIUM = 6,  & ! Equilibrium simulation with Ring Polymer MD
!                          POTENTIALPRINT  = 10    ! Static analysis of the potential 
! 
!    !> Variable to set the print level of the calculation
!    INTEGER :: PrintType
!    INTEGER, PARAMETER :: DEBUG       = 3, &   ! fully detailed information about the trajs
!                          FULL        = 2, &   ! files to plot the make animations, averages for each traj
!                          MINIMAL     = 1      ! minimal level of output, only final averages
! 
! 


!    !> Subroutine to check the availability of a given runtype option
!    SUBROUTINE CheckRunType( IntNr )
!       IMPLICIT NONE
!       INTEGER, INTENT(IN) :: IntNr
!       LOGICAL :: Check 
! 
!       Check = ( IntNr /= EQUILIBRIUM .AND. &
!                 IntNr /= HARMONICMODEL .AND. &
!                 IntNr /= RELAXATION .AND. &
!                 IntNr /= SCATTERING .AND. &
!                 IntNr /= RPMD_RELAXATION .AND. &
!                 IntNr /= RPMD_EQUILIBRIUM .AND. &
!                 IntNr /= POTENTIALPRINT )
!       CALL ERROR( Check, " SharedData.CheckRunType: Invalid RunType option " )
!    END SUBROUTINE CheckRunType
! 
!    !> Subroutine to check the availability of a given printtype option
!    SUBROUTINE CheckPrintType( IntNr )
!       IMPLICIT NONE
!       INTEGER, INTENT(IN) :: IntNr
!       LOGICAL :: Check 
! 
!       Check = ( IntNr /= DEBUG .AND. &
!                 IntNr /= FULL .AND. &
!                 IntNr /= MINIMAL  )
!       CALL ERROR( Check, " SharedData.CheckPrintType: Invalid PrintType option " )
!    END SUBROUTINE CheckPrintType
! 
!    !> Subroutine to check the availability of a given printtype option
!    SUBROUTINE CheckBathType( IntNr )
!       IMPLICIT NONE
!       INTEGER, INTENT(IN) :: IntNr
!       LOGICAL :: Check 
! 
!       Check = ( IntNr /= SLAB_POTENTIAL .AND. &
!                 IntNr /= NORMAL_BATH .AND. &
!                 IntNr /= CHAIN_BATH .AND. &
!                 IntNr /= DOUBLE_CHAIN .AND. &
!                 IntNr /= LANGEVIN_DYN  )
!       CALL ERROR( Check, " SharedData.CheckBathType: Invalid BathType option " )
!    END SUBROUTINE CheckBathType

END MODULE SharedData
