!***************************************************************************************
!*                           PROGRAM PolyOnTheFly
!***************************************************************************************
!
!>  \brief     On-the-fly Ring Polymer Molecular Dynamics Code
!>  \details   The code implements Ring Polymer Molecular Dynamics on a given  \n
!>             system, and is thought to be interfaced with SIESTA to perform  \n
!>             on-the-fly computation of the atomic forces.
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             20 August 2014
!
!***************************************************************************************
!
!>  \par Updates
!>  \arg  6 November 2014 : new MPI version with master reading input files and 
!>                          broadcasting data across the nodes
!
!>  \todo     \arg Set normal termination of the MPI execution in case of error
!>            \arg Initial conditions for the RPMD should be improved      
!>                 
!***************************************************************************************

PROGRAM PolyOnTheFly
#include "preprocessoptions.cpp"
   USE InputField
   USE UnitConversion
   USE SharedData
   USE ClassicalEqMotion
   USE PotentialModule
   USE RandomNumberGenerator
   USE FFTWrapper
   USE MyLinearAlgebra
   USE OutputModule
   USE PeriodicBoundary
   
   IMPLICIT NONE

   ! FEW VARIABLES HERE, DEFINE MOST VARIABLES IN SharedData MODULE

   ! Variable to handle the command line
   INTEGER :: NArgs
   LOGICAL :: Help = .FALSE.
   
   ! Input file name, set from command line arguments
   CHARACTER(120) :: InputFileName
   ! Derived type to handle input data
   TYPE(InputFile) :: InputData
   ! Units of input data, defined from the input file
   INTEGER     :: InputLength, InputEnergy, InputMass, InputTime, InputTemp, InputFreq
      
   ! keep track of initial time and final time
   INTEGER, DIMENSION(8)    :: Time1, Time2

   ! MPI data
   INTEGER, DIMENSION(:), ALLOCATABLE :: NodeOfTraj
   INTEGER, DIMENSION(:), ALLOCATABLE :: NrOfTrajPerNode

   ! General integer counters
   INTEGER :: iCoord, iNode
   
   ! Output data
   INTEGER, DIMENSION(2) :: OutputSet0
   REAL, DIMENSION(3)    :: OutputSet1
   REAL, DIMENSION(4)    :: OutputSet2
   REAL, DIMENSION(4)    :: OutputSet3
   
   ! =========================================================================
  !                         CODE STRUCTURE
   ! =========================================================================
   
   ! (*) = ONLY THE MASTER NODE
   
   ! 1) INPUT SECTION
   !    read input data from file (*)
   !    process data and convert to internal units (*)
   !    syncronizing data across all the nodes
   
   ! NOTE: data about the potential are read and processed at the same time
   !       but the code is included in the potential module
   
   ! 2) SETUP SECTION
   !    set random number seed generator ( PARALLEL: 1 + MPI_ID )
   !    setup potential
   !    memory allocation for the trajectory
   !    setup propagation data for equilibration and microcanonical dynamics
   !    setup of istantaneous and global averages
   !    setup partitioning of the trajectories over the nodes
  
   ! CYCLE OVER THE TRAJECTORIES
   
      ! 3-i) EQUILIBRATION SECTION
   
      ! 4-i) RUN SECTION

   ! 5) CONCLUSION AND AVERAGES
   
   ! =========================================================================
   
#if defined(WITH_MPI)
   CALL MyMPI_Init()
   __MPI_OnlyMasterBEGIN 
   WRITE(*,"(/,A)") " Running application with MPI parallelization "
   WRITE(*,*)        "Executable compilation ", __DATE__, " ", __TIME__
   __MPI_OnlyMasterEND
   CALL MyMPI_Barrier()
   WRITE(*,"(A,A,A,A,A)") ' Process ',TRIM((NumberToString(__MPI_CurrentNrOfProc))), &
                            ' of ', TRIM(NumberToString(__MPI_TotalNrOfProcs)), ' is alive '
   CALL MyMPI_Barrier()
#else
   WRITE(*,"(/,A)") " Running serial application "
   WRITE(*,*)        "Executable compilation ", __DATE__, " ",__TIME__
#endif

   __MPI_OnlyMasterBEGIN
   PRINT "(/,     '                    ==============================')"
   PRINT "(       '                              PolyOnTheFly        ')"
   PRINT "(       '                    ==============================',2/)"
   PRINT "(       '                  Author: Matteo Bonfanti, August 2014  ')"
   PRINT "(       '         On-The-Fly Ring Polymer Molecular Dynamics code in FORTRAN 90/95   ',2/)"
   __MPI_OnlyMasterEND
      
   !*************************************************************
   !         COMMAND LINE ARGUMENT
   !*************************************************************

   ! Check and read from command line the input file name
   NArgs = COMMAND_ARGUMENT_COUNT()
   IF (NArgs<1) THEN
      Help = .TRUE.
   ELSE
      CALL GET_COMMAND_ARGUMENT( 1, InputFileName )
      IF ( trim(adjustl(InputFileName)) == "help" ) Help = .TRUE.
   ENDIF
   IF (Help) THEN ! Call help
      __MPI_OnlyMasterBEGIN
      WRITE(*,"(A)")   '   Launch this program as:'
      WRITE(*,"(A,/)") '   % PolyOnTheFly "InputFileName" '
      __MPI_OnlyMasterEND
      CALL MyMPI_Finalize()
      STOP
   ENDIF

   __MPI_OnlyMasterBEGIN
   CALL date_and_time (values=Time1)
   __MPI_OnlyMasterEND

   !*************************************************************
   !         INPUT SECTION
   !*************************************************************
   !            ### note for parallelization ###
   !     this section should be executed by the master
   !     process only when data is read and processed,
   !     it is distributed to all the slaves
   !*************************************************************

   __MPI_OnlyMasterBEGIN
   ! Open and read from input file the input parameters of the calculation
   CALL OpenFile( InputData, InputFileName )

   ! Read the Input data
   CALL Input( )
   ! Convert to internal units
   CALL ConvertUnitsAndProcessData( )
   
   ! Read the potential data from input file
   CALL InputPotential( InputData, InputUnits )
   
   ! close input file
   CALL CloseFile( InputData )

   WRITE(*, 902) NrTrajs, NBeads, Temperature*TemperatureConversion(InternalUnits, InputUnits), TemperUnit(InputUnits),  &
                 BeadsFrequency*FreqConversion(InternalUnits, InputUnits), FreqUnit(InputUnits)

   WRITE(*, 903) EquilTotalTime*TimeConversion(InternalUnits, InputUnits), TimeUnit(InputUnits),    &
                 EquilTimeStep*TimeConversion(InternalUnits, InputUnits), TimeUnit(InputUnits),     &
                 PrintTimeStep*TimeConversion(InternalUnits, InputUnits), TimeUnit(InputUnits),     &
                 (1./LangevinGamma)*TimeConversion(InternalUnits, InputUnits), TimeUnit(InputUnits) 

   WRITE(*, 904) DynamicsTotalTime*TimeConversion(InternalUnits, InputUnits), TimeUnit(InputUnits), &
                 TimeStep*TimeConversion(InternalUnits, InputUnits), TimeUnit(InputUnits),          &
                 PrintTimeStep*TimeConversion(InternalUnits, InputUnits), TimeUnit(InputUnits) 
   __MPI_OnlyMasterEND
         
   902 FORMAT(" * General Information about the simulation     ",             /, &
              " * Nr of trajectories:                          ", I10,        /, & 
              " * Nr of replicas in the ring polymer dynamics: ", I10,        /, &
              " * Temperature:                                 ", F10.1,1X,A, /, &
              " * Frequency of the force between the beads     ", F10.4,1X,A, /   )

   903 FORMAT(" * Relaxation variables                         ",             /, &
              " * Equilibration time:                          ", F10.4,1X,A, /, & 
              " * Equilibration time step:                     ", F10.4,1X,A, /, & 
              " * Equilibration print step:                    ", F10.4,1X,A, /, & 
              " * Langevin relaxation time of the system       ", F10.4,1X,A, /   ) 

   904 FORMAT(" * Microcanonica propagation variables          ",             /, &
              " * Propagation time:                            ", F10.4,1X,A, /, & 
              " * Propagation time step:                       ", F10.4,1X,A, /, & 
              " * Propagation print step:                      ", F10.4,1X,A, /   )
   
   !*************************************************************
   !         SETUP SECTION 
   !*************************************************************

#if defined(WITH_MPI)
   ! >>> SEND DATA TO SLAVES
   CALL SyncroDataAcrossNodes( )
   CALL SyncroPotentialDataAcrossNodes( )
#endif
   
   ! Initialize random number seed
   CALL SetSeed( RandomNr, - (1 + __MPI_CurrentNrOfProc ) )
   ! Setup the potential
   CALL SetupPotential( )
   ! Setup periodic boundary conditions
   IF ( PotentialIsPeriodic() ) THEN
      CALL PBC_Setup(  GetPotentialUnitCell() )
   END IF
   ! Setup data for propagation and allocate memory
   CALL SetupPropagation( )

#if defined(WITH_MPI)
   ! Define et of Indices of iTraj to run on each slave
   allocate( NodeOfTraj(NrTrajs), NrOfTrajPerNode(0:__MPI_TotalNrOfProcs-1) )  
   DO iTraj = 1, NrTrajs
         NodeOfTraj(iTraj) = MOD(iTraj,__MPI_TotalNrOfProcs)
   END DO
   ! Compute the number of trajectories that are executed by each node
   DO iNode = 0, __MPI_TotalNrOfProcs-1
      NrOfTrajPerNode(iNode) = 0
      DO iTraj = 1, NrTrajs 
         IF ( NodeOfTraj(iTraj) == iNode ) NrOfTrajPerNode(iNode) = NrOfTrajPerNode(iNode) + 1 
      END DO
   END DO
#endif

   ! Print info on trajectory parallelization 
   __MPI_OnlyMasterBEGIN
   WRITE(*,"(A,A,A)")  " Running ", TRIM(NumberToString(NrTrajs)), " trajectories ... "
#if defined(WITH_MPI)
   DO iNode = 0, __MPI_TotalNrOfProcs-1
      WRITE(*,"(A,A,A,A)")     " Trajectories on node # ", TRIM(NumberToString(iNode)),   &
                                 " : ",  TRIM(NumberToString(NrOfTrajPerNode(iNode)))
   END DO
#endif
   __MPI_OnlyMasterEND
   
!    Compare numerical and analytical forces, for debugging purpose
!    CALL CheckForces( )

   !*************************************************************
   !         RUN SECTION 
   !*************************************************************

   DO iTraj = 1, NrTrajs

#if defined(WITH_MPI)
      IF ( NodeOfTraj(iTraj) == __MPI_CurrentNrOfProc ) then
#endif

         ! Store internal state of the random number generator at the beginning of the traj
         OutputSet0 = GetInternalState( RandomNr )
      
         ! Initialize coordinates of the trajectory
         CALL InitializeTrajectory()

         ! Evolve trajectory: equilbration and then microcanonical dynamics
         CALL Thermalization( )
         CALL DynamicsRun( )

         ! Print to standard output info about the trajectory
         WRITE(*, "(2/,A)") "***************************************************"
         WRITE(*, "(A,A)")  "            Trajectory Nr. ", TRIM(NumberToString(iTraj))
         WRITE(*, "(A,A)")  "         propagated by node # ", TRIM(NumberToString(NodeOfTraj(iTraj)))
         WRITE(*, "(A,/)" ) "***************************************************"
         WRITE(*, "(A)") " RNG internal state: " // TRIM(NumberToString(OutputSet0(1))) // " " // &
                                                         TRIM(NumberToString(OutputSet0(2)))
         IF ( NBeads > 1 ) THEN
            WRITE(*,"(/,A)") " Propagating system in the canonical ensamble with PILE ... " 
         ELSE
            WRITE(*,"(/,A)") " Propagating system in the canonical ensamble with symplectic propagator ... " 
         END IF
         WRITE(*,"(A)") " Thermalization completed! "
         WRITE(*,600) OutputSet1(1)*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits),  &
                      OutputSet1(2)*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits),  &
                      OutputSet1(3)*TemperatureConversion(InternalUnits, InputUnits), TemperUnit(InputUnits) 
         WRITE(*,"(/,A)") " Propagating system in the microcanonical ensamble... " 
         WRITE(*,601) OutputSet2(1)*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits),  &
                      OutputSet2(2)*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits),  &
                      OutputSet2(3)*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits),  &
                      OutputSet2(4)*TemperatureConversion(InternalUnits, InputUnits), TemperUnit(InputUnits) 
         WRITE(*,"(A)") " Time propagation completed! "
         WRITE(*,602) OutputSet3(1)*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits),  &
                      OutputSet3(2)*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits),  &
                      OutputSet3(3)*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits),  &
                      OutputSet3(4)*TemperatureConversion(InternalUnits, InputUnits), TemperUnit(InputUnits) 
         WRITE(*,"(A)") " "

#if defined(WITH_MPI)
      END IF
#endif
      
   END DO

   !*************************************************************
   !         AVERAGES AND OUTPUT SECTION 
   !*************************************************************
   !            ### note for parallelization ###
   !     this section should be executed by the master
   !     process only with data collected from all the 
   !     slaves 
   !*************************************************************
   
   ! normalize averages
   CALL DynAveragesOutput( FINALIZE_AVERAGES )
   CALL EquilAveragesOutput( FINALIZE_AVERAGES )
   ! Print averages to outptu files
   CALL DynAveragesOutput( PRINT_AVERAGES_AND_DISPOSE )
   CALL EquilAveragesOutput( PRINT_AVERAGES_AND_DISPOSE )
   
   ! DISPOSE SUBROUTINE AND DATA
   CALL DisposePotential()
   CALL PBC_Dispose()
   CALL DisposeEvolutionData( MolecularDynamics )
   CALL DisposeEvolutionData( InitialConditions )
   CALL DisposeFFT( RingNormalModes )
      
   ! Deallocate MEMORY
   DEALLOCATE( MassVector, X, V, A )
   DEALLOCATE( KinPerCoord, CentroidPos, CentroidVel )
   
   __MPI_OnlyMasterBEGIN
   CALL date_and_time (values=Time2)
   
   WRITE(*,*)
   WRITE(*,"(A,F10.1,A)") " Execution Time : ",TimeDifference( Time2, Time1 )/1000.0, " / s "
   __MPI_OnlyMasterEND

#if defined(WITH_MPI)
   CALL MyMPI_Finalize()
#endif

      600 FORMAT (/, " Average values of the thermalization dynamics "   ,/,  &
                     " * Kinetic Energy "                                ,/,  &
                     "     average                          ",1F10.4,1X,A,/,  &
                     "     standard deviation               ",1F10.4,1X,A,/,  &
                     " * Temperature (time average)         ",1F10.4,1X,A,/) 
                     
      601 FORMAT (/, " Initial condition of the MD trajectory        "   ,/,  &
                     " (for RP, full ring polymer hamiltonian)       "   ,/,  &
                     " * Kinetic Energy                     ",1F10.4,1X,A,/,  &
                     " * Potential Energy                   ",1F10.4,1X,A,/,  &
                     " * Total Energy                       ",1F10.4,1X,A,/,  &
                     " * Istantaneous Temperature           ",1F10.4,1X,A,/) 

      602 FORMAT (/, " Final condition of the MD trajectory          "   ,/,  &
                     " (for RP, full ring polymer hamiltonian)       "   ,/,  &
                     " * Kinetic Energy                     ",1F10.4,1X,A,/,  &
                     " * Potential Energy                   ",1F10.4,1X,A,/,  &
                     " * Total Energy                       ",1F10.4,1X,A,/,  &
                     " * Istantaneous Temperature           ",1F10.4,1X,A,/)                      
  
!============================================================================================
                                  CONTAINS
!============================================================================================


!*******************************************************************************
!                                  Input
!*******************************************************************************
!> Read from input file the necessary variables to run the code.
!*******************************************************************************
   SUBROUTINE Input( )
      IMPLICIT NONE

      ! read input units ( or set them to default value )
      CALL SetFieldFromInput( InputData, "InputLength", InputLength,  UNITS_ANGSTROM )
      CALL SetFieldFromInput( InputData, "InputEnergy", InputEnergy,  UNITS_EV       )
      CALL SetFieldFromInput( InputData, "InputMass",   InputMass,    UNITS_AMU      )
      CALL SetFieldFromInput( InputData, "InputTime",   InputTime,    UNITS_FEMTOS   )
      CALL SetFieldFromInput( InputData, "InputTemp",   InputTemp,    UNITS_KELVIN   )
      CALL SetFieldFromInput( InputData, "InputFreq",   InputFreq,    UNITS_CMMINUS1 )
   
      ! READ THE INFORMATION ABOUT THE SYSTEM
   
      ! Nr of replicas of the system for the ring polymer dynamics
      CALL SetFieldFromInput( InputData, "NBeads", NBeads )
      ! Temperature of the system
      CALL SetFieldFromInput( InputData, "Temperature", Temperature )
      ! Nr of the trajectories of the simulation
      CALL SetFieldFromInput( InputData, "NrTrajs", NrTrajs )
      
      ! READ INFO ABOUT THE OUTPUT
      
      ! Set time between writing the output
      CALL SetFieldFromInput( InputData, "PrintTimeStep", PrintTimeStep )
      
      ! READ THE INFORMATION ABOUT THE EQUILIBRATION
   
      ! Set gamma of the equilibration Langevin dynamics
      CALL SetFieldFromInput( InputData, "EquilRelaxTime", LangevinGamma )
      ! Set equilibration time
      CALL SetFieldFromInput( InputData, "EquilTotalTime", EquilTotalTime )
      ! Set equilibration time step
      CALL SetFieldFromInput( InputData, "EquilTimeStep", EquilTimeStep )
      
      ! READ THE VARIABLES FOR THE TIME PROPAGATION

      ! Timestep of the propagation
      CALL SetFieldFromInput( InputData, "TimeStep", TimeStep )
      ! Set propagation time
      CALL SetFieldFromInput( InputData, "DynamicsTotalTime", DynamicsTotalTime )

      ! READ THE VARIABLES FOR THE ANALYSIS

      ! Kinetic energy distribution
      CALL SetFieldFromInput( InputData, "Out_KinDistrib", Out_KinDistrib, .FALSE. )
      IF ( Out_KinDistrib ) THEN
         CALL SetFieldFromInput( InputData, "Out_KinDistrib_nE", Out_KinDistrib_nE )
         CALL SetFieldFromInput( InputData, "Out_KinDistrib_DE", Out_KinDistrib_DE )
      END IF

   END SUBROUTINE Input

   ! ===================================================================================================

   SUBROUTINE ConvertUnitsAndProcessData( )
      IMPLICIT NONE

      ! Setup conversion factors from input to internal units
      CALL Initialize_UnitConversion( InputUnits, InputLength, InputEnergy, InputMass, 11, InputTime, InputTemp, InputFreq )

      ! Convert input data to internal units
      Temperature = Temperature * TemperatureConversion(InputUnits, InternalUnits)
      PrintTimeStep = PrintTimeStep * TimeConversion(InputUnits, InternalUnits)
      EquilTotalTime = EquilTotalTime * TimeConversion(InputUnits, InternalUnits)
      EquilTimeStep = EquilTimeStep * TimeConversion(InputUnits, InternalUnits)
      TimeStep = TimeStep * TimeConversion(InputUnits, InternalUnits)
      DynamicsTotalTime = DynamicsTotalTime * TimeConversion(InputUnits, InternalUnits)

      ! Friction coefficient of the Langevin dynamics for thermalization
      LangevinGamma = 1. / ( LangevinGamma * TimeConversion(InputUnits, InternalUnits) )

      ! Compute relevant step numbers
      NrSteps = CEILING( DynamicsTotalTime / TimeStep )
      EquilNrSteps = CEILING(  EquilTotalTime / EquilTimeStep )
      ! Set the step interval between each printing step
      PrintStepInterval = NINT( PrintTimeStep / TimeStep )
      EquilPrintStepInterval = NINT( PrintTimeStep / EquilTimeStep )
      
      ! Frequency and force constant of the harmonic force between the beads
      BeadsFrequency = NBeads * Temperature
      BeadsForceConst = ( NBeads * Temperature )**2
           
      ! Variables for the kinetic energy binning
      IF ( Out_KinDistrib ) THEN
         Out_KinDistrib_DE = Out_KinDistrib_DE * EnergyConversion(InputUnits, InternalUnits)
      END IF

   END SUBROUTINE ConvertUnitsAndProcessData

   ! ===================================================================================================

   SUBROUTINE SyncroDataAcrossNodes( )
      IMPLICIT NONE

      ! Data read from input file
      CALL MyMPI_BroadcastToSlaves( NBeads )
      CALL MyMPI_BroadcastToSlaves( Temperature )
      CALL MyMPI_BroadcastToSlaves( NrTrajs )
      CALL MyMPI_BroadcastToSlaves( PrintTimeStep )
      CALL MyMPI_BroadcastToSlaves( EquilTotalTime )
      CALL MyMPI_BroadcastToSlaves( LangevinGamma )
      CALL MyMPI_BroadcastToSlaves( EquilTimeStep )
      CALL MyMPI_BroadcastToSlaves( TimeStep )
      CALL MyMPI_BroadcastToSlaves( DynamicsTotalTime )
      
      ! Data trivially computed from input data 
      CALL MyMPI_BroadcastToSlaves( NrSteps )
      CALL MyMPI_BroadcastToSlaves( EquilNrSteps )
      CALL MyMPI_BroadcastToSlaves( PrintStepInterval )
      CALL MyMPI_BroadcastToSlaves( EquilPrintStepInterval )
      CALL MyMPI_BroadcastToSlaves( BeadsFrequency )
      CALL MyMPI_BroadcastToSlaves( BeadsForceConst )

      ! Variables for the kinetic energy binning
      CALL MyMPI_BroadcastToSlaves( Out_KinDistrib )
      IF ( Out_KinDistrib ) THEN
         CALL MyMPI_BroadcastToSlaves( Out_KinDistrib_nE )
         CALL MyMPI_BroadcastToSlaves( Out_KinDistrib_DE )
      END IF

   END SUBROUTINE SyncroDataAcrossNodes

   ! ===================================================================================================

   SUBROUTINE SetupPropagation( )
      IMPLICIT NONE
      REAL :: WellDepth, EquilDist

      ! Get the dimension of the system to allocate memory
      NAtoms = GetAtomsNumber( )
      NDim = GetDoFNumber( )
      ! Allocate memory to store the mass vector
      ALLOCATE( MassVector(NDim) )
      ! Define masses of the system
      MassVector = GetMasses( )
      
      ! Allocate memory and initialize vectors for trajectory, acceleration and masses
      ALLOCATE( X(NDim*NBeads), V(NDim*NBeads), A(NDim*NBeads) )
      
      ! Set variables for EOM integration in the microcanonical ensamble 
      CALL EvolutionSetup( MolecularDynamics, NDim, MassVector, TimeStep )
      ! Set ring polymer molecular dynamics parameter
      IF ( NBeads > 1 ) CALL SetupRingPolymer( MolecularDynamics, NBeads, BeadsFrequency ) 

      ! Set transform from ring coordinates to normal modes
      CALL SetupFFT( RingNormalModes, NBeads ) 

      ! Set variables for EOM integration of the system in the canonical ensamble with Langevin dynamics
      CALL EvolutionSetup( InitialConditions, NDim, MassVector, EquilTimeStep )
      ! Set ring polymer molecular dynamics parameter
      IF ( NBeads > 1 ) CALL SetupRingPolymer( InitialConditions, NBeads, BeadsFrequency ) 
      ! Set Langevin thermostat for initial equilibration of the system
      CALL SetupThermostat( InitialConditions, LangevinGamma, Temperature )

      ! Allocate memory for istantaneous averages
      ALLOCATE( KinPerCoord(NDim), CentroidPos(NDim), CentroidVel(NDim) )

      ! Setup output of averages over the set of trajectories
      CALL DynAveragesOutput( SETUP_AVERAGES )
      CALL EquilAveragesOutput( SETUP_AVERAGES )

   END SUBROUTINE SetupPropagation
      
   ! ===================================================================================================

   SUBROUTINE InitializeTrajectory( )
      IMPLICIT NONE

      ! *************  Initial conditions of the system *****************

      ! Define appropriate initial conditions depending on the potential
      X(:) = 0.0;  V(:) = 0.0
      CALL GetInitialPositions( X( 1 : NDim ), RandomNr )
      
      ! Copy the coordinates on all the beads of the ring polymer
      IF ( NBeads > 1 ) THEN
         DO iCoord = 2, NBeads
            X( NDim*(iCoord-1)+1 : NDim*iCoord ) = X( 1 : NDim ) 
         END DO
      END IF

      !>>>>>>> A BETTER INITIAL COND FOR RPMD SHOULD BE IMPLEMENTED <<<<<<<
      ! E.g. : ring polymer with gyration radius as in free particle 
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      ! Initialize output for the single trajectory
      CALL SingleTrajectoryOutput( SETUP_OUTPUT )

   END SUBROUTINE InitializeTrajectory

   ! ===================================================================================================

   SUBROUTINE Thermalization( )
      IMPLICIT NONE
      INTEGER :: iStep
      REAL    :: KinTimeAverage, KinTimeStDev


      ! Bring the atomic coordinates to the first unit cell
      CALL PBC_BringToFirstCell( X, NBeads )

      ! Initialize forces
      IF ( NBeads > 1 ) THEN
         CALL EOM_RPMSymplectic( InitialConditions, X, V, A, GetPotential, PotEnergy, RandomNr, 1 )
      ELSE
         PotEnergy = GetPotential( X, A )
         A(:) = A(:) / MassVector(:)
      END IF

      KinTimeAverage = 0.0
      KinTimeStDev   = 0.0

      kStep = 0
      ! Cycle over the equilibration steps
      DO iStep = 0, EquilNrSteps

         IF ( iStep > 0 ) THEN
            ! Propagate for one timestep with Velocity-Verlet
            IF ( NBeads > 1 ) THEN
               CALL EOM_RPMSymplectic( InitialConditions, X, V, A, GetPotential, PotEnergy, RandomNr )
            ELSE
               CALL EOM_LangevinSecondOrder( InitialConditions, X, V, A, GetPotential, PotEnergy, RandomNr )
            END IF

            ! Bring the atomic coordinates to the first unit cell
            CALL PBC_BringToFirstCell( X, NBeads )
         ENDIF

         IF ( MOD(iStep, EquilPrintStepInterval) == 0.0 ) THEN

            ! New print step of the equilibration
            kStep = kStep + 1

            ! Update time value
            Time = REAL(iStep)*EquilTimeStep

            ! Compute average energies
            CALL KineticAverages( X, V, A, MassVector, KinEnergy, KinPerCoord ) 
            TotEnergy = PotEnergy + KinEnergy
            IF ( NBeads > 1 ) THEN
               ! compute mechanical energy of the ring polymer
               RPKinEnergy = EOM_KineticEnergy( InitialConditions, V )
               RPPotEnergy = EOM_InterBeadsPotential( InitialConditions, X ) + NBeads * PotEnergy
               RPTotEnergy = RPKinEnergy + RPPotEnergy
            END IF
            ! Compute coordinate and velocity centroid
            CentroidPos = CentroidCoord( X )
            CentroidVel = CentroidCoord( V )
            ! Compute time average of kinetic energy
            IF ( NBeads > 1 ) THEN
               KinTimeAverage = KinTimeAverage + RPKinEnergy
               KinTimeStDev   = KinTimeStDev   + RPKinEnergy**2
            ELSE
               KinTimeAverage = KinTimeAverage + KinEnergy
               KinTimeStDev   = KinTimeStDev   + KinEnergy**2
            END IF

            ! Print istantaneous values of the trajectory
            CALL SingleTrajectoryOutput( PRINT_OUTPUT )

            ! Update averages over the set of trajectories
            CALL EquilAveragesOutput( UPDATE_AVERAGES )

         END IF
      
      END DO
      
      KinTimeAverage = KinTimeAverage / kStep
      KinTimeStDev = SQRT( KinTimeStDev / kStep - KinTimeAverage**2 )

      OutputSet1(:) = (/ KinTimeAverage, KinTimeStDev, 2.0*KinTimeAverage/NDim/NBeads /)

      ! Write to the output file a separation between equilibration and dynamics
      CALL SingleTrajectoryOutput( DIVIDE_EQ_DYN )


   END SUBROUTINE Thermalization

   ! ===================================================================================================

   SUBROUTINE DynamicsRun( )
      IMPLICIT NONE
      INTEGER :: iStep

      ! Bring the atomic coordinates to the first unit cell
      CALL PBC_BringToFirstCell( X, NBeads )
 
      ! Compute starting potential and forces
      IF ( NBeads > 1 ) THEN
         CALL EOM_RPMSymplectic( MolecularDynamics, X, V, A, GetPotential, PotEnergy, RandomNr, 1 )
      ELSE
         PotEnergy = GetPotential( X, A )
         A(:) = A(:) / MassVector(:)
      END IF

      kStep = 0

      ! cycle over nstep velocity verlet iterations
      DO iStep = 0,NrSteps

         IF ( iStep > 0 ) THEN 
            ! Propagate for one timestep with Velocity-Verlet
            IF ( NBeads > 1 ) THEN
               CALL EOM_RPMSymplectic( MolecularDynamics, X, V, A, GetPotential, PotEnergy, RandomNr )
            ELSE
               CALL EOM_LangevinSecondOrder( MolecularDynamics, X, V, A, GetPotential, PotEnergy, RandomNr )
            END IF

            ! Bring the atomic coordinates to the first unit cell
            CALL PBC_BringToFirstCell( X, NBeads )
         END IF

         ! output to write every nprint steps 
         IF ( mod(iStep,PrintStepInterval) == 0 ) THEN

            ! increment counter for printing steps
            kStep = kStep+1

            ! Update time value
            Time = REAL(iStep)*TimeStep

            ! Compute average energies
            CALL KineticAverages( X, V, A, MassVector, KinEnergy, KinPerCoord ) 
            TotEnergy = PotEnergy + KinEnergy
            IF ( NBeads > 1 ) THEN
               ! compute mechanical energy of the ring polymer
               RPKinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
               RPPotEnergy = EOM_InterBeadsPotential( MolecularDynamics, X ) + NBeads * PotEnergy
               RPTotEnergy = RPKinEnergy + RPPotEnergy
            END IF

            ! Compute coordinate and velocity centroid
            CentroidPos = CentroidCoord( X )
            CentroidVel = CentroidCoord( V )               

            ! For the first step, store initial energy values to print trajectory info to st out
            IF ( iStep == 0 ) THEN
               IF ( NBeads > 1 ) THEN
                  OutputSet2(:) = (/ RPKinEnergy, RPPotEnergy, RPTotEnergy, 2.0*RPKinEnergy/NDim/NBeads /)
               ELSE
                  OutputSet2(:) = (/ KinEnergy, PotEnergy, TotEnergy, 2.0*KinEnergy/NDim /)
               END IF
            END IF

            ! Print istantaneous values of the trajectory
            CALL SingleTrajectoryOutput( PRINT_OUTPUT )

            ! Update averages over the set of trajectories
            CALL DynAveragesOutput( UPDATE_AVERAGES )

         END IF 

      END DO

      ! Store final energy values to print trajectory info to st out
      IF ( NBeads > 1 ) THEN
         OutputSet3(:) = (/ RPKinEnergy, RPPotEnergy, RPTotEnergy, 2.0*RPKinEnergy/NDim/NBeads /)
      ELSE
         OutputSet3(:) = (/ KinEnergy, PotEnergy, TotEnergy, 2.0*KinEnergy/NDim /)
      END IF

      ! Close output file
      CALL SingleTrajectoryOutput( CLOSE_OUTPUT )

   END SUBROUTINE DynamicsRun

   ! ===================================================================================================
   
   SUBROUTINE KineticAverages( X, V, A, Mass, KinTot, KinCoords ) 
      IMPLICIT NONE
      REAL, DIMENSION(NDim*NBeads), INTENT(IN) :: X, V, A
      REAL, DIMENSION(NDim), INTENT(IN)        :: Mass
      REAL, INTENT(OUT)                        :: KinTot
      REAL, DIMENSION(NDim), INTENT(OUT)       :: KinCoords

      REAL, DIMENSION(NDim*NBeads)  :: SngForce
      REAL, DIMENSION(NDim)         :: CentroidV
      INTEGER                       :: iCoord, iBead

      IF ( NBeads > 1 ) THEN

         ! ***************** VIRIAL AVERAGE PER ATOM **************************
         ! centroid of the virial product (virial as average)

         KinCoords(:) = 0.0
         DO iBead = 1, NBeads
            DO iCoord = 1, NDim
               KinCoords(iCoord) = KinCoords(iCoord) - X((iBead-1)*NDim+iCoord) * Mass(iCoord) * A((iBead-1)*NDim+iCoord)
            END DO
         END DO         
         KinCoords(:) = 0.5 * KinCoords(:) / real(NBeads)

         ! ******************* VIRIAL TOTAL ENERGY ****************************
         
         KinTot = 0.0
         DO iCoord = 1, NDim
            KinTot = KinTot + KinCoords(iCoord)
         END DO         

      ELSE IF ( NBeads == 1 ) THEN

         ! ***************** MECHANICAL KIN PER ATOM **************************
         ! regular definition of kinetic energy

         DO iCoord = 1, NDim
            KinCoords(iCoord) = 0.5 * Mass(iCoord) * V(iCoord)**2
         END DO

         ! ******************* MECHANICAL TOTAL KIN ENERGY ********************
         
         KinTot = 0.0
         DO iCoord = 1, NDim
            KinTot = KinTot + KinCoords(iCoord)
         END DO         
         
      END IF

!       ! 2) THERMODYNAMICS TOTAL ENERGY
!       ! sum up energy of the ring oscillators and subtract it to average energy
!       EnergyAverages(2) = 0.0 
!       DO iCoord = 1, NDim
!          IF ( NBeads == 2 ) THEN
!             EnergyAverages(2) = EnergyAverages(2) + Mass(iCoord) * ( X(NDim+iCoord) - X(iCoord) )**2 
!          ELSE IF ( NBeads > 2 ) THEN
!             DO iBead = 1, NBeads-1
!                EnergyAverages(2) = EnergyAverages(2) + Mass(iCoord) * ( X(iBead*NDim+iCoord) - X((iBead-1)*NDim+iCoord) )**2 
!             END DO
!             EnergyAverages(2) = EnergyAverages(2) + Mass(iCoord) * ( X(iCoord) - X((NBeads-1)*NDim+iCoord) )**2 
!          END IF
!       END DO
!       EnergyAverages(2) = 0.5 * ( NBeads * NDim * Temperature - BeadsForceConst * EnergyAverages(2) )
!       ! add potential energy to get full energy
!       EnergyAverages(2) = EnergyAverages(2) + EnergyAverages(4)
! 
!       ! 3) POTENTIAL PLUS KINETIC ENERGY OF THE CENTROID
!       ! Compute centroid of the velocity
!       CentroidV = CentroidCoord( V )
!       EnergyAverages(3) = 0.0
!       DO iCoord = 1, NDim
!          EnergyAverages(3) = EnergyAverages(3) + Mass(iCoord) * CentroidV(iCoord)**2
!       END DO
!       EnergyAverages(3) =  0.5 * EnergyAverages(3)
!       ! add potential energy to get full energy
!       EnergyAverages(3) = EnergyAverages(3) + EnergyAverages(4)

   END SUBROUTINE KineticAverages
      
   ! ===================================================================================================

   FUNCTION CentroidCoord( X )
      IMPLICIT NONE
      REAL, DIMENSION( NDim )   :: CentroidCoord
      REAL, DIMENSION( NDim*NBeads ), INTENT(IN) :: X
      INTEGER :: i, iCoord

      CentroidCoord( 1:NDim ) = 0.0
      DO i = 1, NBeads
         DO iCoord = 1, NDim
            CentroidCoord( iCoord ) = CentroidCoord( iCoord ) + X( (i-1)*NDim + iCoord )
         END DO
      END DO
      CentroidCoord( 1:NDim ) = CentroidCoord( 1:NDim ) / NBeads

   END FUNCTION CentroidCoord

   ! ===================================================================================================

   FUNCTION GyrationRadius( X ) RESULT( Gyration )
      IMPLICIT NONE
      REAL, DIMENSION( NDim*NBeads ), INTENT(IN) :: X
      REAL, DIMENSION( NDim )   :: Centroid, Gyration
      INTEGER :: iBead, iCoord

      Centroid(:) = 0.0
      Gyration(:) = 0.0
      DO iBead = 1, NBeads
         DO iCoord = 1, NDim
            Centroid( iCoord ) = Centroid( iCoord ) + X( (iBead-1)*NDim + iCoord )
            Gyration( iCoord ) = Gyration( iCoord ) + X( (iBead-1)*NDim + iCoord )**2
         END DO
      END DO
      Centroid( 1:NDim ) = Centroid( 1:NDim ) / NBeads
      Gyration( 1:NDim ) = Gyration( 1:NDim ) / NBeads - Centroid( 1:NDim )**2

   END FUNCTION GyrationRadius

   ! ===================================================================================================
   
   REAL FUNCTION CorrelationFunction( X, X0 )
      IMPLICIT NONE
      REAL, DIMENSION( NDim*NBeads ), INTENT(IN) :: X, X0
      REAL, DIMENSION( NDim )                    :: XCentroid, X0Centroid

      XCentroid = CentroidCoord( X )
      X0Centroid = CentroidCoord( X0 )
      CorrelationFunction =  TheOneWithVectorDotVector( XCentroid(:), X0Centroid(:) )

   END FUNCTION CorrelationFunction
   
   ! ===================================================================================================

   FUNCTION TimeDifference( Time1, Time2 )
      IMPLICIT NONE
      INTEGER, DIMENSION(8), INTENT(IN)    :: Time1, Time2         
      REAL :: TimeDifference
   
      TimeDifference =  (Time1(3)-Time2(3))*24.0*60.0*60.0*1000. + &
                        (Time1(5)-Time2(5))*60.0*60.0*1000. + &
                        (Time1(6)-Time2(6))*60.0*1000. + &
                        (Time1(7)-Time2(7))*1000. + &
                        (Time1(8)-Time2(8))
         
   END FUNCTION TimeDifference
      
   ! ===================================================================================================

   SUBROUTINE CheckForces( )
      IMPLICIT NONE

      REAL, DIMENSION(:), ALLOCATABLE :: ForceAnalytic, ForceNumeric, XpD, XmD
      REAL :: VpD, VmD
      INTEGER :: i

      ALLOCATE( ForceAnalytic(NDim), ForceNumeric(NDim), XpD(NDim), XmD(NDim) )
      CALL GetInitialPositions( X,  RandomNr ) 
      VpD =  GetPotential( X, ForceAnalytic )

      PRINT*, " ANALYTIC FORCE "
      PRINT*, ForceAnalytic
      PRINT*, " "
      DO i = 1, size(X)
         XpD = X
         XpD(i) = X(i) + 0.00001
         XmD = X
         XmD(i) = X(i) - 0.00001
         VpD = GetPotential( XpD, ForceAnalytic )
         VmD = GetPotential( XmD, ForceAnalytic )
         ForceNumeric(i) = - ( VpD - VmD ) / 2.0 / 0.00001
      END DO
      
      PRINT*, " NUMERIC FORCE "
      PRINT*, ForceNumeric
      PRINT*, " "

      PRINT*, " % DIFFERENCE "
      DO i = 1, size(X)
         WRITE(*,"(1F20.6)") ABS(ForceNumeric(i) - ForceAnalytic(i))/ABS(ForceNumeric(i))*100.
      END DO

      STOP

   END SUBROUTINE CheckForces

   END PROGRAM PolyOnTheFly
      
      