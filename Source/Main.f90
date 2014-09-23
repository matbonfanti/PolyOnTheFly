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
!>  \arg 
!
!>  \todo          ____________________________
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

   ! FEW VARIABLES HERE
   ! DEFINE MOST VARIABLES IN SharedData MODULE

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
   INTEGER :: CurrentMPITask = 0
   
   INTEGER :: iCoord
   
   ! 1) INPUT SECTION
   !    read info about the molecular system
   !    define parameters of the method
   
   ! 2) SETUP SECTION
   !    memory allocation
   !    setup potential
   !    setup propagation
   !    setup random number generator ( parallel )
   !    (other)
  
   ! CYCLE OVER THE TRAJECTORIES
   
      ! 3-i) EQUILIBRATION SECTION
   
      ! 4-i) RUN SECTION

   ! 5) CONCLUSION AND AVERAGES
   
   PRINT "(/,     '                    ==============================')"
   PRINT "(       '                              PolyOnTheFly        ')"
   PRINT "(       '                    ==============================',2/)"
   PRINT "(       '                  Author: Matteo Bonfanti, August 2014  ')"
   PRINT "(       '         On-The-Fly Ring Polymer Molecular Dynamics code in FORTRAN 90/95   ',2/)"

   
   !*************************************************************
   !         COMMAND LINE ARGUMENT
   !*************************************************************

   ! Check and read from command line the input file name
   NArgs = COMMAND_ARGUMENT_COUNT()
   IF (NArgs<1) THEN
      Help = .TRUE.
   ELSE
      CALL GET_COMMAND_ARGUMENT( 1, InputFileName )
      IF ( trim(InputFileName) == "help" ) Help = .TRUE.
   ENDIF
   IF (Help) THEN ! Call help
      PRINT*, ' Launch this program as:'
      PRINT*, ' % PolyOnTheFly "InputFileName" '
      STOP
   ENDIF

   CALL date_and_time (values=Time1)

   !*************************************************************
   !         INPUT SECTION
   !*************************************************************
   !            ### note for parallelization ###
   !     this section should be executed by the master
   !     process only when data is read and processed,
   !     it is distributed to all the slaves
   !*************************************************************

   ! Open and read from input file the input parameters of the calculation
   CALL OpenFile( InputData, InputFileName )

   ! Read the Input data
   CALL Input( )
   ! Convert to internal units
   CALL ConvertUnitsAndProcessData( )
   
   ! close input file
   CALL CloseFile( InputData )

   !       WRITE(*, 902) NBeads, BeadsFrequency/MyConsts_cmmin1toAU
! 
!       IF ( .NOT. MorsePotential ) THEN
!            WRITE(*, 903) InitEnergy*MyConsts_Hartree2eV, (InitEnergy-MinimumEnergy)*MyConsts_Hartree2eV
!       ELSE 
!            WRITE(*, 903) InitEnergy*MyConsts_Hartree2eV, (InitEnergy+MorseDe)*MyConsts_Hartree2eV
!       END IF
! 
!       WRITE(*, 904) NrTrajs, TimeStep/MyConsts_fs2AU, NrOfSteps, NrOfPrintSteps
! 
!       WRITE(*, 905) Temperature/MyConsts_K2AU
! 
! 
!    902 FORMAT(" * Nr of replicas in the ring polymer dynamics: ", I10,/, &
!               " * Frequency of the force between the beads     ", F10.4,/ )
! 
!    903 FORMAT(" * Initial conditions of the atom-surface system ", /,&
!               " * Absolute initial energy (eV):                ",F10.4,/,& 
!               " *  - w.r.t. the bottom of the well (eV):       ",F10.4,/ )
! 
!    904 FORMAT(" * Dynamical simulation variables               ", /,&
!               " * Nr of trajectories:                          ",I10,  /,& 
!               " * Propagation time step (fs):                  ",F10.4,/,& 
!               " * Nr of time steps of each trajectory:         ",I10,  /,& 
!               " * Nr of print steps of each trajectory:        ",I10,  / )
! 
!    905 FORMAT(" * Bath equilibration variables                 ", /,&
!               " * Temperature of the surface:                  ",F10.4 )


   
   !*************************************************************
   !         SETUP SECTION 
   !*************************************************************

   ! >>>>>>>>>>>>>>>>>>>> HERE THE MPI FORK >>>>>>>>>>>>>>>>>>>>>>>>

   CurrentMPITask = 1      ! SERIAL EXECUTION

   ! >>> SEND DATA TO SLAVES
   ! >>> DEFINE SET OF INDICES of iTraj TO RUN ON EACH SLAVE
   
   CALL Setup( )

!    CALL CheckForces( )

   !*************************************************************
   !         RUN SECTION 
   !*************************************************************

   PRINT "(A,I5,A)"," Running ", NrTrajs, " trajectories ... "

   DO iTraj = 1, NrTrajs

      PRINT "(2/,A)", "***************************************************"
      PRINT "(A,I4)", "            Trajectory Nr. ", iTraj
      PRINT "(A,/)" , "***************************************************"

      CALL InitializeTrajectory()

      CALL Thermalization( )
      CALL DynamicsRun( )

   END DO
   
   !*************************************************************
   !         AVERAGES AND OUTPUT SECTION 
   !*************************************************************
   !            ### note for parallelization ###
   !     this section should be executed by the master
   !     process only with data collected from all the 
   !     slaves 
   !*************************************************************
   
   CALL Averages( )
   
   ! DISPOSE SUBROUTINE
  
   CALL date_and_time (values=Time2)
   
   WRITE(*,*)
   WRITE(*,"(A,F10.1,A)") " Execution Time : ",TimeDifference( Time2, Time1 )/1000.0, " / s "
   
  
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
      ! ********* DEFAULT VALUES *********** !
      !      distance    - Angstrom          !
      !      energy      - electronVolt      !
      !      mass        - AMU               !
      !      time        - femtosecond       !
      !      temperature - Kelvin            !
      !      frequency   - cm-1              !
      ! *************************************!
      CALL SetFieldFromInput( InputData, "InputLength", InputLength,  1 )
      CALL SetFieldFromInput( InputData, "InputEnergy", InputEnergy,  3 )
      CALL SetFieldFromInput( InputData, "InputMass",   InputMass,    8 )
      CALL SetFieldFromInput( InputData, "InputTime",   InputTime,   13 )
      CALL SetFieldFromInput( InputData, "InputTemp",   InputTemp,   16 )
      CALL SetFieldFromInput( InputData, "InputFreq",   InputFreq,   18 )
   
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
           
   END SUBROUTINE ConvertUnitsAndProcessData

   ! ===================================================================================================

   SUBROUTINE Setup( )
      IMPLICIT NONE
      REAL :: WellDepth, EquilDist

      ! Initialize random number seed
      CALL SetSeed( RandomNr, -CurrentMPITask )

      ! Initialize system parameters and potential related subroutines
      CALL SetupPotential( "potential.dat" )

      ! Get the dimension of the system to allocate memory
      NAtoms = GetAtomsNumber( )
      NDim = 3*Natoms
      
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

      ! Set variables for EOM integration of the system only in the microcanonical ensamble 
      ! this will be done in a serial way, so no replicated data
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

   END SUBROUTINE Setup
      
   ! ===================================================================================================

   SUBROUTINE InitializeTrajectory( )
      IMPLICIT NONE

      ! *************  Initial conditions of the system *****************
      ! DEFINE APPROPRIATE INITIAL CONDITIONS 
      X(:) = 0.0;  V(:) = 0.0
      CALL GetInitialPositions( X( 1 : NDim ), RandomNr )
      DO iCoord = 2, NBeads
         X( NDim*(iCoord-1)+1 : NDim*iCoord ) = X( 1 : NDim ) 
      END DO

!       DO iCoord = 1, NDim*NBeads
!          X(iCoord) = X(iCoord) + UniformRandomNr( RandomNr ) * 0.5
!       END DO

      ! >>> ONLY MAIN SHOULD EXECUTE THE FOLLOWING CALL
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
         PRINT "(/,A)", " Propagating system in the canonical ensamble with PILE ... " 
         CALL EOM_RPMSymplectic( InitialConditions, X, V, A, GetPotential, PotEnergy, RandomNr, 1 )
      ELSE
         PRINT "(/,A)", " Propagating system in the canonical ensamble with symplectic propagator ... " 
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
      
      PRINT "(A)", " Thermalization completed! "

      KinTimeAverage = KinTimeAverage / kStep
      KinTimeStDev = SQRT( KinTimeStDev / kStep - KinTimeAverage**2 )

      PRINT 600, KinTimeAverage*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits),  &
                 KinTimeStDev*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits),  &
                 2.0*KinTimeAverage/NDim/NBeads*TemperatureConversion(InternalUnits, InputUnits), TemperUnit(InputUnits) 

      ! Write to the output file a separation between equilibration and dynamics
      CALL SingleTrajectoryOutput( DIVIDE_EQ_DYN )

      600 FORMAT (/, " Average values of the thermalization dynamics "   ,/,  &
                     " * Kinetic Energy "                                ,/,  &
                     "     average                          ",1F10.4,1X,A,/,  &
                     "     standard deviation               ",1F10.4,1X,A,/,  &
                     " * Temperature (time average)         ",1F10.4,1X,A,/) 

   END SUBROUTINE Thermalization

   ! ===================================================================================================

   SUBROUTINE DynamicsRun( )
      IMPLICIT NONE
      INTEGER :: iStep

      PRINT "(/,A)", " Propagating system in the microcanonical ensamble... " 

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

            IF ( iStep == 0 ) THEN
               IF ( NBeads > 1 ) THEN
                  PRINT 601, RPKinEnergy*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits),  &
                           RPPotEnergy*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits),  &
                           RPTotEnergy*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits),  &
                           2.0*RPKinEnergy/NDim/NBeads*TemperatureConversion(InternalUnits, InputUnits), TemperUnit(InputUnits) 
               ELSE
                  PRINT 601, KinEnergy*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits), &
                           PotEnergy*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits), &
                           TotEnergy*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits), &
                           2.0*KinEnergy/NDim*TemperatureConversion(InternalUnits, InputUnits), TemperUnit(InputUnits) 
               END IF
            END IF

            ! Print istantaneous values of the trajectory
            CALL SingleTrajectoryOutput( PRINT_OUTPUT )

            ! Update averages over the set of trajectories
            CALL DynAveragesOutput( UPDATE_AVERAGES )

         END IF 

      END DO

      PRINT "(A)", " Time propagation completed! "

      IF ( NBeads > 1 ) THEN
         PRINT 602, RPKinEnergy*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits),  &
                    RPPotEnergy*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits),  &
                    RPTotEnergy*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits),  &
                    2.0*RPKinEnergy/NDim/NBeads*TemperatureConversion(InternalUnits, InputUnits), TemperUnit(InputUnits) 
      ELSE
         PRINT 602, KinEnergy*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits), &
                    PotEnergy*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits), &
                    TotEnergy*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits), &
                    2.0*KinEnergy/NDim*TemperatureConversion(InternalUnits, InputUnits), TemperUnit(InputUnits) 
      END IF

      ! Close output file
      CALL SingleTrajectoryOutput( CLOSE_OUTPUT )

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

   END SUBROUTINE DynamicsRun

   ! ===================================================================================================

   SUBROUTINE Averages( )
      IMPLICIT NONE

      ! normalize averages
      CALL DynAveragesOutput( FINALIZE_AVERAGES )
      CALL EquilAveragesOutput( FINALIZE_AVERAGES )
      ! Print averages to outptu files
      CALL DynAveragesOutput( PRINT_AVERAGES_AND_DISPOSE )
      CALL EquilAveragesOutput( PRINT_AVERAGES_AND_DISPOSE )


   END SUBROUTINE Averages

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
      
      
      
      
      
      
      
      
      
      
      


   




!       ! ALLOCATE AND INITIALIZE DATA WITH THE AVERAGES OVER THE SET OF TRAJECTORIES
! 
!       ! PRINTTYPE
!       ! MINIMAL - print the average energies over time only
!       ! DEBUG - print also the average coordinates
!       ! FULL - print detailed information about the initial conditions and about the trajectories
! 
!       ! Allocate and initialize the variable for the trajectory average of the position correlation
!       ALLOCATE( PositionCorrelation(0:NrOfPrintSteps) )
!       PositionCorrelation(:)  = 0.0
!       
!       ! Allocate and initialize the variable for the trajectory averages of the energy
!       ALLOCATE( AverageE(NrOfEnergyAverages,0:NrOfPrintSteps) ) 
!       AverageE(:,:) = 0.0
!       
!       ! Average coordinates over time 
!       IF ( PrintType >= FULL ) THEN
!          ALLOCATE( AverageCoord(NDim,0:NrOfPrintSteps) )
!          AverageCoord(1:NDim,0:NrOfPrintSteps) = 0.0
!       END IF

! 
!    SUBROUTINE PolymerVibrationalRelax_Run()
!       IMPLICIT NONE
!       INTEGER :: CurrentThread, NrOfThreads
!       !> Output units: minimal and standard output
!       INTEGER  ::  AvEnergyOutputUnit, AvCoordOutputUnit, AvBathCoordUnit, TLogUnit, PosCorrelationUnit
!       !> Output units: debug output
!       INTEGER  ::  DebugUnitEn, DebugUnitCoord, DebugUnitVel
!       !> Initial energies of the bath
!       REAL     ::  InitKinAverage, InitKinVariance, InitPotAverage, InitPotVariance
!       !> Energy averages
!       REAL     ::  TotEnergy, PotEnergy, KinEnergy
!       !> integer Counters
!       INTEGER  ::  iTraj, iBead, iStep, kStep, iCoord
!       !> Filename for output files
!       CHARACTER(100) :: OutFileName
!       !> Centroid of positions and velocities
!       REAL, DIMENSION(NDim*NBeads) :: X0
!       REAL, DIMENSION(NSystem*NBeads) :: InitX, InitV
!       REAL, DIMENSION(NrOfInitSnapshots,NSystem*2) :: SystemInitConditions
! 
!       REAL, DIMENSION(NBath,NBeads) :: InitQBath, InitVBath
! 
!       REAL, DIMENSION(NDim,0:NrOfPrintSteps) :: GyrationAverage
! 
!       PosCorrelationUnit = LookForFreeUnit()
!       OPEN( FILE="PositionCorrelation.dat", UNIT=PosCorrelationUnit )
!       WRITE(PosCorrelationUnit, "(3A,I6,A,/)") "# <X(0)X(t)> vs time (",TimeUnit(InputUnits)," | au) - ", NrTrajs, " trajectories "
!                                                                                   

! 
! 
! 
!          DO iBead = 1, NBeads
!             X((iBead-1)*NDim+1:(iBead-1)*NDim+NSystem) = InitX((iBead-1)*NSystem+1: iBead*NSystem)
!             V((iBead-1)*NDim+1:(iBead-1)*NDim+NSystem) = InitV((iBead-1)*NSystem+1: iBead*NSystem)
!          END DO
! ! 
! !          DO iBead = 1, NBeads
! !             X((iBead-1)*NDim+1:(iBead-1)*NDim+NSystem) = SystemInitConditions( kStep, 1:NSystem )
! !             V((iBead-1)*NDim+1:(iBead-1)*NDim+NSystem) = SystemInitConditions( kStep, NSystem+1:NSystem*2 )
! !          END DO
! 
!          ! *************  Initial conditions of the bath *****************
! 
!          ! Use harmonic oscillator sampling for the normal modes of the ring polymer bath
!          CALL BathOfRingsThermalConditions( Bath, NBeads, BeadsFrequency, Temperature*NBeads, InitQBath, InitVBath, & 
!                      PotEnergy, KinEnergy, RandomNr, RingNormalModes(CurrentThread) )
! 
!          DO iBead = 1, NBeads
!             X( (iBead-1)*NDim+NSystem+1 : iBead*NDim ) = InitQBath(:,iBead)
!             V( (iBead-1)*NDim+NSystem+1 : iBead*NDim ) = InitVBath(:,iBead)
!          END DO
! 
!          !*************************************************************
!          ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
!          !*************************************************************
!          
!          ! Increment averages of the initial conditions of the bath
!          InitKinAverage  = InitKinAverage + KinEnergy/real(NBath*NBeads) 
!          InitKinVariance = InitKinVariance + (KinEnergy/real(NBath*NBeads))**2
!          InitPotAverage  = InitPotAverage + PotEnergy/real(NBath*NBeads)
!          InitPotVariance = InitPotVariance + (PotEnergy/real(NBath*NBeads))**2
!          
!          ! PRINT INFORMATION ON INITIAL CONDITIONS of THE BATH and THE SYSTEM
! !         __OMP_OnlyMasterBEGIN
! !          WRITE(*,600)  (KSys+VSys)*MyConsts_Hartree2eV, KinEnergy*MyConsts_Hartree2eV, 2.0*KinEnergy/MyConsts_K2AU, &
! !                                                         PotEnergy*MyConsts_Hartree2eV, 2.0*PotEnergy/MyConsts_K2AU
! !         __OMP_OnlyMasterEND
! 
!          ! Store initial coordinates
!          X0 = X
!          ! Compute position correlation function
!          PositionCorrelation(0) = PositionCorrelation(0) + CorrelationFunction( X, X0 )
!          
!          ! Various expectation values of the energy
!          AverageE(:,0) = AverageE(:,0) + EnergyAverages( X, V, MassVector ) 
! 
!          ! Average values of the centroid coordinates
!          IF ( PrintType >= FULL ) AverageCoord(1:NDim,0)    = AverageCoord(1:NDim,0)    + CentroidCoord( X )
!          IF ( PrintType >= FULL ) GyrationAverage(1:NDim,0) = GyrationAverage(1:NDim,0) + GyrationRadius( X )
! 
!         ! Open unit for massive output, with detailed info on trajectories
!          IF ( PrintType == DEBUG ) THEN
!             WRITE(OutFileName,"(A,I4.4,A)") "Traj_",iTraj,"_Energy.dat"
!             DebugUnitEn = LookForFreeUnit()
!             OPEN( Unit=DebugUnitEn, File=OutFileName )
! 
!             WRITE(OutFileName,"(A,I4.4,A)") "Traj_",iTraj,"_Coord.dat"
!             DebugUnitCoord = LookForFreeUnit()
!             OPEN( Unit=DebugUnitCoord, File=OutFileName )
! 
!             WRITE(OutFileName,"(A,I4.4,A)") "Traj_",iTraj,"_Vel.dat"
!             DebugUnitVel = LookForFreeUnit()
!             OPEN( Unit=DebugUnitVel, File=OutFileName )
! 
!             ! Write initial values
!             WRITE( DebugUnitEn, "(/,A)" ) "# TRAJECTORY ENERGY: time / fs | Ekin, Epot, Etot / Eh "
!             WRITE(DebugUnitEn,800) 0.0,  KinEnergy, PotEnergy, KinEnergy+PotEnergy
! 
!             WRITE( DebugUnitCoord, "(/,A)" ) "# TRAJECTORY COORD: time / fs | X(1) X(2) ... X(N) / bohr "
!             WRITE(DebugUnitCoord,800) 0.0, X(:)
! 
!             WRITE( DebugUnitVel, "(/,A)" ) "# TRAJECTORY VELOCITIES: time / fs | X(1) X(2) ... X(N) / au "
!             WRITE(DebugUnitVel,800) 0.0, V(:)
! 
!           END IF
! 
!          !*************************************************************
!          !         TIME EVOLUTION OF THE TRAJECTORY
!          !*************************************************************
!  
!          ! initialize counter for printing steps
!          kStep = 0
! 
!          __OMP_OnlyMasterBEGIN 
!          PRINT "(/,A)", " Propagating system and bath in the microcanonical ensamble... " __OMP_OnlyMasterEND
!          
!          ! Compute starting potential and forces
!          CALL EOM_RPMSymplectic( MolecularDynamics(CurrentThread), X, V, A,  GetPotential, PotEnergy, RandomNr, 1 )
! 
!          ! cycle over nstep velocity verlet iterations
!          DO iStep = 1,NrOfSteps
! 
!             ! Propagate for one timestep
!             CALL EOM_RPMSymplectic( MolecularDynamics(CurrentThread), X, V, A, GetPotential, PotEnergy, RandomNr )
! 
!             ! output to write every nprint steps 
!             IF ( mod(iStep,PrintStepInterval) == 0 ) THEN
! 
!                ! increment counter for printing steps
!                kStep = kStep+1
!                IF ( kStep > NrOfPrintSteps ) CYCLE 
! 
!                ! Various expectation values of the energy
!                AverageE(:,kStep) = AverageE(:,kStep) + EnergyAverages( X, V, MassVector ) 
! 
!                ! Compute position correlation function
!                PositionCorrelation(kStep) = PositionCorrelation(kStep) + CorrelationFunction( X, X0 )
!          
!                ! Average values of the centroid coordinates
!                IF ( PrintType >= FULL )    AverageCoord(1:NDim,kStep) = AverageCoord(1:NDim,kStep) + CentroidCoord( X )
!                IF ( PrintType >= FULL ) GyrationAverage(1:NDim,kStep) = GyrationAverage(1:NDim,kStep) + GyrationRadius( X )
! 
!                ! If massive level of output, print traj information to std out
!                IF ( PrintType == DEBUG ) THEN
!                   WRITE(DebugUnitEn,800) TimeStep*real(iStep)/MyConsts_fs2AU, KinEnergy, PotEnergy, KinEnergy+PotEnergy
!                   WRITE(DebugUnitCoord,800) TimeStep*real(iStep)/MyConsts_fs2AU, X(:)
!                   WRITE(DebugUnitVel,800) TimeStep*real(iStep)/MyConsts_fs2AU, V(:)
!                END IF
! 
!             END IF 
! 
!          END DO
! 
!          __OMP_OnlyMasterBEGIN PRINT "(A)", " Time propagation completed! " __OMP_OnlyMasterEND
! 
!          IF ( PrintType == DEBUG ) THEN
!             CLOSE( Unit=DebugUnitEn )
!             CLOSE( Unit=DebugUnitCoord )
!             CLOSE( Unit=DebugUnitVel )
!          END IF
! 
! !          __OMP_OnlyMasterBEGIN
! !          ! Print log information about the final condition of the trajectory
! !          WRITE(*,601)  (KSys+VSys)*MyConsts_Hartree2eV, KBath*MyConsts_Hartree2eV, IstTemperature
! !          __OMP_OnlyMasterEND
! 
!       END DO
!       !$OMP END DO
!       !$OMP END PARALLEL
!       
!       PRINT "(A)"," Done! "
! 
!       InitKinAverage = InitKinAverage / real(NrTrajs) 
!       InitKinVariance = SQRT( InitKinVariance/real(NrTrajs) - InitKinAverage**2 )
!       InitPotAverage = InitPotAverage / real(NrTrajs) 
!       InitPotVariance = SQRT( InitPotVariance/real(NrTrajs) - InitPotAverage**2 ) 
!       
!       TLogUnit = LookForFreeUnit()
!       OPEN( UNIT = TLogUnit, FILE = "FinalAverTemperature.log" )
!       WRITE(TLogUnit,500)  NrTrajs, 2.*InitKinAverage*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits),  &
!                                     2.*InitKinVariance*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits), &
!                                     2.*InitPotAverage*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits),  &
!                                     2.*InitPotVariance*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits)
!       CLOSE( UNIT = TLogUnit )
! 
!       !*************************************************************
!       !         OUTPUT OF THE RELEVANT AVERAGES 
!       !*************************************************************
! 
!       ! Normalize averages 
!       PositionCorrelation(:)  =  PositionCorrelation(:)  / real(NrTrajs)
!       AverageE(:,:)           = AverageE(:,:)            / real(NrTrajs)
!       IF ( PrintType >= FULL )   AverageCoord(1:NDim,:) = AverageCoord(1:NDim,:) / real(NrTrajs) 
!       IF ( PrintType >= FULL )   GyrationAverage(1:NDim,:) = GyrationAverage(1:NDim,:) / real(NrTrajs) 
! 
!       ! PRINT average energy of the system, of the coupling, of the bath
!       DO iStep = 0, NrOfPrintSteps
!          WRITE(AvEnergyOutputUnit,"(F14.8,100F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, &
!                                                      AverageE(:,iStep)*MyConsts_Hartree2eV
!       END DO
! 
!       DO iStep = 0, NrOfPrintSteps
!          WRITE(PosCorrelationUnit,"(F14.8,100F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, &
!                                                      PositionCorrelation(iStep)*MyConsts_Hartree2eV
!       END DO
! 
!       IF ( PrintType >= FULL ) THEN
! 
!          ! PRINT average coordinates
!          DO iStep = 0, NrOfPrintSteps
!             WRITE(AvCoordOutputUnit,"(F14.8,4F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, &
!                                        AverageCoord(1:NSystem,iStep)*MyConsts_Bohr2Ang 
!          END DO
! 
!          ! PRINT gyration radius
!          DO iStep = 0, NrOfPrintSteps
!             WRITE(987,"(F14.8,10000F14.5)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, &
!                                        GyrationAverage(1:NSystem,iStep)*(MyConsts_Bohr2Ang**2) 
!          END DO
! 
!          ! PRINT average bath coordinates
!          DO iCoord = 1, NBath
!             WRITE(AvBathCoordUnit,"(/,A,I5,/)") "#  Bath Coord # ", iCoord
!             DO iStep = 0, NrOfPrintSteps
!                WRITE(AvBathCoordUnit,"(F20.8,3F20.8)") real(iCoord)+AverageCoord(NSystem+iCoord,iStep)*10, &
!                         TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU
!             END DO
!          END DO
! 
!       END IF
! 
!       ! Close output files
!       CLOSE(PosCorrelationUnit)
!       CLOSE( AvEnergyOutputUnit )
!       IF ( PrintType >= FULL ) THEN
!          CLOSE( AvCoordOutputUnit )
!          CLOSE( AvBathCoordUnit )
!       END IF
! 
!       800 FORMAT()
!       600 FORMAT (/, " Initial condition of the MD trajectory ",/   &
!                      " * Energy of the system (eV)          ",1F10.4,/    &
!                      " * Kinetic Energy of the bath (eV)    ",1F10.4,/    &
!                      " * Bath translational temperature (K) ",1F10.4,/    &
!                      " * Potential Energy of the bath (eV)  ",1F10.4,/    &
!                      " * Bath vibrational temperature (K)   ",1F10.4,/ ) 
!       601 FORMAT (/, " Final condition of the MD trajectory ",/   &
!                      " * Energy of the system (eV)         ",1F10.4,/    &
!                      " * Kinetic Energy of the bath (eV)   ",1F10.4,/    &
!                      " * Istantaneous temperature (K)      ",1F10.4,/ ) 
!       500 FORMAT (/, " Nr of trajectories :                               ",1I10,       2/, &
!                      " Average Translational Temperature :                ",1F15.6,1X,A, /, &
!                      " Standard Deviation of Translational Temperature :  ",1F15.6,1X,A,2/, &
!                      " Average Vibrational Temperature :                  ",1F15.6,1X,A, /, &
!                      " Standard Deviation of Vibrational Temperature :    ",1F15.6,1X,A, / ) 
! 
!    END SUBROUTINE PolymerVibrationalRelax_Run
! 
! !*************************************************************************************************
! 
! !*******************************************************************************
! !> Free the memory which has been used for the vibrational relaxation simulation.
! !>
! !*******************************************************************************
!    SUBROUTINE PolymerVibrationalRelax_Dispose()
!       IMPLICIT NONE
!       INTEGER :: i 
! 
!       ! Deallocate memory 
!       DEALLOCATE( X, V, A, MassVector )
!       DEALLOCATE( AverageE )
!       IF ( PrintType >= FULL ) DEALLOCATE( AverageCoord )
! 
!       ! Unset propagators 
!       DO i = 1, size(MolecularDynamics)
!          CALL DisposeEvolutionData( MolecularDynamics(i) )
!          CALL DisposeFFT( RingNormalModes(i) )
!          CALL DisposeEvolutionData( InitialConditions(i) )
!       END DO
! 
!       DEALLOCATE( MolecularDynamics, RingNormalModes, InitialConditions )
! 
!    END SUBROUTINE PolymerVibrationalRelax_Dispose


