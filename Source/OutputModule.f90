!***************************************************************************************
!*                           MODULE OutputModule
!***************************************************************************************
!
!>  \brief     Writing output files
!>  \details   This module controls the level of output desired and handles the 
!>             files in which the results are written.
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
!>  \arg 31 March 2015 :  implemented for each trajectory velocity distribution analysis
!>                        with a velocity binning
!
!>  \todo          ____________________________
!>                 
!***************************************************************************************
MODULE OutputModule
#include "preprocessoptions.cpp"
   USE SharedData
   USE PotentialModule
   USE VTFFileModule
   USE UnitConversion

   ! LIST OF OUTPUT FILES THAT ARE WRITTEN BY THE FOLLOWING SUBROUTINES:

   ! A) SINGLE TRAJECTORY FILES ( 1 file per each trajectory )
   !    1) VTF trajectory file (trajectory snapshots in VTF format, using VTFFileModule)
   !    2) Total energy file (istantaneous value of kin, pot, total energy and temperature of the trajectory)
   !    3) Velocity distribution of each degree of freedom and of the degrees alltogether

   ! B) AVERAGE VALUE FILES ( 1 file per each simulation )
   !    1) 

   PRIVATE

   PUBLIC :: SingleTrajectoryOutput, DynAveragesOutput, EquilAveragesOutput

   !> \name SINGLE TRAJECTORY OUTPUT ACTIONS
   !> Integers number identifying the kind of action to be performed by
   !> the SingleTrajectoryOutput subroutine
   !> @{
   INTEGER, PARAMETER, PUBLIC  ::  SETUP_OUTPUT  = 1
   INTEGER, PARAMETER, PUBLIC  ::  PRINT_OUTPUT  = 2
   INTEGER, PARAMETER, PUBLIC  ::  CLOSE_OUTPUT  = 3
   INTEGER, PARAMETER, PUBLIC  ::  DIVIDE_EQ_DYN = 4
   !> @}

   !> \name AVERAGES OUTPUT ACTIONS
   !> Integers number identifying the kind of action to be performed by
   !> the EquilAveragesOutput and DynAveragesOutput subroutine
   !> @{
   INTEGER, PARAMETER, PUBLIC  ::  SETUP_AVERAGES             = 11
   INTEGER, PARAMETER, PUBLIC  ::  UPDATE_AVERAGES            = 12
   INTEGER, PARAMETER, PUBLIC  ::  FINALIZE_AVERAGES          = 13
   INTEGER, PARAMETER, PUBLIC  ::  PRINT_AVERAGES_AND_DISPOSE = 14
   !> @}

   !> \name STATUS VARIABLE OF THE AVERAGES OUTPUT
   !> Integers identifying the status of the averages output
   !> @{
   INTEGER, PARAMETER  ::  IS_NOT_SET_UP = 0
   INTEGER, PARAMETER  ::  IS_SET_UP     = 1
   INTEGER, PARAMETER  ::  HAS_DATA      = 2
   INTEGER, PARAMETER  ::  IS_FINALIZED  = 3
   !> @}

   !> Setup variable for the output of the current single trajectory values
   LOGICAL, SAVE :: WritingCurrentTrajectory = .FALSE.

   !> Setup variable for the output of equilibration averages
   INTEGER, SAVE :: EquilibrationAveragesStatus = 0

   !> Setup variable for the output of dynamics averages
   INTEGER, SAVE :: DynamicsAveragesStatus = 0

   !> Trajectory output section: equilibration or dynamics
   INTEGER, SAVE :: nEquilDyn = 0

   !> Nr of steps counter
   INTEGER, DIMENSION(2), SAVE :: NrOfSteps

   ! OUTPUT UNITS for single trajectory output
   INTEGER :: TrajTotEnergyUnit
   INTEGER :: TrajRingPolymerEnergyUnit
   INTEGER :: TrajCentroidXUnit
   INTEGER :: TrajCentroidVUnit
   INTEGER :: TrajVelDistrUnit

   ! OUTPUT UNITS for equilibration averages output
   INTEGER :: EquilTotEnergyUnit
   INTEGER :: EquilRingPolymerEnergyUnit
   INTEGER :: EquilCentroidXUnit
   INTEGER :: EquilCentroidVUnit

   ! OUTPUT UNITS for equilibration averages output
   INTEGER :: DynTotEnergyUnit
   INTEGER :: DynRingPolymerEnergyUnit
   INTEGER :: DynCentroidXUnit
   INTEGER :: DynCentroidVUnit

   !> Number of print steps of the thermalization
   INTEGER :: NrEquilPrintSteps
   !> Number of print steps of the dynamics
   INTEGER :: NrDynPrintSteps

   !> Object to write VTF trajectory file
   TYPE( VTFFile ), SAVE :: TrajectoryVTF

   !> Memory to store single trajectory output
   INTEGER, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: ParticleVelBinning

   !> Memory to store equilibration averages
   REAL, DIMENSION(:), ALLOCATABLE :: EquilKinEnergy, EquilPotEnergy
   REAL, DIMENSION(:), ALLOCATABLE :: EquilRPKinEnergy, EquilRPPotEnergy

   !> Memory to store dynamics averages
   REAL, DIMENSION(:), ALLOCATABLE :: DynKinEnergy, DynPotEnergy
   REAL, DIMENSION(:), ALLOCATABLE :: DynRPKinEnergy, DynRPPotEnergy

   ! Data formats
   !> time vs averages, decimal format
   CHARACTER(20), PARAMETER, PRIVATE :: FF = "(1F12.5,4(1F15.8,1X))"
   !> time vs averages, exponential format
   CHARACTER(20), PARAMETER, PRIVATE :: FE = "(1X,1E12.5,100(1E15.8,1X))"

   
!============================================================================================
                                       CONTAINS
!============================================================================================


   SUBROUTINE SingleTrajectoryOutput( Action )
      IMPLICIT NONE
      INTEGER, INTENT(IN)  ::  Action
      CHARACTER(2), DIMENSION(:), ALLOCATABLE :: AtomsLabels
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: BondsLogical
      INTEGER :: i, j, jStart, jEnd
      CHARACTER(50) :: OutFileName
      REAL, DIMENSION(6) :: UnitCellDim

      SELECT CASE( Action )

!******************************************************************************************************
         CASE(SETUP_OUTPUT)
!******************************************************************************************************

            ! Check if current trajectory output has correct status
            CALL ERROR( WritingCurrentTrajectory, &
                     " OutputModule.SingleTrajectoryOutput: already writing output for current trajectory" )

            ! Initialize object to print the RPMD trajectory in VTF format for VMD
            WRITE(OutFileName,"(A,I4.4,A)") "Traj_",iTraj,"_RPMDTrajectory"
            CALL VTFFile_Setup( TrajectoryVTF, OutFileName )
            ! Write header section of the VTF file 
            ALLOCATE( AtomsLabels(NAtoms*NBeads), BondsLogical(NAtoms*NBeads,NAtoms*NBeads) )
            AtomsLabels(:) = "H "
            BondsLogical(:,:) = .FALSE.
            DO i = 1, NAtoms
               DO j = 1, NBeads-1
                  BondsLogical( i+j*NAtoms, i+(j-1)*NAtoms ) = .TRUE.
               END DO
               BondsLogical( i+(NBeads-1)*NAtoms, i ) = .TRUE.
            END DO
            CALL VTFFile_WriteGeneralData( TrajectoryVTF, AtomsLabels, BondsLogical )
            DEALLOCATE( AtomsLabels, BondsLogical )

            ! Open unit to write trajectory energy
            TrajTotEnergyUnit = LookForFreeUnit()
            WRITE(OutFileName,"(A,I4.4,A)") "Traj_",iTraj,"_TotEnergy.dat"
            OPEN( FILE=OutFileName, UNIT=TrajTotEnergyUnit )
            WRITE(TrajTotEnergyUnit, "(A,I6,/)") "# E/T vs time (" // trim(TimeUnit(InputUnits)) // " "    &
                  // trim(TemperUnit(InputUnits)) // " vs " // trim(EnergyUnit(InputUnits)) // ") - trajectory # ", iTraj
            WRITE(TrajTotEnergyUnit, "(A)") "# Langevin equilibration "

            IF ( NBeads > 1 ) THEN
               ! Open unit to write trajectory energy
               TrajRingPolymerEnergyUnit = LookForFreeUnit()
               WRITE(OutFileName,"(A,I4.4,A)") "Traj_",iTraj,"_RPEnergy.dat"
               OPEN( FILE=OutFileName, UNIT=TrajRingPolymerEnergyUnit )
               WRITE(TrajRingPolymerEnergyUnit, "(A,I6,/)") "# E/T vs time (" // trim(TimeUnit(InputUnits)) // " "    &
                     // trim(TemperUnit(InputUnits)) // " vs " // trim(EnergyUnit(InputUnits)) // ") - trajectory # ", iTraj
               WRITE(TrajRingPolymerEnergyUnit, "(A)") "# Langevin equilibration "
            END IF

            ! Integer variable to check actual section of the calculation (equilibration / dynamics )
            nEquilDyn = 1
            NrOfSteps(:) = 0

            ! Velocity binning
            IF ( Out_VelDistrib ) THEN
               ! Allocate memory and initialize array
               ALLOCATE( ParticleVelBinning( NDim, Out_VelDistrib_nV, 2 ) )
               ParticleVelBinning = 0
            END IF

            ! Now update status variable
            WritingCurrentTrajectory = .TRUE.

!******************************************************************************************************
         CASE(PRINT_OUTPUT)
!******************************************************************************************************

            ! Check if current trajectory output has correct status
            CALL ERROR( .NOT. WritingCurrentTrajectory, &
                     " OutputModule.SingleTrajectoryOutput: output for current trajectory not initialized (print)" )

            ! Write trajectory snapshot to VTF output file
            UnitCellDim = GetUnitCellDimensions( )
            UnitCellDim(1:3) = UnitCellDim(1:3) * LengthConversion(InternalUnits, InputUnits)
            UnitCellDim(4:6) = UnitCellDim(4:6) * AngleConversion(InternalUnits, InputUnits)
            CALL VTFFile_WriteTimeStep( TrajectoryVTF, X*LengthConversion(InternalUnits, InputUnits), UnitCellDim )

            ! Write energy values to the total energy file
            WRITE(TrajTotEnergyUnit,800) Time*TimeConversion(InternalUnits, InputUnits),           &
                         KinEnergy*EnergyConversion(InternalUnits, InputUnits),                    &
                         PotEnergy*EnergyConversion(InternalUnits, InputUnits),                    &
                         TotEnergy*EnergyConversion(InternalUnits, InputUnits),                    &
                         2.0*KinEnergy/NDim*TemperatureConversion(InternalUnits, InputUnits)

            IF ( NBeads > 1 ) THEN
               ! Write energy values to the total energy file
               WRITE(TrajRingPolymerEnergyUnit,800) Time*TimeConversion(InternalUnits, InputUnits),  &
                           RPKinEnergy*EnergyConversion(InternalUnits, InputUnits),                  &
                           RPPotEnergy*EnergyConversion(InternalUnits, InputUnits),                  &
                           RPTotEnergy*EnergyConversion(InternalUnits, InputUnits),                  &
                           2.0*RPKinEnergy/NDim/NBeads*TemperatureConversion(InternalUnits, InputUnits)
            END IF

            ! Increment array for the kinetic energy binning
            IF ( Out_VelDistrib ) THEN
               DO i = 1, NDim
                  IF ( ABS(CentroidVel(i)) < Out_VelDistrib_nV*Out_VelDistrib_DV ) THEN
                     j = FLOOR( ABS(CentroidVel(i))/Out_VelDistrib_DV )+1
                     ParticleVelBinning( i, j, nEquilDyn ) = ParticleVelBinning( i, j, nEquilDyn ) + 1
                  END IF
               END DO
            END IF

            ! Increment nr of steps of current section
            NrOfSteps(nEquilDyn) = NrOfSteps(nEquilDyn) + 1

!******************************************************************************************************
         CASE(CLOSE_OUTPUT)
!******************************************************************************************************

            ! Check if current trajectory output has correct status
            CALL ERROR( .NOT. WritingCurrentTrajectory, &
                     " OutputModule.SingleTrajectoryOutput: output for current trajectory not initialized (close)" )

            ! Close trajectory VTF file
            CALL VTFFile_Dispose( TrajectoryVTF )

            ! Close total energy file
            CLOSE( UNIT=TrajTotEnergyUnit )
            IF ( NBeads > 1 )  CLOSE( UNIT=TrajRingPolymerEnergyUnit )

            ! Write kinetic energy binning to output file
            IF ( Out_VelDistrib ) THEN
               ! Open unit to write kin energy distribution
               TrajVelDistrUnit = LookForFreeUnit()
               WRITE(OutFileName,"(A,I4.4,A)") "Traj_",iTraj,"_VelDistrib.dat"
               OPEN( FILE=OutFileName, UNIT=TrajVelDistrUnit )
               ! Write particle distribution during equilibration
               WRITE(TrajVelDistrUnit, "(/,A,I6)") "# Pdof vs v (" // trim(VelocityUnit(InputUnits)) // ") - equil # ", iTraj
               DO i = 1, Out_VelDistrib_nV
                  WRITE(TrajVelDistrUnit,801) Out_VelDistrib_DV*REAL(i-0.5)*VelocityConversion(InternalUnits,InputUnits),  &
                                                           REAL(ParticleVelBinning(:,i,1))/REAL(NrOfSteps(1))
               END DO
               WRITE(TrajVelDistrUnit, "(/,A,I6)") "# Pdof vs v (" // trim(VelocityUnit(InputUnits)) // ") - dyn # ", iTraj         
               DO i = 1, Out_VelDistrib_nV
                  WRITE(TrajVelDistrUnit,801) Out_VelDistrib_DV*REAL(i-0.5)*VelocityConversion(InternalUnits,InputUnits),  &
                                                           REAL(ParticleVelBinning(:,i,2))/REAL(NrOfSteps(2))
               END DO
               WRITE(TrajVelDistrUnit, "(/,A,I6)") "# Pfull vs v (" // trim(VelocityUnit(InputUnits)) // ") - equil # ", iTraj       
               DO i = 1, Out_VelDistrib_nV
                  WRITE(TrajVelDistrUnit,801) Out_VelDistrib_DV*REAL(i-0.5)*VelocityConversion(InternalUnits,InputUnits),  &
                                                           REAL(SUM(ParticleVelBinning(:,i,1)))/REAL(NrOfSteps(1)*NDim)
               END DO
               WRITE(TrajVelDistrUnit, "(/,A,I6)") "# Pfull vs v (" // trim(VelocityUnit(InputUnits)) // ") - dyn # ", iTraj        
               DO i = 1, Out_VelDistrib_nV
                  WRITE(TrajVelDistrUnit,801) Out_VelDistrib_DV*REAL(i-0.5)*VelocityConversion(InternalUnits,InputUnits),  &
                                                           REAL(SUM(ParticleVelBinning(:,i,2)))/REAL(NrOfSteps(2)*NDim)
               END DO
            END IF
            CLOSE( TrajVelDistrUnit )
            DEALLOCATE( ParticleVelBinning )

            ! Now update status variable
            WritingCurrentTrajectory = .FALSE.

!******************************************************************************************************
         CASE(DIVIDE_EQ_DYN)
!******************************************************************************************************

            ! Go from the equil section to the dynamics one
            nEquilDyn = 2

            ! Check if current trajectory output has correct status
            CALL ERROR( .NOT. WritingCurrentTrajectory, &
                     " OutputModule.SingleTrajectoryOutput: output for current trajectory not initialized (divide)" )

            WRITE(TrajTotEnergyUnit, "(/,A)") "# Microcanonical dynamics "
            IF ( NBeads > 1 )  WRITE(TrajRingPolymerEnergyUnit, "(/,A)") "# Microcanonical dynamics " 

         CASE DEFAULT
            CALL AbortWithError( " OutputModule.SingleTrajectoryOutput: given action is not defined ")
      END SELECT

      ! Format for output printing
      800 FORMAT( 1F12.5,4(1F15.8,1X) )
      801 FORMAT( 1F12.5,1000(1F8.4,1X) )

   END SUBROUTINE SingleTrajectoryOutput

!============================================================================================

   SUBROUTINE DynAveragesOutput( Action )
      IMPLICIT NONE
      INTEGER, INTENT(IN)  ::  Action
      INTEGER :: i
      REAL    :: Time

      SELECT CASE( Action )

!******************************************************************************************************
         CASE( SETUP_AVERAGES )
!******************************************************************************************************

            ! Check if current trajectory output has correct status
            CALL ERROR( DynamicsAveragesStatus /= IS_NOT_SET_UP, &
                     " OutputModule.DynAveragesOutput: already writing output for dynamical averages" )

            NrDynPrintSteps   = NrSteps / PrintStepInterval + 1

            ! Allocate memory
            ALLOCATE( DynKinEnergy(NrDynPrintSteps), DynPotEnergy(NrDynPrintSteps))
            IF (NBeads > 1) ALLOCATE( DynRPKinEnergy(NrDynPrintSteps), DynRPPotEnergy(NrDynPrintSteps)) 

            ! Initialize arrays
            DynKinEnergy(:) = 0.0
            DynPotEnergy(:) = 0.0
            IF ( NBeads > 1 ) THEN
               DynRPKinEnergy(:) = 0.0
               DynRPPotEnergy(:) = 0.0
            END IF

            ! Now update status variable
            DynamicsAveragesStatus = IS_SET_UP

!******************************************************************************************************
         CASE( UPDATE_AVERAGES )
!******************************************************************************************************

            ! Check if current trajectory output has correct status
            CALL ERROR( DynamicsAveragesStatus /= IS_SET_UP .AND. DynamicsAveragesStatus /= HAS_DATA, &
                     " OutputModule.DynAveragesOutput: output for dynamical averages not initialized" )

            ! Energy averages
            DynKinEnergy(kStep) = DynKinEnergy(kStep) + KinEnergy
            DynPotEnergy(kStep) = DynPotEnergy(kStep) + PotEnergy
            IF ( NBeads > 1 ) THEN
               DynRPKinEnergy(kStep) = DynRPKinEnergy(kStep) + RPKinEnergy
               DynRPPotEnergy(kStep) = DynRPPotEnergy(kStep) + RPPotEnergy
            END IF

            ! Now update status variable
            DynamicsAveragesStatus = HAS_DATA


!******************************************************************************************************
         CASE( FINALIZE_AVERAGES )
!******************************************************************************************************

#if defined(WITH_MPI)
            CALL MyMPI_AllReduceMaxValue(DynamicsAveragesStatus)
#endif

            ! Check if current trajectory output has correct status
            CALL ERROR( DynamicsAveragesStatus /= HAS_DATA, &
                     " OutputModule.DynAveragesOutput: no data to print" )

#if defined(WITH_MPI)
	    CALL MyMPI_SumToMaster( DynKinEnergy )
	    CALL MyMPI_SumToMaster( DynPotEnergy )
            IF ( NBeads > 1 ) THEN
	       CALL MyMPI_SumToMaster( DynRPKinEnergy )
	       CALL MyMPI_SumToMaster( DynRPPotEnergy )
	    END IF
#endif

	    __MPI_OnlyMasterBEGIN
            ! Normalize by number of trajectories
            DynKinEnergy(:) = DynKinEnergy(:) / NrTrajs
            DynPotEnergy(:) = DynPotEnergy(:) / NrTrajs
            IF ( NBeads > 1 ) THEN
               DynRPKinEnergy(:) = DynRPKinEnergy(:) / NrTrajs
               DynRPPotEnergy(:) = DynRPPotEnergy(:) / NrTrajs
            END IF            
	    __MPI_OnlyMasterEND	

            ! Now update status variable
            DynamicsAveragesStatus = IS_FINALIZED


!******************************************************************************************************
         CASE( PRINT_AVERAGES_AND_DISPOSE )
!******************************************************************************************************

            ! Check if current trajectory output has correct status
            CALL ERROR( DynamicsAveragesStatus /= IS_FINALIZED, &
                     " OutputModule.DynAveragesOutput: data has not been finalized" )

	    __MPI_OnlyMasterBEGIN
            ! Open unit to write average energy
            DynTotEnergyUnit = LookForFreeUnit()
            OPEN( FILE="Dyn_TotEnergy.dat", UNIT=DynTotEnergyUnit )
            WRITE(DynTotEnergyUnit, "(A,I6,/)") "# E/T vs time (" // trim(TimeUnit(InputUnits)) // " "    &
                  // trim(TemperUnit(InputUnits)) // " vs " // trim(EnergyUnit(InputUnits)) // ") - trajectory # ", iTraj

            IF ( NBeads > 1 ) THEN
               ! Open unit to write ring polymer average energy
               DynRingPolymerEnergyUnit = LookForFreeUnit()
               OPEN( FILE="Dyn_RPEnergy.dat", UNIT=DynRingPolymerEnergyUnit )
               WRITE(DynRingPolymerEnergyUnit, "(A,I6,/)") "# E/T vs time (" // trim(TimeUnit(InputUnits)) // " "    &
                     // trim(TemperUnit(InputUnits)) // " vs " // trim(EnergyUnit(InputUnits)) // ") - trajectory # ", iTraj
            END IF

            DO i = 1, NrDynPrintSteps
               
               Time = REAL((i-1)*PrintStepInterval)*TimeStep

               ! Write energy values to the total energy file
               WRITE(DynTotEnergyUnit,800) Time*TimeConversion(InternalUnits, InputUnits),                        &
                           DynKinEnergy(i)*EnergyConversion(InternalUnits, InputUnits),                         &
                           DynPotEnergy(i)*EnergyConversion(InternalUnits, InputUnits),                         &
                           (DynKinEnergy(i)+DynPotEnergy(i))*EnergyConversion(InternalUnits, InputUnits),     &
                           2.0*DynKinEnergy(i)/NDim*TemperatureConversion(InternalUnits, InputUnits)

               IF ( NBeads > 1 ) THEN
                  ! Write energy values to the ring polymer energy file
                  WRITE(DynRingPolymerEnergyUnit,800) Time*TimeConversion(InternalUnits, InputUnits),                 &
                              DynRPKinEnergy(i)*EnergyConversion(InternalUnits, InputUnits),                        &
                              DynRPPotEnergy(i)*EnergyConversion(InternalUnits, InputUnits),                        &
                              (DynRPKinEnergy(i)+DynRPPotEnergy(i))*EnergyConversion(InternalUnits, InputUnits),  &
                              2.0*DynRPKinEnergy(i)/NDim/NBeads*TemperatureConversion(InternalUnits, InputUnits)
               END IF
            END DO

            ! Close files
            CLOSE( DynTotEnergyUnit )
            CLOSE( DynRingPolymerEnergyUnit )

            ! Deallocate memory
            DEALLOCATE( DynKinEnergy, DynPotEnergy )
            IF (NBeads > 1) DEALLOCATE( DynRPKinEnergy, DynRPPotEnergy ) 
	    __MPI_OnlyMasterEND

            ! Now update status variable
            DynamicsAveragesStatus = IS_NOT_SET_UP

         CASE DEFAULT
            CALL AbortWithError( " OutputModule.DynAveragesOutput: given action is not defined ")
      END SELECT

      ! Format for output printing
      800 FORMAT( 1F12.5,4(1F15.8,1X) )

   END SUBROUTINE DynAveragesOutput

!============================================================================================

   SUBROUTINE EquilAveragesOutput( Action )
      IMPLICIT NONE
      INTEGER, INTENT(IN)  ::  Action
      INTEGER :: i
      REAL    :: Time

      SELECT CASE( Action )

!******************************************************************************************************
         CASE( SETUP_AVERAGES )
!******************************************************************************************************

            ! Check if current trajectory output has correct status
            CALL ERROR( EquilibrationAveragesStatus /= IS_NOT_SET_UP, &
                     " OutputModule.EquilAveragesOutput: already writing output for dynamical averages" )

            NrEquilPrintSteps = EquilNrSteps / EquilPrintStepInterval + 1

            ! Allocate memory
            ALLOCATE( EquilKinEnergy(NrEquilPrintSteps), EquilPotEnergy(NrEquilPrintSteps))
            IF (NBeads > 1) ALLOCATE( EquilRPKinEnergy(NrEquilPrintSteps), EquilRPPotEnergy(NrEquilPrintSteps)) 

            ! Initialize arrays
            EquilKinEnergy(:) = 0.0
            EquilPotEnergy(:) = 0.0
            IF ( NBeads > 1 ) THEN
               EquilRPKinEnergy(:) = 0.0
               EquilRPPotEnergy(:) = 0.0
            END IF

            ! Now update status variable
            EquilibrationAveragesStatus = IS_SET_UP

!******************************************************************************************************
         CASE( UPDATE_AVERAGES )
!******************************************************************************************************

            ! Check if current trajectory output has correct status
            CALL ERROR( EquilibrationAveragesStatus /= IS_SET_UP .AND. EquilibrationAveragesStatus /= HAS_DATA, &
                     " OutputModule.EquilAveragesOutput: output for dynamical averages not initialized" )

            ! Energy averages
            EquilKinEnergy(kStep) = EquilKinEnergy(kStep) + KinEnergy
            EquilPotEnergy(kStep) = EquilPotEnergy(kStep) + PotEnergy
            IF ( NBeads > 1 ) THEN
               EquilRPKinEnergy(kStep) = EquilRPKinEnergy(kStep) + RPKinEnergy
               EquilRPPotEnergy(kStep) = EquilRPPotEnergy(kStep) + RPPotEnergy
            END IF

            ! Now update status variable
            EquilibrationAveragesStatus = HAS_DATA


!******************************************************************************************************
         CASE( FINALIZE_AVERAGES )
!******************************************************************************************************

#if defined(WITH_MPI)
            CALL MyMPI_AllReduceMaxValue(EquilibrationAveragesStatus)
#endif

            ! Check if current trajectory output has correct status
            CALL ERROR( EquilibrationAveragesStatus /= HAS_DATA, &
                     " OutputModule.EquilAveragesOutput: no data to print" )

#if defined(WITH_MPI)
	    CALL MyMPI_SumToMaster( EquilKinEnergy )
	    CALL MyMPI_SumToMaster( EquilPotEnergy )
            IF ( NBeads > 1 ) THEN
	       CALL MyMPI_SumToMaster( EquilRPKinEnergy )
	       CALL MyMPI_SumToMaster( EquilRPPotEnergy )
	    END IF
#endif

	    __MPI_OnlyMasterBEGIN
            ! Normalize by number of trajectories
            EquilKinEnergy(:) = EquilKinEnergy(:) / NrTrajs
            EquilPotEnergy(:) = EquilPotEnergy(:) / NrTrajs
            IF ( NBeads > 1 ) THEN
               EquilRPKinEnergy(:) = EquilRPKinEnergy(:) / NrTrajs
               EquilRPPotEnergy(:) = EquilRPPotEnergy(:) / NrTrajs
            END IF            
	    __MPI_OnlyMasterEND

            ! Now update status variable
            EquilibrationAveragesStatus = IS_FINALIZED


!******************************************************************************************************
         CASE( PRINT_AVERAGES_AND_DISPOSE )
!******************************************************************************************************

            ! Check if current trajectory output has correct status
            CALL ERROR( EquilibrationAveragesStatus /= IS_FINALIZED, &
                     " OutputModule.EquilAveragesOutput: data has not been finalized" )

	    __MPI_OnlyMasterBEGIN
            ! Open unit to write average energy
            EquilTotEnergyUnit = LookForFreeUnit()
            OPEN( FILE="Equil_TotEnergy.dat", UNIT=EquilTotEnergyUnit )
            WRITE(EquilTotEnergyUnit, "(A,I6,/)") "# E/T vs time (" // trim(TimeUnit(InputUnits)) // " "    &
                  // trim(TemperUnit(InputUnits)) // " vs " // trim(EnergyUnit(InputUnits)) // ") - trajectory # ", iTraj

            IF ( NBeads > 1 ) THEN
               ! Open unit to write ring polymer average energy
               EquilRingPolymerEnergyUnit = LookForFreeUnit()
               OPEN( FILE="Equil_RPEnergy.dat", UNIT=EquilRingPolymerEnergyUnit )
               WRITE(EquilRingPolymerEnergyUnit, "(A,I6,/)") "# E/T vs time (" // trim(TimeUnit(InputUnits)) // " "    &
                     // trim(TemperUnit(InputUnits)) // " vs " // trim(EnergyUnit(InputUnits)) // ") - trajectory # ", iTraj
            END IF

            DO i = 1, NrEquilPrintSteps
               
               Time = REAL((i-1)*EquilPrintStepInterval)*EquilTimeStep

               ! Write energy values to the total energy file
               WRITE(EquilTotEnergyUnit,800) Time*TimeConversion(InternalUnits, InputUnits),                        &
                           EquilKinEnergy(i)*EnergyConversion(InternalUnits, InputUnits),                         &
                           EquilPotEnergy(i)*EnergyConversion(InternalUnits, InputUnits),                         &
                           (EquilKinEnergy(i)+EquilPotEnergy(i))*EnergyConversion(InternalUnits, InputUnits),     &
                           2.0*EquilKinEnergy(i)/NDim*TemperatureConversion(InternalUnits, InputUnits)

               IF ( NBeads > 1 ) THEN
                  ! Write energy values to the ring polymer energy file
                  WRITE(EquilRingPolymerEnergyUnit,800) Time*TimeConversion(InternalUnits, InputUnits),                 &
                              EquilRPKinEnergy(i)*EnergyConversion(InternalUnits, InputUnits),                        &
                              EquilRPPotEnergy(i)*EnergyConversion(InternalUnits, InputUnits),                        &
                              (EquilRPKinEnergy(i)+EquilRPPotEnergy(i))*EnergyConversion(InternalUnits, InputUnits),  &
                              2.0*EquilRPKinEnergy(i)/NDim/NBeads*TemperatureConversion(InternalUnits, InputUnits)
               END IF
            END DO

            ! Close files
            CLOSE( EquilTotEnergyUnit )
            CLOSE( EquilRingPolymerEnergyUnit )

            ! Deallocate memory
            DEALLOCATE( EquilKinEnergy, EquilPotEnergy )
            IF (NBeads > 1) DEALLOCATE( EquilRPKinEnergy, EquilRPPotEnergy ) 
	    __MPI_OnlyMasterEND

            ! Now update status variable
            EquilibrationAveragesStatus = IS_NOT_SET_UP

         CASE DEFAULT
            CALL AbortWithError( " OutputModule.EquilAveragesOutput: given action is not defined ")
      END SELECT

      ! Format for output printing
      800 FORMAT( 1F12.5,4(1F15.8,1X) )

   END SUBROUTINE EquilAveragesOutput

!============================================================================================

END MODULE OutputModule
 
