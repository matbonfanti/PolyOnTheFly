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
!>  \arg 
!
!>  \todo          ____________________________
!>                 
!***************************************************************************************
MODULE OutputModule
#include "preprocessoptions.cpp"
   USE SharedData
   USE VTFFileModule
   USE UnitConversion

   ! LIST OF OUTPUT FILES THAT ARE WRITTEN BY THE FOLLOWING SUBROUTINES:

   ! A) SINGLE TRAJECTORY FILES ( 1 file per each trajectory )
   !    1) VTF trajectory file (trajectory snapshots in VTF format, using VTFFileModule)
   !    2) Total energy file (istantaneous value of kin, pot, total energy and temperature of the trajectory)

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

   ! OUTPUT UNITS for single trajectory output
   INTEGER :: TrajTotEnergyUnit
   INTEGER :: TrajRingPolymerEnergyUnit
   INTEGER :: TrajCentroidXUnit
   INTEGER :: TrajCentroidVUnit

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

            ! Now update status variable
            WritingCurrentTrajectory = .TRUE.

!******************************************************************************************************
         CASE(PRINT_OUTPUT)
!******************************************************************************************************

            ! Check if current trajectory output has correct status
            CALL ERROR( .NOT. WritingCurrentTrajectory, &
                     " OutputModule.SingleTrajectoryOutput: output for current trajectory not initialized (print)" )

            ! Write trajectory snapshot to VTF output file
            CALL VTFFile_WriteTimeStep( TrajectoryVTF, X, (/ 10., 10., 10., 90., 90., 90. /) )

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

            ! Now update status variable
            WritingCurrentTrajectory = .FALSE.

!******************************************************************************************************
         CASE(DIVIDE_EQ_DYN)
!******************************************************************************************************

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

            ! Check if current trajectory output has correct status
            CALL ERROR( DynamicsAveragesStatus /= HAS_DATA, &
                     " OutputModule.DynAveragesOutput: no data to print" )

            ! Normalize by number of trajectories
            DynKinEnergy(:) = DynKinEnergy(:) / NrTrajs
            DynPotEnergy(:) = DynPotEnergy(:) / NrTrajs
            IF ( NBeads > 1 ) THEN
               DynRPKinEnergy(:) = DynRPKinEnergy(:) / NrTrajs
               DynRPPotEnergy(:) = DynRPPotEnergy(:) / NrTrajs
            END IF            

            ! Now update status variable
            DynamicsAveragesStatus = IS_FINALIZED


!******************************************************************************************************
         CASE( PRINT_AVERAGES_AND_DISPOSE )
!******************************************************************************************************

            ! Check if current trajectory output has correct status
            CALL ERROR( DynamicsAveragesStatus /= IS_FINALIZED, &
                     " OutputModule.DynAveragesOutput: data has not been finalized" )

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

            ! Check if current trajectory output has correct status
            CALL ERROR( EquilibrationAveragesStatus /= HAS_DATA, &
                     " OutputModule.EquilAveragesOutput: no data to print" )

            ! Normalize by number of trajectories
            EquilKinEnergy(:) = EquilKinEnergy(:) / NrTrajs
            EquilPotEnergy(:) = EquilPotEnergy(:) / NrTrajs
            IF ( NBeads > 1 ) THEN
               EquilRPKinEnergy(:) = EquilRPKinEnergy(:) / NrTrajs
               EquilRPPotEnergy(:) = EquilRPPotEnergy(:) / NrTrajs
            END IF            

            ! Now update status variable
            EquilibrationAveragesStatus = IS_FINALIZED


!******************************************************************************************************
         CASE( PRINT_AVERAGES_AND_DISPOSE )
!******************************************************************************************************

            ! Check if current trajectory output has correct status
            CALL ERROR( EquilibrationAveragesStatus /= IS_FINALIZED, &
                     " OutputModule.EquilAveragesOutput: data has not been finalized" )

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
 
