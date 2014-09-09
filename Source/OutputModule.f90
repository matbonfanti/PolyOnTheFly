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

   PUBLIC :: SingleTrajectoryOutput

   !> \name ACTIONS
   !> Integers number identifying the kind of action to be performed by the subroutines
   !> @{
   INTEGER, PARAMETER, PUBLIC  ::  SETUP_OUTPUT  = 1
   INTEGER, PARAMETER, PUBLIC  ::  PRINT_OUTPUT  = 2
   INTEGER, PARAMETER, PUBLIC  ::  CLOSE_OUTPUT  = 3
   INTEGER, PARAMETER, PUBLIC  ::  DIVIDE_EQ_DYN = 4
   !> @}

   !> Setup variable for the output of the current single trajectory values
   LOGICAL, SAVE :: WritingCurrentTrajectory = .FALSE.

   ! OUTPUT UNITS
   INTEGER :: TrajTotEnergyUnit
   INTEGER :: TrajRingPolymerEnergyUnit
   INTEGER :: TrajCentroidXUnit
   INTEGER :: TrajCentroidVUnit
   INTEGER :: TrajCoordEnergyUnit

   !> Object to write VTF trajectory file
   TYPE( VTFFile ), SAVE :: TrajectoryVTF

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
                           PotEnergy*NBeads*EnergyConversion(InternalUnits, InputUnits),             &
                           RPTotEnergy*EnergyConversion(InternalUnits, InputUnits),                  &
                           2.0*RPKinEnergy/NDim*TemperatureConversion(InternalUnits, InputUnits)
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

   SUBROUTINE AverageOutput( Action )
      IMPLICIT NONE
      INTEGER, INTENT(IN)  ::  Action

      SELECT CASE( Action )

         CASE(SETUP_OUTPUT)


         CASE(PRINT_OUTPUT)

         CASE(CLOSE_OUTPUT)

         CASE(DIVIDE_EQ_DYN)

         CASE DEFAULT
            CALL AbortWithError( " OutputModule.SingleTrajectoryOutput: given action is not defined ")
      END SELECT

   END SUBROUTINE AverageOutput

!============================================================================================

END MODULE OutputModule
 
