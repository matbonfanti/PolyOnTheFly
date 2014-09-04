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

   PRIVATE

   PUBLIC :: SetupOutput, PrintOutput, DisposeOutput

   !> Setup variable for the module
   LOGICAL, SAVE :: OutputModuleIsSetup = .FALSE.

   ! OUTPUT UNITS

   INTEGER :: TrajTotEnergyUnit
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

   SUBROUTINE SetupOutput(   )
      IMPLICIT NONE
      CHARACTER(2), DIMENSION(:), ALLOCATABLE :: AtomsLabels
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: BondsLogical
      INTEGER :: i, j, jStart, jEnd
      CHARACTER(50) :: OutFileName

      ! exit if module is setup
      IF ( OutputModuleIsSetup ) RETURN

      ! Open output files
      
      ! Initialize object to print the RPMD trajectory in VTF format for VMD
      CALL VTFFile_Setup( TrajectoryVTF, "RPMDTrajectory" )
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
      WRITE(OutFileName,"(A,I4.4,A)") "Traj_",1,"_TotEnergy.dat"
      OPEN( FILE=OutFileName, UNIT=TrajTotEnergyUnit )
      WRITE(TrajTotEnergyUnit, "(A,I6,/)") "# E/T vs time (" // trim(TimeUnit(InputUnits)) // " "    &
                  // trim(TemperUnit(InputUnits)) // " vs " // trim(EnergyUnit(InputUnits)) // ") - trajectory # ", 1

      ! Module is now ready
      OutputModuleIsSetup = .TRUE.
      
   END SUBROUTINE SetupOutput

!============================================================================================

   SUBROUTINE PrintOutput( Time )
      IMPLICIT NONE
      REAL, INTENT(IN) :: Time

      ! Error if module not have been setup yet
      CALL ERROR( .NOT. OutputModuleIsSetup, " OutputModule.PrintOutput : Module not Setup" )
      
      CALL VTFFile_WriteTimeStep( TrajectoryVTF, X, (/ 10., 10., 10., 90., 90., 90. /) )

       WRITE(TrajTotEnergyUnit,800) Time*TimeConversion(InternalUnits, InputUnits),           &
                                      KinEnergy*EnergyConversion(InternalUnits, InputUnits),                    &
                                      PotEnergy*EnergyConversion(InternalUnits, InputUnits),                    &
                                      TotEnergy*EnergyConversion(InternalUnits, InputUnits),                    &
                                      2.0*KinEnergy/NDim*TemperatureConversion(InternalUnits, InputUnits)

      800 FORMAT( 1F12.5,4(1F15.8,1X) )

   END SUBROUTINE PrintOutput

!============================================================================================

   SUBROUTINE DisposeOutput(  )
      IMPLICIT NONE

      ! exit if module is not setup
      IF ( .NOT. OutputModuleIsSetup ) RETURN

      CALL VTFFile_Dispose( TrajectoryVTF )

      CLOSE( UNIT=TrajTotEnergyUnit )


      OutputModuleIsSetup = .FALSE.
      
   END SUBROUTINE DisposeOutput

!============================================================================================

END MODULE OutputModule
 
