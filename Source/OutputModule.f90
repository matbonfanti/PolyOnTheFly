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

   PRIVATE

   PUBLIC :: SetupOutput, PrintOutput, DisposeOutput

   !> Setup variable for the module
   LOGICAL, SAVE :: OutputModuleIsSetup = .FALSE.

   !> Object to write VTF trajectory file
   TYPE( VTFFile ), SAVE :: TrajectoryVTF

   
!============================================================================================
                                       CONTAINS
!============================================================================================

   SUBROUTINE SetupOutput(   )
      IMPLICIT NONE
      CHARACTER(2), DIMENSION(:), ALLOCATABLE :: AtomsLabels
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: BondsLogical
      INTEGER :: i, j, jStart, jEnd

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

      ! Module is now ready
      OutputModuleIsSetup = .TRUE.
      
   END SUBROUTINE SetupOutput

!============================================================================================

   SUBROUTINE PrintOutput(  )
      IMPLICIT NONE

      ! Error if module not have been setup yet
      CALL ERROR( .NOT. OutputModuleIsSetup, " OutputModule.PrintOutput : Module not Setup" )
      
      CALL VTFFile_WriteTimeStep( TrajectoryVTF, X, (/ 10., 10., 10., 90., 90., 90. /) )

   END SUBROUTINE PrintOutput

!============================================================================================

   SUBROUTINE DisposeOutput(  )
      IMPLICIT NONE

      ! exit if module is not setup
      IF ( .NOT. OutputModuleIsSetup ) RETURN

      CALL VTFFile_Dispose( TrajectoryVTF )

      OutputModuleIsSetup = .FALSE.
      
   END SUBROUTINE DisposeOutput

!============================================================================================

END MODULE OutputModule
 
