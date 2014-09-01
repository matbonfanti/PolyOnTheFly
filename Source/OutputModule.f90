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

   PRIVATE

   PUBLIC :: SetupOutput, PrintOutput, DisposeOutput

   !> Setup variable for the module
   LOGICAL, SAVE :: OutputModuleIsSetup = .FALSE.

   
!============================================================================================
                                       CONTAINS
!============================================================================================

   SUBROUTINE SetupOutput(   )
      IMPLICIT NONE

      ! exit if module is setup
      IF ( OutputModuleIsSetup ) RETURN

      ! Open output files
      
      
      ! Module is now ready
      OutputModuleIsSetup = .TRUE.
      
   END SUBROUTINE SetupOutput

!============================================================================================

   SUBROUTINE PrintOutput(  )
      IMPLICIT NONE

      ! Error if module not have been setup yet
      CALL ERROR( .NOT. OutputModuleIsSetup, " OutputModule.PrintOutput : Module not Setup" )
      
   END SUBROUTINE PrintOutput

!============================================================================================

   SUBROUTINE DisposeOutput(  )
      IMPLICIT NONE

      ! exit if module is not setup
      IF ( .NOT. OutputModuleIsSetup ) RETURN

      OutputModuleIsSetup = .FALSE.
      
   END SUBROUTINE DisposeOutput

!============================================================================================

END MODULE OutputModule
 
