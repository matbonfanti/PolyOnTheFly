!***************************************************************************************
!*                           MODULE VTFFileModule
!***************************************************************************************
!
!>  \brief     Write VTF trajectory file
!>  \details   This module defines an object to write a molecular trajectory to 
!>             an output file in the VTF format, which is suited to 
!>             visualize the trajectory with VMD.
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
MODULE VTFFileModule
#include "preprocessoptions.cpp"

   PRIVATE

   PUBLIC :: VTFFile_Setup, VTFFile_WriteGeneralData, VTFFile_WriteTimeStep, VTFFile_Dispose
   PUBLIC :: VTFFile

   !> VTF file data
   TYPE :: VTFFile
      PRIVATE
      LOGICAL              ::  isSetup        !< Setup variable for the object
      CHARACTER(len=30)    ::  FileName       !< The name of the file
      INTEGER              ::  Unit           !< The unit to which it is being written
      INTEGER              ::  NrAtoms        !< Number of atoms in the unit cell       
      INTEGER              ::  StepCounter    !< Counter for the steps of the trajectory
   END TYPE VTFFile

   
!============================================================================================
                                       CONTAINS
!============================================================================================

   SUBROUTINE VTFFile_Setup( VTFFileObj, Name  )
      IMPLICIT NONE
      TYPE(VTFFile), INTENT(INOUT) :: VTFFileObj
      CHARACTER(*), INTENT(IN)     :: Name

      ! exit if module is setup
      IF ( VTFFileObj%isSetup ) THEN
#if defined(LOG_FILE)
         __OPEN_LOG_FILE; 
         WRITE(__LOG_UNIT,*) " VTFFileModule: error writing data for file ", TRIM(Name)
         WRITE(__LOG_UNIT,*) "                current object is is already initialized";
         __CLOSE_LOG_FILE
#endif
         RETURN
      END IF

      ! Initialize variables
      VTFFileObj%FileName = TRIM(ADJUSTL(Name)) // ".vtf"
      VTFFileObj%Unit = LookForFreeUnit()
      VTFFileObj%StepCounter = 0

      ! Open output files
      OPEN( UNIT = VTFFileObj%Unit, FILE = VTFFileObj%FileName )
      
      ! Module is now ready
      VTFFileObj%isSetup = .TRUE.
      
#if defined(LOG_FILE)
      __OPEN_LOG_FILE
      WRITE(__LOG_UNIT,*) " VTFFileModule: object is now initialized"
      WRITE(__LOG_UNIT,*) "                writing file ",TRIM(VTFFileObj%FileName)," to unit ", VTFFileObj%Unit
      __CLOSE_LOG_FILE
#endif

   END SUBROUTINE VTFFile_Setup

!============================================================================================

   SUBROUTINE VTFFile_WriteGeneralData( VTFFileObj, Atoms, Bonds )
      IMPLICIT NONE
      TYPE(VTFFile), INTENT(INOUT)                  :: VTFFileObj
      CHARACTER(2), DIMENSION(:), INTENT(IN)        :: Atoms 
      LOGICAL, DIMENSION(:,:), INTENT(IN), OPTIONAL :: Bonds
      CHARACTER(20) :: iString, jString
      INTEGER :: i, j

      ! Error if module not have been setup yet
      CALL ERROR( .NOT. VTFFileObj%isSetup, " VTFFile_WriteGeneralData : current object is not set up " )
      
      ! Store the number of atoms for later consistency check
      VTFFileObj%NrAtoms = SIZE( Atoms )

      ! If bond matrix is given, check its dimension
      IF ( PRESENT(Bonds) ) THEN
         CALL ERROR( SIZE( Bonds,1) /= VTFFileObj%NrAtoms, " VTFFile_WriteGeneralData: inconsistent 1st dimension of bond matrix" )
         CALL ERROR( SIZE( Bonds,2) /= VTFFileObj%NrAtoms, " VTFFile_WriteGeneralData: inconsistent 2nd dimension of bond matrix" )
      END IF

      ! Write comment line
      WRITE( VTFFileObj%Unit, "(A,/)") " # VTF file written by FORTRAN module 'VTFFileModule' "
      ! Write list of atoms 
      WRITE( VTFFileObj%Unit, "(A)" )  " # STRUCTURE BLOCK "
      DO i = 1,  VTFFileObj%NrAtoms
         WRITE( VTFFileObj%Unit, "(A,I3,A,A)" ) " atom ",i-1,"   name ",Atoms(i)
      END DO

#if defined(LOG_FILE)
      __OPEN_LOG_FILE
      WRITE(__LOG_UNIT,*) " VTFFileModule: atoms list has been written to file ",TRIM(VTFFileObj%FileName)
      __CLOSE_LOG_FILE
#endif

      IF ( PRESENT(Bonds) ) THEN

         ! Write list of bonds (read from lower triangle of the matrix)
         DO j = 1,  VTFFileObj%NrAtoms
            DO i = j+1,  VTFFileObj%NrAtoms
               WRITE(iString,*) i-1
               WRITE(jString,*) j-1
               IF ( Bonds(i,j) ) WRITE( VTFFileObj%Unit, *) "bond ",TRIM(ADJUSTL(iString)),":",TRIM(ADJUSTL(jString))
            END DO
         END DO

#if defined(LOG_FILE)
         __OPEN_LOG_FILE
         WRITE(__LOG_UNIT,*) " VTFFileModule: bonds list has been written to file ",TRIM(VTFFileObj%FileName)
         __CLOSE_LOG_FILE
      END IF
#endif

      WRITE( VTFFileObj%Unit, "(/,A,/)" )  " # TIME STEP BLOCK "

   END SUBROUTINE VTFFile_WriteGeneralData

!============================================================================================

   SUBROUTINE VTFFile_WriteTimeStep( VTFFileObj, Coordinates, CellSize  )
      IMPLICIT NONE
      TYPE(VTFFile), INTENT(INOUT)       :: VTFFileObj
      REAL, DIMENSION(:)               :: Coordinates
      REAL, DIMENSION(6), OPTIONAL       :: CellSize
      INTEGER :: i

      ! Error if module not have been setup yet
      CALL ERROR( .NOT. VTFFileObj%isSetup, " VTFFile_WriteTimeStep : current object is not set up " )
      
      ! check dimensions of the coordinates array
      CALL ERROR( SIZE(Coordinates) /= 3*VTFFileObj%NrAtoms, " VTFFile_WriteTimeStep: inconsistent dim of coords matrix" )

      ! Increment the counter of the timesteps and write comment to VTF file
      VTFFileObj%StepCounter = VTFFileObj%StepCounter + 1
      WRITE( VTFFileObj%Unit, "(/,A,I6.6)" )  " # time step number ", VTFFileObj%StepCounter
      ! Write timestep directive
      WRITE( VTFFileObj%Unit, "(A)" )  " timestep "

      IF ( PRESENT(CellSize) ) THEN
         ! Write unit cell size in ( a, b, c, alpha, beta, gamma) format
         WRITE( VTFFileObj%Unit, "(A,3F16.8,3F12.4)") " pbc   ", CellSize(:)
      END IF

      DO i = 1, VTFFileObj%NrAtoms
         ! WRite coordinates of the i-th atoms
         WRITE( VTFFileObj%Unit, "(3F15.6)") Coordinates( 3*(i-1)+1 : 3*i )
      END DO

#if defined(LOG_FILE)
      __OPEN_LOG_FILE
      WRITE(__LOG_UNIT,*) " VTFFileModule: atomic coordinates have been written to file ",TRIM(VTFFileObj%FileName)
      WRITE(__LOG_UNIT,*) "                for time step number ", VTFFileObj%StepCounter
      __CLOSE_LOG_FILE
#endif

   END SUBROUTINE VTFFile_WriteTimeStep

!============================================================================================

   SUBROUTINE VTFFile_Dispose( VTFFileObj  )
      IMPLICIT NONE
      TYPE(VTFFile), INTENT(INOUT)       :: VTFFileObj

      ! exit if module is not setup
      IF ( .NOT. VTFFileObj%isSetup ) THEN
#if defined(LOG_FILE)
         __OPEN_LOG_FILE
         WRITE(__LOG_UNIT,*) " VTFFileModule: error disposing data, current object is not yet initialized"
         __CLOSE_LOG_FILE
#endif
         RETURN
      END IF

      CLOSE( UNIT = VTFFileObj%Unit )

      VTFFileObj%isSetup = .FALSE.
      
   END SUBROUTINE VTFFile_Dispose

!============================================================================================

END MODULE VTFFileModule
 
