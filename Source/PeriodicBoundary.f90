!***************************************************************************************
!*                           MODULE PeriodicBoundary
!***************************************************************************************
!
!>  \brief     Periodic boundary conditions
!>  \details   This module defines all the relevant subroutines to implement 
!>             periodic boundary conditions. When using the module without
!>             setting it up, a non periodic system is assumed.
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
!>  \arg 7 November 2014 : setup is modified, now only one matrix is given
!>                         as input, with columns corresponding to the 
!>                         unit cell vectors
!>  \arg 7 November 2014 : added logical function to check if PBC are defined
!>  \arg 7 November 2014 : added datatype and subroutine to define the 
!>                         neighbour cells withing a given cutoff value
!
!>  \todo          ____________________________
!>                 
!***************************************************************************************

MODULE PeriodicBoundary
#include "preprocessoptions.cpp"
   USE FFTWrapper

   PRIVATE

   PUBLIC :: PBC_Setup, PBC_Dispose, PBC_SystemIsPeriodic
   PUBLIC :: PBC_BringToFirstCell, PBC_SetNearTranslations

   PUBLIC :: PBC_NearCellTranslations

   ! Parameters to define maximum number of near cells to include in a summation
   INTEGER, PARAMETER  ::  CellsMaxNr = 5
   
   ! Datatype to define the cell translations within a given cutoff
   TYPE PBC_NearCellTranslations
      INTEGER                       :: Nr
      REAL, DIMENSION(:,:), POINTER :: TranslVectors
      LOGICAL                       :: IsSetup
   END TYPE PBC_NearCellTranslations
   
   !> Vectors of the unit cell
   REAL, DIMENSION(3,3) :: UnitCellVector

   !> Matrix transformation for coordinates
   REAL, DIMENSION(3,3) :: Fractional2Cartesian, Cartesian2Fractional

   !> Setup variable for the module
   LOGICAL, SAVE :: ModuleisSetup = .FALSE.

   
!============================================================================================
                                       CONTAINS
!============================================================================================


   !*******************************************************************************
   !                                   PBC_Setup
   !*******************************************************************************
   !> Setup module: get from input values of the unit cell vectors, and define
   !> the transformations from cartesian coordinates to fractional coordinates
   !> and viceversa. Output is written to the standard log file.
   !>
   !> @param      CellVectors, 3x3 array with unit cell vectors as columns
   !*******************************************************************************
   
   SUBROUTINE PBC_Setup( CellVectors  )
      IMPLICIT NONE
      REAL, DIMENSION(3,3), INTENT(IN) :: CellVectors

      ! warning if module is setup
      IF ( ModuleisSetup ) THEN
#if defined(LOG_FILE)
         __OPEN_LOG_FILE
         WRITE(__LOG_UNIT,*) " PeriodicBoundary: module is already setup. Overwriting data "
         __CLOSE_LOG_FILE
#endif
      END IF

      ! Store Unit Cell Vectors
      UnitCellVector(:,1) = CellVectors(:,1)
      UnitCellVector(:,2) = CellVectors(:,2)
      UnitCellVector(:,3) = CellVectors(:,3)

      ! Define matrix transformations to and from fractional coordinates
      Fractional2Cartesian = UnitCellVector(:,:)
      Cartesian2Fractional = TheOneWithInverseMatrix( Fractional2Cartesian, 3 ) 

      ! Module is now ready
      ModuleisSetup = .TRUE.
      
#if defined(LOG_FILE)
      __OPEN_LOG_FILE
      WRITE(__LOG_UNIT,*) " PeriodicBoundary: module is now initialized"
      WRITE(__LOG_UNIT,*) " PeriodicBoundary: vector A - ", UnitCellVector(:,1)
      WRITE(__LOG_UNIT,*) " PeriodicBoundary: vector B - ", UnitCellVector(:,2)
      WRITE(__LOG_UNIT,*) " PeriodicBoundary: vector C - ", UnitCellVector(:,3)
#if defined(PASS_FRACTIONAL_COORDS) 
      WRITE(__LOG_UNIT,"(/,A)") " PeriodicBoundary: in input and output fractional coordinates are assumed"
#else
      WRITE(__LOG_UNIT,"(/,A)") " PeriodicBoundary: in input and output cartesian coordinates are assumed"
#endif
      __CLOSE_LOG_FILE
#endif

   END SUBROUTINE PBC_Setup

   
!============================================================================================


   !*******************************************************************************
   !          PBC_SystemIsPeriodic
   !*******************************************************************************
   !> Gives logical variable ModuleisSetup to check whether PBCs are defined.
   !>
   !> @returns    Logical variable, to see whether PBC are defined or not
   !*******************************************************************************
   
   LOGICAL FUNCTION PBC_SystemIsPeriodic( )
      IMPLICIT NONE
      
      PBC_SystemIsPeriodic = ModuleisSetup
   END FUNCTION PBC_SystemIsPeriodic
   
   
!============================================================================================


   !*******************************************************************************
   !                            PBC_BringToFirstCell
   !*******************************************************************************
   !> Translate coordinates to the symmetric point in the first unit cell. 
   !> If the external program use cartesian coordinates, the corresponding 
   !> fractional coordinate are computed and then the fractional coordinates 
   !> are taken and trasformed back to cartesian coordinates.
   !> The number of beads in RPMD is required to make sure that the degree of 
   !> freedom is moved to another unit cell only when the centroid of the ring
   !> polymer is outside the first unit cell. 
   !>
   !> @param      Array of dimension 3*n with the arbitrary coordinates
   !> @param      Number of replicas of each degree of freedom
   !> @returns    Array of dimension 3*n with the first unit cell coordinates.
   !*******************************************************************************
   
   SUBROUTINE PBC_BringToFirstCell( X, InputNBeads )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(INOUT) :: X
      INTEGER, INTENT(IN), OPTIONAL     :: InputNBeads
      REAL, DIMENSION(3) :: Vector
      REAL, DIMENSION(:), ALLOCATABLE :: Centroid
      INTEGER :: i, j, NBeads, NCoord

      ! Return if module not have been setup yet ( NON PERIODIC SYSTEM )
      IF ( .NOT. ModuleisSetup )  RETURN 
      
      ! Check the dimension of the vector
      CALL ERROR( MOD( size(X), 3 ) /= 0, " PeriodicBoundary.PBC_BringToFirstCell: X size is not multiple of 3 " )

      ! Check the consistency of the number of beads and store the number
      IF (PRESENT( InputNBeads )) THEN
         CALL ERROR( MOD(size(X),InputNBeads) /= 0, " PeriodicBoundary.PBC_BringToFirstCell: X size is not multiple of NBeads ")
         NBeads = InputNBeads
      ELSE
         NBeads = 1
      END IF
      NCoord = size(X)/NBeads

#if !defined(PASS_FRACTIONAL_COORDS) 
      ! Trasform to fractional coordinates
      CALL InPlaceCartesianToFractional( X )
#endif

      IF ( NBeads > 1 ) THEN

         ! Allocate memory to store the centroid of the ring polymer
         ALLOCATE( Centroid(NCoord) )

         ! Compute centroid of the RP
         Centroid(:) = 0.0
         DO i = 1, NBeads
            Centroid(:) = Centroid(:) + X( (i-1)*NCoord+1 : i*NCoord )
         END DO
         Centroid(:) = Centroid(:) / NBeads

         ! Translate ring polymer to the first unit cell
         DO i = 1, NBeads
            X( (i-1)*NCoord+1 : i*NCoord ) = X( (i-1)*NCoord+1 : i*NCoord ) - FLOOR( Centroid(:) )
         END DO

      ELSE IF ( NBeads == 1 ) THEN

         ! Translate coordinates to the first unit cell
         X( : ) = X( : ) - FLOOR( X(:) )

      END IF

#if !defined(PASS_FRACTIONAL_COORDS) 
      ! Trasform back to cartesian coordinates
      CALL InPlaceFractionalToCartesian( X )
#endif

   END SUBROUTINE PBC_BringToFirstCell


!============================================================================================


   !*******************************************************************************
   !                         PBC_SetNearTranslations
   !*******************************************************************************
   !> Given the PBC conditions previously set (no PBC assumed if setup is skipped)
   !> and given a cutoff value, the subroutine computes the set of unit cell
   !> translations within the cutoff. The nearest neighbour cells are always 
   !> included. Data is stored in the PBC_NearCellTranslations datatype.
   !>
   !> @param      NeighbourCells   TYPE(PBC_NearCellTranslations) to store translations 
   !> @param      CutOffValue      Input real with the cutoff of the translations
   !*******************************************************************************

   SUBROUTINE PBC_SetNearTranslations( NeighbourCells, CutOffValue )
      IMPLICIT NONE
      TYPE(PBC_NearCellTranslations) :: NeighbourCells
      REAL, INTENT(IN)               :: CutOffValue
      
      REAL, DIMENSION(3, (2*CellsMaxNr+1)**3 ) :: TmpNearTranslations
      INTEGER                :: i, j, k
      REAL, DIMENSION(3)     :: Vector
      REAL                   :: Distance

      ! In case data has already been used, deallocate memory
      IF ( NeighbourCells%IsSetup ) THEN
         DEALLOCATE( NeighbourCells%TranslVectors )
      END IF
      
      ! Initialize nr of cells 
      NeighbourCells%Nr = 0
      
      ! If module is not setup, it is assumed that no PBCs are present 
      IF ( .NOT. ModuleisSetup ) THEN 
      
         ! Only the 0,0,0 translation is included
         NeighbourCells%Nr = 1
         ALLOCATE( NeighbourCells%TranslVectors(3,NeighbourCells%Nr) )
         NeighbourCells%TranslVectors(:,1) =  (/ 0., 0., 0. /) 
                    
      ! Otherwise, normally define the translations, always including the nearest neighbour cells
      ELSE

         ! Cycle over many neighbour cells
         DO i = -CellsMaxNr, +CellsMaxNr
            DO j = -CellsMaxNr, +CellsMaxNr
               DO k = -CellsMaxNr, +CellsMaxNr
               
                  ! Compute the translation vectors in cartensian coordinates
                  Vector = FractionalToCartesian( (/ REAL(i), REAL(j), REAL(k) /) )
                  ! Compute the distance between the cells
                  Distance = SQRT( TheOneWithVectorDotVector( Vector, Vector ) )
                  
                  ! If cell is within cutoff or if cell is nearest neighbour, include it in the tranlation set
                  IF ( Distance < CutOffValue .OR. ( abs(i) <= 1 .AND. abs(j) <= 1 .AND. abs(k) <= 1 )) THEN
                     NeighbourCells%Nr = NeighbourCells%Nr + 1
                     TmpNearTranslations(:,NeighbourCells%Nr) =  Vector
                  END IF
                  
               END DO
            END DO
         END DO
         
         ! Allocate memory and store translation vectors in the array
         ALLOCATE( NeighbourCells%TranslVectors(3,NeighbourCells%Nr) )
         NeighbourCells%TranslVectors(:,:) =  TmpNearTranslations(:,1:NeighbourCells%Nr )
            
      END IF
      
      ! Now data is setup
      NeighbourCells%IsSetup = .TRUE.

#if defined(LOG_FILE)
         __OPEN_LOG_FILE
!          PRINT*, " Number of cells included in the summation: ", NearPeriodicImages
!          PRINT*, " "
!          DO i = 1, NearPeriodicImages
!             PRINT*, " Translation # ", i, "  Vector: ", NearTranslations(:,i)
!          END DO
!          WRITE(__LOG_UNIT,*) " PeriodicBoundary: module is already setup. Overwriting data "
         __CLOSE_LOG_FILE
#endif

   END SUBROUTINE PBC_SetNearTranslations


!============================================================================================


   !*******************************************************************************
   !                          FractionalToCartesian
   !*******************************************************************************
   !> Convert fractional coordinates to cartesian coordinates.
   !> Works only if the slab geometry has been already setup.
   !>
   !> @param      Array of dimension 3 with the fractional coordinates.
   !> @returns    Array of dimension 3 with the cartesian coordinates
   !*******************************************************************************
   
   FUNCTION FractionalToCartesian( Fractional ) RESULT( Cartesian )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN) :: Fractional
      REAL, DIMENSION(size(Fractional)) :: Cartesian
      INTEGER :: i

      ! Return if module not have been setup yet ( NON PERIODIC SYSTEM )
      IF ( .NOT. ModuleisSetup )  RETURN 

      ! Cycle over the atoms
      DO i = 1, size(Fractional)/3
         ! trasform fractional coordinates to cartesian coordinates
         Cartesian( (i-1)*3+1 : i*3 ) = TheOneWithMatrixVectorProduct( Fractional2Cartesian, Fractional( (i-1)*3+1 : i*3 ) )
      END DO

   END FUNCTION FractionalToCartesian

   !*******************************************************************************
   !                     InPlaceFractionalToCartesian
   !*******************************************************************************
   !> Modification of the previous function to perform the transformation
   !> Fractional -> Cartesian coords in place.
   !>
   !> @param X     Array of dimension 3 with the fractional coordinates.
   !*******************************************************************************
   
   SUBROUTINE InPlaceFractionalToCartesian( X ) 
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(INOUT) :: X
      INTEGER :: i

      ! Return if module not have been setup yet ( NON PERIODIC SYSTEM )
      IF ( .NOT. ModuleisSetup )  RETURN 

      ! Cycle over the atoms
      DO i = 1, size(X)/3
         ! trasform fractional coordinates to cartesian coordinates
         X( (i-1)*3+1 : i*3 ) = TheOneWithMatrixVectorProduct( Fractional2Cartesian, X( (i-1)*3+1 : i*3 ) )
      END DO

   END SUBROUTINE InPlaceFractionalToCartesian


   !*******************************************************************************
   !                             CartesianToFractional
   !*******************************************************************************
   !> Convert cartesian coordinates to fractional coordinates.
   !> Works only if the slab geometry has been already setup.
   !>
   !> @param      Array of dimension 3 with the cartesian coordinates
   !> @returns    Array of dimension 3 with the fractional coordinates.
   !*******************************************************************************
   
   FUNCTION CartesianToFractional( Cartesian ) RESULT( Fractional )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN) :: Cartesian
      REAL, DIMENSION(size(Cartesian)) :: Fractional
      INTEGER :: i

      ! Return if module not have been setup yet ( NON PERIODIC SYSTEM )
      IF ( .NOT. ModuleisSetup )  RETURN 

      ! Cycle over the atoms
      DO i = 1, size(Cartesian)/3
         ! trasform cartesian coordinates to fractional coordinates
         Fractional( (i-1)*3+1 : i*3 ) = TheOneWithMatrixVectorProduct( Cartesian2Fractional, Cartesian( (i-1)*3+1 : i*3 ) )
      END DO

   END FUNCTION CartesianToFractional

   !*******************************************************************************
   !                     InPlaceCartesianToFractional
   !*******************************************************************************
   !> Modification of the previous function to perform the transformation
   !> Cartesian -> Fractional coords in place.
   !>
   !> @param X     Array of dimension 3 with the cartesian coordinates.
   !*******************************************************************************
   
   SUBROUTINE InPlaceCartesianToFractional( X ) 
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(INOUT) :: X
      INTEGER :: i

      ! Return if module not have been setup yet ( NON PERIODIC SYSTEM )
      IF ( .NOT. ModuleisSetup )  RETURN 

      ! Cycle over the atoms
      DO i = 1, size(X)/3
         ! trasform cartesian coordinates to fractional coordinates
         X( (i-1)*3+1 : i*3 ) = TheOneWithMatrixVectorProduct( Cartesian2Fractional, X( (i-1)*3+1 : i*3 ) )
      END DO

   END SUBROUTINE InPlaceCartesianToFractional

   
!============================================================================================


   SUBROUTINE PBC_Dispose(  )
      IMPLICIT NONE

      ! exit if module is not setup
      IF ( .NOT. ModuleisSetup ) THEN
#if defined(LOG_FILE)
         __OPEN_LOG_FILE
         WRITE(__LOG_UNIT,*) " PeriodicBoundary: error disposing data, module is not yet initialized"
         __CLOSE_LOG_FILE
#endif
         RETURN
      END IF

      ModuleisSetup = .FALSE.
      
   END SUBROUTINE PBC_Dispose

   
!============================================================================================

END MODULE PeriodicBoundary
 
