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
!>  \arg 
!
!>  \todo          ____________________________
!>                 
!***************************************************************************************

MODULE PeriodicBoundary
#include "preprocessoptions.cpp"
   USE FFTWrapper

   PRIVATE

   PUBLIC :: PBC_Setup, PBC_Dispose, PBC_BringToFirstCell
   PUBLIC :: FractionalToCartesian, CartesianToFractional

   !> Vectors of the unit cell
   REAL, DIMENSION(3,3) :: UnitCellVector

   !> Matrix transformation for coordinates
   REAL, DIMENSION(3,3) :: Fractional2Cartesian, Cartesian2Fractional

   !> Setup variable for the module
   LOGICAL, SAVE :: ModuleisSetup = .FALSE.

   
!============================================================================================
                                       CONTAINS
!============================================================================================

   SUBROUTINE PBC_Setup( VectorA, VectorB, VectorC  )
      IMPLICIT NONE
      REAL, DIMENSION(3), INTENT(IN) :: VectorA, VectorB, VectorC

      ! warning if module is setup
      IF ( ModuleisSetup ) THEN
#if defined(LOG_FILE)
         __OPEN_LOG_FILE
         WRITE(__LOG_UNIT,*) " PeriodicBoundary: module is already setup. Overwriting data "
         __CLOSE_LOG_FILE
#endif
      END IF

      ! Store Unit Cell Vectors
      UnitCellVector(:,1) = VectorA
      UnitCellVector(:,2) = VectorB
      UnitCellVector(:,3) = VectorC

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
!          PBC_BringToFirstCell
!*******************************************************************************
!> Translate coordinates to the symmetric point in the first unit cell. 
!> If the external program use cartesian coordinates, the corresponding fractional
!> coordinate are computed and then the fractional coordinates are taken and trasformed
!> back to cartesian coordinates.
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
!          FractionalToCartesian
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
!          InPlaceFractionalToCartesian
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
!          CartesianToFractional
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
!          InPlaceCartesianToFractional
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
 
