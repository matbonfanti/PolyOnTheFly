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

   !> Coordinates which are used in the external program, fractional or cartesian
   LOGICAL, SAVE :: PassFractionalCoords = .TRUE.

   !> Matrix transformation for coordinates
   REAL, DIMENSION(3,3) :: Fractional2Cartesian, Cartesian2Fractional

   !> Setup variable for the module
   LOGICAL, SAVE :: ModuleisSetup = .FALSE.

   
!============================================================================================
                                       CONTAINS
!============================================================================================

   SUBROUTINE PBC_Setup( VectorA, VectorB, VectorC, FractionalCoords  )
      IMPLICIT NONE
      REAL, DIMENSION(3), INTENT(IN) :: VectorA, VectorB, VectorC
      LOGICAL, INTENT(IN), OPTIONAL  :: FractionalCoords

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

      ! Store external coordinates definition
      IF ( PRESENT( FractionalCoords) )   PassFractionalCoords = FractionalCoords

      ! Module is now ready
      ModuleisSetup = .TRUE.
      
#if defined(LOG_FILE)
      __OPEN_LOG_FILE
      WRITE(__LOG_UNIT,*) " PeriodicBoundary: module is now initialized"
      WRITE(__LOG_UNIT,*) " PeriodicBoundary: vector A - ", UnitCellVector(:,1)
      WRITE(__LOG_UNIT,*) " PeriodicBoundary: vector B - ", UnitCellVector(:,2)
      WRITE(__LOG_UNIT,*) " PeriodicBoundary: vector C - ", UnitCellVector(:,3)
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
      INTEGER :: i, j, NBeads

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

      ! Cycle over the atoms
      DO i = 1, size(X)/3

         ! trasform cartesian coordinates to fractional coordinates
         IF ( PassFractionalCoords ) THEN
            Vector = X( (i-1)*3+1 : i*3 )
         ELSE
            Vector = TheOneWithMatrixVectorProduct( Cartesian2Fractional, X( (i-1)*3+1 : i*3 ) )
         END IF

         ! move in first unit cell
         Vector(1) = Vector(1) - FLOOR( Vector(1) )
         Vector(2) = Vector(2) - FLOOR( Vector(2) )
         Vector(3) = Vector(3) - FLOOR( Vector(3) )

         ! trasform back coordinate to cartesian 
         IF ( PassFractionalCoords ) THEN
            X( (i-1)*3+1 : i*3 ) = Vector
         ELSE
             X( (i-1)*3+1 : i*3 ) = TheOneWithMatrixVectorProduct( Fractional2Cartesian, Vector )
         END IF

      END DO

   END SUBROUTINE PBC_BringToFirstCell


!============================================================================================

!*******************************************************************************
!          FractionalToCartesian
!*******************************************************************************
!> Convert fractional coordinates to cartesian coordinates for X and Y. \n
!> The coordinate along Z is left unchanged. Works only if the slab geometry has
!> been already setup.
!>
!> @param      Array of dimension 3 with the fractional coordinates.
!> @returns    Array of dimension 3 with the cartesian coordinates
!*******************************************************************************
   FUNCTION FractionalToCartesian( Fractional ) RESULT( Cartesian )
      IMPLICIT NONE
      REAL, DIMENSION(3) :: Cartesian
      REAL, DIMENSION(3), INTENT(IN) :: Fractional

      ! Return if module not have been setup yet ( NON PERIODIC SYSTEM )
      IF ( .NOT. ModuleisSetup )  RETURN 

      ! trasform coordinate for X and Y
      Cartesian(1:3) = TheOneWithMatrixVectorProduct( Fractional2Cartesian, Fractional(1:3) )
   END FUNCTION FractionalToCartesian


!*******************************************************************************
!          CartesianToFractional
!*******************************************************************************
!> Convert cartesian  coordinates to fractional coordinates for X and Y. \n
!> The coordinate along Z is left unchanged. Works only if the slab geometry has
!> been already setup.
!>
!> @param      Array of dimension 3 with the cartesian coordinates
!> @returns    Array of dimension 3 with the fractional coordinates.
!*******************************************************************************
   FUNCTION CartesianToFractional( Cartesian ) RESULT( Fractional )
      IMPLICIT NONE
      REAL, DIMENSION(3) :: Fractional
      REAL, DIMENSION(3), INTENT(IN) :: Cartesian

      ! Return if module not have been setup yet ( NON PERIODIC SYSTEM )
      IF ( .NOT. ModuleisSetup )  RETURN 

      ! trasform coordinate for X and Y
      Fractional(1:3) = TheOneWithMatrixVectorProduct( Cartesian2Fractional, Cartesian(1:3) )
   END FUNCTION CartesianToFractional


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
 
