!***************************************************************************************
!*                              MODULE FFTWrapper
!***************************************************************************************
!
!>  \brief     Wrapper of fast fourier transform
!>  \details   This module contains subroutines 
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             28 November 2013
!>
!***************************************************************************************
!
!>  \par Updates
!>  \arg 
!
!>  \todo         
!>                 
!***************************************************************************************

#if !defined(WITH_FFTW3) 
#warning "FFTWrapper: FFTW3 not available: using matrix-vector product DFT instead."
#define WITH_DFT_MATRIX
#else
#warning "FFTWrapper: compiling with FFTW3.3 libraries."
#undef WITH_DFT_MATRIX
#endif

#if defined(WITH_DFT_MATRIX) 
!#error "FFTWrapper: DFT has been implemented but not tested yet "
#endif

MODULE FFTWrapper
#include "preprocessoptions.cpp"
#if defined(WITH_FFTW3)
   USE, INTRINSIC :: ISO_C_BINDING
#endif
#if defined(WITH_DFT_MATRIX)
   USE MyLinearAlgebra
#endif
   IMPLICIT NONE

#if defined(WITH_FFTW3)
   INCLUDE 'fftw3.f03'
#endif

   PRIVATE
   PUBLIC ::  FFTComplexType, FFTHalfComplexType
   PUBLIC ::  SetupFFT, ExecuteFFT, DisposeFFT

   INTEGER, PARAMETER, PUBLIC :: DIRECT_FFT  = 0
   INTEGER, PARAMETER, PUBLIC :: INVERSE_FFT = 1

   TYPE FFTComplexType
#if defined(WITH_FFTW3)
      INTEGER(C_SIZE_T) :: NData
      COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:), POINTER :: In, Out
      TYPE(C_PTR) :: WorkIn, WorkOut
      TYPE(C_PTR) :: DirectPlan, InversePlan
#endif
#if defined(WITH_DFT_MATRIX)
      INTEGER :: NData
      COMPLEX, DIMENSION(:,:), POINTER :: DirectDFT, InverseDFT
#endif
      LOGICAL     :: isSetup = .FALSE.
   END TYPE FFTComplexType

   TYPE FFTHalfComplexType
#if defined(WITH_FFTW3)
      INTEGER(C_SIZE_T) :: NData
      REAL(C_DOUBLE), DIMENSION(:), POINTER :: In, Out
      TYPE(C_PTR) :: WorkIn, WorkOut
      TYPE(C_PTR) :: DirectPlan, InversePlan
      REAL, DIMENSION(:), POINTER :: NormalizationDir, NormalizationInv
#endif
#if defined(WITH_DFT_MATRIX)
      INTEGER :: NData
      REAL, DIMENSION(:,:), POINTER :: DirectDFT, InverseDFT
#endif
      LOGICAL     :: isSetup = .FALSE.
   END TYPE FFTHalfComplexType

   INTERFACE SetupFFT
      MODULE PROCEDURE SetupFFTComplex, SetupFFTHalfComplex
   END INTERFACE 

   INTERFACE ExecuteFFT
      MODULE PROCEDURE ExecuteFFTComplex, ExecuteFFTHalfComplex
   END INTERFACE 

   INTERFACE DisposeFFT
      MODULE PROCEDURE DisposeFFTComplex, DisposeFFTHalfComplex
   END INTERFACE 

!*******************************************************************************************************************************
                                                         CONTAINS
!*******************************************************************************************************************************

   SUBROUTINE SetupFFTComplex( PlanData, N )
      IMPLICIT NONE
      TYPE( FFTComplexType ) :: PlanData
      INTEGER, INTENT(IN)    :: N
      INTEGER :: i, j
      COMPLEX :: i2pioverN

      CALL ERROR( PlanData%isSetup, " SetupFFTComplex: PlanData is already set " )

      PlanData%NData = N

#if defined(WITH_FFTW3)
      PlanData%WorkIn = FFTW_ALLOC_COMPLEX( PlanData%NData )
      PlanData%WorkOut = FFTW_ALLOC_COMPLEX( PlanData%NData )

      CALL C_F_POINTER( PlanData%WorkIn, PlanData%In, [ PlanData%NData ] )
      CALL C_F_POINTER( PlanData%WorkOut, PlanData%Out, [ PlanData%NData ] )

      PlanData%DirectPlan  = FFTW_PLAN_DFT_1D( N, PlanData%In, PlanData%Out, FFTW_FORWARD, FFTW_MEASURE )
      PlanData%InversePlan = FFTW_PLAN_DFT_1D( N, PlanData%In, PlanData%Out, FFTW_BACKWARD, FFTW_MEASURE )
#endif
#if defined(WITH_DFT_MATRIX)
      ALLOCATE( PlanData%DirectDFT( PlanData%NData, PlanData%NData ), PlanData%InverseDFT( PlanData%NData, PlanData%NData ) )
      i2pioverN = 2.0*MyConsts_PI*MyConsts_I / real(PlanData%NData)
      DO i = 0, PlanData%NData-1
         DO j = 0, PlanData%NData-1
            PlanData%DirectDFT(i+1,j+1) = EXP( - i2pioverN * real(i*j) )
            PlanData%InverseDFT(i+1,j+1) = EXP( i2pioverN * real(i*j) ) / real(PlanData%NData)
         END DO
      END DO
#endif

      PlanData%isSetup = .TRUE. 

   END SUBROUTINE SetupFFTComplex

   SUBROUTINE SetupFFTHalfComplex( PlanData, N )
      IMPLICIT NONE
      TYPE( FFTHalfComplexType ) :: PlanData
      INTEGER, INTENT(IN)        :: N
      INTEGER :: i, j
      REAL :: pi2overN

      CALL ERROR( PlanData%isSetup, " SetupFFTHalfComplex: PlanData is already set " )

      PlanData%NData = N

#if defined(WITH_FFTW3)
      PlanData%WorkIn = FFTW_ALLOC_REAL( PlanData%NData )
      PlanData%WorkOut = FFTW_ALLOC_REAL( PlanData%NData )

      CALL C_F_POINTER( PlanData%WorkIn, PlanData%In, [ PlanData%NData ] )
      CALL C_F_POINTER( PlanData%WorkOut, PlanData%Out, [ PlanData%NData ] )

      PlanData%DirectPlan  = FFTW_PLAN_R2R_1D ( N, PlanData%In, PlanData%Out, FFTW_R2HC, FFTW_MEASURE )
      PlanData%InversePlan = FFTW_PLAN_R2R_1D ( N, PlanData%In, PlanData%Out, FFTW_HC2R, FFTW_MEASURE )

      ALLOCATE( PlanData%NormalizationDir( N ), PlanData%NormalizationInv( N ) )
      DO i = 0, N-1
         IF ( i == 0 .OR. real(i) == real(N)/2.0 ) THEN
            PlanData%NormalizationDir(i+1) = SQRT( 1.0 / real(N) )
            PlanData%NormalizationInv(i+1) = SQRT( 1.0 / real(N) )
         ELSE 
            PlanData%NormalizationDir(i+1) = SQRT( 2.0 / real(N) ) 
            PlanData%NormalizationInv(i+1) = SQRT( 1.0 / ( real(N) * 2.0 ) )
         END IF
      END DO
#endif
#if defined(WITH_DFT_MATRIX)
      ALLOCATE( PlanData%DirectDFT( PlanData%NData, PlanData%NData ), PlanData%InverseDFT( PlanData%NData, PlanData%NData ) )
      pi2overN = 2.0*MyConsts_PI / real(PlanData%NData)
      DO i = 0, PlanData%NData-1
         DO j = 0, PlanData%NData-1
            IF ( i == 0 ) THEN
               PlanData%DirectDFT(i+1,j+1)  = SQRT( 1.0 /  real(PlanData%NData) )
               PlanData%InverseDFT(j+1,i+1) = SQRT( 1.0 /  real(PlanData%NData) )
            ELSE IF ( i > 0 .AND. real(i) < real(N)/2.0 ) THEN
               PlanData%DirectDFT(i+1,j+1)  = COS( pi2overN * real(i*j) ) * SQRT( 2.0 /  real(PlanData%NData) )
               PlanData%InverseDFT(j+1,i+1) = COS( pi2overN * real(i*j) ) * SQRT( 2.0 /  real(PlanData%NData) )
            ELSE IF ( real(i) == real(N)/2.0 ) THEN
               PlanData%DirectDFT(i+1,j+1)  = (-1.0)**j * SQRT( 1.0 /  real(PlanData%NData) )
               PlanData%InverseDFT(j+1,i+1) = (-1.0)**j * SQRT( 1.0 /  real(PlanData%NData) )
            ELSE IF ( i > N/2 .AND. i < N ) THEN
               PlanData%DirectDFT(i+1,j+1)  = SIN( pi2overN * real(i*j) ) * SQRT( 2.0 /  real(PlanData%NData) )
               PlanData%InverseDFT(j+1,i+1) = SIN( pi2overN * real(i*j) ) * SQRT( 2.0 /  real(PlanData%NData) )
            END IF
         END DO
      END DO
#endif

      PlanData%isSetup = .TRUE. 

   END SUBROUTINE SetupFFTHalfComplex

!*******************************************************************************************************************************

   SUBROUTINE ExecuteFFTComplex( PlanData, Vector, Direction )
      IMPLICIT NONE
      TYPE( FFTComplexType ), INTENT(INOUT) :: PlanData
      COMPLEX, DIMENSION(:), INTENT(INOUT)  :: Vector
      INTEGER, INTENT(IN)                   :: Direction

      CALL ERROR( size(Vector) /= PlanData%NData, " ExecuteFFTComplex: wrong array size " )
      CALL ERROR( .NOT. PlanData%isSetup, " ExecuteFFTComplex: PlanData is not set " )

#if defined(WITH_FFTW3)
      PlanData%In = Vector
      IF ( Direction == DIRECT_FFT ) THEN
         CALL FFTW_EXECUTE_DFT( PlanData%DirectPlan,  PlanData%In, PlanData%Out)
         Vector = PlanData%Out
      ELSE IF ( Direction == INVERSE_FFT ) THEN
         CALL FFTW_EXECUTE_DFT( PlanData%InversePlan, PlanData%In, PlanData%Out)
         Vector = PlanData%Out / real( PlanData%NData )
      END IF
#endif
#if defined(WITH_DFT_MATRIX)
      IF ( Direction == DIRECT_FFT ) THEN
         Vector = TheOneWithMatrixVectorProduct( PlanData%DirectDFT, Vector )
      ELSE IF ( Direction == INVERSE_FFT ) THEN
         Vector = TheOneWithMatrixVectorProduct( PlanData%InverseDFT, Vector )
      END IF
#endif

   END SUBROUTINE ExecuteFFTComplex

   SUBROUTINE ExecuteFFTHalfComplex( PlanData, Vector, Direction )
      IMPLICIT NONE
      TYPE( FFTHalfComplexType ), INTENT(INOUT) :: PlanData
      REAL, DIMENSION(:), INTENT(INOUT)         :: Vector
      INTEGER, INTENT(IN)                       :: Direction

      CALL ERROR( size(Vector) /= PlanData%NData, " ExecuteFFTHalfComplex: wrong array size " )
      CALL ERROR( .NOT. PlanData%isSetup, " ExecuteFFTHalfComplex: PlanData is not set " )

#if defined(WITH_FFTW3)
      IF ( Direction == DIRECT_FFT ) THEN
         PlanData%In = Vector
         CALL FFTW_EXECUTE_R2R( PlanData%DirectPlan,  PlanData%In, PlanData%Out)
         Vector(:) = PlanData%Out(:) * PlanData%NormalizationDir(:)
      ELSE IF ( Direction == INVERSE_FFT ) THEN
         PlanData%In = Vector * PlanData%NormalizationInv(:)
         CALL FFTW_EXECUTE_R2R( PlanData%InversePlan, PlanData%In, PlanData%Out)
         Vector(:) = PlanData%Out(:)
      END IF
#endif
#if defined(WITH_DFT_MATRIX)
      IF ( Direction == DIRECT_FFT ) THEN
         Vector = TheOneWithMatrixVectorProduct( PlanData%DirectDFT, Vector )
      ELSE IF ( Direction == INVERSE_FFT ) THEN
         Vector = TheOneWithMatrixVectorProduct( PlanData%InverseDFT, Vector )
      END IF
#endif

   END SUBROUTINE ExecuteFFTHalfComplex

!*******************************************************************************************************************************

   SUBROUTINE DisposeFFTComplex( PlanData )
      IMPLICIT NONE
      TYPE( FFTComplexType ) :: PlanData

#if defined(WITH_FFTW3)
      CALL FFTW_DESTROY_PLAN( PlanData%DirectPlan )
      CALL FFTW_DESTROY_PLAN( PlanData%InversePlan )
      CALL FFTW_FREE( PlanData%WorkIn )
      CALL FFTW_FREE( PlanData%WorkOut )
#endif
#if defined(WITH_DFT_MATRIX)
      DEALLOCATE( PlanData%DirectDFT, PlanData%InverseDFT )
#endif
      PlanData%isSetup = .FALSE. 

   END SUBROUTINE DisposeFFTComplex

   SUBROUTINE DisposeFFTHalfComplex( PlanData )
      IMPLICIT NONE
      TYPE( FFTHalfComplexType ) :: PlanData

#if defined(WITH_FFTW3)
      CALL FFTW_DESTROY_PLAN( PlanData%DirectPlan )
      CALL FFTW_DESTROY_PLAN( PlanData%InversePlan )
      CALL FFTW_FREE( PlanData%WorkIn )
      CALL FFTW_FREE( PlanData%WorkOut )
#endif
#if defined(WITH_DFT_MATRIX)
      DEALLOCATE( PlanData%DirectDFT, PlanData%InverseDFT )
#endif
      PlanData%isSetup = .FALSE. 

   END SUBROUTINE DisposeFFTHalfComplex

END MODULE FFTWrapper
