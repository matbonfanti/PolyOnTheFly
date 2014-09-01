!***************************************************************************************
!*                              MODULE RandomNumberGenerator
!***************************************************************************************
!
!>  \brief     Random number generators
!>  \details   This class defines ...
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             15 January 2013
!>
!***************************************************************************************
!
!>  \par Updates
!>  \arg 7 December 2013 : the subroutines UniformRandomNr and  GaussianRandomNr 
!>       are now thread safe, and different initial seeds can be set by using 
!>       different RNGInternalState variables and initializing them with different
!>       calls of SetSeed. Thus the module is ready to be used for parallel applications
!
!>  \todo          ....
!>                 
!***************************************************************************************
!
!>   \remark     The function for uniform random number generation is taken from  \n 
!>               (*) Press, Teukolsky, Vetterling, Flannery \n
!>                   Numerical Recipes, the art of scientific computing, \n
!>                   Pag. 1142 Vol. 2 (FORTRAN 90) \n
!>                   Cambridge University Press \n
!>                   http://apps.nrbook.com/fortran/index.html
!
!***************************************************************************************
MODULE RandomNumberGenerator
#include "preprocessoptions.cpp"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: RNGInternalState
   PUBLIC :: SetSeed, UniformRandomNr, GaussianRandomNr
   PUBLIC :: TestGaussianDistribution, TestGaussianDistribution2, TestCorrelations

   ! Integer type for the random number generation
   INTEGER, PARAMETER, PUBLIC :: K4B=selected_int_kind(9)

   TYPE RNGInternalState
      ! These variables define the internal state of the random number generator
      REAL         :: am
      INTEGER(K4B) :: ix, iy, k
      INTEGER(K4B) :: seed
      ! Since the gaussian random nr are generated in couples, the following variable
      ! store the non used gaussian number for later calls of the subroutine
      REAL         :: TempGaussian
      LOGICAL      :: GaussianAvail = .FALSE.
   END TYPE

 CONTAINS   

!****************************************************************************

! Set seed for random number generation and initialize the internal variables
! of the RNGInternalState data type
   
   SUBROUTINE SetSeed( IntState, Seed )
      IMPLICIT NONE
      TYPE( RNGInternalState), INTENT(INOUT) :: IntState
      INTEGER, INTENT(IN)                    :: Seed

      CALL ERROR( Seed > 0, " RandomNumberGenerator.SetSeed: negative seed required " )

      ! Initialize seed
      IntState%seed = INT( Seed, K4B )
      
      ! Initialize internal state of the RNG 
      IntState%am   = 0.0
      IntState%ix   = -1
      IntState%iy   = -1
      IntState%k    = 0
      
      ! Initialize variables to store gaussian random number
      IntState%TempGaussian = 0.0
      IntState%GaussianAvail = .FALSE.
      
   END SUBROUTINE SetSeed
 
  
!****************************************************************************

! The subroutine computes a preudo-random real number in the 0-1 interval
! This subroutine has been adapted from Numerical Recipes for Fortran 90
! to take into account the possibility of using the RNG for parallel 
! applications (for details on the algorithm, see pag 1142 NR for FORTRAN)
   
   REAL FUNCTION UniformRandomNr( IntState ) RESULT(Ran)
      IMPLICIT NONE
      TYPE( RNGInternalState ), INTENT(INOUT) :: IntState 
      INTEGER(K4B), PARAMETER :: IA=16807, IM=2147483647, IQ=127773, IR=2836

      if (IntState%seed <= 0 .or. IntState%iy < 0) then 
         IntState%am = nearest(1.0,-1.0)/IM
         IntState%iy = ior(ieor(888889999,abs(IntState%seed)),1)
         IntState%ix = ieor(777755555,abs(IntState%seed))
         IntState%seed = abs(IntState%seed)+1
      end if
      
      IntState%ix = ieor(IntState%ix, ishft(IntState%ix,13))
      IntState%ix = ieor(IntState%ix, ishft(IntState%ix,-17))
      IntState%ix = ieor(IntState%ix, ishft(IntState%ix,5))
      IntState%k=IntState%iy/IQ
      IntState%iy=IA*(IntState%iy-IntState%k*IQ)-IR*IntState%k
      if (IntState%iy < 0) IntState%iy=IntState%iy+IM
      ran=IntState%am*ior(iand(IM,ieor(IntState%ix,IntState%iy)),1)

   END FUNCTION UniformRandomNr
   
!****************************************************************************

! Generate random numbers distribuited according to a gaussian distribution
! with standard deviation sigma ( Box-Muller algorithm ) 

   REAL FUNCTION GaussianRandomNr( IntState ) RESULT( RandNr )
      IMPLICIT NONE
      TYPE( RNGInternalState ), INTENT(INOUT) :: IntState 
      REAL :: Theta, R, Random

      IF ( .NOT. IntState%GaussianAvail ) THEN

            ! Generate random number for 2D gaussian function
            Theta = 2. * MyConsts_PI * UniformRandomNr(IntState)
            R     = SQRT( -2. * LOG( UniformRandomNr(IntState) ) )

            ! TRansform in cartesian coordinate, Return one number and store the other
            RandNr = R * sin( Theta )
            IntState%TempGaussian = R * cos( Theta )
            IntState%GaussianAvail = .TRUE.

      ELSE IF ( IntState%GaussianAvail ) THEN

            RandNr = IntState%TempGaussian
            IntState%GaussianAvail = .FALSE.

      END IF
  
   END FUNCTION GaussianRandomNr

   
   SUBROUTINE TestCorrelations( MaxN, NrThreads )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MaxN, NrThreads
      TYPE( RNGInternalState ), DIMENSION(NrThreads) :: RandomNrGen
      INTEGER :: PrintStep
      REAL :: Xcorr, Xmean, Random, Previous
      INTEGER :: iN, iRNG, NTot
      
      ! Set the output steps
      PrintStep = 1000
      
      ! Print intestation of the table
      WRITE(125,*) "# Nr of random numbers, correlation between consecutive nr: <x_i x_i+1> - <x_i>**2 "

      DO iRNG = 1, NrThreads
         CALL SetSeed( RandomNrGen(iRNG), -iRNG )
      END DO      
      
      ! Cycle over nr of number to generate
      NTot = 0
      XCorr = 0.0
      Xmean = 0.0
      DO iN = 1, MaxN
         DO iRNG = 1, NrThreads

            ! Generate uniform rnd number 
            Previous = Random
            Random = UniformRandomNr( RandomNrGen(iRNG) ) 
            NTot = NTot + 1

            ! increment sum and correlation product
            IF ( NTot > 1 ) XCorr = XCorr + Previous * Random
            Xmean = Xmean + Random
            
            ! IF it's a printing step, print average and st dev
            IF ( mod( NTot, PrintStep ) == 0 ) THEN
               WRITE(125,*) NTot, XCorr/(NTot-1) - (Xmean/NTot)**2
            ENDIF
         END DO
      END DO

   END SUBROUTINE TestCorrelations


   SUBROUTINE TestGaussianDistribution( MaxN, NrThreads )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MaxN, NrThreads
      TYPE( RNGInternalState ), DIMENSION(NrThreads) :: RandomNrGen
      REAL :: AverageEst, SigmaEst, Random
      INTEGER :: NrOfThreads, CurrentThread, i, j, N
      INTEGER :: NrOfStepEachPring, PrintStep

      ! Print intestation of the table
      WRITE(123,*) "# Nr of random numbers, error in Average value, error in Standard deviation "
      
      ! Set the output steps
      NrOfStepEachPring = 10000
      PrintStep = MaxN*NrThreads / NrOfStepEachPring

      DO i = 1, NrThreads
         CALL SetSeed( RandomNrGen(i), -i )
      END DO

      ! Initialize variables
      AverageEst = 0.0
      SigmaEst = 0.0
      N = 0
      
      DO i = 1, MaxN
         DO j = 1, NrThreads
            Random = GaussianRandomNr( RandomNrGen(j) )
            AverageEst = AverageEst + Random
            SigmaEst   = SigmaEst   + Random**2
            N = N + 1 
            IF  ( MOD(N,NrOfStepEachPring) == 0 ) THEN
               WRITE(123,*) N, AverageEst/N, sqrt(SigmaEst/N-(AverageEst/N)**2) 
            END IF
        END DO
      END DO
      
   END SUBROUTINE TestGaussianDistribution
 
   SUBROUTINE TestGaussianDistribution2( MaxN )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MaxN
      INTEGER, DIMENSION(61) :: IntCounter
      INTEGER :: i, NrOfThreads, CurrentThread, RanIndex
      TYPE( RNGInternalState ) :: RandomNrGen
      REAL :: Random

      WRITE(122,*) "# Interval of x, Nr of random numbers "

      !$OMP PARALLEL PRIVATE(CurrentThread, Random, RandomNrGen, RanIndex )
      !$OMP MASTER
      NrOfThreads = __OMP_TotalNrOfThreads
      IntCounter(:) = 0.0
      !$OMP END MASTER
   
      CurrentThread = __OMP_CurrentThreadNum
      CALL SetSeed( RandomNrGen, -1-CurrentThread+1 )
   
      !$OMP DO REDUCTION(+:IntCounter)
      DO i = 1, MaxN
         Random = GaussianRandomNr( RandomNrGen ) 
         RanIndex = CEILING( 30.0 + 10.*Random )
         IF ( RanIndex >= 1 .AND. RanIndex <= 61 ) THEN
            IntCounter(RanIndex) =  IntCounter(RanIndex) + 1
         END IF
      END DO
      !$OMP END DO 
      !$OMP END PARALLEL

      DO i = 1, 61
         WRITE(122,*)  (real(i-31)/10.), IntCounter(i)
      END DO

   END SUBROUTINE TestGaussianDistribution2


   !****************************************************************************

 
END MODULE RandomNumberGenerator
