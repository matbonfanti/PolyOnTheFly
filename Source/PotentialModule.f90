!***************************************************************************************
!*                           MODULE PotentialModule
!***************************************************************************************
!
!>  \brief     Potential Energy Module
!>  \details   This module setup and compute the potential energy for the
!>             moleucular dynamics propagation.
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
MODULE PotentialModule
#include "preprocessoptions.cpp"
   USE PeriodicBoundary

   PRIVATE

   PUBLIC :: SetupPairPotential, DisposePotential, GetPotential


   !> Setup variable for the potential
   LOGICAL, SAVE :: PotentialModuleIsSetup = .FALSE.

   !> Use model pair potential instead of on-the-fly computation of the forces
   LOGICAL, SAVE :: ModelPotential = .FALSE.

   !> How many neighbour atoms are included in the pair potential summation
   INTEGER, SAVE :: NearPeriodicImages = 0
   !> Lattice translations of the near neighbour cells
   REAL, DIMENSION(:,:), ALLOCATABLE :: NearTranslations
   !> Cutoff distance of the pair potential
   REAL, SAVE   :: CutOff = 100.0


   !> Number of atoms of the system
   INTEGER, SAVE :: AtomNo
   
   !> \name LENNARD-JONES POTENTIAL PARAMETERS
   !> Parameters of the Lennard-Jones pair potential
   !> @{
   REAL, SAVE :: LJ_WellDepth = 1.0       !< Potential well depth (energy)
   REAL, SAVE :: LJ_EquilDist = 1.0       !< Potential equilibrium distance (lenght)
   !> @}
   
   
!============================================================================================
                                       CONTAINS
!============================================================================================

   SUBROUTINE SetupPairPotential( N, WellDepth, EquilDist, PBC  )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      REAL, INTENT(IN)    :: WellDepth, EquilDist
      LOGICAL, INTENT(IN) :: PBC
      INTEGER :: i, j, k
      REAL, DIMENSION(3,125) :: TmpNearTranslations
      REAL, DIMENSION(3) :: Vector
      REAL :: Distance

      ! exit if module is setup
      IF ( PotentialModuleIsSetup ) RETURN

      ! use model potential
      ModelPotential = .TRUE.
      
      ! Set the number of atoms of the system
      AtomNo = N
      
      ! Store the pair potential parameters
      LJ_WellDepth = WellDepth
      LJ_EquilDist = EquilDist

      ! Set the cutoff distance of the pair potential
      CutOff = EquilDist * 6.0

      IF ( PBC ) THEN
         ! Define how many neighbour cells are checked in the pair potential summation
         DO i = -5, +5
            DO j = -5, +5
               DO k = -5, +5
                  Vector = FractionalToCartesian( (/ REAL(i), REAL(j), REAL(k) /) )
                  Distance = SQRT( TheOneWithVectorDotVector( Vector, Vector ) )
                  IF (Distance < CutOff) THEN
                     NearPeriodicImages = NearPeriodicImages + 1
                     TmpNearTranslations(:,NearPeriodicImages) =  Vector
                  END IF
               END DO
            END DO
         END DO
         ALLOCATE( NearTranslations(3,NearPeriodicImages) )
         NearTranslations(:,:) =  TmpNearTranslations(:,1:NearPeriodicImages )
      ELSE 
         NearPeriodicImages = 1
         ALLOCATE( NearTranslations(3,NearPeriodicImages) )
         NearTranslations(:,1) =  (/ 0, 0, 0 /) 
      ENDIF

      ! Module is now ready
      PotentialModuleIsSetup = .TRUE.
      
   END SUBROUTINE SetupPairPotential

!============================================================================================

   REAL FUNCTION GetPotential( X, Force )
      IMPLICIT NONE
      REAL, DIMENSION(:), TARGET, INTENT(IN)  :: X
      REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Force

      ! Error if module not have been setup yet
      CALL ERROR( .NOT. PotentialModuleIsSetup, " PotentialModule.GetPotential : Module not Setup" )

      ! Check dimension of function arguments
      CALL ERROR( SIZE(X) /= 3*AtomNo,    " PotentialModule.GetPotential : wrong dimension of coordinate array " )
      CALL ERROR( SIZE(X) /= SIZE(Force), " PotentialModule.GetPotential : wrong dimension of forces array " )
      
      ! If model potential, use lennard-jones pair potential
      IF ( ModelPotential ) THEN
         CALL PairPotential( X(:), GetPotential, Force(:) )

      ELSE IF ( .NOT.  ModelPotential ) THEN
         CALL AbortWithError( " GetPotential: true potential is not yet available " )
      END IF
      
   END FUNCTION GetPotential

!============================================================================================

   SUBROUTINE DisposePotential(  )
      IMPLICIT NONE

      ! exit if module is not setup
      IF ( .NOT. PotentialModuleIsSetup ) RETURN

      PotentialModuleIsSetup = .FALSE.
      
   END SUBROUTINE DisposePotential

!============================================================================================

   SUBROUTINE PairPotential( Positions, V, Forces )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN)  :: Positions(:)
      REAL, INTENT(OUT)               :: V
      REAL, DIMENSION(:), INTENT(OUT) :: Forces(:)

      INTEGER :: iAtom, jAtom, iTrasl
      REAL, DIMENSION(3) :: iCoord, FirstDist, TranslatedDist
      REAL    :: Distance, LJDerivative

      ! Initialize output variables
      V = 0.0
      Forces(:) = 0.0
      
      DO iAtom = 1, AtomNo
         ! Extract coordinates of the i-th atom
         iCoord = Positions( (iAtom-1)*3+1 : iAtom*3 )
         
         DO jAtom = iAtom+1, AtomNo
            ! Extract distance vector between the i-th and j-th atoms
            FirstDist = iCoord( : ) - Positions( (jAtom-1)*3+1 : jAtom*3 )

            DO iTrasl = 1, NearPeriodicImages

               ! Compute periodic image of the distance
               TranslatedDist = FirstDist + NearTranslations(:,iTrasl) 
               Distance = SQRT( TheOneWithVectorDotVector( TranslatedDist , TranslatedDist ) )

               IF ( Distance > CutOff ) CYCLE

               ! Compute the pair potential
               V = V + LennardJones( Distance, LJDerivative )

               ! Update forces
               Forces( (iAtom-1)*3+1 ) = Forces( (iAtom-1)*3+1 ) + LJDerivative * ( TranslatedDist(1) ) / Distance
               Forces( (iAtom-1)*3+2 ) = Forces( (iAtom-1)*3+3 ) + LJDerivative * ( TranslatedDist(2) ) / Distance
               Forces(  iAtom*3      ) = Forces(  iAtom*3      ) + LJDerivative * ( TranslatedDist(3) ) / Distance
               Forces( (jAtom-1)*3+1 ) = Forces( (jAtom-1)*3+1 ) - LJDerivative * ( TranslatedDist(1) ) / Distance
               Forces( (jAtom-1)*3+2 ) = Forces( (jAtom-1)*3+3 ) - LJDerivative * ( TranslatedDist(2) ) / Distance
               Forces(  jAtom*3      ) = Forces(  jAtom*3      ) - LJDerivative * ( TranslatedDist(3) ) / Distance

            END DO
         END DO
      END DO
      
   END SUBROUTINE PairPotential

!============================================================================================

!    REAL FUNCTION MorseV( Positions, Forces ) RESULT(V) 
!       IMPLICIT NONE
!       REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Positions
!       REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Forces 
! 
!       V = MorseDe * ( exp(-2.0*MorseAlpha*Positions(1)) - 2.0 * exp(-MorseAlpha*Positions(1)) )  
!       Forces(1) = 2.0 * MorseAlpha * MorseDe * (  exp(-2.0*MorseAlpha*Positions(1)) - exp(-MorseAlpha*Positions(1)) )  
! 
!    END FUNCTION MorseV
      
   REAL FUNCTION LennardJones( Distance, Force ) RESULT(V) 
      IMPLICIT NONE
      REAL, INTENT(IN)  :: Distance
      REAL, INTENT(OUT) :: Force 

      V = LJ_WellDepth *( (LJ_EquilDist/Distance)**12 - 2.0*(LJ_EquilDist/Distance)**6 )
      Force = - 12.0 * LJ_WellDepth / Distance * (LJ_EquilDist/Distance)**7 *( (LJ_EquilDist/Distance)**6 - 1.0 )
   END FUNCTION LennardJones

END MODULE PotentialModule


! 
!    REAL, PARAMETER :: SmallDelta = 1.E-04
! 
!    !> Max nr of iterations for potential optimization
!    INTEGER, PARAMETER :: MaxIter = 10000
!    !> Threshold for conjugate gradient convergence
!    REAL, PARAMETER :: GradEps = 1.0E-6
!    !> Parameter for finite difference computation
!    REAL, PARAMETER :: Delta = 1.0E-4
! 
!    CONTAINS
! 
! 
! ! ************************************************************************************
! 
! 
!       SUBROUTINE SetupPotential( MassHydro, MassCarb, Collinear )
!          IMPLICIT NONE
!          REAL, INTENT(IN)     :: MassHydro, MassCarb
!          LOGICAL, OPTIONAL    :: Collinear
!          REAL, DIMENSION(124) :: Positions
!          INTEGER :: iCoord
!          REAL    :: Value
!          REAL, DIMENSION(4,4)       :: HessianSystem
! 
!          ! exit if module is setup
!          IF ( PotentialModuleIsSetup ) RETURN
! 
!          ! setup force constant ( in eV per Ang^2 )
!          rkc = (36.0*gam2+6.0*delt)/(bndprm**2)
! 
!          ! Setup if potential is collinear or not
!          IF ( PRESENT( Collinear ) ) THEN
!             CollinearPES =  Collinear
!          ELSE
!             CollinearPES = .FALSE.
!          END IF
! 
! 
!          ! Define the reference 4D potential for the graphite in minimum E
! 
!          ! Set guess starting coordinate for the minimization of the slab potential
!          Positions(1:124) = 0.0
!          Positions(3) = HZEquilibrium   ! reasonable guess for H Z coordinate
!          Positions(4) = C1Puckering
! 
!          ! Minimize potential
!          MinimumEnergy =  MinimizePotential( Positions, (/ (.TRUE., iCoord=1,124)  /) )      
! 
!          ! Translate to bring C3,C4,C5 in the Z=0 plane
!          Value = (Positions(5)+Positions(6)+Positions(7))/3.0
!          DO iCoord= 3,124
!             Positions(iCoord) = Positions(iCoord) - Value
!          END DO
! 
!          ! Store the coordinate of the slab
!          MinSlab(:) = Positions(5:124)
! 
!          ! Store the carbon puckering and the H Z at equilibrium
!          C1Puckering = Positions(4)
!          HZEquilibrium = Positions(3)
! 
!          ! Set the normal modes of the 4D potential
! 
!          ! compute the hessian
!          HessianSystem = HessianOfTheSystem( Positions, MassHydro, MassCarb )
! 
!          ! Diagonalize the hessian
!          CALL TheOneWithDiagonalization( HessianSystem, NormalModes4D_Vecs, NormalModes4D_Freq )
! 
! 
! 
! 
! #if defined(VERBOSE_OUTPUT)
!          WRITE(*,502) SQRT(NormalModes4D_Freq(1)), "au", "au", NormalModes4D_Vecs(:,1), &
!                       SQRT(NormalModes4D_Freq(2)), "au", "au", NormalModes4D_Vecs(:,2), &
!                       SQRT(NormalModes4D_Freq(3)), "au", "au", NormalModes4D_Vecs(:,3), &
!                       SQRT(NormalModes4D_Freq(4)), "au", "au", NormalModes4D_Vecs(:,4)
! 
!          502 FORMAT (/, " 4D system potential normal modes             ",/, &
!                         " 1) Frequency:                   ",1F15.2,1X,A, /, &
!                         "    Mass-scaled coords of the normal mode / ",A," : ",4F12.6, /, &
!                         " 2) Frequency:                   ",1F15.2,1X,A, /, &
!                         "    Mass-scaled coords of the normal mode / ",A," : ",4F12.6, /, &
!                         " 3) Frequency:                   ",1F15.2,1X,A, /, &
!                         "    Mass-scaled coords of the normal mode / ",A," : ",4F12.6, /, &
!                         " 4) Frequency:                   ",1F15.2,1X,A, /, &
!                         "    Mass-scaled coords of the normal mode / ",A," : ",4F12.6, 2/)
! 
!          WRITE(*,*) " Potential has been setup"
!          WRITE(*,*) " "
! #endif
!       END SUBROUTINE
! 
! 
! ! ************************************************************************************
! 
!       ! Setup initial conditions for the H atom + C slab for 
!       ! a simulation of vibrational relaxation
!       ! The slab is fixed in the equilibrium position with no momentum ( classical 0K )
!       ! The initial position and momenta of C and H are randomly chosen among a set of 
!       ! conditions which are given as input
!       ! data are initialized in ATOMIC UNITS
!       SUBROUTINE ZeroKelvinSlabConditions( Positions, Velocities, CHInitialConditions, RandomNr )
!          IMPLICIT NONE
!          REAL, DIMENSION(:), INTENT(OUT) :: Positions, Velocities
!          REAL, DIMENSION(:,:), INTENT(IN) :: CHInitialConditions
!          TYPE(RNGInternalState), INTENT(INOUT) :: RandomNr
!          INTEGER :: NDoF, iBath, NRandom, NInit
!          REAL :: Value
! 
!          ! Check the number of non frozen degree of freedom
!          NDoF = size( Positions )
!          CALL ERROR( size(Velocities) /= NDoF, "PotentialModule.ZeroKelvinSlabConditions: array dimension mismatch" )
! 
!          ! Check if the nr of dimension is compatible with the slab maximum size
!          CALL ERROR( (NDoF > 124) .OR. (NDoF < 4), "PotentialModule.ZeroKelvinSlabConditions: wrong number of DoFs" )
! 
!          ! Check the nr of starting conditions given ( there should be 8 coordinates: 4 positions and 4 momenta )
!          NRandom = size( CHInitialConditions, 1 )
!          CALL ERROR( size( CHInitialConditions, 2 ) /= 8, "PotentialModule.ZeroKelvinSlabConditions: wrong number of coords " )
!       
!          ! Set the velocities to zero
!          Velocities(:) = 0.0
! 
!          ! Set the slab in the equilibrium geometry
!          Positions(5:NDoF) = MinSlab(1:NDoF-4)
! 
!          ! Choose a random initial set of coordinates
!          NInit = CEILING( UniformRandomNr(RandomNr)*real(NRandom)  )
! 
!          ! Accordingly set position and velocity
!          Positions(1:4) = CHInitialConditions( NInit, 1:4 )
!          Velocities(1:4) = CHInitialConditions( NInit, 5:8 )
! 
!       END SUBROUTINE ZeroKelvinSlabConditions
! 
! 
!       ! Setup initial conditions for the H atom + C slab
!       ! data are initialized in ATOMIC UNITS
!       SUBROUTINE ThermalEquilibriumConditions( Positions, Velocities, Temperature, MassHydro, MassCarb, RandomNr )
!          IMPLICIT NONE
! 
!          REAL, DIMENSION(:), INTENT(OUT) :: Positions, Velocities
!          REAL, INTENT(IN)  :: Temperature, MassCarb, MassHydro
!          TYPE(RNGInternalState), INTENT(INOUT) :: RandomNr
!          INTEGER           :: nCarbon, NDoF
!          REAL              :: SigmaCarbonVelocity, SigmaHydroVelocity
! 
!          ! All the atoms are initially at the equilibrium position for stable chemisorption 
!          ! Value for the puckering are taken from J. Phys. Chem. B, 2006, 110, 18811-18817
!          ! Equilibrium position of zH obtained instead from plot of the PES
!          ! Velocities are sampled according to a Maxwell-Boltzmann distribution at temperature T
! 
!          ! Check the number of non frozen degree of freedom
!          NDoF = size( Positions )
!          CALL ERROR( size(Velocities) /= NDoF, "PotentialModule.ThermalEquilibriumConditions: array dimension mismatch" )
! 
!          ! Check if the nr of dimension is compatible with the slab maximum size
!          CALL ERROR( (NDoF > 124) .OR. (NDoF < 4), "PotentialModule.ThermalEquilibriumConditions: wrong number of DoFs" )
!             
!          ! Equilibrium position of H atom
!          Positions(1) = 0.0000
!          Positions(2) = 0.0000
!          Positions(3) = 1.483 / MyConsts_Bohr2Ang
! 
!          ! Equilibrium position of C1 atom
!          Positions(4) = C1Puckering
! 
!          ! Equilibrium position of the other carbon atoms 
!          DO nCarbon = 5,NDoF
!             Positions(nCarbon)   = 0.0
!          END DO
! 
!          ! Compute st deviation of Maxwell-Boltzmann distribution ( for the VELOCITY, not momenta!)
!          SigmaCarbonVelocity = sqrt( Temperature / MassCarb )
!          SigmaHydroVelocity  = sqrt( Temperature / MassHydro )
! 
!          ! Random velocities according to Maxwell-Boltzmann
!          IF ( CollinearPES ) THEN
!                Velocities(1) = 0.0
!                Velocities(2) = 0.0
!          ELSE 
!                Velocities(1) = GaussianRandomNr( RandomNr ) * SigmaHydroVelocity
!                Velocities(2) = GaussianRandomNr( RandomNr ) * SigmaHydroVelocity
!          END IF
!          Velocities(3) = GaussianRandomNr( RandomNr ) * SigmaHydroVelocity
!          DO nCarbon = 4,NDoF
!             Velocities(nCarbon) = GaussianRandomNr( RandomNr ) * SigmaCarbonVelocity 
!          END DO
! !         Velocities(4) = 0.0      ! TO FIX EVEN C1 ATOM
! 
! !          DO nCarbon = 1, Size( EdgeCarbons ) 
! !             Velocities(EdgeCarbons(nCarbon)+3) = 0.0
! !          END DO
! 
!       END SUBROUTINE ThermalEquilibriumConditions
! 
! 
! ! ************************************************************************************
! 
!       ! Setup initial conditions for the scattering of H atom on a thermalized C slab
!       ! data are initialized in ATOMIC UNITS
!       SUBROUTINE ScatteringConditions( Positions, Velocities, ImpactParam, InitZ, IncEnergy, Temperature, MassHydro, MassCarb )
!          IMPLICIT NONE
! 
!          REAL, DIMENSION(:), INTENT(OUT) :: Positions, Velocities
!          REAL, INTENT(IN)  :: Temperature, MassCarb, MassHydro
!          REAL, INTENT(IN)  :: ImpactParam, InitZ, IncEnergy
!          INTEGER           :: nCarbon, NDoF
!          REAL              :: SigmaMomentum, SigmaPosition
!          REAL              :: gaus1, gaus2, gvar1, gvar2, gr1, gr2, gs1, gs2
! 
!          ! Error if module not have been setup yet
!          CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.CarbonForceConstant : Module not Setup" )
! 
!          ! Check the number of non frozen degree of freedom
!          NDoF = size( Positions )
!          CALL ERROR( size(Velocities) /= NDoF, "PotentialModule.ScatteringConditions: array dimension mismatch" )
! 
!          ! Check if the nr of dimension is compatible with the slab maximum size
!          CALL ERROR( (NDoF > 124) .OR. (NDoF < 4), "PotentialModule.ScatteringConditions: wrong number of DoFs" )
! 
!          ! Scattering position of H atom
!          Positions(1) = ImpactParam
!          Positions(2) = 0.00001
!          Positions(3) = InitZ
! 
!          ! Velocity of the H atom
!          Velocities(1) = 0.0
!          Velocities(2) = 0.0
!          Velocities(3) = - sqrt( 2.0* IncEnergy / MassHydro )
! 
!          ! standard deviation of the position distribution (force constant needs to be in AU)
!          SigmaPosition = sqrt( Temperature / (rkc*(MyConsts_Bohr2Ang)**2/MyConsts_Hartree2eV) )
!          ! standard deviation of the momentum distribution
!          SigmaMomentum = sqrt( MassCarb * Temperature )
! 
!          ! Cycle over carbon atoms in the slab
!          DO nCarbon = 4,NDoF
! 
!             ! Initialization
!             Positions(nCarbon)   = 0.0
!             Velocities(nCarbon)   = 0.0
! 
!             ! Generate gaussian random numbers for position and velocity
!             DO 
!                call random_number(gaus1)
!                call random_number(gaus2)
!                gvar1=2.0*gaus1-1
!                gvar2=2.0*gaus2-1
!                gr1=gvar1**2+gvar2**2
!                IF (gr1 < 1.0) EXIT
!             END DO
!             gr2=sqrt(-2.0*alog(gr1)/gr1)
!             gs1=gvar1*gr2
!             gs2=gvar2*gr2
! 
!             ! Set position and velocity of the carbon atom
!             Positions(nCarbon)  = SigmaPosition*gs1
!             Velocities(nCarbon) = SigmaMomentum*gs2/MassCarb
! 
!          END DO
! 
!       END SUBROUTINE ScatteringConditions
! 
! 
! ! ******************************************************************************************      
! 
!       REAL FUNCTION MinimizePotential( Coords, Mask ) RESULT( Pot )
!          IMPLICIT NONE
!          REAL, INTENT(INOUT), DIMENSION(:)               :: Coords
!          LOGICAL, INTENT(IN), DIMENSION(size(Coords)) :: Mask
! 
!          INTEGER :: NrDimension, NrOptimization
!          INTEGER :: iIter, iCoord
!          REAL, DIMENSION(size(Coords)) :: Gradient
!          REAL :: Norm
! 
!          ! Set dimension number
!          NrDimension = size(Coords)
!          ! Set optimization coordinates nr
!          NrOptimization = count( Mask )
!          ! Check if the nr of dimension is compatible with the slab maximum size
!          CALL ERROR( (NrDimension > 124) .OR. (NrDimension < 4), "PotentialModule.MinimizePotential: wrong number of DoFs" )
! 
!          ! Cycle over steepest descent iterations
!          DO iIter = 1, MaxIter
! 
!             ! compute negative of the gradient
!             Pot = VHSticking( Coords, Gradient )
! 
!             ! compute norm of the gradient
!             Norm = 0.0
!             DO iCoord = 1, NrDimension
!                IF ( Mask( iCoord ) ) THEN
!                   Norm = Norm + Gradient(iCoord)**2
!                END IF
!             END DO
!             Norm = SQRT( Norm / NrOptimization )
! 
!             ! check convergence
!             IF (Norm < GradEps) EXIT
!       
!             ! move geometry along gradient
!             DO iCoord = 1, NrDimension
!                IF ( Mask( iCoord ) ) THEN
!                   Coords(iCoord) = Coords(iCoord) + Gradient(iCoord)
!                END IF
!             END DO
! 
!          END DO
! 
! #if defined(VERBOSE_OUTPUT)
!          WRITE(*,"(/,A,I6,A)") "Convergence in ", iIter, " steps"
! #endif
!          CALL WARN( iIter == MaxIter, "PotentialModule. MinimizePotential: convergence not reached" )
! 
!       END FUNCTION MinimizePotential
! 
!    ! ************************************************************************************************
! 
!    FUNCTION HessianOfTheSystem( AtPoint, MassHydro, MassCarb ) RESULT( Hessian )
!       IMPLICIT NONE
!       REAL, DIMENSION(4,4) :: Hessian
!       REAL, INTENT(IN)     :: MassHydro, MassCarb
!       REAL, DIMENSION(4), INTENT(IN) :: AtPoint
!       REAL, DIMENSION(4) :: Coordinates, FirstDerivative, MassVector
!       REAL :: Potential
!       INTEGER :: i, k
! 
!       REAL, DIMENSION(4), PARAMETER :: Deltas = (/ -2.0,    -1.0,    +1.0,    +2.0    /)
!       REAL, DIMENSION(4), PARAMETER :: Coeffs = (/ +1./12., -8./12., +8./12., -1./12. /) 
! 
!       REAL, DIMENSION(3), PARAMETER :: ForwardDeltas = (/  0.0,   +1.0,  +2.0   /)
!       REAL, DIMENSION(3), PARAMETER :: ForwardCoeffs = (/ -3./2., +2.0,  -1./2. /) 
! 
!       MassVector(1:3) = MassHydro
!       MassVector(4)   = MassCarb
! 
!       Hessian(:,:) = 0.0
! 
!       ! Compute the second derivatives for displacements of x and y
!       ! IMPORTANT!!! since rho = 0 is a singular value of the function, 
!       ! the derivative is computed slightly off the minimum, and is computed for x,y > 0
!       DO i = 1, 2
!          DO k = 1, size(ForwardDeltas)
! 
!             ! Define small displacement from the point where compute the derivative
!             Coordinates(:) = AtPoint(:)
!             IF ( Coordinates(i) < 0.0 ) THEN
!                Coordinates(i) = - Coordinates(i)
!             END IF
! 
!             IF ( Coordinates(i) < 0.001 ) THEN
!                Coordinates(i) = Coordinates(i) + 0.001 + ForwardDeltas(k)*SmallDelta
!             ELSE
!                Coordinates(i) = Coordinates(i) + ForwardDeltas(k)*SmallDelta
!             END IF
! 
!             ! Compute potential and forces in the displaced coordinate
!             Potential = VHFourDimensional( Coordinates, FirstDerivative )
!             FirstDerivative = - FirstDerivative
! 
!             ! Increment numerical derivative of the analytical derivative
!             Hessian(i,:) = Hessian(i,:) + ForwardCoeffs(k)*FirstDerivative(:)/SmallDelta
! 
!          END DO
!       END DO
! 
!       DO i = 3, 4
!          DO k = 1, size(Deltas)
! 
!             ! Define small displacement from the point where compute the derivative
!             Coordinates(:) = AtPoint(:)
!             Coordinates(i) = Coordinates(i) + Deltas(k)*SmallDelta
! 
!             ! Compute potential and forces in the displaced coordinate
!             Potential = VHFourDimensional( Coordinates, FirstDerivative )
!             FirstDerivative = - FirstDerivative
! 
!             ! Increment numerical derivative of the analytical derivative
!             Hessian(i,:) = Hessian(i,:) + Coeffs(k)*FirstDerivative(:)/SmallDelta
! 
!          END DO
!       END DO
! 
!       DO k = 1, 4
!          DO i = 1, 4
!             Hessian(i,k) = Hessian(i,k) / SQRT( MassVector(i)*MassVector(k) )
!          END DO
!       END DO
! 
! !       CALL TheOneWithMatrixPrintedLineAfterLine( Hessian )
! 
!    END FUNCTION HessianOfTheSystem
! 
! ! ******************************************************************************************      

!    
!       
! END MODULE PotentialModule
! 
