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
   USE InputField
   USE UnitConversion
   USE RandomNumberGenerator
   USE fsiesta

   PRIVATE

   PUBLIC :: SetupPotential, GetPotential, DisposePotential
   PUBLIC :: GetUnitCellDimensions, GetInitialPositions, GetAtomsNumber, GetMasses

   !> \name SYSTEM ID
   !> Integers number identifying the kind of potential adopted
   !> @{
   INTEGER, PARAMETER  ::  FREE_PARTICLES    = 0
   INTEGER, PARAMETER  ::  LJ_PAIR_POTENTIAL = 1
   INTEGER, PARAMETER  ::  SIESTA_ONTHEFLY   = 2
   !> @}

   !> Setup variable for the potential
   LOGICAL, SAVE :: PotentialModuleIsSetup = .FALSE.

   !> Identificative number of the potential currently used
   INTEGER, SAVE :: SystemNumber

   !> Number of atoms of the system
   INTEGER, SAVE :: AtomNo

   !> Periodic boundary conditions Flag
   LOGICAL, SAVE :: PBC = .FALSE.

   !> Size for cubic cells
   REAL, SAVE :: BoxSize
   !> Vectors of the unit cell
   REAL, DIMENSION(3,3), SAVE :: UnitVectors

   !> Label of the system for SIESTA
   CHARACTER(100) :: SystemLabel
   !> Atomic numbers of the atoms 
   INTEGER, DIMENSION(:), SAVE, ALLOCATABLE  :: AtomicNumbers
     !> Positions of the atoms
   REAL, DIMENSION(:,:), SAVE, ALLOCATABLE :: AtomicPositions

   !> Cutoff distance for a pair potential
   REAL, SAVE   :: CutOff = 100.0
   !> How many neighbour atoms are included in the pair potential summation
   INTEGER, SAVE :: NearPeriodicImages = 0
   !> Lattice translations of the near neighbour cells
   REAL, DIMENSION(:,:), ALLOCATABLE :: NearTranslations

   !> \name LENNARD-JONES POTENTIAL PARAMETERS
   !> Parameters of the Lennard-Jones pair potential
   !> @{
   REAL, SAVE :: LJ_WellDepth = 1.0       !< Potential well depth (energy)
   REAL, SAVE :: LJ_EquilDist = 1.0       !< Potential equilibrium distance (lenght)
   !> @}


!============================================================================================
                                       CONTAINS
!============================================================================================

   SUBROUTINE SetupPotential( PotentialFileName )
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN)   :: PotentialFileName

      CHARACTER(100)    :: ErrorMsg, String
      TYPE(InputFile)   :: PotentialData
      INTEGER           :: InputLength, InputEnergy, InputMass       ! Units of input data, defined from the input file
      INTEGER           :: iAtom
      LOGICAL           :: FileExists

      ! exit if module is setup
      IF ( PotentialModuleIsSetup ) RETURN

      __MPI_OnlyMasterBEGIN
      ! Open potential input file
      CALL OpenFile( PotentialData, PotentialFileName )

      ! Read input unit
      CALL SetFieldFromInput( PotentialData, "InputLength", InputLength,  1 )
      CALL SetFieldFromInput( PotentialData, "InputEnergy", InputEnergy,  3 )
      CALL SetFieldFromInput( PotentialData, "InputMass", InputMass,  8 )

      ! Read system number
      CALL SetFieldFromInput( PotentialData, "SystemNumber", SystemNumber )
      __MPI_OnlyMasterEND

      CALL MyMPI_BroadcastToSlaves( InputLength ) 
      CALL MyMPI_BroadcastToSlaves( InputEnergy ) 
      CALL MyMPI_BroadcastToSlaves( InputMass ) 
      CALL MyMPI_BroadcastToSlaves( SystemNumber ) 

      ! Check if the potential id number is valid
      WRITE(ErrorMsg,*) " PotentialModule.SetupPotential : no potential corresponds to Id = ", SystemNumber
      CALL ERROR ( .NOT. PotentialIdExists(SystemNumber),  ErrorMsg )

      ! Setup conversion factors from potential to internal units
      CALL Initialize_UnitConversion( PotentialUnits, InputLength, InputEnergy, InputMass, 11, 12, 15, 17 )

      ! Depending on system number, different setup operations
      SELECT CASE( SystemNumber )

         ! Gas of free particles
         CASE( FREE_PARTICLES )

            ! number of particles
            CALL SetFieldFromInput( PotentialData, "AtomNo", AtomNo, 1 )
            ! periodic system
            CALL SetFieldFromInput( PotentialData, "PBC", PBC, .FALSE. )
            ! Size of the simulation box
            IF ( PBC ) THEN
               CALL SetFieldFromInput( PotentialData, "BoxSize", BoxSize )
               CALL PBC_Setup( (/ BoxSize, 0.0, 0.0 /), (/ 0.0, BoxSize, 0.0 /), (/ 0.0, 0.0, BoxSize /) )
            END IF

         ! Gas of LJ particles
         CASE( LJ_PAIR_POTENTIAL )

            __MPI_OnlyMasterBEGIN
            ! number of particles
            CALL SetFieldFromInput( PotentialData, "AtomNo", AtomNo, 1 )
            ! periodic system
            CALL SetFieldFromInput( PotentialData, "PBC", PBC, .FALSE. )
            ! Size of the simulation box
            IF ( PBC ) THEN
               CALL SetFieldFromInput( PotentialData, "BoxSize", BoxSize )
            END IF

            ! Store the pair potential parameters
            CALL SetFieldFromInput( PotentialData, "WellDepth", LJ_WellDepth )
            CALL SetFieldFromInput( PotentialData, "EquilDist", LJ_EquilDist )
            LJ_WellDepth = LJ_WellDepth * EnergyConversion( PotentialUnits, InternalUnits )
            LJ_EquilDist = LJ_EquilDist * LengthConversion( PotentialUnits, InternalUnits )

            ! Set the cutoff distance of the pair potential
            CALL SetFieldFromInput(PotentialData,"CutOff",CutOff,LJ_EquilDist*LengthConversion(InternalUnits,PotentialUnits)*6.0)
            CutOff = CutOff * LengthConversion( PotentialUnits, InternalUnits )
            __MPI_OnlyMasterEND

            CALL MyMPI_BroadcastToSlaves( AtomNo )
            CALL MyMPI_BroadcastToSlaves( PBC )
            CALL MyMPI_BroadcastToSlaves( BoxSize )
            CALL MyMPI_BroadcastToSlaves( LJ_WellDepth )
            CALL MyMPI_BroadcastToSlaves( LJ_EquilDist )
            CALL MyMPI_BroadcastToSlaves( CutOff )

            IF ( PBC ) THEN
              CALL PBC_Setup( (/ BoxSize, 0.0, 0.0 /), (/ 0.0, BoxSize, 0.0 /), (/ 0.0, 0.0, BoxSize /) )
            END IF

            ! Check how many image atoms should be included in the pair potential summation
            CALL SetNearTranslations( )

         ! Gas of LJ particles
         CASE( SIESTA_ONTHEFLY )

            ! SystemLabel
            WRITE(SystemLabel,"(A,I0.3)") "PolyOnTheFly", 0
            ! HERE THE NUMBER ZERO SHOULD BE THE SLAVE ID NUMBER!

            ! number of particles
            CALL SetFieldFromInput( PotentialData, "AtomNo", AtomNo, 1 )
            ALLOCATE( AtomicNumbers(AtomNo), AtomicPositions(3,AtomNo) )

            ! periodic system
            CALL SetFieldFromInput( PotentialData, "PBC", PBC, .FALSE. )
            ! Vectors of the unit cell
            IF ( PBC ) THEN
               CALL SetFieldFromInput( PotentialData, "AUnitVector", UnitVectors(:,1) )
               CALL SetFieldFromInput( PotentialData, "BUnitVector", UnitVectors(:,2) )
               CALL SetFieldFromInput( PotentialData, "CUnitVector", UnitVectors(:,3) )
               UnitVectors = UnitVectors * LengthConversion( PotentialUnits, InternalUnits )
               CALL PBC_Setup( UnitVectors(:,1), UnitVectors(:,2), UnitVectors(:,3) )
            END IF

            ! Read labels of the atoms 
            CALL SetFieldFromInput( PotentialData, "AtomicNumbers", AtomicNumbers )
            ! Read Atomic positions
            DO iAtom = 1, AtomNo
               WRITE(String,"(A,I0.3)") "Atom",iAtom
               CALL SetFieldFromInput( PotentialData, trim(adjustl(String)), AtomicPositions(:,iAtom) )
            END DO
            AtomicPositions = AtomicPositions * LengthConversion( PotentialUnits, InternalUnits )

            ! Check if pseudopotentials are present
            DO iAtom = 1, AtomNo
               INQUIRE(FILE=TRIM(AtomicLabel(AtomicNumbers(iAtom)))//".psf" , EXIST=FileExists)
               CALL ERROR( .NOT. FileExists, " SetupPotential: missing pseudopotential "//TRIM(AtomicLabel(iAtom))//".psf" )
            END DO

            CALL WriteSIESTAInput()

            ! Initialize SIESTA input units
             call siesta_units( 'Ang', 'eV' )
            ! Initialize siesta processes
            CALL siesta_launch( TRIM(ADJUSTL(SystemLabel)), 1 )

      END SELECT

      __MPI_OnlyMasterBEGIN
      ! close input file
      CALL CloseFile( PotentialData )
      __MPI_OnlyMasterEND

      ! Module is now ready
      PotentialModuleIsSetup = .TRUE.

   END SUBROUTINE SetupPotential

!============================================================================================

   LOGICAL FUNCTION PotentialIdExists( IdNr )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IdNr

      PotentialIdExists =  ( ( IdNr == FREE_PARTICLES    ) .OR. &
                             ( IdNr == LJ_PAIR_POTENTIAL ) .OR. &
                             ( IdNr == SIESTA_ONTHEFLY   )  )

   END FUNCTION PotentialIdExists

!============================================================================================

   SUBROUTINE GetInitialPositions( X, RandomNr )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(OUT) :: X
      TYPE(RNGInternalState)          :: RandomNr
      INTEGER :: i, j
      REAL :: MinDistance, Distance
      REAL, DIMENSION(3) :: VecDist

      ! Error if module not have been setup yet
      CALL ERROR( .NOT. PotentialModuleIsSetup, " PotentialModule.InitialPositions : Module not Setup" )

      ! Check dimension of subroutine argument
      CALL ERROR( SIZE(X) /= 3*AtomNo,    " PotentialModule.InitialPositions : wrong dimension of coordinate array " )

      ! Depending on system number, different initialization operations
      SELECT CASE( SystemNumber )

         CASE( FREE_PARTICLES, LJ_PAIR_POTENTIAL )

            X( 1 ) = UniformRandomNr( RandomNr ) * BoxSize
            X( 2 ) = UniformRandomNr( RandomNr ) * BoxSize
            X( 3 ) = UniformRandomNr( RandomNr ) * BoxSize

            DO i = 2, AtomNo
               X( (i-1)*3+1 ) = UniformRandomNr( RandomNr ) * BoxSize
               X( (i-1)*3+2 ) = UniformRandomNr( RandomNr ) * BoxSize
               X( (i-1)*3+3 ) = UniformRandomNr( RandomNr ) * BoxSize

               DO
                  MinDistance = 1000.0
                  DO j = 1, i-1
                     VecDist = X( (i-1)*3+1 : i*3 ) - X( (j-1)*3+1 : j*3 )
                     Distance = SQRT( TheOneWithVectorDotVector( VecDist, VecDist ) )
                     MinDistance = MIN( MinDistance, Distance )
                  END DO
                  IF  ( MinDistance > 5.0 ) EXIT
                  X( (i-1)*3+1 ) = UniformRandomNr( RandomNr ) * BoxSize
                  X( (i-1)*3+2 ) = UniformRandomNr( RandomNr ) * BoxSize
                  X( (i-1)*3+3 ) = UniformRandomNr( RandomNr ) * BoxSize
               END DO
            END DO

         CASE( SIESTA_ONTHEFLY )

            DO i = 1, AtomNo
               X( (i-1)*3+1 : (i-1)*3+3 ) = AtomicPositions( :,i )
            END DO

      END SELECT

   END SUBROUTINE GetInitialPositions

!============================================================================================

   INTEGER FUNCTION GetAtomsNumber( )
      IMPLICIT NONE

      ! Error if module not have been setup yet
      CALL ERROR( .NOT. PotentialModuleIsSetup, " PotentialModule.GetDOFNumber : Module not Setup" )

      GetAtomsNumber = AtomNo

   END FUNCTION GetAtomsNumber

!============================================================================================

   FUNCTION GetMasses( )
      IMPLICIT NONE
      REAL, DIMENSION( 3*AtomNo ) :: GetMasses
      INTEGER :: i

      ! Error if module not have been setup yet
      CALL ERROR( .NOT. PotentialModuleIsSetup, " PotentialModule.GetMasses : Module not Setup" )

      ! Depending on system number, different operations
      SELECT CASE( SystemNumber )

         CASE( FREE_PARTICLES, LJ_PAIR_POTENTIAL )
            GetMasses(:) = 1.0 * MassConversion(PotentialUnits, InternalUnits)

         CASE( SIESTA_ONTHEFLY )
            DO i = 1, AtomNo
               GetMasses( (i-1)*3+1 : (i-1)*3+3 ) = AtomicMass( AtomicNumbers(i) )
            END DO

      END SELECT

   END FUNCTION GetMasses

!============================================================================================

   FUNCTION GetUnitCellDimensions( ) RESULT(UnitCell)
      IMPLICIT NONE
      REAL, DIMENSION(6) :: UnitCell

      ! Error if module not have been setup yet
      CALL ERROR( .NOT. PotentialModuleIsSetup, " PotentialModule.GetUnitCellDimensions : Module not Setup" )

      IF ( PBC ) THEN
         UnitCell(1) = SQRT( TheOneWithVectorDotVector( UnitVectors(:,1), UnitVectors(:,1) ) )
         UnitCell(2) = SQRT( TheOneWithVectorDotVector( UnitVectors(:,2), UnitVectors(:,2) ) )
         UnitCell(3) = SQRT( TheOneWithVectorDotVector( UnitVectors(:,3), UnitVectors(:,3) ) )
         UnitCell(4) = ACOS( TheOneWithVectorDotVector( UnitVectors(:,2), UnitVectors(:,3) ) / UnitCell(2) / UnitCell(3) )
         UnitCell(5) = ACOS( TheOneWithVectorDotVector( UnitVectors(:,1), UnitVectors(:,3) ) / UnitCell(1) / UnitCell(3) )
         UnitCell(6) = ACOS( TheOneWithVectorDotVector( UnitVectors(:,1), UnitVectors(:,2) ) / UnitCell(1) / UnitCell(2) )
         UnitCell(4:6) = UnitCell(4:6) 
      ELSE
         UnitCell = (/ REAL(1.E+100), REAL(1.E+100), REAL(1.E+100), REAL(90.), REAL(90.), REAL(90.) /)
      ENDIF

   END FUNCTION GetUnitCellDimensions

!============================================================================================

   REAL FUNCTION GetPotential( X, Force )
      IMPLICIT NONE
      REAL, DIMENSION(:), TARGET, INTENT(IN)  :: X
      REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Force
      REAL, DIMENSION(3,AtomNo) :: TmpForces
      REAL, DIMENSION(3,3)      :: Stress

      ! Error if module not have been setup yet
      CALL ERROR( .NOT. PotentialModuleIsSetup, " PotentialModule.GetPotential : Module not Setup" )

      ! Check dimension of function arguments
      CALL ERROR( SIZE(X) /= 3*AtomNo,    " PotentialModule.GetPotential : wrong dimension of coordinate array " )
      CALL ERROR( SIZE(X) /= SIZE(Force), " PotentialModule.GetPotential : wrong dimension of forces array " )
      
      ! Depending on system number, different operations
      SELECT CASE( SystemNumber )

         CASE( FREE_PARTICLES )
            GetPotential = 0.0
            Force(:) = 0.0

         CASE( LJ_PAIR_POTENTIAL )
            CALL PairPotential( X(:), GetPotential, Force(:) )

         CASE( SIESTA_ONTHEFLY ) 
            ! Copy coordianates in x,y,z format
            AtomicPositions = RESHAPE( X, (/ 3, AtomNo /) )
            AtomicPositions = AtomicPositions * MyConsts_Bohr2Ang 
            ! compute forces with siesta
            IF ( PBC ) THEN
               CALL siesta_forces( SystemLabel, AtomNo, AtomicPositions, UnitVectors*MyConsts_Bohr2Ang, &
                             GetPotential, TmpForces, Stress )
            ELSE
               CALL siesta_forces( label=SystemLabel, na=AtomNo, xa=AtomicPositions, energy=GetPotential, fa=TmpForces )
            END IF
            Force = RESHAPE( TmpForces, (/ 3*AtomNo /) )
            Force = Force / MyConsts_Hartree2eV * MyConsts_Bohr2Ang 
            GetPotential = GetPotential / MyConsts_Hartree2eV
      END SELECT

   END FUNCTION GetPotential

!============================================================================================

   SUBROUTINE DisposePotential(  )
      IMPLICIT NONE

      ! exit if module is not setup
      IF ( .NOT. PotentialModuleIsSetup ) RETURN

      PotentialModuleIsSetup = .FALSE.
      
   END SUBROUTINE DisposePotential

!============================================================================================

   SUBROUTINE SetNearTranslations( )
      IMPLICIT NONE
      INTEGER :: i, j, k
      REAL, DIMENSION(3,216) :: TmpNearTranslations
      REAL, DIMENSION(3) :: Vector
      REAL :: Distance

      NearPeriodicImages = 0
      IF ( PBC ) THEN
         ! Define how many neighbour cells are checked in the pair potential summation
         DO i = -5, +5
            DO j = -5, +5
               DO k = -5, +5
                  Vector = FractionalToCartesian( (/ REAL(i), REAL(j), REAL(k) /) )
                  Distance = SQRT( TheOneWithVectorDotVector( Vector, Vector ) )
                  IF ( Distance < CutOff .OR. ( abs(i) <= 1 .AND. abs(j) <= 1 .AND. abs(k) <= 1 )) THEN
                     NearPeriodicImages = NearPeriodicImages + 1
                     TmpNearTranslations(:,NearPeriodicImages) =  Vector
                  END IF
               END DO
            END DO
         END DO
!          PRINT*, " Number of cells included in the summation: ", NearPeriodicImages
         ALLOCATE( NearTranslations(3,NearPeriodicImages) )
         NearTranslations(:,:) =  TmpNearTranslations(:,1:NearPeriodicImages )

!          PRINT*, " "
!          DO i = 1, NearPeriodicImages
!             PRINT*, " Translation # ", i, "  Vector: ", NearTranslations(:,i)
!          END DO

      ELSE 
         NearPeriodicImages = 1
         ALLOCATE( NearTranslations(3,NearPeriodicImages) )
         NearTranslations(:,1) =  (/ 0., 0., 0. /) 
      ENDIF

   END SUBROUTINE SetNearTranslations

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
         
         DO jAtom = 1, AtomNo
            ! Extract distance vector between the i-th and j-th atoms
            FirstDist = iCoord( : ) - Positions( (jAtom-1)*3+1 : jAtom*3 )

            DO iTrasl = 1, NearPeriodicImages

               ! Compute periodic image of the distance
               TranslatedDist = FirstDist + NearTranslations(:,iTrasl) 
               Distance = SQRT( TheOneWithVectorDotVector( TranslatedDist , TranslatedDist ) )

               IF (( iAtom == jAtom ) .AND. ( ALL(NearTranslations(:,iTrasl) == (/ 0., 0., 0. /)) )) CYCLE

               IF ( Distance > CutOff ) CYCLE

               ! Compute the pair potential
               V = V + LennardJones( Distance, LJDerivative )

               ! Update forces
               Forces( (iAtom-1)*3+1 ) = Forces( (iAtom-1)*3+1 ) - LJDerivative * ( TranslatedDist(1) ) / Distance
               Forces( (iAtom-1)*3+2 ) = Forces( (iAtom-1)*3+2 ) - LJDerivative * ( TranslatedDist(2) ) / Distance
               Forces(  iAtom*3      ) = Forces(  iAtom*3      ) - LJDerivative * ( TranslatedDist(3) ) / Distance
               Forces( (jAtom-1)*3+1 ) = Forces( (jAtom-1)*3+1 ) + LJDerivative * ( TranslatedDist(1) ) / Distance
               Forces( (jAtom-1)*3+2 ) = Forces( (jAtom-1)*3+2 ) + LJDerivative * ( TranslatedDist(2) ) / Distance
               Forces(  jAtom*3      ) = Forces(  jAtom*3      ) + LJDerivative * ( TranslatedDist(3) ) / Distance

            END DO
         END DO
      END DO

      V = 0.5 * V
      Forces(:) = 0.5 * Forces(:)
      
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
      
   REAL FUNCTION LennardJones( Distance, Derivative ) RESULT(V) 
      IMPLICIT NONE
      REAL, INTENT(IN)  :: Distance
      REAL, INTENT(OUT) :: Derivative 

      V = LJ_WellDepth *( (LJ_EquilDist/Distance)**12 - 2.0*(LJ_EquilDist/Distance)**6 )
      Derivative = - 12.0 * LJ_WellDepth / Distance * ( (LJ_EquilDist/Distance)**12 - (LJ_EquilDist/Distance)**6 )
   END FUNCTION LennardJones


   SUBROUTINE WriteSIESTAInput()
      IMPLICIT NONE
      INTEGER :: OutputUnit 
      INTEGER, DIMENSION(AtomNo) :: UniqueKindOfAtoms
      INTEGER :: KindOfAtomsNr, i, j

      ! Extract the unique kind of atoms from the list of atomic numbers
      UniqueKindOfAtoms = AtomicNumbers
      CALL RemoveDups(UniqueKindOfAtoms, KindOfAtomsNr) 

      ! Open output unit
      OutputUnit = LookForFreeUnit()
      OPEN( UNIT=OutputUnit, FILE=TRIM(ADJUSTL(SystemLabel))//".fdf" )

      WRITE( OutputUnit, 300 ) "PolyOnTheFly - Ab Initio Forces", TRIM(ADJUSTL(SystemLabel)), AtomNo, KindOfAtomsNr, &
                               "forces", "100 Ry", "T"

      WRITE( OutputUnit, * ) " OccupationFunction         MP      "
      WRITE( OutputUnit, * ) " OccupationMPOrder          5       "
      WRITE( OutputUnit, * ) " ElectronicTemperature      50.0 K  "
      WRITE( OutputUnit, * ) " PAO.BasisSize              DZ      "
      WRITE( OutputUnit, * ) " XC.functional              GGA     "
      WRITE( OutputUnit, * ) " XC.authors                 RPBE    "
      WRITE( OutputUnit, * ) " DM.MixingWeight            0.500   "
      WRITE( OutputUnit, * ) " DM.NumberPulay             5       "
      WRITE( OutputUnit, * ) " DM.Tolerance               1.0d-4  "

      WRITE( OutputUnit, 400 ) 
      DO i = 1, KindOfAtomsNr
         WRITE( OutputUnit, 401 ) i,  UniqueKindOfAtoms(i), ADJUSTR(AtomicLabel(UniqueKindOfAtoms(i)))
      END DO
      WRITE( OutputUnit, 402 )

      WRITE( OutputUnit, 500 ) "Bohr"
      DO i = 1, AtomNo
         DO j = 1, KindOfAtomsNr
            IF ( UniqueKindOfAtoms(j) == AtomicNumbers(i) ) EXIT
         END DO
         WRITE( OutputUnit, 501 ) AtomicPositions(:,i), j
      END DO
      WRITE( OutputUnit, 502 ) 

   300 FORMAT( 3X,"SystemName",     T30,A40, /,  &
               3X,"SystemLabel",    T30,A40, /,  &
               3X,"NumberOfAtoms",  T30,I40, /,  &
               3X,"NumberOfSpecies",T30,I40, 2/, &
               3X,"MD.TypeOfRun",   T30,A40, 2/, &
               3X,"MeshCutoff",     T30,A40, /,  &
               3X,"DM.UseSaveDM",   T30,A40, /   )

   400 FORMAT( 3X,"%block ChemicalSpeciesLabel" )
   401 FORMAT( 3X,I4,I4,A4 )
   402 FORMAT( 3X,"%endblock ChemicalSpeciesLabel", / )

   500 FORMAT( 3X,"AtomicCoordinatesFormat", T30,A40, /, &
               3X,"%block AtomicCoordinatesAndAtomicSpecies" )
   501 FORMAT( 3X, 3F15.6,I10 )
   502 FORMAT( 3X,"%endblock AtomicCoordinatesAndAtomicSpecies", / )

   END SUBROUTINE WriteSIESTAInput


   FUNCTION AtomicLabel( AtomicNumber )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: AtomicNumber
      CHARACTER(3)        :: AtomicLabel
      CHARACTER(100)      :: ErrorMsg

      SELECT CASE( AtomicNumber )
         CASE ( 1 )
            AtomicLabel = "H"
         CASE ( 2 )
            AtomicLabel = "He"
         CASE ( 3 )
            AtomicLabel = "Li"
         CASE ( 4 )
            AtomicLabel = "Be"
         CASE ( 5 )
            AtomicLabel = "B"
         CASE ( 6 )
            AtomicLabel = "C"
         CASE ( 7 )
            AtomicLabel = "N"
         CASE ( 8 )
            AtomicLabel = "O"
         CASE ( 9 )
            AtomicLabel = "F"
         CASE ( 10 )
            AtomicLabel = "Ne"
         CASE ( 11 )
            AtomicLabel = "Na"
         CASE ( 12 )
            AtomicLabel = "Mg"
         CASE ( 13 )
            AtomicLabel = "Al"
         CASE ( 14 )
            AtomicLabel = "Si"
         CASE ( 15 )
            AtomicLabel = "P"
         CASE ( 16 )
            AtomicLabel = "S"
         CASE ( 17 )
            AtomicLabel = "Cl"
         CASE ( 18 )
            AtomicLabel = "Ar"
         CASE DEFAULT
            WRITE( ErrorMsg, * ) " AtomicLabel: atomic number ",AtomicNumber," is not available "
            CALL AbortWithError( ErrorMsg )
      END SELECT

   END FUNCTION AtomicLabel


   FUNCTION AtomicMass( AtomicNumber )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: AtomicNumber
      REAL                :: AtomicMass
      CHARACTER(100)      :: ErrorMsg

      SELECT CASE( AtomicNumber )
         CASE ( 1 )
            AtomicMass = 1.00782503207
         CASE ( 2 )
            AtomicMass = 4.00260325415
         CASE ( 3 )
            AtomicMass = 7.01600455
         CASE ( 4 )
            AtomicMass = 9.0121822
         CASE ( 5 )
            AtomicMass = 11.0093054
         CASE ( 6 )
            AtomicMass = 12.0000000
         CASE ( 7 )
            AtomicMass = 14.0030740048
         CASE ( 8 )
            AtomicMass = 15.99491461956
         CASE ( 9 )
            AtomicMass = 18.99840322
         CASE ( 10 )
            AtomicMass = 19.9924401754
         CASE ( 11 )
            AtomicMass = 22.9897692809
         CASE ( 12 )
            AtomicMass = 23.985041700
         CASE ( 13 )
            AtomicMass = 26.98153863
         CASE ( 14 )
            AtomicMass = 27.9769265325
         CASE ( 15 )
            AtomicMass = 30.97376163
         CASE ( 16 )
            AtomicMass = 31.97207100
         CASE ( 17 )
            AtomicMass = 34.96885268
         CASE ( 18 )
            AtomicMass = 39.9623831225
         CASE DEFAULT
            WRITE( ErrorMsg, * ) " AtomicMass: atomic number ",AtomicNumber," is not available "
            CALL AbortWithError( ErrorMsg )
      END SELECT
      AtomicMass = AtomicMass * MyConsts_Uma2Au

   END FUNCTION AtomicMass


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
