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
!>  \arg 8 November 2014: input and setup are now separated, data synchronization
!>                        for MPI application has been added
!>  \arg 8 November 2014: SetNearTranslations has been moved to PBC Module,
!>                        new subs to query the PBC defined in the potential module
!
!>  \todo          \arg logging of input and setup subroutines
!>                 \arg SIESTA: run calculation in different dir
!>                 \arg SIESTA: improve definition of input directives
!>                 \arg define initial velocities according to maxwell-boltzmann
!>                 \arg include minimization of the potential and hessian computation
!>                 \arg move atomic dictionaries to another module
!>                 
!***************************************************************************************
MODULE PotentialModule
#include "preprocessoptions.cpp"
   USE PeriodicBoundary
   USE InputField
   USE UnitConversion
   USE RandomNumberGenerator
   USE fsiesta
   USE DFTBWrapper

   PRIVATE

   ! Setup subroutines
   PUBLIC  ::  InputPotential, SyncroPotentialDataAcrossNodes, SetupPotential
   ! Queries on the periodic boundary conditions
   PUBLIC  ::  PotentialIsPeriodic, GetPotentialUnitCell, GetUnitCellDimensions
   ! Queries on the number of degrees of freedom
   PUBLIC  ::  GetAtomsNumber, GetDoFNumber, GetMasses
   ! Compute initial positions for the system in agreement with the defined potential
   PUBLIC  ::  GetInitialPositions
   ! Get the potential and the forces at given coordinates
   PUBLIC  ::  GetPotential
   ! Dispose subroutine
   PUBLIC  ::  DisposePotential

   
   !> \name SYSTEM ID
   !> Integers number identifying the kind of potential adopted
   !> @{
   INTEGER, PARAMETER  ::  FREE_PARTICLES    = 0
   INTEGER, PARAMETER  ::  LJ_PAIR_POTENTIAL = 1
   INTEGER, PARAMETER  ::  SIESTA_ONTHEFLY   = 2
   INTEGER, PARAMETER  ::  DFTB_ONTHEFLY     = 3
   !> @}

   !> Setup variable for the potential
   LOGICAL, SAVE :: PotentialModuleIsSetup = .FALSE.

   !> Identificative number of the potential currently used
   INTEGER, SAVE :: SystemNumber

   !> Number of atoms of the system
   INTEGER, SAVE :: AtomNo
   !> Number of degrees of freedom of the system
   INTEGER, SAVE :: DoFNo

   !> Periodic boundary conditions Flag
   LOGICAL, SAVE :: PBC = .FALSE.
   !> Vectors of the unit cell
   REAL, DIMENSION(3,3), SAVE :: UnitVectors = RESHAPE( (/1.E+100, 0.0, 0.0 , 0.0, 1.E+100, 0.0, 0.0, 0.0, 1.E+100/), (/3,3/) )
   !> Atomic postions read from input file
   REAL, DIMENSION(:,:), SAVE, ALLOCATABLE  :: AtomicPositions 

   !> \name PAIR POTENTIALS 
   !> Parameters for a generic pair potential 
   !> @{
   REAL, SAVE                      :: CutOff = 100.0   !< Cutoff distance for a pair potential
   TYPE(PBC_NearCellTranslations)  :: NeighbourCells   !< Translations to the neighbour cells within the cutoff radius 
   !> @}
   
   !> \name LENNARD-JONES POTENTIAL PARAMETERS
   !> Parameters of the Lennard-Jones pair potential
   !> @{
   REAL, SAVE :: LJ_WellDepth = 1.0                            !< Potential well depth (energy)
   REAL, SAVE :: LJ_EquilDist = 1.0                            !< Potential equilibrium distance (length)
   !> @}
   
   !> \name MORSE POTENTIAL PARAMETERS
   !> Parameters of the Morse pair potential
   !> @{
   REAL, SAVE :: MorseDe        = 1.0                            !< Potential well depth (energy)
   REAL, SAVE :: MorseEquilDist = 1.0                            !< Potential equilibrium distance (length)
   REAL, SAVE :: MorseAlpha     = 1.0                            !< Potential curvature (1/length)
   !> @}

   !> \name SIESTA ON-THE-FLY FORCES
   !> Variables for SIESTA
   !> @{
   CHARACTER(100), SAVE                     :: SystemLabel     !< Label of the system for SIESTA
   INTEGER, DIMENSION(:), SAVE, ALLOCATABLE :: AtomicNumbers   !< Atomic numbers of the atoms 
   TYPE(Units), SAVE                        :: SIESTAUnits     !< Units definition for SIESTA
   !> @}

   !> \name DFTB+ ON-THE-FLY FORCES
   !> Variables for DFTB+
   !> @{
   CHARACTER(100), SAVE                          :: LocalDir        !< local directory where to run DFTB+
   CHARACTER(100), SAVE                          :: HSDInputFile    !< name of the HSD file
   CHARACTER(3), DIMENSION(:), SAVE, ALLOCATABLE :: AtomicLabels    !< Atomic labels of the atoms 
   CHARACTER(120), DIMENSION(200), SAVE          :: ReadFileContent !< Lines of the hsd file other than the geometry
   CHARACTER(120), SAVE                          :: AtomTypeString  !< string defining the DFTB+ atom types
   INTEGER, DIMENSION(:), SAVE, ALLOCATABLE      :: AtomTypes       !< Atom types of the atoms 
   !> @}


!============================================================================================
                                       CONTAINS
!============================================================================================


   !*******************************************************************************
   !                               InputPotential
   !*******************************************************************************
   !>  Read from input file the variables to define the potential used in the
   !>  dynamical simulation, process and convert them to internal units.
   !>  Input file and unit definitions are given to the subroutines as input
   !>  arguments, via the specific data type TYPE(InputFile) [module InputField]
   !>  and TYPE(Units) [module UnitConversion]. Then the input file definition
   !>  and units definition are elsewhere, and they can be the very same of the 
   !>  other part of the input mechanism.
   !> 
   !> @param    PotentialData   InputFile datatype - I/O unit for reading
   !> @param    PotentialUnits  Units datatype - conversion factors
   !*******************************************************************************

   SUBROUTINE InputPotential( PotentialData, PotentialUnits )
      IMPLICIT NONE
      TYPE(InputFile)               :: PotentialData
      TYPE (Units), INTENT(IN)      :: PotentialUnits
      
      CHARACTER(100)    :: ErrorMsg, String
      INTEGER           :: iAtom
      
      ! If the potential has been already set up, print warning and skip the rest of the subroutine
      CALL WARN( PotentialModuleIsSetup, " InputPotential: Potential module has been already set up " ) 
      IF ( PotentialModuleIsSetup ) RETURN
      
      ! Read system id number
      CALL SetFieldFromInput( PotentialData, "SystemNumber", SystemNumber )
      
      ! Check if the potential id number is valid
      WRITE(ErrorMsg,*) " PotentialModule.SetupPotential : no potential corresponds to Id = ", SystemNumber
      CALL ERROR ( .NOT. PotentialIdExists(SystemNumber),  ErrorMsg )

     
      ! SYSTEM SPECIFIC VARIABLES 
      
      SELECT CASE( SystemNumber )

         ! ===========================================================================================
         ! Gas of free particles
         CASE( FREE_PARTICLES )
         
            ! number of particles
            CALL SetFieldFromInput( PotentialData, "AtomNo", AtomNo, 1 )
            
            ! periodic system and size of the simulation cell
            CALL SetFieldFromInput( PotentialData, "PBC", PBC, .FALSE. )
            IF ( PBC ) THEN
               CALL SetFieldFromInput( PotentialData, "AUnitVector", UnitVectors(:,1) )
               CALL SetFieldFromInput( PotentialData, "BUnitVector", UnitVectors(:,2) )
               CALL SetFieldFromInput( PotentialData, "CUnitVector", UnitVectors(:,3) )
               UnitVectors = UnitVectors * LengthConversion( PotentialUnits, InternalUnits )
            END IF

         ! ===========================================================================================
         ! Gas of Lennard-Jones particles
         CASE( LJ_PAIR_POTENTIAL )

            ! number of particles
            CALL SetFieldFromInput( PotentialData, "AtomNo", AtomNo, 1 )
            
            ! periodic system and size of the simulation cell
            CALL SetFieldFromInput( PotentialData, "PBC", PBC, .FALSE. )
            IF ( PBC ) THEN
               CALL SetFieldFromInput( PotentialData, "AUnitVector", UnitVectors(:,1) )
               CALL SetFieldFromInput( PotentialData, "BUnitVector", UnitVectors(:,2) )
               CALL SetFieldFromInput( PotentialData, "CUnitVector", UnitVectors(:,3) )
               UnitVectors = UnitVectors * LengthConversion( PotentialUnits, InternalUnits )
            END IF

            ! Store the pair potential parameters
            CALL SetFieldFromInput( PotentialData, "WellDepth", LJ_WellDepth )
            CALL SetFieldFromInput( PotentialData, "EquilDist", LJ_EquilDist )
            ! Set the cutoff distance of the pair potential
            CALL SetFieldFromInput( PotentialData, "CutOff", CutOff, LJ_EquilDist*6.0 )

            ! Convert them to internal units
            LJ_WellDepth = LJ_WellDepth * EnergyConversion( PotentialUnits, InternalUnits )
            LJ_EquilDist = LJ_EquilDist * LengthConversion( PotentialUnits, InternalUnits )
            CutOff       = CutOff       * LengthConversion( PotentialUnits, InternalUnits )
 
         ! ===========================================================================================
         ! SIESTA on-the-fly computation of the forces
         CASE( SIESTA_ONTHEFLY )

            ! number of particles
            CALL SetFieldFromInput( PotentialData, "AtomNo", AtomNo, 1 )
            
            ! periodic system and size of the simulation cell
            CALL SetFieldFromInput( PotentialData, "PBC", PBC, .FALSE. )
            IF ( PBC ) THEN
               CALL SetFieldFromInput( PotentialData, "AUnitVector", UnitVectors(:,1) )
               CALL SetFieldFromInput( PotentialData, "BUnitVector", UnitVectors(:,2) )
               CALL SetFieldFromInput( PotentialData, "CUnitVector", UnitVectors(:,3) )
               UnitVectors = UnitVectors * LengthConversion( PotentialUnits, InternalUnits )
            END IF

            ! Allocate memory to store initial definition of atomic info
            ALLOCATE( AtomicNumbers(AtomNo), AtomicPositions(3,AtomNo) )

            ! Read labels of the atoms 
            CALL SetFieldFromInput( PotentialData, "AtomicNumbers", AtomicNumbers )
            
            ! Read Atomic positions and convert them to internal units
            DO iAtom = 1, AtomNo
               WRITE(String,"(A,I0.3)") "Atom",iAtom
               CALL SetFieldFromInput( PotentialData, trim(adjustl(String)), AtomicPositions(:,iAtom) )
            END DO
            AtomicPositions = AtomicPositions * LengthConversion( PotentialUnits, InternalUnits )
 
        ! ===========================================================================================
        ! DFTB+ on-the-fly computation of the forces
         CASE( DFTB_ONTHEFLY )

            ! name of the hsd input file
            CALL SetFieldFromInput( PotentialData, "HSDInputFile", HSDInputFile, "dftb_in.hsd" )

            ! Initialize ReadFileContent with comment lines
            ReadFileContent(:) = "#"

            ! Read from DFTB+ input file the parameters for the calculation and the initial geometry
            CALL DFTBInputParser( HSDInputFile, ReadFileContent, PBC, UnitVectors, AtomNo, AtomicLabels, AtomicPositions )

            ! When XYZ file is given, read initial set of snapshots from XYZ
!             CALL ReadXYZSnapshots()

      END SELECT
      
   END SUBROUTINE InputPotential

   
!============================================================================================


   !*******************************************************************************
   !                               PotentialIdExists
   !*******************************************************************************
   !>  Check if the input integer number corresponds to a system potential
   !>  definition which is coded in the current module. Should be updated 
   !>  as well if an additional potential is included in the module
   !> 
   !> @param       IdNr            Integer number, possible system id
   !*******************************************************************************

   LOGICAL FUNCTION PotentialIdExists( IdNr )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IdNr

      PotentialIdExists =  ( ( IdNr == FREE_PARTICLES    ) .OR. &
                             ( IdNr == LJ_PAIR_POTENTIAL ) .OR. &
                             ( IdNr == DFTB_ONTHEFLY     ) .OR. &
                             ( IdNr == SIESTA_ONTHEFLY   )  )

   END FUNCTION PotentialIdExists

   
!============================================================================================


   !*******************************************************************************
   !                    SyncroPotentialDataAcrossNodes
   !*******************************************************************************
   !>  Syncronize potential data across all the nodes in a MPI run.
   !>  Allocatable arrays may be allocated on the master node, while
   !>  on the slaves they might be non allocated since the slaves where 
   !>  idle during the input procedure. Pay attantion to allocation! 
   !*******************************************************************************

   SUBROUTINE SyncroPotentialDataAcrossNodes( )
      IMPLICIT NONE
      INTEGER           :: iAtom
      
      ! If the potential has been already set up, print warning and skip the rest of the subroutine
      CALL WARN( PotentialModuleIsSetup, " SyncroPotentialDataAcrossNodes: Potential module has been already set up " ) 
      IF ( PotentialModuleIsSetup ) RETURN
    
      ! Syncro GENERAL VARIABLES that are relevant for all the systems
      CALL MyMPI_BroadcastToSlaves( SystemNumber )
      CALL MyMPI_BroadcastToSlaves( AtomNo )
      CALL MyMPI_BroadcastToSlaves( PBC )
      IF ( PBC ) THEN
         CALL MyMPI_BroadcastToSlaves( UnitVectors(:,1) )
         CALL MyMPI_BroadcastToSlaves( UnitVectors(:,2) )
         CALL MyMPI_BroadcastToSlaves( UnitVectors(:,3) )
      END IF
      
      ! synchro SYSTEM SPECIFIC VARIABLES 
      SELECT CASE( SystemNumber )

         ! Gas of free particles
         CASE( FREE_PARTICLES )
         
            ! --------------------------------------|
            ! NO ADDITIONAL VARIABLE TO SYNCHRONIZE |
            ! --------------------------------------|
            
         ! Gas of Lennard-Jones particles
         CASE( LJ_PAIR_POTENTIAL )

            CALL MyMPI_BroadcastToSlaves( LJ_WellDepth )
            CALL MyMPI_BroadcastToSlaves( LJ_EquilDist )
            CALL MyMPI_BroadcastToSlaves( CutOff )

         ! SIESTA on-the-fly computation of the forces
         CASE( SIESTA_ONTHEFLY )

            ! Allocate memory to store initial definition of atomic info
            IF ( .NOT. ALLOCATED(AtomicNumbers) )    ALLOCATE( AtomicNumbers(AtomNo) )
            IF ( .NOT. ALLOCATED(AtomicPositions) )  ALLOCATE( AtomicPositions(3,AtomNo) )

            ! Synchro atomic labes and atomic positions
            CALL MyMPI_BroadcastToSlaves( AtomicNumbers )
            DO iAtom = 1, AtomNo
               CALL MyMPI_BroadcastToSlaves( AtomicPositions(:,iAtom) )
            END DO

        ! DFTB+ on-the-fly computation of the forces
         CASE( DFTB_ONTHEFLY )

            ! Allocate memory to store initial definition of atomic info
            IF ( .NOT. ALLOCATED(AtomicLabels) )    ALLOCATE( AtomicLabels(AtomNo) )
            IF ( .NOT. ALLOCATED(AtomicPositions) )  ALLOCATE( AtomicPositions(3,AtomNo) )

            ! Synchro atomic labes and atomic positions
            CALL MyMPI_BroadcastToSlaves( AtomicLabels )
            DO iAtom = 1, AtomNo
               CALL MyMPI_BroadcastToSlaves( AtomicPositions(:,iAtom) )
            END DO

            ! Synchro the DFTB+ hamiltonian input lines
            CALL MyMPI_BroadcastToSlaves( ReadFileContent )

      END SELECT
      
   END SUBROUTINE SyncroPotentialDataAcrossNodes


!============================================================================================

   
   !*******************************************************************************
   !                    SetupPotential
   !*******************************************************************************
   !>  Once data is read and synchronized on all nodes, setup data
   !>  which is needed to compute the potential.
   !>  E.g.: setup the communication interface between the code and 
   !>  the external electronic structure program, define the variables 
   !>  for pair potential summation, etc etc ...
   !*******************************************************************************

   SUBROUTINE SetupPotential(  )
      IMPLICIT NONE
      LOGICAL           :: FileExists
      INTEGER           :: iAtom, nUnique, jUnique
      CHARACTER(3), DIMENSION(100) :: UniqueAtomicLabels

      ! exit if module is setup
      CALL WARN( PotentialModuleIsSetup, " SetupPotential: Potential module has been already set up " ) 
      IF ( PotentialModuleIsSetup ) RETURN

      ! Depending on system number, different setup operations
      SELECT CASE( SystemNumber )

         ! Gas of free particles
         CASE( FREE_PARTICLES )
         
            ! In a 3D gas each particle has 3 dof
            DoFNo = 3 * AtomNo

         ! Gas of Lennard-Jones particles
         CASE( LJ_PAIR_POTENTIAL )

            ! In a 3D gas each particle has 3 dof
            DoFNo = 3 * AtomNo

            ! Check how many image atoms should be included in the pair potential summation
            CALL PBC_SetNearTranslations( NeighbourCells, CutOff )

         ! SIESTA on-the-fly computation of the forces
         CASE( SIESTA_ONTHEFLY )

            ! In a 3D system, each particle has 3 dof
            DoFNo = 3 * AtomNo

            ! Set units for SIESTA
            CALL Initialize_UnitConversion( SIESTAUnits, UNITS_ANGSTROM, UNITS_EV, UNITS_AMU, UNITS_DEGREE, UNITS_FEMTOS, &
                                                         UNITS_KELVIN, UNITS_CMMINUS1 )

            ! SystemLabel, defined according to the mpi id
            WRITE(SystemLabel,"(A,I0.3)") "PolyOnTheFly", __MPI_CurrentNrOfProc

            ! Check if pseudopotentials are present
            DO iAtom = 1, AtomNo
               INQUIRE(FILE=TRIM(AtomicLabel(AtomicNumbers(iAtom)))//".psf" , EXIST=FileExists)
               CALL ERROR( .NOT. FileExists, " SetupPotential: missing pseudopotential "//TRIM(AtomicLabel(iAtom))//".psf" )
            END DO

            ! Write SIESTA input file
            CALL WriteSIESTAInput()

            ! Initialize SIESTA input units
            call siesta_units( 'Ang', 'eV' )
            ! Initialize siesta processes
            CALL siesta_launch( TRIM(ADJUSTL(SystemLabel)), 1 )

        ! DFTB+ on-the-fly computation of the forces
         CASE( DFTB_ONTHEFLY )

            ! In a 3D system, each particle has 3 dof
            DoFNo = 3 * AtomNo

            ! Define string with atomic labels and the kind definition of the atoms

            ! Find unique labels
            UniqueAtomicLabels(:) = " " 
            UniqueAtomicLabels(1) = AtomicLabels(1)
            nUnique = 1
            extcycle: DO iAtom = 2, AtomNo
               DO jUnique = 1, nUnique
                  IF ( TRIM(ADJUSTL(AtomicLabels(iAtom))) == TRIM(ADJUSTL(UniqueAtomicLabels(jUnique))) ) THEN
                     CYCLE extcycle
                  END IF
               END DO
               nUnique = nUnique + 1
               UniqueAtomicLabels(nUnique) = AtomicLabels(iAtom)
            END DO extcycle

            ! join labels in the atom type string 
            AtomTypeString = " "
            DO jUnique = 1, nUnique
               AtomTypeString = TRIM(ADJUSTL(AtomTypeString)) // " " // UniqueAtomicLabels(jUnique)
            END DO

            ! define the type of each atom
            ALLOCATE( AtomTypes( AtomNo ) )
            DO iAtom = 1, AtomNo
               DO jUnique = 1, nUnique
                  IF ( TRIM(ADJUSTL(AtomicLabels(iAtom))) == TRIM(ADJUSTL(UniqueAtomicLabels(jUnique))) ) THEN
                     AtomTypes(iAtom) = jUnique
                     EXIT
                  END IF
               END DO
            END DO

            ! Set scratch directory to run DFTB+
            WRITE(LocalDir,"(A,I0.3)") "POTF_DFTB_", __MPI_CurrentNrOfProc

            ! Create scratch directory and move there
            ! TODO: ADD ERROR MESSAGES!!!
            CALL EXECUTE_COMMAND_LINE( "mkdir -p "//TRIM(ADJUSTL(LocalDir)) )

      END SELECT

      ! Module is now ready
      PotentialModuleIsSetup = .TRUE.

   END SUBROUTINE SetupPotential


!============================================================================================


   !*******************************************************************************
   !                    DisposePotential
   !*******************************************************************************
   !>  Deallocate memory and set status variable to false
   !*******************************************************************************

   SUBROUTINE DisposePotential(  )
      IMPLICIT NONE

      ! exit if module is not setup
      IF ( .NOT. PotentialModuleIsSetup ) RETURN

      ! Deallocate allocated memory
      IF ( ALLOCATED(AtomicNumbers) )    DEALLOCATE( AtomicNumbers )
      IF ( ALLOCATED(AtomicLabels) )     DEALLOCATE( AtomicLabels )
      IF ( ALLOCATED(AtomicPositions) )  DEALLOCATE( AtomicPositions )
      IF ( ALLOCATED(AtomTypes) )        DEALLOCATE( AtomTypes )

      ! Set status variable
      PotentialModuleIsSetup = .FALSE.
      
   END SUBROUTINE DisposePotential

   
!============================================================================================


   !*******************************************************************************
   !             PotentialIsPeriodic and GetPotentialUnitCell
   !*******************************************************************************
   !>  Functions to query the PBC definitions within the potential subroutines.
   !>  PotentialIsPeriodic gives a logical value, which is TRUE if PBC are defined.
   !>  GetPotentialUnitCell gives the array with the unit cell definitions. A cubic 
   !>  huge cell is given if PBC are not defined.
   !>  GetUnitCellDimensions gives the unit cell as (|a|,|b|,|c|,alpha,beta,gamma)  
   !*******************************************************************************

   LOGICAL FUNCTION PotentialIsPeriodic(  )
      IMPLICIT NONE

      ! Error if module not have been setup yet
      CALL ERROR( .NOT. PotentialModuleIsSetup, " PotentialModule.PotentialIsPeriodic : Module not Setup" )

      PotentialIsPeriodic = PBC
      
   END FUNCTION PotentialIsPeriodic

   FUNCTION GetPotentialUnitCell( )
      IMPLICIT NONE
      REAL, DIMENSION(3,3) :: GetPotentialUnitCell
      
      ! Error if module not have been setup yet
      CALL ERROR( .NOT. PotentialModuleIsSetup, " PotentialModule.GetPotentialUnitCell : Module not Setup" )

      GetPotentialUnitCell = UnitVectors
      
   END FUNCTION GetPotentialUnitCell

   FUNCTION GetUnitCellDimensions( ) RESULT(UnitCell)
      IMPLICIT NONE
      REAL, DIMENSION(6) :: UnitCell

      ! Error if module not have been setup yet
      CALL ERROR( .NOT. PotentialModuleIsSetup, " PotentialModule.GetUnitCellDimensions : Module not Setup" )

      UnitCell(1) = SQRT( TheOneWithVectorDotVector( UnitVectors(:,1), UnitVectors(:,1) ) )
      UnitCell(2) = SQRT( TheOneWithVectorDotVector( UnitVectors(:,2), UnitVectors(:,2) ) )
      UnitCell(3) = SQRT( TheOneWithVectorDotVector( UnitVectors(:,3), UnitVectors(:,3) ) )
      UnitCell(4) = ACOS( TheOneWithVectorDotVector( UnitVectors(:,2), UnitVectors(:,3) ) / UnitCell(2) / UnitCell(3) )
      UnitCell(5) = ACOS( TheOneWithVectorDotVector( UnitVectors(:,1), UnitVectors(:,3) ) / UnitCell(1) / UnitCell(3) )
      UnitCell(6) = ACOS( TheOneWithVectorDotVector( UnitVectors(:,1), UnitVectors(:,2) ) / UnitCell(1) / UnitCell(2) )
      UnitCell(4:6) = UnitCell(4:6) 

   END FUNCTION GetUnitCellDimensions

   
!============================================================================================


   !*******************************************************************************
   !             GetAtomsNumber, GetDoFNumber and GetMasses
   !*******************************************************************************
   !>  Functions to query the degrees of freedom defined in the system.
   !>  GetAtomsNumber returns the number of particles of the system
   !>  GetDoFNumber returns the number of DoF defined in the potential.
   !>  GetMasses returns an array with the masses associated to the degrees of freedom.
   !*******************************************************************************

   INTEGER FUNCTION GetAtomsNumber( )
      IMPLICIT NONE

      ! Error if module not have been setup yet
      CALL ERROR( .NOT. PotentialModuleIsSetup, " PotentialModule.GetAtomsNumber : Module not Setup" )

      GetAtomsNumber = AtomNo

   END FUNCTION GetAtomsNumber
   
   INTEGER FUNCTION GetDoFNumber( )
      IMPLICIT NONE

      ! Error if module not have been setup yet
      CALL ERROR( .NOT. PotentialModuleIsSetup, " PotentialModule.GetDoFNumber : Module not Setup" )

      GetDoFNumber = DoFNo

   END FUNCTION GetDoFNumber
   
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
         CASE( DFTB_ONTHEFLY )
            DO i = 1, AtomNo
               GetMasses( (i-1)*3+1 : (i-1)*3+3 ) = AtomicMass( AtomicLabelToNumber( AtomicLabels(i)) )
            END DO
         CASE DEFAULT
               CALL AbortWithError( " PotentialModule.GetMasses : undefined SystemNumber " )
      END SELECT

   END FUNCTION GetMasses
   
   
!============================================================================================


   !*******************************************************************************
   !                    GetInitialPositions
   !*******************************************************************************
   !>  When the system potential is setup, use this subroutine to get the
   !>  initial coordinates of the particles.
   !> 
   !> @param X         real array of size DoFNo, on exit has the initial coords
   !> @param RandonNr  variable with the internal status of the RandomNumberGenerator
   !*******************************************************************************

   SUBROUTINE GetInitialPositions( X, RandomNr )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(OUT) :: X
      TYPE(RNGInternalState)          :: RandomNr
      INTEGER            :: i, j
      REAL               :: MinDistance, Distance
      REAL, DIMENSION(3) :: VecDist

      ! Error if module not have been setup yet
      CALL ERROR( .NOT. PotentialModuleIsSetup, " PotentialModule.InitialPositions : Module not Setup" )
      
      ! Check dimension of subroutine argument
      CALL ERROR( SIZE(X) /= DoFNo,    " PotentialModule.InitialPositions : wrong dimension of coordinate array " )

      ! Depending on system number, different initialization operations
      SELECT CASE( SystemNumber )

         ! For a 3D gas, distribute the initial coordinates randomly in the simulation box
         CASE( FREE_PARTICLES, LJ_PAIR_POTENTIAL )

            ! Set the coordinates of the first particle
            X( 1 ) = UniformRandomNr( RandomNr )
            X( 2 ) = UniformRandomNr( RandomNr )
            X( 3 ) = UniformRandomNr( RandomNr )

            ! Generate the following particles at a minimum distance from the others
            DO i = 2, AtomNo
               X( (i-1)*3+1 ) = UniformRandomNr( RandomNr )
               X( (i-1)*3+2 ) = UniformRandomNr( RandomNr )
               X( (i-1)*3+3 ) = UniformRandomNr( RandomNr )

               DO
                  ! Compute the minimum distance between each couples of particles
                  MinDistance = 1000.0
                  DO j = 1, i-1
                     VecDist = X( (i-1)*3+1 : i*3 ) - X( (j-1)*3+1 : j*3 )
                     Distance = SQRT( TheOneWithVectorDotVector( VecDist, VecDist ) )
                     MinDistance = MIN( MinDistance, Distance )
                  END DO
                  
                  ! Exit from the cycle if all the particles have a minimum distance of 2% (in internal coords)
                  IF  ( MinDistance > 0.02 ) EXIT
                  
                  ! Otherwise, generate new coordinates for the current particle
                  X( (i-1)*3+1 ) = UniformRandomNr( RandomNr )
                  X( (i-1)*3+2 ) = UniformRandomNr( RandomNr )
                  X( (i-1)*3+3 ) = UniformRandomNr( RandomNr )
               END DO
            END DO
            
            ! Transform internal coordinates to cartesian coordinates
            DO i = 1, AtomNo
               X( (i-1)*3+1:(i-1)*3+3 ) =  TheOneWithMatrixVectorProduct( UnitVectors ,  X( (i-1)*3+1:(i-1)*3+3 ) )
            END DO

         ! For SIESTA, use initial coordinates defined from input file
         CASE( SIESTA_ONTHEFLY )

            DO i = 1, AtomNo
               X( (i-1)*3+1 : (i-1)*3+3 ) = AtomicPositions( :,i )
            END DO

         ! For DFTB+, use initial coordinates defined from input file
         CASE( DFTB_ONTHEFLY )

            DO i = 1, AtomNo
               X( (i-1)*3+1 : (i-1)*3+3 ) = AtomicPositions( :,i )
            END DO

      END SELECT

   END SUBROUTINE GetInitialPositions

   
!============================================================================================


   !*******************************************************************************
   !                    GetPotential
   !*******************************************************************************
   !>  Compute potential and forces at given particle coordinates using the
   !>  potential defined. The function directly returns the potential value
   !>  while the forces are given as function argument.
   !> 
   !> @param   X       Array of size DoFNo, with the input coords of the particles
   !> @param   Force   Array of size DoFNo, on exit has the forces at X
   !> @returns GetPotential   Real value of the potential at X 
   !*******************************************************************************
   
   REAL FUNCTION GetPotential( X, Force )
      IMPLICIT NONE
      REAL, DIMENSION(:), TARGET, INTENT(IN)  :: X
      REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Force
      REAL, DIMENSION(3,AtomNo) :: TmpForces
      REAL, DIMENSION(3,3)      :: Stress

      ! Error if module not have been setup yet
      CALL ERROR( .NOT. PotentialModuleIsSetup, " PotentialModule.GetPotential : Module not Setup" )

      ! Check dimension of function arguments
      CALL ERROR( SIZE(X) /= DoFNo,       " PotentialModule.GetPotential : wrong dimension of coordinate array " )
      CALL ERROR( SIZE(X) /= SIZE(Force), " PotentialModule.GetPotential : wrong dimension of forces array " )
      
      ! Depending on system number, different operations
      SELECT CASE( SystemNumber )

         ! Gas of free particles, potential and forces are null
         CASE( FREE_PARTICLES )
            GetPotential = 0.0
            Force(:) = 0.0

         ! Gas of Lennard-Jones particles, compute pair potential forces
         CASE( LJ_PAIR_POTENTIAL )
            CALL PairPotential( X(:), GetPotential, Force(:) )

         ! Call SIESTA to compute the forces
         CASE( SIESTA_ONTHEFLY ) 
         
            ! Copy coordianates in x,y,z format, convert them to 
            AtomicPositions = RESHAPE( X, (/ 3, AtomNo /) )
            AtomicPositions = AtomicPositions * LengthConversion( InternalUnits, SIESTAUnits ) 
            
            ! compute forces with siesta
            IF ( PBC ) THEN
               CALL siesta_forces( SystemLabel, AtomNo, AtomicPositions, UnitVectors*LengthConversion(InternalUnits,SIESTAUnits), &
                             GetPotential, TmpForces, Stress )
            ELSE
               CALL siesta_forces( label=SystemLabel, na=AtomNo, xa=AtomicPositions, energy=GetPotential, fa=TmpForces )
            END IF
            
            ! Copy forces in vector format and convert them to internal units 
            Force = RESHAPE( TmpForces, (/ 3*AtomNo /) )
            Force = Force * ForceConversion(SIESTAUnits,InternalUnits)
            ! Convert potential to internal units
            GetPotential = GetPotential * EnergyConversion(SIESTAUnits,InternalUnits)
          
         ! Call DFTB+ to compute the forces
         CASE( DFTB_ONTHEFLY )   

            ! Copy coordianates in x,y,z format
            AtomicPositions = RESHAPE( X, (/ 3, AtomNo /) )
            
            ! Write input file for DFTB+
            CALL WriteDFTBInput( LocalDir, PBC, UnitVectors, AtomNo, AtomTypeString, AtomicPositions, AtomTypes, ReadFileContent  )

            ! Run calculation with DFTB+
            CALL EXECUTE_COMMAND_LINE( "cd "//TRIM(ADJUSTL(LocalDir))//"; dftb+ > log" )

            ! Read forces from DFTB+ output file
            CALL ReadDFTBForces( LocalDir, AtomNo, PBC, GetPotential, TmpForces, Stress )

            ! Copy forces in vector format 
            Force = RESHAPE( TmpForces, (/ 3*AtomNo /) )

      END SELECT

   END FUNCTION GetPotential

   
!============================================================================================
!                  PRIVATE SUBROUTINES
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

            DO iTrasl = 1, NeighbourCells%Nr

               ! Compute periodic image of the distance
               TranslatedDist = FirstDist + NeighbourCells%TranslVectors(:,iTrasl) 
               Distance = SQRT( TheOneWithVectorDotVector( TranslatedDist , TranslatedDist ) )

               IF (( iAtom == jAtom ) .AND. ( ALL(NeighbourCells%TranslVectors(:,iTrasl) == (/ 0., 0., 0. /)) )) CYCLE

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

    REAL FUNCTION MorseV( Position, Force ) RESULT(V) 
       IMPLICIT NONE
       REAL, INTENT(IN)  :: Position
       REAL, INTENT(OUT) :: Force 
 
       V = MorseDe * ( exp(-2.0*MorseAlpha*(Position-MorseEquilDist)) - 2.0 * exp(-MorseAlpha*(Position-MorseEquilDist)) )  
       Force = 2.0 * MorseAlpha * MorseDe * (  exp(-2.0*MorseAlpha*(Position-MorseEquilDist)) -         &
                     exp(-MorseAlpha*(Position-MorseEquilDist)) )  
 
    END FUNCTION MorseV
      
   REAL FUNCTION LennardJones( Distance, Derivative ) RESULT(V) 
      IMPLICIT NONE
      REAL, INTENT(IN)  :: Distance
      REAL, INTENT(OUT) :: Derivative 

      V = LJ_WellDepth *( (LJ_EquilDist/Distance)**12 - 2.0*(LJ_EquilDist/Distance)**6 )
      Derivative = - 12.0 * LJ_WellDepth / Distance * ( (LJ_EquilDist/Distance)**12 - (LJ_EquilDist/Distance)**6 )
   END FUNCTION LennardJones

!============================================================================================

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


   FUNCTION AtomicLabelToNumber( AtomicLabel ) RESULT( AtomicNumber )
      IMPLICIT NONE
      INTEGER                  :: AtomicNumber
      CHARACTER(3), INTENT(IN) :: AtomicLabel
      CHARACTER(100)           :: ErrorMsg
      INTEGER   ::  i 

      CHARACTER(3), PARAMETER, DIMENSION(18) :: &
      Labels =  (/ "H  ", "He ", "Li ", "Be ", "B  ", "C  ", "N  ", "O  ", "F  ", "Ne ", "Na ", "Mg ", "Al ", &
                   "Si ", "P  ", "S  ", "Cl ", "Ar " /)
      INTEGER, PARAMETER, DIMENSION(18) :: &
      Numbers = (/ (i, i=1, 18 ) /)

      i = 0
      DO 
         i = i + 1
         IF ( i > SIZE(Labels) ) THEN
            WRITE( ErrorMsg, * ) " AtomicLabel: atomic label ",AtomicLabel," is not available "
            CALL AbortWithError( ErrorMsg )
         END IF
         IF ( TRIM(ADJUSTL(AtomicLabel)) == TRIM(ADJUSTL(Labels(i))) ) THEN
            AtomicNumber = Numbers(i)
            EXIT
         END IF
      END DO

   END FUNCTION AtomicLabelToNumber



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
