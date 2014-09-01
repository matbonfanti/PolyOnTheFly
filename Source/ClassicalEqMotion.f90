!***************************************************************************************
!*                              MODULE ClassicalEqMotion
!***************************************************************************************
!
!>  \brief     Subroutines for integration of classical equations of motion
!>  \details   This module contains subroutines to evolve a classical system \n
!>             for one timestep. Some parameters (like timestep, thermostat options ...) \n
!>             are setup in a derived data type and setup as initialization of the module. \n
!>             Variables changing at each  timestep (positions, velocities, accelerations) \n
!>             instead are given as input of the propagation subroutines. \n
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             20 January 2013
!>
!***************************************************************************************
!
!>  \par Updates
!>  \arg 8 Novembre 2013: thermoswitch is now an optional argument of the
!>                        SetupThermostat subroutine (default: TRUE for all entries)
!>  \arg 28 November 2013: ring polymer propagation implemented, with symplectic
!>                         integrator
!>  \arg 9 December 2013: old integrators have been commented and will be removed
!
!>  \todo   clean up the code: leave only 1 propagator for RPMD and 1 propagator for
!>          normal MD, with or without langevin friction ( in case of RPMD, Parrinello
!>          algorithm, in case of regular MD, stochastic integrated algorith - both have
!>          Velocity Verlet as limit for gamma = 0 ) then rename subroutines for propagation
!>          and clean the datatype with the only data that are necessary finally fix 
!>          the setup subroutines including all the necessary checks and setup options
!>                 
!***************************************************************************************

MODULE ClassicalEqMotion
#include "preprocessoptions.cpp"
   USE RandomNumberGenerator
   USE FFTWrapper
   IMPLICIT NONE
   
      PRIVATE
      PUBLIC :: Evolution
      PUBLIC :: EvolutionSetup, SetupThermostat, SetupRingPolymer
      PUBLIC :: DisposeEvolutionData, DisposeThermostat, DisposeRingPolymer
      PUBLIC :: EOM_KineticEnergy, EOM_LangevinSecondOrder, EOM_RPMSymplectic, EOM_VelocityVerlet

      REAL, PARAMETER :: Over2Sqrt3 = 1.0 / ( 2.0 * SQRT(3.0) )

      !> Evolution datatype, storing all the general data required
      !> for integration of the EoM, with or without Langevin thermostat and with or without Ring Polymer MD
      TYPE Evolution
         INTEGER   :: NDoF                               !< Nr of degrees of freedom
         REAL, DIMENSION(:), POINTER :: Mass             !< Vector with the masses of the system
         REAL  :: dt                                     !< Time step of integration
         REAL  :: Gamma                                  !< Integrated Langevin friction for one timestep
         REAL, DIMENSION(:), POINTER :: ThermalNoise     !< Vector of the thermal noise sigma
         LOGICAL, DIMENSION(:), POINTER :: ThermoSwitch  !< Vector to set the thermostat on or off for the dof

         INTEGER :: NBeads                                 !< Nr of system replicas 
         REAL  :: RPFreq                                   !< Frequency of the harmonic force between the beads
         REAL, DIMENSION(:), POINTER :: NormModesFreq      !< Vector with the frequencies of the normal modes 
         REAL, DIMENSION(:,:), POINTER :: NormModesPropag  !< Propagation coefficients of normal coordinates and velocities
         TYPE(FFTHalfComplexType) :: RingNormalModes       !< FFT transform to compute normal modes of the RP 
         REAL, DIMENSION(:), POINTER :: GammaLang          !< Friction coeffs of PILE
         REAL, DIMENSION(:), POINTER :: AlphaLang          !< Set of coeffs for PILE integration
         REAL, DIMENSION(:), POINTER :: BetaLang           !< Set of coeffs for PILE integration

         LOGICAL :: HasThermostat = .FALSE.              !< Thermostat data has been setup
         LOGICAL :: HasRingPolymer = .FALSE.             !< Ring polymer data has been setup
         LOGICAL :: IsSetup = .FALSE.                    !< Evolution data has been setup
      END TYPE Evolution

   CONTAINS
   

!================================================================================================================================
!                                  SETUP SUBROUTINES
!================================================================================================================================


!*******************************************************************************
!                        EvolutionSetup
!*******************************************************************************
!> Setup evolution general data, store them in the Evolution datatype.
!>
!> @param EvolutionData    Evolution data type to setup
!> @param NDoF             Integer var with the number of degree of freedom
!> @param MassVector       NDoF-dim real vector with the masses
!*******************************************************************************
   SUBROUTINE EvolutionSetup( EvolutionData, NDoF, MassVector, TimeStep )
      IMPLICIT NONE
      TYPE( Evolution ), INTENT(INOUT)  :: EvolutionData
      INTEGER, INTENT(IN)               :: NDoF
      REAL, DIMENSION(:), INTENT(IN)    :: MassVector
      REAL, INTENT(IN)                  :: TimeStep

      ! Check dimension of mass vector
      CALL ERROR( size(MassVector) /= NDoF, "ClassicalEqMotion.EvolutionSetup: MassVector array dimension mismatch" )

      ! warn user if overwriting previously setup data
      CALL WARN( EvolutionData%IsSetup, "ClassicalEqMotion.EvolutionSetup: overwriting evolution data" ) 

      ! Store number of degree of freedom
      EvolutionData%NDoF = NDoF
      
      ! If necessary, deallocate array and then store masses
      IF (EvolutionData%IsSetup)  DEALLOCATE( EvolutionData%Mass )
      ALLOCATE( EvolutionData%Mass( NDoF ) )
      EvolutionData%Mass = MassVector
      
      ! Store timestep
      EvolutionData%dt = TimeStep
      
      ! evolutiondata is now setup
      EvolutionData%IsSetup = .TRUE.
      EvolutionData%HasRingPolymer = .FALSE.
      EvolutionData%HasThermostat = .FALSE.
      
#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A,I5,A,1F8.3)") "Evolution data type is setup: Nr DoF is ",NDoF," and TimeStep ", TimeStep
#endif
      
   END SUBROUTINE EvolutionSetup

   
!*******************************************************************************
!                       SetupThermostat
!*******************************************************************************
!> Setup thermostat data, store them in the Evolution datatype.
!>
!> @param EvolData     Evolution data type to setup
!> @param Gamma        Langevin friction coefficient
!> @param Temperature  Temperature of thermostat
!*******************************************************************************
   SUBROUTINE SetupThermostat( EvolData, Gamma, Temperature, ThermoSwitch )
      IMPLICIT NONE
      TYPE( Evolution ), INTENT(INOUT)            :: EvolData
      REAL, INTENT(IN)                            :: Gamma, Temperature
      LOGICAL, DIMENSION(:), INTENT(IN), OPTIONAL :: ThermoSwitch
      INTEGER :: iDoF
      
      ! error if trying to setup thermostat of a non-setup evolution type
      CALL ERROR( .NOT. EvolData%IsSetup, "ClassicalEqMotion.SetupThermostat: evolution data not setup" )
      ! warn user if overwriting previously setup data
      CALL WARN( EvolData%HasThermostat, "ClassicalEqMotion.SetupThermostat: overwriting thermostat data" ) 
      ! error if the thermostat switch has wrong dimension
      IF ( PRESENT(ThermoSwitch) )  CALL ERROR( size(ThermoSwitch) /= EvolData%NDoF , &
                                          "ClassicalEqMotion.SetupThermostat: thermostat switch mismatch" )

      ! Store gamma value
      EvolData%Gamma = Gamma
      
      ! Store the Thermostat switch array
      ALLOCATE( EvolData%ThermoSwitch(EvolData%NDoF) )
      IF ( PRESENT(ThermoSwitch) ) THEN
         EvolData%ThermoSwitch(:) = ThermoSwitch(:)
      ELSE
         EvolData%ThermoSwitch(:) = .TRUE.
      END IF

      ! Set standard deviations of the thermal noise
      IF ( .NOT. EvolData%HasThermostat )  ALLOCATE( EvolData%ThermalNoise( EvolData%NDoF )  )
      EvolData%ThermalNoise(:) = 0.0
      DO iDoF = 1, EvolData%NDoF
         IF ( EvolData%ThermoSwitch(iDoF) ) THEN
!             EvolData%ThermalNoise(iDoF) = sqrt( 2.0*Temperature*EvolData%Mass(iDoF)*Gamma/EvolData%dt )
            EvolData%ThermalNoise(iDoF) = sqrt( 2.0*Temperature*Gamma/EvolData%Mass(iDoF) )
         END IF
      END DO

      ! Set arrays for PILE integration when setting thermostat for a ring polymer propagator
      IF ( EvolData%HasRingPolymer ) THEN
         IF ( .NOT. EvolData%HasThermostat )     ALLOCATE( EvolData%GammaLang(EvolData%NBeads), &
                                                EvolData%AlphaLang(EvolData%NBeads), EvolData%BetaLang(EvolData%NBeads) )

         ! Set friction parameters for the ring normal modes
         EvolData%GammaLang(1) = Gamma
         DO iDoF = 2, EvolData%NBeads
            EvolData%GammaLang(iDoF) = 2.0 * EvolData%NormModesFreq(iDoF)
         END DO

         ! Set integration coefficients
         DO iDoF = 1, EvolData%NBeads
            EvolData%AlphaLang(iDoF) = EXP( - 0.5 * EvolData%dt * EvolData%GammaLang(iDoF) )
            EvolData%BetaLang(iDoF)  = SQRT( (1.0 - EvolData%AlphaLang(iDoF)**2) * Temperature * EvolData%NBeads )
         END DO
      END IF

      ! Themostat data is now setup
      EvolData%HasThermostat = .TRUE.

#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A,1F8.3,A,1F8.3)") "Thermostat is setup with Gamma = ",Gamma," and Temperature = ", Temperature
      WRITE(*,*) " Langevin DoFs: ", EvolData%ThermoSwitch(:)
#endif
      
   END SUBROUTINE SetupThermostat
   

!*******************************************************************************
!                          SetupRingPolymer
!*******************************************************************************
!> Setup data for ring polymer progation, store them in the Evolution datatype.
!>
!> @param EvolData     Evolution data type to setup
!> @param NBeads       Nr of replicas of the system for the RPMD
!*******************************************************************************
   SUBROUTINE SetupRingPolymer( EvolData, NBeads, PolymerFreq )
      IMPLICIT NONE
      TYPE( Evolution ), INTENT(INOUT)    :: EvolData
      INTEGER, INTENT(IN)                 :: NBeads
      REAL, INTENT(IN)                    :: PolymerFreq
      INTEGER :: i

      ! error if trying to setup thermostat of a non-setup evolution type
      CALL ERROR( .NOT. EvolData%IsSetup, "ClassicalEqMotion.SetupRingPolymer: evolution data not setup" )
      ! warn user if overwriting previously setup data
      CALL WARN( EvolData%HasRingPolymer, "ClassicalEqMotion.SetupRingPolymer: overwriting ring polymer data" ) 

      ! Store nr of beads and ring polymer harmonic frequency
      EvolData%NBeads = NBeads
      EvolData%RPFreq = PolymerFreq

      ! Setup FFT to compute normal modes of the free ring polymer
      CALL SetupFFT( EvolData%RingNormalModes, EvolData%NBeads ) 

      ! Set frequencies of the normal modes
      ALLOCATE( EvolData%NormModesFreq( EvolData%NBeads ) )
      DO i = 1, EvolData%NBeads
         EvolData%NormModesFreq(i) = 2.0 * EvolData%RPFreq * SIN( MyConsts_PI * real(i-1) / real(EvolData%NBeads) )
      END DO

      ! Set free propagator of the normal modes
      ALLOCATE( EvolData%NormModesPropag( 4,EvolData%NBeads ) )
      EvolData%NormModesPropag(:,1) = (/ 1.0, EvolData%dt, 0.0, 1.0 /)  ! Free evolution of the w=0 normal mode 
      DO i = 2, EvolData%NBeads                                         ! Harmonic oscillator evolution of other modes
         EvolData%NormModesPropag(1,i) = COS( EvolData%NormModesFreq(i)*EvolData%dt )
         EvolData%NormModesPropag(2,i) = SIN( EvolData%NormModesFreq(i)*EvolData%dt ) / EvolData%NormModesFreq(i)
         EvolData%NormModesPropag(3,i) = - SIN( EvolData%NormModesFreq(i)*EvolData%dt ) * EvolData%NormModesFreq(i)
         EvolData%NormModesPropag(4,i) = COS( EvolData%NormModesFreq(i)*EvolData%dt )
      END DO

      ! RPMD data is now setup
      EvolData%HasRingPolymer = .TRUE.

#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A,1I8,A)") "RPMD is setup with ",EvolData%NBeads," replicas "
#endif
      
   END SUBROUTINE SetupRingPolymer



!================================================================================================================================
!                              DISPOSE SUBROUTINES
!================================================================================================================================


!*******************************************************************************
!                          DisposeEvolutionData
!*******************************************************************************
!> Dispose evolution data.
!>
!> @param EvolData     Evolution data type  to dispose
!*******************************************************************************
   SUBROUTINE DisposeEvolutionData( EvolData )
      IMPLICIT NONE
      TYPE( Evolution ), INTENT(INOUT)  :: EvolData
      INTEGER :: iDoF

      ! continue if trying to dispose data that is not setup
      IF (.NOT. EvolData%IsSetup)  RETURN

      ! dispose thermostat and ringpolymer if it is the case
      IF ( EvolData%HasThermostat )  CALL DisposeThermostat( EvolData )
      IF ( EvolData%HasRingPolymer ) CALL DisposeRingPolymer( EvolData )

      ! Deallocate standard deviations of the thermal noise
      DEALLOCATE( EvolData%Mass )

      ! Themostat data is now disposed
      EvolData%HasThermostat = .FALSE.

#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A)") "Themostat has been disposed"
#endif
   END SUBROUTINE DisposeEvolutionData


!*******************************************************************************
!                          DisposeThermostat
!*******************************************************************************
!> Dispose thermostat data.
!>
!> @param EvolData     Evolution data type with the thermostat to dispose
!*******************************************************************************
   SUBROUTINE DisposeThermostat( EvolData )
      IMPLICIT NONE
      TYPE( Evolution ), INTENT(INOUT)  :: EvolData
      INTEGER :: iDoF
      
      ! continue if trying to dispose thermostat that is not setup
      IF ((.NOT. EvolData%IsSetup) .OR. (.NOT. EvolData%HasThermostat))  RETURN

      ! Set coefficient with zero friction
      EvolData%Gamma = 0.0
      
      ! Deallocate standard deviations of the thermal noise
      DEALLOCATE( EvolData%ThermalNoise, EvolData%ThermoSwitch )
      
      ! Themostat data is now disposed
      EvolData%HasThermostat = .FALSE.

#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A)") "Themostat has been disposed"
#endif
   END SUBROUTINE DisposeThermostat


!*******************************************************************************
!                          DisposeRingPolymer
!*******************************************************************************
!> Dispose ring polymer data.
!>
!> @param EvolData     Evolution data type with the thermostat to dispose
!*******************************************************************************
   SUBROUTINE DisposeRingPolymer( EvolData )
      IMPLICIT NONE
      TYPE( Evolution ), INTENT(INOUT)  :: EvolData
      INTEGER :: iDoF
      
      ! continue if trying to dispose thermostat that is not setup
      IF ((.NOT. EvolData%IsSetup) .OR. (.NOT. EvolData%HasRingPolymer))  RETURN

      ! Dispose data for FFT
      CALL DisposeFFT( EvolData%RingNormalModes )
   
      ! Ring polymer data is now disposed
      EvolData%HasRingPolymer = .FALSE.

#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A)") "Ring polymer has been disposed"
#endif
   END SUBROUTINE DisposeRingPolymer


!================================================================================================================================
!                         PROPAGATION SUBROUTINES
!================================================================================================================================


!*******************************************************************************
!                          EOM_VelocityVerlet
!*******************************************************************************
!> Propagate trajectory with Velocity-Verlet algorith.
!> If the Langevin parameters are setup, propagation is done in the
!> canonical ensamble with a Langevin thermostat
!> NOTE: THIS INTEGRATOR IS BETTER SUITED FOR MICROCANONICAL DYNAMICS 
!>       in case of Langevin MD, use Beeman's algorithm!
!> @ref http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
!>
!> @param EvolData     Evolution data type
!*******************************************************************************
   SUBROUTINE EOM_VelocityVerlet( EvolData, Pos, Vel, Acc, GetPotential, V, RandomNr )
      IMPLICIT NONE

      TYPE( Evolution ), INTENT(INOUT)                 :: EvolData
      REAL, DIMENSION( EvolData%NDoF ), INTENT(INOUT)  :: Pos, Vel, Acc
      REAL, INTENT(OUT)                                :: V
      TYPE(RNGInternalState), INTENT(INOUT) :: RandomNr

      INTERFACE
         REAL FUNCTION GetPotential( X, Force )
            REAL, DIMENSION(:), TARGET, INTENT(IN)  :: X
            REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Force
         END FUNCTION GetPotential
      END INTERFACE
      
      INTEGER :: iDoF

      ! (1) FULL TIME STEP FOR THE POSITIONS
      Pos(:) = Pos(:) + Vel(:)*EvolData%dt + 0.5*Acc(:)*(EvolData%dt**2)
 
      ! (2) HALF TIME STEP FOR THE VELOCITIES
      Vel(:) = Vel(:) + 0.5*Acc(:)*EvolData%dt

      ! (3) NEW FORCES AND ACCELERATIONS 
      V = GetPotential( Pos, Acc )       ! Compute new forces and store the potential value
      Acc(:) = Acc(:)/EvolData%Mass(:)   ! only potential forces

      ! (4) HALF TIME STEP AGAIN FOR THE VELOCITIES
      Vel(:) = Vel(:) + 0.5*Acc(:)*EvolData%dt

   END SUBROUTINE EOM_VelocityVerlet   
   


! !*******************************************************************************
! !> Propagate trajectory with Beeman's algorith.
! !> If the Langevin parameters are setup, propagation is done in the
! !> canonical ensamble with a Langevin thermostat
! !> NOTE: THIS INTEGRATOR IS BETTER SUITED FOR LANGEVIN DYNAMICS 
! !>       in case of microcanonical MD, use Velocity-Verlet!
! !> @ref http://en.wikipedia.org/wiki/Beeman%27s_algorithm
! !>
! !> @param EvolData     Evolution data type
! !*******************************************************************************
!    SUBROUTINE EOM_Beeman( EvolData, Pos, Vel, Acc, PreAcc, GetPotential, V, RandomNr )
!       IMPLICIT NONE
! 
!       TYPE( Evolution ), INTENT(INOUT)                 :: EvolData
!       REAL, DIMENSION( EvolData%NDoF ), INTENT(INOUT)  :: Pos, Vel, Acc, PreAcc
!       REAL, INTENT(OUT)                                :: V
!       TYPE(RNGInternalState), INTENT(INOUT) :: RandomNr
! 
!       INTERFACE
!          REAL FUNCTION GetPotential( X, Force )
!             REAL, DIMENSION(:), TARGET, INTENT(IN)  :: X
!             REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Force
!          END FUNCTION GetPotential
!       END INTERFACE
!       
!       INTEGER :: iDoF
!       ! Temporary array for predicted velocity and new accelerations
!       REAL, DIMENSION( EvolData%NDoF ) :: NewPos, NewVel, NewAcc
!  
!       ! (1) PREDICTED POSITIONS
!       NewPos(:) = Pos(:) + Vel(:)*EvolData%dt + (4.*Acc(:)-PreAcc(:))*(EvolData%dt**2)/6.0
! 
!       ! (2) PREDICTED VELOCITY
!       NewVel(:) = Vel(:) + (3.*Acc(:)-PreAcc(:))*EvolData%dt/2.0
!  
!       IF ( .NOT. EvolData%HasThermostat ) THEN        ! Integration without Langevin thermostat
! 
!          ! (2) NEW ACCELERATION
!          V = GetPotential( NewPos, NewAcc )         ! Compute new forces and store the potential value
!          NewAcc(:) =  NewAcc(:) / EvolData%Mass(:)   ! Devide by the mass
! 
!       ELSE IF ( ( EvolData%HasThermostat ) ) THEN
!             
!          ! (3) NEW ACCELERATION
!          V = GetPotential( NewPos, NewAcc )         ! Compute new forces and store the potential value
! 
!          DO iDoF = 1, EvolData%NDoF
!             IF ( EvolData%ThermoSwitch(iDoF) ) THEN
!                NewAcc(iDoF) = ( NewAcc(iDoF) + GaussianRandomNr(RandomNr)*EvolData%ThermalNoise(iDoF) ) &
!                                                         / EvolData%Mass(iDoF) - EvolData%Gamma*NewVel(iDoF)
!             ELSE IF ( .NOT. EvolData%ThermoSwitch(iDoF) ) THEN
!                NewAcc(iDoF) = NewAcc(iDoF)  / EvolData%Mass(iDoF)
!             END IF
!          END DO
! 
!       END IF
! 
!       ! CORRECTED POSITIONS
!       Pos(:) = Pos(:) + Vel(:)*EvolData%dt + (NewAcc(:)+2*Acc(:))*(EvolData%dt**2)/6.0
! ! 
!       ! (4) CORRECTED VELOCITIES
!       Vel(:) = Vel(:) + (Acc(:) + NewAcc(:))*EvolData%dt/2.0    
! 
!       ! Store new acceleration
!       PreAcc(:) = Acc(:)
!       Acc(:) = NewAcc(:)
! 
!    END SUBROUTINE EOM_Beeman   


!*******************************************************************************
!>                   EOM_LangevinSecondOrder
!*******************************************************************************
!> Propagate trajectory with Vanden-Eijnden and Ciccoti algorith, when
!> the Langevin parameters are setup. When microcanonical propagation 
!> is assumed, the algorithm is equivalent to Velocity-Verlet
!> THIS INTEGRATOR IS A GENERAL SYMPLECTIC PROPAGATOR FOR
!> MICROCANONICAL MD AND CANONICAL LANGEVIN MD
!> @ref http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
!> @ref M.E. Tuckerman "Statistical Mechanics: Theory and Molecular Simulation"
!>
!> @param EvolData      Evolution data type
!> @param Pos           In/Out coordinates vector
!> @param Vel           In/Out velocities vector
!> @param Acc           In/Out acceleration vector (need to be properly computed as input)
!> @param GetPotential  Function to evaluete potential and forces
!> @param RandomNr      Internal state of the random number generator
!*******************************************************************************
   SUBROUTINE EOM_LangevinSecondOrder( EvolData, Pos, Vel, Acc, GetPotential, V, RandomNr )
      IMPLICIT NONE

      TYPE( Evolution ), INTENT(INOUT)                 :: EvolData
      REAL, DIMENSION( EvolData%NDoF ), INTENT(INOUT)  :: Pos, Vel, Acc
      REAL, INTENT(OUT)                                :: V
      TYPE(RNGInternalState), INTENT(INOUT)            :: RandomNr

      INTERFACE
         REAL FUNCTION GetPotential( X, Force )
            REAL, DIMENSION(:), TARGET, INTENT(IN)  :: X
            REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Force
         END FUNCTION GetPotential
      END INTERFACE

      INTEGER :: iDoF
      REAL, DIMENSION( EvolData%NDoF ) :: A
      REAL, DIMENSION( EvolData%NDoF ) :: Xi, Eta 

      IF ( .NOT. EvolData%HasThermostat ) THEN        ! Integration without Langevin thermostat

         ! (1) FULL TIME STEP FOR THE POSITIONS
         Pos(:) = Pos(:) + Vel(:)*EvolData%dt + 0.5*Acc(:)*(EvolData%dt**2)
   
         ! (2) HALF TIME STEP FOR THE VELOCITIES
         Vel(:) = Vel(:) + 0.5*Acc(:)*EvolData%dt

         ! (3) NEW FORCES AND ACCELERATIONS 
         V = GetPotential( Pos, Acc )       ! Compute new forces and store the potential value
         Acc(:) = Acc(:)/EvolData%Mass(:)   ! only potential forces

         ! (4) HALF TIME STEP AGAIN FOR THE VELOCITIES
         Vel(:) = Vel(:) + 0.5*Acc(:)*EvolData%dt

      ELSE IF ( ( EvolData%HasThermostat ) ) THEN     ! Integration with Langevin thermostat

         ! (0) COMPUTE NECESSARY RANDOM VALUES
         DO iDoF = 1, EvolData%NDoF
            IF ( EvolData%ThermoSwitch(iDoF) ) THEN
               Xi(iDoF)  = GaussianRandomNr(RandomNr)
               Eta(iDoF) = GaussianRandomNr(RandomNr)
               A(iDoF) = 0.5 * EvolData%dt**2 * ( Acc(iDoF) - EvolData%Gamma*Vel(iDoF) ) + &
                         EvolData%ThermalNoise(iDoF) * EvolData%dt**(1.5) * ( 0.5 * Xi(iDoF) + Over2Sqrt3 * Eta(iDoF)  )
            ELSE IF ( .NOT. EvolData%ThermoSwitch(iDoF) ) THEN
                A(iDoF) = 0.5*Acc(iDoF)*(EvolData%dt**2)
            END IF
         END DO

         ! (1) UPDATE POSITION
         Pos(:) = Pos(:) + Vel(:)*EvolData%dt + A(:)

         ! (2) PARTIAL UPDATE OF THE VELOCITIES
         DO iDoF = 1, EvolData%NDoF
            IF ( EvolData%ThermoSwitch(iDoF) ) THEN
               Vel(iDoF) = (1.0 - EvolData%Gamma*EvolData%dt) * Vel(iDoF) + 0.5*Acc(iDoF)*EvolData%dt
            ELSE IF ( .NOT. EvolData%ThermoSwitch(iDoF) ) THEN
               Vel(iDoF) = Vel(iDoF) + 0.5*Acc(iDoF)*EvolData%dt
            END IF
         END DO

         ! (3) NEW FORCES AND ACCELERATIONS 
         V = GetPotential( Pos, Acc )       ! Compute new forces and store the potential value
         Acc(:) = Acc(:)  / EvolData%Mass(:)

         ! (4) FINAL UPDATE OF THE VELOCITIES
         DO iDoF = 1, EvolData%NDoF
            IF ( EvolData%ThermoSwitch(iDoF) ) THEN
               Vel(iDoF) = Vel(iDoF) + 0.5*Acc(iDoF)*EvolData%dt + SQRT(EvolData%dt) * EvolData%ThermalNoise(iDoF) * Xi(iDoF) &
                                    - EvolData%Gamma * A(iDoF)
            ELSE IF ( .NOT. EvolData%ThermoSwitch(iDoF) ) THEN
               Vel(iDoF) = Vel(iDoF) + 0.5*Acc(iDoF)*EvolData%dt
            END IF
         END DO
            
      END IF

   END SUBROUTINE EOM_LangevinSecondOrder   


!*******************************************************************************
!>                   EOM_RPMSymplectic
!*******************************************************************************
!> Propagate trajectory with symplectic algorithm for Ring-Polymer MD. 
!> SpecialOptions : 1 - initialize acceleration and do not propagate
!>                  2 - equilibrate only the internal mode of the ring polymer 
!>
!> @param EvolData        Evolution data type
!> @param Pos             In/Out coordinates vector
!> @param Vel             In/Out velocities vector
!> @param Acc             In/Out acceleration vector (need to be properly computed as input)
!> @param GetPotential    Function to evaluete potential and forces
!> @param RandomNr        Internal state of the random number generator
!> @param SpecialOptions  Define special behaviour of the subroutine
!*******************************************************************************
   SUBROUTINE EOM_RPMSymplectic( EvolData, Pos, Vel, Acc, GetPotential, V, RandomNr, SpecialOptions )
      IMPLICIT NONE
      TYPE( Evolution ), INTENT(INOUT)                                   :: EvolData
      REAL, DIMENSION( EvolData%NDoF * EvolData%NBeads ), INTENT(INOUT)  :: Pos, Vel, Acc
      REAL, INTENT(OUT)                                                  :: V
      TYPE(RNGInternalState), INTENT(INOUT)                              :: RandomNr
      INTEGER, INTENT(IN), OPTIONAL                                      :: SpecialOptions
   
      INTERFACE
         REAL FUNCTION GetPotential( X, Force )
            REAL, DIMENSION(:), TARGET, INTENT(IN)  :: X
            REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Force
         END FUNCTION GetPotential
      END INTERFACE

      REAL, DIMENSION(EvolData%NBeads)  :: BeadQAt0, BeadVAt0, BeadQAtT, BeadVAtT
      INTEGER :: iDoF, iBead, iStart, iEnd, PropagOption
      REAL    :: VBead
      REAL, DIMENSION(EvolData%NBeads, EvolData%NDoF) :: StoreAccel

      CALL ERROR( .NOT. EvolData%HasRingPolymer, " EOM_RPMSymplectic: data for RPMD are needed "  ) 

      ! the subroutine can be used for special kind of propagations
      IF ( PRESENT(SpecialOptions) ) THEN
         PropagOption = SpecialOptions      ! in such case, a special option is enabled
      ELSE
         PropagOption = 0                   ! 0 means normal propagation
      END IF

      IF ( PropagOption /= 1 ) THEN

      ! (1)  HALF TIME STEP FOR THERMOSTATTING THE SYSTEM (only when langevin dynamics)
         IF ( EvolData%HasThermostat ) THEN
            DO iDoF = 1, EvolData%NDoF

               DO iBead = 1, EvolData%NBeads                        ! extract single bead velocities
                  BeadVAt0(iBead) = Vel( (iBead-1) * EvolData%NDoF + iDoF )
               END DO
               CALL ExecuteFFT( EvolData%RingNormalModes, BeadVAt0, DIRECT_FFT ) ! transform to normal modes velocities

               IF ( PropagOption == 2 ) THEN
                  BeadVAtT(1) = BeadVAt0(1)
                  DO iBead = 2, EvolData%NBeads                                   ! evolve normal modes
                     BeadVAtT(iBead) = BeadVAt0(iBead) * EvolData%AlphaLang(iBead) + &
                                       GaussianRandomNr( RandomNr ) * EvolData%BetaLang(iBead) / SQRT( EvolData%Mass(iDoF) )
                  END DO
               ELSE
                  DO iBead = 1, EvolData%NBeads                                   ! evolve normal modes
                     BeadVAtT(iBead) = BeadVAt0(iBead) * EvolData%AlphaLang(iBead) + &
                                       GaussianRandomNr( RandomNr ) * EvolData%BetaLang(iBead) / SQRT( EvolData%Mass(iDoF) )
                  END DO
               END IF

               CALL ExecuteFFT( EvolData%RingNormalModes, BeadVAtT, INVERSE_FFT ) ! transform back to original coords
               DO iBead = 1, EvolData%NBeads               ! copy single bead velocities to input arrays
                  Vel( (iBead-1) * EvolData%NDoF + iDoF ) = BeadVAtT(iBead)
               END DO

            END DO
         END IF

      ! (2) HALF TIME STEP FOR THE VELOCITIES
         DO iBead = 1, EvolData%NBeads
            iStart = (iBead-1) * EvolData%NDoF + 1
            iEnd   = iBead * EvolData%NDoF
            Vel( iStart:iEnd ) = Vel( iStart:iEnd ) + 0.5*Acc( iStart:iEnd )*EvolData%dt
         END DO
      END IF

      ! (3) EXACT PROPAGATION WITH THE INTERBEADS POTENTIAL
      DO iDoF = 1, EvolData%NDoF

         DO iBead = 1, EvolData%NBeads      ! extract single bead positions and velocities
            BeadQAt0(iBead) = Pos( (iBead-1) * EvolData%NDoF + iDoF )
            BeadVAt0(iBead) = Vel( (iBead-1) * EvolData%NDoF + iDoF )
         END DO

         CALL ExecuteFFT( EvolData%RingNormalModes, BeadQAt0, DIRECT_FFT ) ! transform to normal modes coords
         CALL ExecuteFFT( EvolData%RingNormalModes, BeadVAt0, DIRECT_FFT )

         IF ( PropagOption /= 1 ) THEN
            IF ( PropagOption == 2 ) THEN
            BeadQAtT(1) = BeadQAt0(1)
            BeadVAtT(1) = BeadVAt0(1)
            DO iBead = 2, EvolData%NBeads                                   ! evolve normal modes
               BeadQAtT(iBead) = EvolData%NormModesPropag(1,iBead) * BeadQAt0(iBead) + &
                                                                          EvolData%NormModesPropag(2,iBead) * BeadVAt0(iBead)
               BeadVAtT(iBead) = EvolData%NormModesPropag(3,iBead) * BeadQAt0(iBead) + &
                                                                          EvolData%NormModesPropag(4,iBead) * BeadVAt0(iBead)
            END DO
            ELSE
            DO iBead = 1, EvolData%NBeads                                   ! evolve normal modes
               BeadQAtT(iBead) = EvolData%NormModesPropag(1,iBead) * BeadQAt0(iBead) + &
                                                                          EvolData%NormModesPropag(2,iBead) * BeadVAt0(iBead)
               BeadVAtT(iBead) = EvolData%NormModesPropag(3,iBead) * BeadQAt0(iBead) + &
                                                                          EvolData%NormModesPropag(4,iBead) * BeadVAt0(iBead)
            END DO
            END IF
         ELSE
            BeadQAtT(:) = BeadQAt0(:)
            BeadVAtT(:) = BeadVAt0(:)
         END IF

         CALL ExecuteFFT( EvolData%RingNormalModes, BeadQAtT, INVERSE_FFT ) ! transform back to original coords
         CALL ExecuteFFT( EvolData%RingNormalModes, BeadVAtT, INVERSE_FFT )

         DO iBead = 1, EvolData%NBeads      ! copy single bead positions and velocities to input arrays
            Pos( (iBead-1) * EvolData%NDoF + iDoF ) = BeadQAtT(iBead)
            Vel( (iBead-1) * EvolData%NDoF + iDoF ) = BeadVAtT(iBead)
         END DO

      END DO

      ! (4) UPDATE FORCES
      V = 0.0
      DO iBead = 1, EvolData%NBeads
         iStart = (iBead-1) * EvolData%NDoF + 1
         iEnd   = iBead * EvolData%NDoF
         VBead = GetPotential( Pos( iStart:iEnd ), Acc( iStart:iEnd ) )   ! Compute new forces and store the potential value
         Acc( iStart:iEnd ) = Acc( iStart:iEnd ) / EvolData%Mass(:)       ! only potential forces 
         V = V + VBead
      END DO
      V = V / real(EvolData%NBeads)
      
      IF ( PropagOption == 2 ) THEN
         DO iBead = 1, EvolData%NBeads
            StoreAccel(iBead, :) = Acc( (iBead-1) * EvolData%NDoF + 1 : iBead * EvolData%NDoF )
         END DO
         DO iDoF = 1, EvolData%NDoF
            CALL ExecuteFFT( EvolData%RingNormalModes, StoreAccel(:, iDoF), DIRECT_FFT )
            StoreAccel(1, iDoF) = 0.0
            CALL ExecuteFFT( EvolData%RingNormalModes, StoreAccel(:, iDoF), INVERSE_FFT )
         END DO
         DO iBead = 1, EvolData%NBeads
            Acc( (iBead-1) * EvolData%NDoF + 1 : iBead * EvolData%NDoF ) = StoreAccel(iBead, :) 
         END DO
      END IF

      IF ( PropagOption /= 1 ) THEN
      ! (5) HALF TIME STEP FOR THE VELOCITIES
         DO iBead = 1, EvolData%NBeads
            iStart = (iBead-1) * EvolData%NDoF + 1
            iEnd   = iBead * EvolData%NDoF
            Vel( iStart:iEnd ) = Vel( iStart:iEnd ) + 0.5*Acc( iStart:iEnd )*EvolData%dt
         END DO

      ! (6)  HALF TIME STEP FOR THERMOSTATTING THE SYSTEM (only when langevin dynamics)
         IF ( EvolData%HasThermostat ) THEN
            DO iDoF = 1, EvolData%NDoF

               DO iBead = 1, EvolData%NBeads                        ! extract single bead velocities
                  BeadVAt0(iBead) = Vel( (iBead-1) * EvolData%NDoF + iDoF )
               END DO
               CALL ExecuteFFT( EvolData%RingNormalModes, BeadVAt0, DIRECT_FFT ) ! transform to normal modes velocities

               IF ( PropagOption == 2 ) THEN
                  BeadVAtT(1) = BeadVAt0(1)
                  DO iBead = 2, EvolData%NBeads                                   ! evolve normal modes
                     BeadVAtT(iBead) = BeadVAt0(iBead) * EvolData%AlphaLang(iBead) + &
                                       GaussianRandomNr(RandomNr) * EvolData%BetaLang(iBead) / SQRT( EvolData%Mass(iDoF) )
                  END DO
               ELSE
                  DO iBead = 1, EvolData%NBeads                                   ! evolve normal modes
                     BeadVAtT(iBead) = BeadVAt0(iBead) * EvolData%AlphaLang(iBead) + &
                                       GaussianRandomNr(RandomNr) * EvolData%BetaLang(iBead) / SQRT( EvolData%Mass(iDoF) )
                  END DO
               END IF

               CALL ExecuteFFT( EvolData%RingNormalModes, BeadVAtT, INVERSE_FFT ) ! transform back to original coords
               DO iBead = 1, EvolData%NBeads               ! copy single bead velocities to input arrays
                  Vel( (iBead-1) * EvolData%NDoF + iDoF ) = BeadVAtT(iBead)
               END DO

            END DO
         END IF

      END IF

   END SUBROUTINE EOM_RPMSymplectic


!================================================================================================================================
!                              OTHER SUBROUTINES
!================================================================================================================================
   
!*******************************************************************************
!> Compute kinetic energy corresponding to a given velocity vector.
!>
!> @param EvolData     Evolution data type
!> @param Velocity     Array containing the velocity at given time step
!*******************************************************************************
   REAL FUNCTION EOM_KineticEnergy( EvolData, Vel, NMax ) RESULT( KinEnergy )
      IMPLICIT NONE
      TYPE( Evolution ), INTENT(INOUT)    :: EvolData
      REAL, DIMENSION(:), INTENT(INOUT)   :: Vel
      INTEGER, INTENT(IN), OPTIONAL       :: NMax
      INTEGER :: iDoF, iBead, N, NDoF
      
      IF (.NOT. PRESENT( NMax )) THEN
         NDoF = EvolData%NDoF
      ELSE 
         NDoF = MIN(NMax,EvolData%NDoF)
      END IF

      KinEnergy = 0.0
      IF ( EvolData%HasRingPolymer ) THEN
         N = 0
         DO iBead = 1, EvolData%NBeads
            DO iDoF = 1, NDoF
               N = N + 1
               KinEnergy = KinEnergy + 0.5 * EvolData%Mass(iDoF) * Vel(N)**2
            END DO
         END DO

      ELSE IF ( .NOT. EvolData%HasRingPolymer ) THEN
         DO iDoF = 1, NDoF
            KinEnergy = KinEnergy + 0.5 * EvolData%Mass(iDoF) * Vel(iDoF)**2
         END DO

      ENDIF

   END FUNCTION EOM_KineticEnergy


!================================================================================================================================
!                              END OF MODULE
!================================================================================================================================

END MODULE ClassicalEqMotion
