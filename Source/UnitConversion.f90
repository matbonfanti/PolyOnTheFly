!***************************************************************************************
!*                              MODULE UnitConversion
!***************************************************************************************
!
!>  \brief         Unit conversion module.
!>         
!>  \details       Define datatype to store information about units and 
!>                 gives conversion factors from given units dataset. Each unit dataset
!>                 consists of length, energy, mass, angle, time, temperature
!>                 and frequency units. 
!
!***************************************************************************************
!
!>  \author        Mark Wijzenbroek, Matteo Bonfanti
!>  \version       2.0
!>  \date          7 May 2013
!>
!***************************************************************************************
!
!>  \par          7 May 2013, major revision \n
!>                the subroutines now do not convert the numbers but just give
!>                conversion factors. the errors are now handled through ErrorTrap.
!>  \par          27 November 2013, revision \n
!>                included new units: time, temperature, frequency \n
!>                new functions to print unit label for given unit in as set
!
!***************************************************************************************

MODULE UnitConversion
#include "preprocessoptions.cpp"

   PRIVATE
   PUBLIC :: Units
   PUBLIC :: Initialize_UnitConversion, Print_Unit_Definition
   PUBLIC :: LengthUnit, EnergyUnit, MassUnit, AngleUnit, TimeUnit, TemperUnit, FreqUnit
   PUBLIC :: LengthConversion, EnergyConversion, MassConversion, AngleConversion, ForceConversion
   PUBLIC :: TimeConversion, TemperatureConversion, FreqConversion


! **************** Length units ********************

   INTEGER, PARAMETER, PUBLIC :: UNITS_BOHR     = 0 !< Index for bohr
   INTEGER, PARAMETER, PUBLIC :: UNITS_ANGSTROM = 1 !< Index for Angstrom

   REAL, PARAMETER :: CONV_ANGSTROM = MyConsts_Bohr2Ang  !< Conversion factor from bohr to Angstrom

! **************** Energy units ********************

   INTEGER, PARAMETER, PUBLIC :: UNITS_HARTREE = 2 !< Index for Hartree
   INTEGER, PARAMETER, PUBLIC :: UNITS_EV      = 3 !< Index for eV
   INTEGER, PARAMETER, PUBLIC :: UNITS_KCALMOL = 4 !< Index for kcal/mol
   INTEGER, PARAMETER, PUBLIC :: UNITS_KJMOL   = 5 !< Index for kJ/mol
   INTEGER, PARAMETER, PUBLIC :: UNITS_JOULE   = 6 !< Index for Joule

   REAL, PARAMETER :: CONV_EV = MyConsts_Hartree2eV            !< Conversion factor from Hartree to eV
   REAL, PARAMETER :: CONV_KCALMOL = MyConsts_Hatree2kcalmol   !< Conversion factor from Hartree to kcal/mol
   REAL, PARAMETER :: CONV_KJMOL = MyConsts_Hatree2kjmol       !< Conversion factor from Hartree to kJ/mol
   REAL, PARAMETER :: CONV_JOULE = MyConsts_Hatree2joule       !< Conversion factor from Hartree to Joule

! **************** Mass units ********************

   INTEGER, PARAMETER, PUBLIC :: UNITS_MELEC = 7 !< Index for electron mass
   INTEGER, PARAMETER, PUBLIC :: UNITS_AMU = 8   !< Index for amu
   INTEGER, PARAMETER, PUBLIC :: UNITS_KG = 9    !< Index for kg

   REAL, PARAMETER :: CONV_AMU = 1/MyConsts_Uma2Au  !< Conversion factor from electron mass to amu
   REAL, PARAMETER :: CONV_KG = MyConsts_mel        !< Conversion factor from electron mass to kg

! **************** Angular units ********************

   INTEGER, PARAMETER, PUBLIC :: UNITS_RAD    = 10 !< Index for radiant
   INTEGER, PARAMETER, PUBLIC :: UNITS_DEGREE = 11 !< Index for sexagesimal degree
 
   REAL, PARAMETER :: CONV_DEGREE = 1./MyConsts_Deg2Rad  !< Conversion factor from radiant to degree

! **************** Time units **********************

   INTEGER, PARAMETER, PUBLIC :: UNITS_AUTIME = 12 !< Index for atomic unit of time
   INTEGER, PARAMETER, PUBLIC :: UNITS_FEMTOS = 13 !< Index for femtosecond
   INTEGER, PARAMETER, PUBLIC :: UNITS_PICOS  = 14 !< Index for picosecond
 
   REAL, PARAMETER :: CONV_FEMTOS = MyConsts_AUofTime        !< Conversion factor from au to femtoseconds
   REAL, PARAMETER :: CONV_PICOS  = 0.001*MyConsts_AUofTime  !< Conversion factor from au to picoseconds

! **************** Temperature **********************

   INTEGER, PARAMETER, PUBLIC :: UNITS_AUTEMP = 15 !< Index for atomic unit of temperature (hartree via boltzmann constant)
   INTEGER, PARAMETER, PUBLIC :: UNITS_KELVIN = 16 !< Index for kelvin
 
   REAL, PARAMETER :: CONV_KELVIN = 1.0/MyConsts_K2AU        !< Conversion factor from au to kelvin

! **************** Frequency ************************

   INTEGER, PARAMETER, PUBLIC :: UNITS_AUFREQ   = 17  !< Index for au of frequency (hartree, since hbar = 1 )
   INTEGER, PARAMETER, PUBLIC :: UNITS_CMMINUS1 = 18  !< Index for cm^-1
   INTEGER, PARAMETER, PUBLIC :: UNITS_FSMINUS1 = 19  !< Index for fs^-1

   REAL, PARAMETER :: CONV_CMMINUS1 = 1.0/MyConsts_cmmin1toAU                          !< Conversion factor from au to cm^-1
   REAL, PARAMETER :: CONV_FSMINUS1 = 1.0/MyConsts_cmmin1toAU*MyConsts_cmmin1tofsmin1  !< Conversion factor from au to fs^-1

! **************** Units labels ********************

   CHARACTER(8), DIMENSION(0:19), PARAMETER :: UNITS_LABELS = &
      (/ "bohr    ", "Ang     ",                                     &    ! labels for the distance
         "Eh      ", "eV      ", "kcal/mol", "kJ/mol  ", "J       ", &    ! labels for the energy
         "m_el    ", "AMU     ", "kg      ",                         &    ! labels for the mass
         "rad     ", "deg     ",                                     &    ! labels for angles        
         "au      ", "fs      ", "ps      ",                         &    ! labels for time        
         "Eh      ", "K       ",                                     &    ! labels for temperature
         "Eh      ", "cm^-1   ", "fs^-1   "                          /)   ! labels for frequency

!>  A set of units. By default, intialized with atomic units and radiant
   TYPE Units
      INTEGER :: LengthUnit  = UNITS_BOHR
      INTEGER :: EnergyUnit  = UNITS_HARTREE
      INTEGER :: MassUnit    = UNITS_MELEC
      INTEGER :: AngleUnit   = UNITS_RAD
      INTEGER :: TimeUnit    = UNITS_AUTIME
      INTEGER :: TemperUnit  = UNITS_AUTEMP
      INTEGER :: FreqUnit    = UNITS_AUFREQ
   END TYPE

!>  Units internally used in the program (set to atomic units).
   TYPE ( Units ), SAVE, PUBLIC :: InternalUnits
  
!>  Units that are used in reading input files, and to write output files.
   TYPE ( Units ), SAVE, PUBLIC :: InputUnits

!> Units which are used in read and write, in the MAIN file only.
   TYPE ( Units ), SAVE, PUBLIC :: MainUnits


!********************************************************************************************************
CONTAINS
!********************************************************************************************************


!*******************************************************************************
! Setup_UnitConversion
!*******************************************************************************
!> Initialize the a set of units.
!*******************************************************************************
   SUBROUTINE Initialize_UnitConversion( ConvUnit, Length, Energy, Mass, Angle, Time, Temp, Freq )
      IMPLICIT NONE
      TYPE (Units)          :: ConvUnit
      INTEGER, INTENT(IN)   :: Length, Energy, Mass, Angle, Time, Temp, Freq

      ! check if length unit exists
      IF ( ( Length == UNITS_BOHR ) .OR. ( Length == UNITS_ANGSTROM ) ) THEN
            ConvUnit%LengthUnit = Length
      ELSE
            CALL AbortWithError( "Initialize_UnitConversion: Unknown length input unit" )
      END IF

      ! check if energy unit exists
      IF ( ( Energy == UNITS_HARTREE ) .OR. ( Energy == UNITS_EV ) .OR.    &
           ( Energy == UNITS_KCALMOL ) .OR. ( Energy == UNITS_KJMOL ) .OR. & 
           ( Energy == UNITS_JOULE )  )  THEN
            ConvUnit%EnergyUnit = Energy
      ELSE
            CALL AbortWithError( "Initialize_UnitConversion: Unknown energy input unit" )
      END IF

      ! check if mass unit exists
      IF ( ( Mass == UNITS_MELEC ) .OR. ( Mass == UNITS_AMU ) .OR. ( Mass == UNITS_KG ) ) THEN
            ConvUnit%MassUnit = Mass
      ELSE
            CALL AbortWithError( "Initialize_UnitConversion: Unknown mass input unit" )
      END IF

      ! check if angle unit exists
      IF ( ( Angle == UNITS_RAD ) .OR. ( Angle == UNITS_DEGREE ) ) THEN
            ConvUnit%AngleUnit = Angle
      ELSE
            CALL AbortWithError( "Initialize_UnitConversion: Unknown angle input unit" )
      END IF

      ! check if time unit exists
      IF ( ( Time == UNITS_AUTIME ) .OR. ( Time == UNITS_FEMTOS ) .OR. ( Time == UNITS_PICOS ) ) THEN
            ConvUnit%TimeUnit = Time
      ELSE
            CALL AbortWithError( "Initialize_UnitConversion: Unknown time input unit" )
      END IF

      ! check if temperature unit exists
      IF ( ( Temp == UNITS_AUTEMP ) .OR. ( Temp == UNITS_KELVIN ) ) THEN
            ConvUnit%TemperUnit = Temp
      ELSE
            CALL AbortWithError( "Initialize_UnitConversion: Unknown temperature input unit" )
      END IF

      ! check if frequency unit exists
      IF ( ( Freq == UNITS_AUFREQ ) .OR. ( Freq == UNITS_CMMINUS1 ) .OR. ( Freq == UNITS_FSMINUS1 ) ) THEN
            ConvUnit%FreqUnit = Freq
      ELSE
            CALL AbortWithError( "Initialize_UnitConversion: Unknown frequency input unit" )
      END IF

  END SUBROUTINE Initialize_UnitConversion


!********************************************************************************************************

!*******************************************************************************
!    Print_Unit_Definition
!*******************************************************************************
!> Print a string describing the units defined in a given Units type.
!>
!> @param Units     The Units that are described by the output string.
!> @returns         A string describing the units set.
!*******************************************************************************
   FUNCTION Print_Unit_Definition( UnitsSet ) RESULT ( String )
      IMPLICIT NONE
      TYPE (Units), INTENT(IN)     :: UnitsSet
      CHARACTER(154)   :: String

      WRITE(String, "(7(1X,A21))")  "Length:      "  // (UNITS_LABELS( UnitsSet%LengthUnit )), &
                                    "Energy:      "  // (UNITS_LABELS( UnitsSet%EnergyUnit )), &
                                    "Mass:        "  // (UNITS_LABELS( UnitsSet%MassUnit )),   &
                                    "Angle:       "  // (UNITS_LABELS( UnitsSet%AngleUnit )),  &
                                    "Time:        "  // (UNITS_LABELS( UnitsSet%TimeUnit )),   &
                                    "Temperature: "  // (UNITS_LABELS( UnitsSet%TemperUnit )),   &
                                    "Frequency:   "  // (UNITS_LABELS( UnitsSet%FreqUnit ))
   END FUNCTION Print_Unit_Definition

   FUNCTION LengthUnit( UnitsSet ) RESULT ( String )
      IMPLICIT NONE
      TYPE (Units), INTENT(IN)     :: UnitsSet
      CHARACTER(8)   :: String

      WRITE(String, "(A)")  UNITS_LABELS( UnitsSet%LengthUnit )
   END FUNCTION LengthUnit

   FUNCTION EnergyUnit( UnitsSet ) RESULT ( String )
      IMPLICIT NONE
      TYPE (Units), INTENT(IN)     :: UnitsSet
      CHARACTER(8)   :: String

      WRITE(String, "(A)")  UNITS_LABELS( UnitsSet%EnergyUnit )
   END FUNCTION EnergyUnit

   FUNCTION MassUnit( UnitsSet ) RESULT ( String )
      IMPLICIT NONE
      TYPE (Units), INTENT(IN)     :: UnitsSet
      CHARACTER(8)   :: String

      WRITE(String, "(A)")  UNITS_LABELS( UnitsSet%MassUnit )
   END FUNCTION MassUnit

   FUNCTION AngleUnit( UnitsSet ) RESULT ( String )
      IMPLICIT NONE
      TYPE (Units), INTENT(IN)     :: UnitsSet
      CHARACTER(8)   :: String

      WRITE(String, "(A)")  UNITS_LABELS( UnitsSet%AngleUnit )
   END FUNCTION AngleUnit

   FUNCTION TimeUnit( UnitsSet ) RESULT ( String )
      IMPLICIT NONE
      TYPE (Units), INTENT(IN)     :: UnitsSet
      CHARACTER(8)   :: String

      WRITE(String, "(A)")  UNITS_LABELS( UnitsSet%TimeUnit )
   END FUNCTION TimeUnit

   FUNCTION TemperUnit( UnitsSet ) RESULT ( String )
      IMPLICIT NONE
      TYPE (Units), INTENT(IN)     :: UnitsSet
      CHARACTER(8)   :: String

      WRITE(String, "(A)")  UNITS_LABELS( UnitsSet%TemperUnit )
   END FUNCTION TemperUnit

   FUNCTION FreqUnit( UnitsSet ) RESULT ( String )
      IMPLICIT NONE
      TYPE (Units), INTENT(IN)     :: UnitsSet
      CHARACTER(8)   :: String

      WRITE(String, "(A)")  UNITS_LABELS( UnitsSet%FreqUnit )
   END FUNCTION FreqUnit

!********************************************************************************************************


!*******************************************************************************
!      LengthConversion
!*******************************************************************************
!> Give length conversion factor from InUnits to OutUnits.
!>
!> @param InUnits    Units data type with the input units.
!> @param OutUnits   Units data type with the output units.
!> @returns          The conversion factor.
!*******************************************************************************
   FUNCTION LengthConversion( InUnits, OutUnits )  RESULT( Conversion )
      IMPLICIT NONE
      TYPE (Units), INTENT(IN)   :: InUnits, OutUnits
      REAL                       :: Conversion

      IF ( InUnits%LengthUnit == OutUnits%LengthUnit ) THEN
         Conversion = 1.0
         RETURN
      END IF

      SELECT CASE ( InUnits%LengthUnit )
         CASE ( UNITS_BOHR )
            Conversion = 1.0
         CASE ( UNITS_ANGSTROM )
            Conversion = 1.0 / CONV_ANGSTROM
         CASE DEFAULT
            CALL AbortWithError( "LengthConversion: Unknown input units" )
      END SELECT

      SELECT CASE ( OutUnits%LengthUnit )
         CASE ( UNITS_BOHR )
            Conversion = Conversion
         CASE ( UNITS_ANGSTROM )
            Conversion = Conversion * CONV_ANGSTROM
         CASE DEFAULT
            CALL AbortWithError( "LengthConversion: Unknown output units" )
      END SELECT

   END FUNCTION LengthConversion


!********************************************************************************************************


!*******************************************************************************
!      EnergyConversion
!*******************************************************************************
!> Give energy conversion factor from InUnits to OutUnits.
!>
!> @param InUnits    Units data type with the input units.
!> @param OutUnits   Units data type with the output units.
!> @returns          The conversion factor.
!*******************************************************************************
   FUNCTION EnergyConversion( InUnits, OutUnits )  RESULT( Conversion )
      IMPLICIT NONE
      TYPE (Units), INTENT(IN)   :: InUnits, OutUnits
      REAL                       :: Conversion

      IF ( InUnits%EnergyUnit == OutUnits%EnergyUnit ) THEN
         Conversion = 1.0
         RETURN
      END IF 

      SELECT CASE ( InUnits%EnergyUnit )
         CASE ( UNITS_HARTREE )
            Conversion = 1.0
         CASE ( UNITS_EV )
            Conversion = 1.0 / CONV_EV
         CASE ( UNITS_KCALMOL )
            Conversion = 1.0 / CONV_KCALMOL
         CASE ( UNITS_KJMOL )
            Conversion = 1.0 / CONV_KJMOL
         CASE ( UNITS_JOULE )
            Conversion = 1.0 / CONV_JOULE
         CASE DEFAULT
            CALL AbortWithError( "EnergyConversion: Unknown input units" )
      END SELECT

      SELECT CASE ( OutUnits%EnergyUnit )
         CASE ( UNITS_HARTREE )
            Conversion = Conversion
         CASE ( UNITS_EV )
            Conversion = Conversion * CONV_EV
         CASE ( UNITS_KCALMOL )
            Conversion = Conversion * CONV_KCALMOL
         CASE ( UNITS_KJMOL )
            Conversion = Conversion * CONV_KJMOL
         CASE ( UNITS_JOULE )
            Conversion = Conversion * CONV_JOULE
         CASE DEFAULT
            CALL AbortWithError( "EnergyConversion: Unknown output units" )
      END SELECT

   END FUNCTION EnergyConversion


!********************************************************************************************************


!*******************************************************************************
!      MassConversion
!*******************************************************************************
!> Give mass conversion factor from InUnits to OutUnits.
!>
!> @param InUnits    Units data type with the input units.
!> @param OutUnits   Units data type with the output units.
!> @returns          The conversion factor.
!*******************************************************************************
   FUNCTION MassConversion( InUnits, OutUnits )  RESULT( Conversion )
      IMPLICIT NONE
      TYPE (Units), INTENT(IN)   :: InUnits, OutUnits
      REAL                       :: Conversion

      IF ( InUnits%MassUnit == OutUnits%MassUnit ) THEN
         Conversion = 1.0
         RETURN
      END IF

      SELECT CASE ( InUnits%MassUnit )
         CASE ( UNITS_MELEC )
            Conversion = 1.0
         CASE ( UNITS_AMU )
            Conversion = 1.0 / CONV_AMU
         CASE ( UNITS_KG )
            Conversion = 1.0 / CONV_KG
         CASE DEFAULT
            CALL AbortWithError( "MassConversion: Unknown input units" )
      END SELECT

      SELECT CASE ( OutUnits%MassUnit )
         CASE ( UNITS_MELEC )
            Conversion = Conversion
         CASE ( UNITS_AMU )
            Conversion = Conversion * CONV_AMU
         CASE ( UNITS_KG )
            Conversion = Conversion * CONV_KG
         CASE DEFAULT
            CALL AbortWithError( "MassConversion: Unknown output units" )
      END SELECT

   END FUNCTION MassConversion


!********************************************************************************************************


!*******************************************************************************
!      AngleConversion
!*******************************************************************************
!> Give angle conversion factor from InUnits to OutUnits.
!>
!> @param InUnits    Units data type with the input units.
!> @param OutUnits   Units data type with the output units.
!> @returns          The conversion factor.
!*******************************************************************************
   FUNCTION AngleConversion( InUnits, OutUnits )  RESULT( Conversion )
      IMPLICIT NONE
      TYPE (Units), INTENT(IN)   :: InUnits, OutUnits
      REAL                       :: Conversion

      IF ( InUnits%AngleUnit == OutUnits%AngleUnit ) THEN
         Conversion = 1.0
         RETURN
      END IF

      SELECT CASE ( InUnits%AngleUnit )
         CASE ( UNITS_RAD )
            Conversion = 1.0
         CASE ( UNITS_DEGREE )
            Conversion = 1.0 / CONV_DEGREE
         CASE DEFAULT
            CALL AbortWithError( "AngleConversion: Unknown input units" )
      END SELECT

      SELECT CASE ( OutUnits%AngleUnit )
         CASE ( UNITS_RAD )
            Conversion = Conversion
         CASE ( UNITS_DEGREE )
            Conversion = Conversion * CONV_DEGREE
         CASE DEFAULT
            CALL AbortWithError( "AngleConversion: Unknown output units" )
      END SELECT

   END FUNCTION AngleConversion


!********************************************************************************************************


!*******************************************************************************
!      TimeConversion
!*******************************************************************************
!> Give time conversion factor from InUnits to OutUnits.
!>
!> @param InUnits    Units data type with the input units.
!> @param OutUnits   Units data type with the output units.
!> @returns          The conversion factor.
!*******************************************************************************
   FUNCTION TimeConversion( InUnits, OutUnits )  RESULT( Conversion )
      IMPLICIT NONE
      TYPE (Units), INTENT(IN)   :: InUnits, OutUnits
      REAL                       :: Conversion

      IF ( InUnits%TimeUnit == OutUnits%TimeUnit ) THEN
         Conversion = 1.0
         RETURN
      END IF

      SELECT CASE ( InUnits%TimeUnit )
         CASE ( UNITS_AUTIME )
            Conversion = 1.0
         CASE ( UNITS_FEMTOS )
            Conversion = 1.0 / CONV_FEMTOS
         CASE ( UNITS_PICOS )
            Conversion = 1.0 / CONV_PICOS
         CASE DEFAULT
            CALL AbortWithError( "TimeConversion: Unknown input units" )
      END SELECT

      SELECT CASE ( OutUnits%TimeUnit )
         CASE ( UNITS_AUTIME )
            Conversion = Conversion
         CASE ( UNITS_FEMTOS )
            Conversion = Conversion * CONV_FEMTOS
         CASE ( UNITS_PICOS )
            Conversion = Conversion * CONV_PICOS
         CASE DEFAULT
            CALL AbortWithError( "TimeConversion: Unknown output units" )
      END SELECT

   END FUNCTION TimeConversion

!********************************************************************************************************

!*******************************************************************************
!      TemperatureConversion
!*******************************************************************************
!> Give temperature conversion factor from InUnits to OutUnits.
!>
!> @param InUnits    Units data type with the input units.
!> @param OutUnits   Units data type with the output units.
!> @returns          The conversion factor.
!*******************************************************************************
   FUNCTION TemperatureConversion( InUnits, OutUnits )  RESULT( Conversion )
      IMPLICIT NONE
      TYPE (Units), INTENT(IN)   :: InUnits, OutUnits
      REAL                       :: Conversion

      IF ( InUnits%TemperUnit == OutUnits%TemperUnit ) THEN
         Conversion = 1.0
         RETURN
      END IF

      SELECT CASE ( InUnits%TemperUnit )
         CASE ( UNITS_AUTEMP )
            Conversion = 1.0
         CASE ( UNITS_KELVIN )
            Conversion = 1.0 / CONV_KELVIN
         CASE DEFAULT
            CALL AbortWithError( "TemperatureConversion: Unknown input units" )
      END SELECT

      SELECT CASE ( OutUnits%TemperUnit )
         CASE ( UNITS_AUTEMP )
            Conversion = Conversion
         CASE ( UNITS_KELVIN )
            Conversion = Conversion * CONV_KELVIN
         CASE DEFAULT
            CALL AbortWithError( "TemperatureConversion: Unknown output units" )
      END SELECT

   END FUNCTION TemperatureConversion


!********************************************************************************************************

!*******************************************************************************
!      FreqConversion
!*******************************************************************************
!> Give frequency conversion factor from InUnits to OutUnits.
!>
!> @param InUnits    Units data type with the input units.
!> @param OutUnits   Units data type with the output units.
!> @returns          The conversion factor.
!*******************************************************************************
   FUNCTION FreqConversion( InUnits, OutUnits )  RESULT( Conversion )
      IMPLICIT NONE
      TYPE (Units), INTENT(IN)   :: InUnits, OutUnits
      REAL                       :: Conversion

      IF ( InUnits%FreqUnit == OutUnits%FreqUnit ) THEN
         Conversion = 1.0
         RETURN
      END IF

      SELECT CASE ( InUnits%FreqUnit )
         CASE ( UNITS_AUFREQ )
            Conversion = 1.0
         CASE ( UNITS_CMMINUS1 )
            Conversion = 1.0 / CONV_CMMINUS1
         CASE ( UNITS_FSMINUS1 )
            Conversion = 1.0 / CONV_FSMINUS1
         CASE DEFAULT
            CALL AbortWithError( "FreqConversion: Unknown input units" )
      END SELECT

      SELECT CASE ( OutUnits%FreqUnit )
         CASE ( UNITS_AUFREQ )
            Conversion = Conversion
         CASE ( UNITS_CMMINUS1 )
            Conversion = Conversion * CONV_CMMINUS1
         CASE ( UNITS_FSMINUS1 )
            Conversion = Conversion * CONV_FSMINUS1
         CASE DEFAULT
            CALL AbortWithError( "FreqConversion: Unknown output units" )
      END SELECT

   END FUNCTION FreqConversion

!********************************************************************************************************

!*******************************************************************************
!   ForceConversion
!*******************************************************************************
!> Give force conversion factor from InUnits to OutUnits.
!>
!> @param InUnits    Units data type with the input units.
!> @param OutUnits   Units data type with the output units.
!> @returns          The conversion factor.
!*******************************************************************************
   FUNCTION ForceConversion( InUnits, OutUnits ) RESULT ( Conversion )
      IMPLICIT NONE
      TYPE (Units), INTENT(IN)   :: InUnits, OutUnits
      REAL                       :: Conversion

      Conversion = EnergyConversion( InUnits, OutUnits ) / LengthConversion( InUnits, OutUnits )
   END FUNCTION ForceConversion


END MODULE UnitConversion

