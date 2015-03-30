!***************************************************************************************
!*                           MODULE DFTBWrapper
!***************************************************************************************
!
!>  \brief     Wrapper for DFTB+ code
!>  \details   This module contains subroutines to use in the molecular
!>             dynamics simulation forces coming from the DFTB+ code.
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             30 March 2015
!
!***************************************************************************************
!
!>  \par Updates
!>  \arg _______________: ____________________
!
!>  \todo          \arg __________________________________________________
!>                 
!***************************************************************************************
MODULE DFTBWrapper
#include "preprocessoptions.cpp"
   USE UnitConversion

   PRIVATE

   PUBLIC :: DFTBInputParser, WriteDFTBInput

   TYPE(Units), SAVE     :: DFTBUnits     !< Units definition for DFTB+

 CONTAINS

   SUBROUTINE DFTBInputParser( HSDFile, ReadFile, PBC, UnitVectors, AtomNo, &
                                                                    AtomicLabels, AtomicPositions )
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: HSDFile
      CHARACTER(120), DIMENSION(200), INTENT(OUT)          :: ReadFile
      LOGICAL, INTENT(OUT)                                 :: PBC
      REAL, DIMENSION(3,3), INTENT(OUT)                    :: UnitVectors
      INTEGER, INTENT(OUT)                                 :: AtomNo
      CHARACTER(3), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: AtomicLabels
      REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)       :: AtomicPositions

      CHARACTER(1) :: Descriptor
      CHARACTER(120) :: Line, AtomLabelsString
      CHARACTER(3), DIMENSION(20) :: LabelsArray
      INTEGER :: UnitNr, NrOfLines, Stat
      INTEGER :: GeomStart, GeomEnd, NrOfLinesWithoutGeom
      INTEGER :: LineCount, i, n, Dummy, AtomKind
      INTEGER :: pos1, pos2

      ! Set units for DFTB+
      CALL Initialize_UnitConversion( DFTBUnits, UNITS_ANGSTROM, UNITS_EV, UNITS_AMU, UNITS_DEGREE, UNITS_FEMTOS, &
                                                            UNITS_KELVIN, UNITS_CMMINUS1 )

      ! Open HSD file
      UnitNr = LookForFreeUnit()
      OPEN( UNIT = UnitNr , FILE = TRIM(ADJUSTL(HSDFile)) ) 

      ! Count nr of lines
      DO NrOfLines = 0,100000
         READ( UnitNr, *, IOSTAT=Stat ) 
         IF (Stat /= 0) EXIT
      END DO
      REWIND(UnitNr)

      ! Identify beginning and end of geometry section
      LineCount = 0
      Ext:DO 
        READ(UnitNr,"(A120)") Line
        LineCount = LineCount + 1
        Line = adjustl(Line)
        IF ( Line(1:8) == "Geometry" ) THEN
          GeomStart = LineCount
          DO 
            READ(UnitNr,"(A120)") Line
            LineCount = LineCount + 1
            Line = ADJUSTL(Line)
            IF ( Line(1:1) == "}" ) THEN
               GeomEnd = LineCount
               EXIT Ext
            END IF    
          END DO
        END IF
      END DO Ext
      NrOfLinesWithoutGeom = NrOfLines - ( GeomEnd-GeomStart+1 ) 
      REWIND(UnitNr)

      ! Read geometry section
      LineCount = 0
      DO 
         READ(UnitNr,"(A120)") Line
         LineCount = LineCount + 1

         ! IF geometry section start is encountered
         IF ( LineCount == GeomStart ) THEN

            ! Read number of atoms and kind of geometry
            READ(UnitNr,*) AtomNo, Descriptor
            ! Allocate memory
            ALLOCATE( AtomicLabels(AtomNo), AtomicPositions(3,AtomNo) )

            ! Read labels and tokenize them
            READ(UnitNr,"(A120)") AtomLabelsString
            pos1 = 1; n = 0
            DO
               pos2 = INDEX((AtomLabelsString(pos1:)), " ")
               IF (pos2 == 0) THEN
                  n = n + 1
                  LabelsArray(n) = AtomLabelsString(pos1:)
                  EXIT
               END IF
               IF ( TRIM(AtomLabelsString(pos1:pos1+pos2-2)) == " " ) THEN
                  pos1 = pos2+pos1
                  CYCLE
               END IF
               n = n + 1
               LabelsArray(n) = AtomLabelsString(pos1:pos1+pos2-2)
               pos1 = pos2+pos1
            END DO

            ! Read atoms and assign correct coordinates and labels
            DO i = 1, AtomNo
               READ(UnitNr,*) Dummy, AtomKind, AtomicPositions(:,i)
               AtomicLabels(i) = TRIM(ADJUSTL( LabelsArray(AtomKind) ))
            END DO

            ! Define periodic boundary conditions
            IF ( TRIM(ADJUSTL(Descriptor)) == "F" .OR. TRIM(ADJUSTL(Descriptor)) == "S" ) THEN
               PBC = .TRUE.
               READ(UnitNr,*)
               READ(UnitNr,*) UnitVectors(1,:)
               READ(UnitNr,*) UnitVectors(2,:)
               READ(UnitNr,*) UnitVectors(3,:)
               UnitVectors = UnitVectors*LengthConversion(DFTBUnits,InternalUnits)
            ELSE IF ( TRIM(ADJUSTL(Descriptor)) == "C" ) THEN
               PBC = .FALSE.
               UnitVectors = RESHAPE( (/1.E+100, 0.0, 0.0 , 0.0, 1.E+100, 0.0, 0.0, 0.0, 1.E+100/), (/3,3/) )
            END IF

            ! Convert fractional coordinates to cartesian or anstrom coordiantes to atomic units
            IF ( TRIM(ADJUSTL(Descriptor)) == "F" ) THEN
               DO i = 1, AtomNo
                  AtomicPositions(:,i) = TheOneWithMatrixVectorProduct( UnitVectors, AtomicPositions(:,i) )
               END DO
            ELSE IF ( TRIM(ADJUSTL(Descriptor)) == "C" .OR. TRIM(ADJUSTL(Descriptor)) == "S" ) THEN
               AtomicPositions = AtomicPositions*LengthConversion(DFTBUnits,InternalUnits)
            END IF

            ! Exit from input file cycle
            EXIT

         END IF
      END DO

      ! Rewing input unit and read lines other than geometry
      REWIND(UnitNr)
      CALL ERROR( NrOfLinesWithoutGeom > SIZE(ReadFile), " DFTBInputParserWithGeometry: too many .HSD input lines... " )

      n = 0
      DO i = 1, NrOfLines
         READ(UnitNr,"(A120)") Line
         IF ( (i<GeomStart) .OR. (i>GeomEnd) ) THEN
            n = n+1
            ReadFile(n) = TRIM(ADJUSTL(Line))
         END IF
      END DO

      ! Close input unit
      CLOSE( UnitNr )
      
   END SUBROUTINE DFTBInputParser

! ================================================================================================================

   SUBROUTINE WriteDFTBInput( Dir, PBC, UnitVectors, AtomNo, AtomTypeString, AtomicPositions, AtomTypes, ReadFileContent )
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN)                            :: Dir
      LOGICAL, INTENT(IN)                                 :: PBC
      REAL, DIMENSION(3,3), INTENT(IN)                    :: UnitVectors
      INTEGER, INTENT(IN)                                 :: AtomNo
      CHARACTER(120)                                      :: AtomTypeString
      REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(IN)       :: AtomicPositions
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN)      :: AtomTypes
      CHARACTER(120), DIMENSION(200), INTENT(IN)          :: ReadFileContent

      INTEGER :: OutputUnit, i

      ! Open output unit
      OutputUnit = LookForFreeUnit()
      OPEN( UNIT=OutputUnit, FILE=TRIM(ADJUSTL(Dir))//"/dftb_in.hsd" )

      ! Write geometry section in GenFormat
      IF ( PBC ) THEN
         WRITE(OutputUnit, 300) AtomNo, "S", AtomTypeString
      ELSE IF ( .NOT. PBC ) THEN
         WRITE(OutputUnit, 300) AtomNo, "C", AtomTypeString
      END IF
      DO i = 1, AtomNo
         WRITE(OutputUnit, 301) i, AtomTypes(i), AtomicPositions(:,i)
      END DO
      IF ( PBC ) WRITE(OutputUnit, 302) 0., 0., 0., UnitVectors(1,:), UnitVectors(2,:), UnitVectors(3,:)
      WRITE(OutputUnit, 303)

      ! Write rest of the output
      DO i = 1, SIZE(ReadFileContent)
         IF ( ReadFileContent(i)(1:1) /= "#" )  WRITE(OutputUnit,*) ReadFileContent(i)
      END DO

      WRITE(OutputUnit, 400)

      ! Close output unit
      CLOSE( OutputUnit )

      300 FORMAT( " Geometry = GenFormat { ", /, &
                              5X, I5, 1X, A1, /, &
                                           A  ) 
      301 FORMAT( I5, 1X, I3, 1X, 3F20.8 )
      302 FORMAT( 3F20.8, /, 3F20.8, /, 3F20.8, /, 3F20.8 )
      303 FORMAT( " } ")

      400 FORMAT( " Options = { ", /, "    CalculateForces = Yes ", /, " } " )


   END SUBROUTINE WriteDFTBInput

! ================================================================================================================


SUBROUTINE ReadDFTBForces()

END SUBROUTINE ReadDFTBForces

END MODULE DFTBWrapper
