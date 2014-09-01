!***************************************************************************************
!*                              MODULE InputField
!***************************************************************************************
!
!>  \brief     Read data from input file
!>  \details   This module is a wrapper for input data procedure. \n
!>             Data related to an input file is store in the InputFile data type. 
!>             An input file can be opened and closed with the subs OpenFile()
!>             and CloseFile(). The fields contained in an open input file can be
!>             read with the subroutine  SetFieldFromInput(), specifying the 
!>             InputFile data, the name of the field to look for in the input file 
!>             and the variable where to store the input data. \n
!>             This module allows a semi-free input format where the data are 
!>             written as " FIELDNAME : VARIABLE "
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             July 2012
!>
!***************************************************************************************
!
!>   \remark         The module is based on Mark Somers' inputfiels.f90
!
!***************************************************************************************
MODULE InputField
   USE ErrorTrap
   USE MyConsts

   PRIVATE
   PUBLIC :: InputFile, OpenFile, CloseFile, SetFieldFromInput

   INTEGER, PARAMETER, PRIVATE :: FILE_IS_OPEN   = 1
   INTEGER, PARAMETER, PRIVATE :: FILE_IS_CLOSED = 0

   !> Data type to store information on the input file.
   TYPE InputFile
      PRIVATE
      CHARACTER(100) :: Name
      INTEGER        :: Unit
      INTEGER        :: Status = FILE_IS_CLOSED
   END TYPE InputFile

   !> A wrapper for different reading different kind of input data.
   INTERFACE SetFieldFromInput
      MODULE PROCEDURE SetRealFieldFromInput, &
               SetIntegerFieldFromInput, SetStringFieldFromInput, SetLogicalFieldFromInput
   END INTERFACE

!********************************************************************************************************
   CONTAINS
!********************************************************************************************************

!*******************************************************************************
!          OpenFile
!*******************************************************************************
!> Wrapper to open an input file with a given name 
!>
!> @param Input     Data type to store the input file related variables
!> @param FileName  Name of the file to open
!*******************************************************************************
   SUBROUTINE OpenFile( Input, FileName )
      IMPLICIT NONE
      TYPE( InputFile ), INTENT(INOUT) :: Input
      CHARACTER(*)                     :: FileName
      INTEGER                          :: InputStatus

      ! check if file of the datatype is closed
      CALL ERROR( Input%Status == FILE_IS_OPEN, " OpenFile: trying to open already open file " )

      ! store input name
      Input%Name = TRIM(ADJUSTL(FileName))
      ! store input unit
      Input%Unit = LookForFreeUnit()

      ! Open file
      OPEN( FILE = Input%Name, UNIT=Input%Unit, IOSTAT=InputStatus )

      ! Check if input file was found on disk
      CALL ERROR( ( InputStatus /= 0 ), " OpenFile: could not find input file " // Input%Name )

      ! Set status of the input file as "open"
      Input%Status = FILE_IS_OPEN

#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A,A,A,I3)") "OpenFile: File ", trim(FileName), " opened as unit ", Input%Unit
#endif

   END SUBROUTINE OpenFile




!*******************************************************************************
!           CloseFile
!*******************************************************************************
!>  Wrapper to close the input file 
!>
!>  @param Input     Data type to store the input file related variables
!*******************************************************************************
   SUBROUTINE CloseFile( Input )
      IMPLICIT NONE
      TYPE( InputFile ), INTENT(INOUT) :: Input
      INTEGER                          :: InputStatus

#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A,A,A,I3)") " CloseFile: Closing file ", trim(Input%Name), " opened as unit ", Input%Unit
#endif

      ! check if file of the datatype is open
      CALL ERROR( Input%Status == FILE_IS_CLOSED, " CloseFile: trying to close an already closed file " )

      ! CLose file
      CLOSE( Input%Unit )

      ! destroy file info
      Input%Name = ""
      ! store input unit
      Input%Unit = 0

      ! Set status of the input file as "closed"
      Input%Status = FILE_IS_CLOSED

   END SUBROUTINE CloseFile


!*******************************************************************************
!           SetRealFieldFromInput
!*******************************************************************************
!>  Look for the field name in input file and store the corresponding
!>  real variable. If the optional argument DefaultV is absent, the subroutine
!>  stops with error in case the Field is not found in the input file. Otherwise,
!>  the default value DefaultV is assigned to the variable.
!>
!>  @param Input       Data type to store the input file related variables
!>  @param FieldName   Field name to find in the input file 
!>  @param Variable    Real variable to store input data
!>  @param DefaultV    Default value to assign to the variable
!*******************************************************************************
   SUBROUTINE SetRealFieldFromInput( Input, FieldName, Variable, DefaultV)
      IMPLICIT NONE
      TYPE( InputFile ), INTENT(INOUT) :: Input
      CHARACTER(*), INTENT(IN)         :: FieldName
      REAL, INTENT(OUT)                :: Variable
      REAL, INTENT(IN), OPTIONAL       :: DefaultV

      CHARACTER(len=1024)              :: Line
      INTEGER                          :: LineLength, ColonPos, ReadingStatus
      CHARACTER(30)                    :: String, RdFormat

      ! check if file of the datatype is open
      CALL ERROR( Input%Status == FILE_IS_CLOSED, " SetRealFieldFromInput: trying to read a closed file input file " )

#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A,A,A,I3)") " SetRealFieldFromInput: Looking for real field ", trim(FieldName), " in unit ", Input%Unit
#endif

      ! Rewind input file
      REWIND( Input%Unit )

      ! cycle over the lines of the input file
      DO 
         ! read input file
         READ( Input%Unit, "(A1024)", IOSTAT=ReadingStatus ) Line

         ! if file is finished without finding the fild, give error or set default value
         IF ( ReadingStatus /= 0 ) THEN
            IF ( .NOT. PRESENT(DefaultV ) ) THEN
               CALL AbortWithError( " SetRealFieldFromInput: could not find field "//TRIM(FieldName) )
            ELSE
               Variable = DefaultV
               WRITE(String,*) DefaultV
               CALL ShowWarning( "Variable "//FieldName//" set to default value of "//String ); EXIT
            END IF
         END IF

         ! remove leading blanks...
         Line = ADJUSTL( Line )
         ! get length of line without the spaces at the end...
         LineLength = LEN_TRIM( Line )

         ! skip line if empty or starts with #
         IF ( LineLength == 0 .OR. Line(1:1) == '#' ) CYCLE

         ! Find position of the first colon
         ColonPos = SCAN( Line, ':' )
         ! skip line if not present
         IF ( ColonPos <= 1 ) CYCLE

         ! check if the field name corresponds
         IF ( TRIM( FieldName ) ==  TRIM( ADJUSTL( Line(1:ColonPos-1) ) )  ) THEN
#if defined(VERBOSE_OUTPUT)
            WRITE(*,*) " Found field ", trim(FieldName)
#endif
            ! rewind one record
            BACKSPACE( Input%Unit )
            ! define reading format
            WRITE( String, * ) ColonPos
            RdFormat = "(A"//TRIM( ADJUSTL( String))//",F200.0)" 
            ! store the value of the variable
            READ(  Input%Unit, RdFormat ) String, Variable
#if defined(VERBOSE_OUTPUT)
            WRITE(*,*) " Set variable equal to ", Variable
#endif
            EXIT
         ENDIF
      END DO

   END SUBROUTINE SetRealFieldFromInput

!*******************************************************************************
!           SetIntegerFieldFromInput
!*******************************************************************************
!>  Look for the field name in input file and store the corresponding
!>  integer variable. If the optional argument DefaultV is absent, the subroutine
!>  stops with error in case the Field is not found in the input file. Otherwise,
!>  the default value DefaultV is assigned to the variable.
!>
!>  @param Input       Data type to store the input file related variables
!>  @param FieldName   Field name to find in the input file 
!>  @param Variable    Integer variable to store input data
!>  @param DefaultV    Default value to assign to the variable
!*******************************************************************************
   SUBROUTINE SetIntegerFieldFromInput( Input, FieldName, Variable, DefaultV)
      IMPLICIT NONE
      TYPE( InputFile ), INTENT(INOUT) :: Input
      CHARACTER(*), INTENT(IN)         :: FieldName
      INTEGER, INTENT(OUT)             :: Variable
      INTEGER, INTENT(IN), OPTIONAL    :: DefaultV

      CHARACTER(len=1024)              :: Line
      INTEGER                          :: LineLength, ColonPos, ReadingStatus
      CHARACTER(30)                    :: String, RdFormat

      ! check if file of the datatype is open
      CALL ERROR( Input%Status == FILE_IS_CLOSED, " SetIntegerFieldFromInput: trying to read a closed file input file " )

#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A,A,A,I3)") " SetIntegerFieldFromInput: Looking for integer field ", trim(FieldName), " in unit ", Input%Unit
#endif

      ! Rewind input file
      REWIND( Input%Unit )

      ! cycle over the lines of the input file
      DO 
         ! read input file
         READ( Input%Unit, "(A1024)", IOSTAT=ReadingStatus ) Line

         ! if file is finished without finding the fild, give error or set default value
         IF ( ReadingStatus /= 0 ) THEN
            IF ( .NOT. PRESENT(DefaultV ) ) THEN
               CALL AbortWithError( " SetIntegerFieldFromInput: could not find field "//TRIM(FieldName) )
            ELSE
               Variable = DefaultV
               WRITE(String,*) DefaultV
               CALL ShowWarning( "Variable "//FieldName//" set to default value of "//String ); EXIT
            END IF
         END IF

         ! remove leading blanks...
         Line = ADJUSTL( Line )
         ! get length of line without the spaces at the end...
         LineLength = LEN_TRIM( Line )

         ! skip line if empty or starts with #
         IF ( LineLength == 0 .OR. Line(1:1) == '#' ) CYCLE

         ! Find position of the first colon
         ColonPos = SCAN( Line, ':' )
         ! skip line if not present
         IF ( ColonPos <= 1 ) CYCLE

         ! check if the field name corresponds
         IF ( TRIM( FieldName ) ==  TRIM( ADJUSTL( Line(1:ColonPos-1) ) )  ) THEN
#if defined(VERBOSE_OUTPUT)
            WRITE(*,*) " Found field ", trim(FieldName)
#endif
            ! rewind one record
            BACKSPACE( Input%Unit )
            ! define reading format
            WRITE( String, * ) ColonPos
            RdFormat = "(A"//TRIM( ADJUSTL( String))//",I200)" 
            ! store the value of the variable
            READ(  Input%Unit, RdFormat ) String, Variable
#if defined(VERBOSE_OUTPUT)
            WRITE(*,*) " Set variable equal to ", Variable
#endif
            EXIT
         ENDIF
      END DO

   END SUBROUTINE SetIntegerFieldFromInput


!*******************************************************************************
!           SetStringFieldFromInput
!*******************************************************************************
!>  Look for the field name in input file and store the corresponding
!>  character variable. If the optional argument DefaultV is absent, the subroutine
!>  stops with error in case the Field is not found in the input file. Otherwise,
!>  the default value DefaultV is assigned to the variable.
!>
!>  @param Input       Data type to store the input file related variables
!>  @param FieldName   Field name to find in the input file 
!>  @param Variable    Character variable to store input data
!>  @param DefaultV    Default value to assign to the variable
!*******************************************************************************
   SUBROUTINE SetStringFieldFromInput( Input, FieldName, Variable, DefaultV)
      IMPLICIT NONE
      TYPE( InputFile ), INTENT(INOUT) :: Input
      CHARACTER(*), INTENT(IN)         :: FieldName
      CHARACTER(*), INTENT(OUT)          :: Variable
      CHARACTER(*), INTENT(IN), OPTIONAL :: DefaultV

      CHARACTER(len=1024)              :: Line
      INTEGER                          :: LineLength, ColonPos, ReadingStatus
      CHARACTER(30)                    :: String, RdFormat

      ! check if file of the datatype is open
      CALL ERROR( Input%Status == FILE_IS_CLOSED, " SetStringFieldFromInput: trying to read a closed file input file " )

#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A,A,A,I3)") " SetStringFieldFromInput: Looking for character field ", trim(FieldName), " in unit ", Input%Unit
#endif

      ! Rewind input file
      REWIND( Input%Unit )

      ! cycle over the lines of the input file
      DO 
         ! read input file
         READ( Input%Unit, "(A1024)", IOSTAT=ReadingStatus ) Line

         ! if file is finished without finding the fild, give error or set default value
         IF ( ReadingStatus /= 0 ) THEN
            IF ( .NOT. PRESENT(DefaultV ) ) THEN
               CALL AbortWithError( " SetStringFieldFromInput: could not find field "//TRIM(FieldName) )
            ELSE
               Variable = trim(DefaultV)
               WRITE(String,*) trim(DefaultV)
               CALL ShowWarning( "Variable "//FieldName//" set to default value of "//String ); EXIT
            END IF
         END IF

         ! remove leading blanks...
         Line = ADJUSTL( Line )
         ! get length of line without the spaces at the end...
         LineLength = LEN_TRIM( Line )

         ! skip line if empty or starts with #
         IF ( LineLength == 0 .OR. Line(1:1) == '#' ) CYCLE

         ! Find position of the first colon
         ColonPos = SCAN( Line, ':' )
         ! skip line if not present
         IF ( ColonPos <= 1 ) CYCLE

         ! check if the field name corresponds
         IF ( TRIM( FieldName ) ==  TRIM( ADJUSTL( Line(1:ColonPos-1) ) )  ) THEN
#if defined(VERBOSE_OUTPUT)
            WRITE(*,*) " Found field ", trim(FieldName)
#endif
            ! rewind one record
            BACKSPACE( Input%Unit )
            ! define reading format
            WRITE( String, * ) ColonPos
            RdFormat = "(A"//TRIM( ADJUSTL( String))//",A)" 
            ! store the value of the variable
            READ(  Input%Unit, RdFormat ) String, Variable
#if defined(VERBOSE_OUTPUT)
            WRITE(*,*) " Set variable equal to ", Variable
#endif
            EXIT
         ENDIF
      END DO

   END SUBROUTINE SetStringFieldFromInput

!*******************************************************************************
!           SetLogicalFieldFromInput
!*******************************************************************************
!>  Look for the field name in input file and store the corresponding
!>  boolean variable. If the optional argument DefaultV is absent, the subroutine
!>  stops with error in case the Field is not found in the input file. Otherwise,
!>  the default value DefaultV is assigned to the variable.
!>
!>  @param Input       Data type to store the input file related variables
!>  @param FieldName   Field name to find in the input file 
!>  @param Variable    Logical variable to store input data
!>  @param DefaultV    Default value to assign to the variable
!*******************************************************************************
   SUBROUTINE SetLogicalFieldFromInput( Input, FieldName, Variable, DefaultV)
      IMPLICIT NONE
      TYPE( InputFile ), INTENT(INOUT) :: Input
      CHARACTER(*), INTENT(IN)         :: FieldName
      LOGICAL, INTENT(OUT)             :: Variable
      LOGICAL, INTENT(IN), OPTIONAL    :: DefaultV

      CHARACTER(len=1024)              :: Line
      INTEGER                          :: LineLength, ColonPos, ReadingStatus
      CHARACTER(30)                    :: String, RdFormat

      ! check if file of the datatype is open
      CALL ERROR( Input%Status == FILE_IS_CLOSED, " SetLogicalFieldFromInput: trying to read a closed file input file " )

#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A,A,A,I3)") " SetLogicalFieldFromInput: Looking for boolean field ", trim(FieldName), " in unit ", Input%Unit
#endif

      ! Rewind input file
      REWIND( Input%Unit )

      ! cycle over the lines of the input file
      DO 
         ! read input file
         READ( Input%Unit, "(A1024)", IOSTAT=ReadingStatus ) Line

         ! if file is finished without finding the fild, give error or set default value
         IF ( ReadingStatus /= 0 ) THEN
            IF ( .NOT. PRESENT(DefaultV ) ) THEN
               CALL AbortWithError( " SetLogicalFieldFromInput: could not find field "//TRIM(FieldName) )
            ELSE
               Variable = DefaultV
               WRITE(String,*) DefaultV
               CALL ShowWarning( "Variable "//FieldName//" set to default value of "//String ); EXIT
            END IF
         END IF

         ! remove leading blanks...
         Line = ADJUSTL( Line )
         ! get length of line without the spaces at the end...
         LineLength = LEN_TRIM( Line )

         ! skip line if empty or starts with #
         IF ( LineLength == 0 .OR. Line(1:1) == '#' ) CYCLE

         ! Find position of the first colon
         ColonPos = SCAN( Line, ':' )
         ! skip line if not present
         IF ( ColonPos <= 1 ) CYCLE

         ! check if the field name corresponds
         IF ( TRIM( FieldName ) ==  TRIM( ADJUSTL( Line(1:ColonPos-1) ) )  ) THEN
#if defined(VERBOSE_OUTPUT)
            WRITE(*,*) " Found field ", trim(FieldName)
#endif
            ! rewind one record
            BACKSPACE( Input%Unit )
            ! define reading format
            WRITE( String, * ) ColonPos
            RdFormat = "(A"//TRIM( ADJUSTL( String))//",L200)" 
            ! store the value of the variable
            READ(  Input%Unit, RdFormat ) String, Variable
#if defined(VERBOSE_OUTPUT)
            WRITE(*,*) " Set variable equal to ", Variable
#endif
            EXIT
         ENDIF
      END DO

   END SUBROUTINE SetLogicalFieldFromInput


END MODULE InputField
