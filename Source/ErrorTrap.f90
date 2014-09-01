!***************************************************************************************
!*                              ErrorTrap
!***************************************************************************************
!*
!*        CONTAINS:           SetupErrorTrapModule()
!*                            DisposeErrorTrapModule()
!*
!*                            ERROR( Test, Message )
!*                            WARN( Test, Message )
!*                            AbortWithError( Message )
!*                            ShowWarning( Message )
!*                            ShowError( Message )
!*
!***************************************************************************************
!*
!*        AUTHOR(S):          1st: M.F.Somers, may 2002.
!*
!***************************************************************************************

MODULE ErrorTrap

INTEGER :: ErrorTrapModuleSetupFlag = 0

CONTAINS

!***************************************************************************************
!* Setup the module ...

SUBROUTINE SetupErrorTrapModule( )
IMPLICIT NONE

    ErrorTrapModuleSetupFlag = ErrorTrapModuleSetupFlag + 1
    IF( ErrorTrapModuleSetupFlag > 1 ) RETURN

END SUBROUTINE SetupErrorTrapModule

!***************************************************************************************
!* Setup the module ...

SUBROUTINE DisposeErrorTrapModule( )
IMPLICIT NONE
    
    IF(ErrorTrapModuleSetupFlag<=0) RETURN
    ErrorTrapModuleSetupFlag = ErrorTrapModuleSetupFlag - 1
    IF( ErrorTrapModuleSetupFlag > 0 ) RETURN

END SUBROUTINE DisposeErrorTrapModule

!***************************************************************************************
!* Displays an error message to stdout...

SUBROUTINE ShowError( Message )
IMPLICIT NONE
CHARACTER(*) :: Message

    IF( ErrorTrapModuleSetupFlag <= 0 ) CALL SetupErrorTrapModule()

    PRINT *
    PRINT '(A,A)', 'ERROR: ', TRIM(Message)
    PRINT *

END SUBROUTINE ShowError

!***************************************************************************************
!* Displays a warning message to stdout...

SUBROUTINE ShowWarning( Message )
IMPLICIT NONE
CHARACTER(*) :: Message

    IF( ErrorTrapModuleSetupFlag <= 0 ) CALL SetupErrorTrapModule()

    PRINT *
    PRINT '(A,A)', 'WARNING: ', TRIM(Message)
    PRINT *

END SUBROUTINE ShowWarning

!***************************************************************************************
!* Displays a message to stdout en stops the program...

SUBROUTINE AbortWithError( String )
IMPLICIT NONE
CHARACTER(*) :: String

         IF( ErrorTrapModuleSetupFlag <= 0 ) CALL SetupErrorTrapModule()

         CALL ShowError( String )

         STOP

END SUBROUTINE AbortWithError

!***************************************************************************************
!* Test the test given and if true, error message is displayed and program aborted...

SUBROUTINE ERROR( Test, Message )
IMPLICIT NONE
LOGICAL       :: Test
CHARACTER(*)  :: Message

    IF( ErrorTrapModuleSetupFlag <= 0 ) CALL SetupErrorTrapModule()

    IF( Test ) CALL AbortWithError(Message)

END SUBROUTINE ERROR

!***************************************************************************************
!* Test the test given and if true the warning is displayed...

SUBROUTINE WARN( Test, Message )
IMPLICIT NONE
LOGICAL       :: Test
CHARACTER(*)  :: Message

    IF( ErrorTrapModuleSetupFlag <= 0 ) CALL SetupErrorTrapModule()

    IF( Test ) CALL ShowWarning(Message)

END SUBROUTINE WARN

!***************************************************************************************

END MODULE ErrorTrap
