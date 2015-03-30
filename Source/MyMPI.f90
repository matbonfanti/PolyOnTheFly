!***************************************************************************************
!*                              MODULE MyMPI
!***************************************************************************************
!
!>  \brief     Wrapper for MPI
!>  \details   
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             October 2014
!>
!***************************************************************************************
!
!>  \pre              The module needs to be initialized by calling the MyMPI_Init 
!>                    subroutine.
!
!***************************************************************************************
!
!>  \par Updates
!>  \arg ________________________________
!
!>  \todo            \arg Put status check at the beginning of all the subroutines
!>                 
!***************************************************************************************
MODULE MyMPI
   IMPLICIT NONE
#if defined(WITH_MPI)
   include 'mpif.h'
#endif

   PRIVATE
   PUBLIC :: MyMPI_Init, MyMPI_Barrier, MyMPI_Finalize
   PUBLIC :: MyMPI_SumToMaster, MyMPI_BroadcastToSlaves, MyMPI_AllReduceMaxValue
   
   !> Setup variable of the module
   LOGICAL :: MyMPISetup = .FALSE.

   INTEGER,PUBLIC :: my_rank, num_procs, err  
   CHARACTER(3),PUBLIC   :: rank


   INTERFACE MyMPI_BroadcastToSlaves
     MODULE PROCEDURE MyMPI_BroadcastToSlaves_RA, MyMPI_BroadcastToSlaves_RS, & 
                      MyMPI_BroadcastToSlaves_IS, MyMPI_BroadcastToSlaves_LS, &
                      MyMPI_BroadcastToSlaves_IA, MyMPI_BroadcastToSlaves_SA
   END INTERFACE


   CONTAINS


!*******************************************************************************
! MyMPI_Init
!*******************************************************************************
!>  Setup MPI and set relevant variables (rank of current proc, total nr of MPI procs)
!*******************************************************************************
   SUBROUTINE MyMPI_Init()
      IMPLICIT NONE

      ! INSERT HERE THE CHECK ON THE MODULE STATUS

#if defined(WITH_MPI)
      call MPI_INIT( err )
      call MPI_COMM_RANK( MPI_COMM_WORLD, my_rank, err )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, num_procs, err )
#else
      my_rank = 0
      num_procs = 1
#endif

      ! Store string of the rank
      WRITE(rank,'(I3.3)') my_rank

      ! MyMPI is now setup
      MyMPISetup = .TRUE.

   END SUBROUTINE MyMPI_Init

   
!*******************************************************************************
! MyMPI_Barrier
!*******************************************************************************
!>  Wrapper for defining a MPI barrier
!*******************************************************************************
   SUBROUTINE MyMPI_Barrier(  )
      IMPLICIT NONE

      ! INSERT HERE THE CHECK ON THE MODULE STATUS
 
      CALL FLUSH()
#if defined(WITH_MPI)
      CALL MPI_Barrier(MPI_COMM_WORLD,err)
#endif

   END SUBROUTINE MyMPI_Barrier
   
   
!*******************************************************************************
! MyMPI_SumToMaster
!*******************************************************************************
!>  Sum the copies of a real array and store the sum on the master node
!*******************************************************************************
   SUBROUTINE MyMPI_SumToMaster( Array )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(INOUT) :: Array
      REAL, DIMENSION(size(Array))      :: Temp

      ! INSERT HERE THE CHECK ON THE MODULE STATUS

#if defined(WITH_MPI)
      CALL MPI_Barrier(MPI_COMM_WORLD,err)
      CALL MPI_Reduce(Array,Temp,size(Array),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,err)
      IF ( my_rank == 0 ) THEN
	 Array = Temp
      END IF
#endif

   END SUBROUTINE MyMPI_SumToMaster
   
!*******************************************************************************
! MyMPI_AllReduceMaxValue
!*******************************************************************************
!>  On all the nodes set the maximum values of integer variables
!*******************************************************************************
   SUBROUTINE MyMPI_AllReduceMaxValue( IntegerVar )
      IMPLICIT NONE
      INTEGER :: IntegerVar

      ! INSERT HERE THE CHECK ON THE MODULE STATUS

#if defined(WITH_MPI)
      CALL MPI_Barrier(MPI_COMM_WORLD,err)
      CALL MPI_AllReduce( IntegerVar, IntegerVar, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, err)
#endif

   END SUBROUTINE MyMPI_AllReduceMaxValue


!*******************************************************************************
! MyMPI_Init
!*******************************************************************************
!>  Setup MPI and set relevant variables (rank of current proc, total nr of MPI procs)
!*******************************************************************************
   SUBROUTINE MyMPI_Finalize()
      IMPLICIT NONE

      ! INSERT HERE THE CHECK ON THE MODULE STATUS

#if defined(WITH_MPI)
      call MPI_FINALIZE(err)
#endif

      ! MyMPI is now not setup
      MyMPISetup = .FALSE.

   END SUBROUTINE MyMPI_Finalize

!*******************************************************************************
! MyMPI_BroadcastToSlaves
!*******************************************************************************
!>  Broadcast the value of a given variable to all the slaves.
!*******************************************************************************
   SUBROUTINE MyMPI_BroadcastToSlaves_RA( Array )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(INOUT) :: Array

      ! INSERT HERE THE CHECK ON THE MODULE STATUS

#if defined(WITH_MPI)
      CALL MPI_Bcast( Array, size(Array), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
#endif

   END SUBROUTINE MyMPI_BroadcastToSlaves_RA

   SUBROUTINE MyMPI_BroadcastToSlaves_RS( RealVar )
      IMPLICIT NONE
      REAL, INTENT(INOUT) :: RealVar

      ! INSERT HERE THE CHECK ON THE MODULE STATUS

#if defined(WITH_MPI)
      CALL MPI_Bcast( RealVar, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
#endif

   END SUBROUTINE MyMPI_BroadcastToSlaves_RS

   SUBROUTINE MyMPI_BroadcastToSlaves_IS( IntVar )
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: IntVar

      ! INSERT HERE THE CHECK ON THE MODULE STATUS

#if defined(WITH_MPI)
      CALL MPI_Bcast( IntVar, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
#endif

   END SUBROUTINE MyMPI_BroadcastToSlaves_IS

   SUBROUTINE MyMPI_BroadcastToSlaves_IA( IntArray )
      IMPLICIT NONE
      INTEGER, DIMENSION(:), INTENT(INOUT) :: IntArray

      ! INSERT HERE THE CHECK ON THE MODULE STATUS

#if defined(WITH_MPI)
      CALL MPI_Bcast( IntArray, size(IntArray), MPI_INTEGER, 0, MPI_COMM_WORLD, err)
#endif

   END SUBROUTINE MyMPI_BroadcastToSlaves_IA

   SUBROUTINE MyMPI_BroadcastToSlaves_LS( LogVar )
      IMPLICIT NONE
      LOGICAL, INTENT(INOUT) :: LogVar

      ! INSERT HERE THE CHECK ON THE MODULE STATUS

#if defined(WITH_MPI)
      CALL MPI_Bcast( LogVar, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, err)
#endif

   END SUBROUTINE MyMPI_BroadcastToSlaves_LS

   SUBROUTINE MyMPI_BroadcastToSlaves_SA( Array )
      IMPLICIT NONE
      CHARACTER(*), DIMENSION(:), INTENT(INOUT) :: Array
      INTEGER :: i
      ! INSERT HERE THE CHECK ON THE MODULE STATUS

#if defined(WITH_MPI)
      DO i = 1, SIZE(Array)
         CALL MPI_BCAST( Array(i), LEN(Array(i)), MPI_CHARACTER, 0, MPI_COMM_WORLD, err)
      END DO
#endif

   END SUBROUTINE MyMPI_BroadcastToSlaves_SA

END MODULE MyMPI



!********************************************* END OF FILE *******************************
