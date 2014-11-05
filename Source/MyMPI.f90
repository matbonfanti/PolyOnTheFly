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
!>  \todo  
!>                 
!***************************************************************************************
MODULE MyMPI
   IMPLICIT NONE
#if defined(WITH_MPI)
   include 'mpif.h'
#endif


   PRIVATE
   PUBLIC :: MyMPI_Init, MyMPI_SumToMaster, MyMPI_Finalize, MyMPI_BroadcastToSlaves

   !> Setup variable of the module
   LOGICAL :: MyMPISetup = .FALSE.

   INTEGER,PUBLIC :: my_rank, num_procs, err  
   CHARACTER(3),PUBLIC   :: rank


   INTERFACE MyMPI_BroadcastToSlaves
     MODULE PROCEDURE MyMPI_BroadcastToSlaves_RA, MyMPI_BroadcastToSlaves_RS, & 
                      MyMPI_BroadcastToSlaves_IS, MyMPI_BroadcastToSlaves_LS
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

   SUBROUTINE MyMPI_BroadcastToSlaves_LS( LogVar )
      IMPLICIT NONE
      LOGICAL, INTENT(INOUT) :: LogVar

      ! INSERT HERE THE CHECK ON THE MODULE STATUS

#if defined(WITH_MPI)
      CALL MPI_Bcast( LogVar, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, err)
#endif

   END SUBROUTINE MyMPI_BroadcastToSlaves_LS


END MODULE MyMPI



!********************************************* END OF FILE *******************************
