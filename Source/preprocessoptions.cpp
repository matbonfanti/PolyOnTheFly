/* This is the general preprocess include file. It should be included in every fortran code at
   the place where the USE statements are located. This ensures the correct working of all the 
   macros and defines located here and supplied by the make files... 
*/

/*==============================================================================*
 * Define alias for log file management.                                        *
 *==============================================================================*/

#if defined(WITH_MPI)
#if defined(LOG_FILE) 
#define __LOG_UNIT              19+my_rank
#define __OPEN_LOG_FILE         OPEN( UNIT=19+my_rank, FILE=LOG_FILE//rank, ACTION="write", POSITION="append" )
#define __CLOSE_LOG_FILE        CLOSE( UNIT=19+my_rank )
#endif
#else
#if defined(LOG_FILE) 
#define __LOG_UNIT              19
#define __OPEN_LOG_FILE         OPEN( UNIT=19, FILE=LOG_FILE, ACTION="write", POSITION="append" )
#define __CLOSE_LOG_FILE        CLOSE( UNIT=19 )
#endif
#endif

/*==============================================================================*
 * Define alias for OMP parallelization.                                        *
 *==============================================================================*/

#if defined(WITH_OPENMP)
#define __OMP_TotalNrOfThreads      OMP_GET_NUM_THREADS()
#define __OMP_CurrentThreadNum      OMP_GET_THREAD_NUM() + 1
#define __OMP_OnlyMasterBEGIN       IF ( OMP_GET_THREAD_NUM() == 0 ) THEN;
#define __OMP_OnlyMasterEND         ; END IF
#else
#define __OMP_TotalNrOfThreads      1
#define __OMP_CurrentThreadNum      1
#define __OMP_OnlyMasterBEGIN   
#define __OMP_OnlyMasterEND     
#endif

/*==============================================================================*
 * Define alias for MPI parallelization.                                        *
 *==============================================================================*/

#if defined(WITH_MPI)
#define __MPI_CurrentNrOfProc       my_rank
#define __MPI_TotalNrOfProcs        num_procs
#define __MPI_OnlyMasterBEGIN       IF (my_rank == 0) THEN;
#define __MPI_OnlyMasterEND         ; END IF
#else
#define __MPI_CurrentNrOfProc       0
#define __MPI_TotalNrOfProcs        1
#define __MPI_OnlyMasterBEGIN      
#define __MPI_OnlyMasterEND        
#endif

/*==============================================================================*
 * Define macro for Siesta-As-A-Subroutine module.                              *
 *==============================================================================*/

#define FC_HAVE_FLUSH


/*==============================================================================*
 * Load common modules.                                                         *
 *==============================================================================*/

   /* These are modules that are used in most parts of the code */
   USE ErrorTrap
   USE MyConsts
   USE MyLinearAlgebra
#if defined(WITH_OPENMP)
   USE omp_lib
#endif
   USE MyMPI


