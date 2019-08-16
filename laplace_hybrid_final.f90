!*************************************************
! Laplace MPI Fortran Version
!
! Temperature is initially 0.0
! Boundaries are as follows:
!
!          T=0.   
!       _________                    ____|____|____|____  
!       |       | 0                 |    |    |    |    | 
!       |       |                   |    |    |    |    | 
! T=0.  | T=0.0 | T                 | 0  | 1  | 2  | 3  | 
!       |       |                   |    |    |    |    | 
!       |_______| 100               |    |    |    |    | 
!       0     100                   |____|____|____|____|
!                                        |    |    |
! Each Processor works on a sub grid and then sends its
! boundaries to neighbours
!
!  Problem designed by: John Urbanic, PSC 2014
!  Problem Solved by  : Masrul Huda 
!
!*************************************************
program mpilaplace

      use mpi
      use openacc
      implicit none

      !Size of plate
      integer, parameter             :: columns_global=10000
      integer, parameter             :: rows=10000

      !these are the new parameters for parallel purposes
      integer, parameter             :: total_pes=16   ! 4 MPI process/node and 4 GPUS/node
      integer, parameter             :: columns=columns_global/total_pes
      integer, parameter             :: left=100, right=101

      !usual mpi variables
      integer                        :: mype, npes, ierr
      integer                        :: status(MPI_STATUS_SIZE)

      double precision, parameter    :: max_temp_error=0.01

      integer                        :: i, j, max_iterations, iteration=1
      double precision               :: dt, dt_global=100.0
      double precision               :: start_time, stop_time

      ! GPU settings and non-blocking mpi settings
      integer                        ::ngpus,mydevice
      integer                        ::request(4)

      double precision, dimension(0:rows+1,0:columns+1) :: temperature, temperature_last
      double precision::checkpoint


      !usual mpi startup routines
      call MPI_Init(ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, npes, ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, mype, ierr)

      !Get GPU/node and bind to 1 GPU/MPI process
      ngpus=acc_get_num_devices(acc_device_nvidia)
      mydevice=mod(mype,ngpus)
      call acc_set_device_num(mydevice,acc_device_nvidia)


      !It is nice to verify that proper number of PEs are running
      if ( npes /= total_pes ) then
         if( mype == 0 ) then
            print *,'This example is hardwired to run only on ', total_pes, ' PEs'
         endif
         !call MPI_Finalize(ierr)
         !stop
      endif

      !Only one PE should prompt user
      if( mype == 0 ) then
         max_iterations = 4000
         print*, 'Maximum iterations = ', max_iterations
      endif
      
      !Other PEs need to recieve this information
      call MPI_Bcast(max_iterations, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      start_time = MPI_Wtime()

      call initialize(temperature_last, npes, mype)

      !$acc data copy(temperature_last) create(temperature)
      !do until global error is minimal or until maximum steps
      do while ( dt_global > max_temp_error .and. iteration <= max_iterations)

         !$acc kernels 
         do i=1,rows
            do j=1,columns
               temperature(i,j)=0.25*(temperature_last(i+1,j)+temperature_last(i-1,j)+ &
                                      temperature_last(i,j+1)+temperature_last(i,j-1) )
            enddo
         enddo
         !$acc end kernels



         !send data
         if (mype < npes-1) then
         !$acc update host(temperature(1:rows,columns))
            call MPI_Isend(temperature(1,columns), rows, MPI_DOUBLE_PRECISION, &
                          mype+1, RIGHT, MPI_COMM_WORLD, request(1),ierr)
         endif
         if (mype /= 0) then
         !$acc update host(temperature(1:rows,1))
            call MPI_ISend(temperature(1,1), rows, MPI_DOUBLE_PRECISION, &
                          mype-1, LEFT, MPI_COMM_WORLD, request(2),ierr)
         endif

         !receive data
         if (mype /= 0) then
            call MPI_IRecv(temperature_last(1,0), rows, MPI_DOUBLE_PRECISION, &
                          MPI_ANY_SOURCE, RIGHT, MPI_COMM_WORLD, request(3), ierr)
            !$acc update device(temperature_last(1:rows,0))
         endif
         if (mype /= npes-1) then
            call MPI_IRecv(temperature_last(1,columns+1), rows, MPI_DOUBLE_PRECISION, &
                          MPI_ANY_SOURCE, LEFT, MPI_COMM_WORLD, request(4), ierr)
            !$acc update device(temperature_last(1:rows,columns+1))
         endif

         ! Wait to complete data communication
         do i=1,4
              call mpi_wait(request(i),status,ierr)
         end do 


         dt=0.0

         !$acc kernels loop reduction(max:dt)
         do i=1,rows
            do j=1,columns
               dt = max( abs(temperature(i,j) - temperature_last(i,j)), dt )
               temperature_last(i,j) = temperature(i,j)
            enddo
         enddo
         !$acc end kernels

         !Need to determine and communicate maximum error
         call MPI_Reduce(dt, dt_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(dt_global, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

         !periodically print test values - only for PE in lower corner
         if( mod(iteration,100).eq.0 ) then
            if( mype == npes-1 ) then
               !$acc update host(temperature(rows-5:rows,columns-5:columns))
               call track_progress(temperature, iteration)
            endif
         endif

         iteration = iteration+1

      enddo

      if (mype==11)then
          checkpoint=temperature(9950,625)    ! (9950,625) at 11 th process corressponds to (9950,7750) at global array
          print*,'Temperature at checkpoint point',checkpoint
      end if 
      !$acc end data


      !Slightly more accurate timing and cleaner output
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      stop_time = MPI_Wtime()

      if( mype == 0 ) then
         print*, 'Max error at iteration ', iteration-1, ' was ',dt_global
         print*, 'Total time was ',stop_time-start_time, ' seconds.'
      endif

      call MPI_Finalize(ierr)

end program mpilaplace

!Parallel version requires more attention to global coordinates
subroutine initialize(temperature_last, npes, mype )
      implicit none

      integer, parameter             :: columns_global=10000
      integer, parameter             :: rows=10000
      integer, parameter             :: total_pes=16
      integer, parameter             :: columns=columns_global/total_pes

      integer                        :: i,j
      integer                        :: npes,mype

      double precision, dimension(0:rows+1,0:columns+1) :: temperature_last
      double precision               :: tmin, tmax

      temperature_last = 0

      !Left and Right Boundaries
      if( mype == 0 ) then
         do i=0,rows+1
            temperature_last(i,0) = 0.0
         enddo
      endif
      if( mype == npes-1 ) then
         do i=0,rows+1
            temperature_last(i,columns+1) = (100.0/rows) * i
         enddo
      endif

      !Top and Bottom Boundaries
      tmin =  mype    * 100.0/npes
      tmax = (mype+1) * 100.0/npes
      do j=0,columns+1
         temperature_last(0,j) = 0.0
         temperature_last(rows+1,j) = tmin + ((tmax-tmin)/columns) * j
      enddo

end subroutine initialize

subroutine track_progress(temperature, iteration)
      implicit none

      integer, parameter             :: columns_global=10000
      integer, parameter             :: rows=10000
      integer, parameter             :: total_pes=16
      integer, parameter             :: columns=columns_global/total_pes

      integer                        :: i,iteration

      double precision, dimension(0:rows+1,0:columns+1) :: temperature

!Parallel version uses global coordinate output so users don't need
!to understand decomposition
      print *, '---------- Iteration number: ', iteration, ' ---------------'
      do i=5,0,-1
         write (*,'("("i4,",",i4,"):",f6.2,"  ")',advance='no'), &
                   rows-i,columns_global-i,temperature(rows-i,columns-i)
      enddo
      print *
end subroutine track_progress
