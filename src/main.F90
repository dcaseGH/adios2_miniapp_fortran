program main

  ! test miniapp for adios2
  ! hard set local dimensions ndx, ndy - global calculated from these (easier than odd numbers per local domain)
  !

  use mpi
  use adios2
  use settings, only : run_settings!, define_local_settings
  use heat_transfer, only : data_object, apply_diffusion, exchange, initialise, &
                            clear_data, apply_heat, apply_column_addition

  implicit none

  integer :: ierr, num_ranks, my_rank, &
             t, iter, &
             num_args, i, j, k, itemp

  ! stuff to run calc
  TYPE(run_settings) :: local_settings
  TYPE(data_object)  :: local_data
  double precision   :: edge_temp

  !in built profile
  integer(8)   :: time_start, time_end, system_count_rate

  !declare adios variables
  type(adios2_adios) :: adios
  type(adios2_io) :: ioPut, ioGet
  type(adios2_variable) :: temperature, temperature_full, var_3d
  type(adios2_engine) :: bpWriter, hdf5Writer!, hdf5Reader !Do I want a reader?
  integer(kind=8), dimension(2) :: ishape, istart, icount
  integer(kind=8), dimension(3) :: ijkshape, ijkstart, ijkcount

  call MPI_INIT(ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, my_rank, ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, num_ranks, ierr)
  if (my_rank .eq. 0) THEN
        write(0,*) 'Running an ADIOS2 miniapp - eg command mpirun -np 2 [exe] 1 2'
        ! Timing
        call system_clock(count_rate=system_count_rate)
        call system_clock(time_start)
  end if

  write(0,*) "Running with rank ", my_rank, " of ", num_ranks
  call local_settings%define_local_settings(my_rank)

  call initialise(local_settings, local_data)
  write(0,*) "x and y distribution = ", local_settings%npx, " ", local_settings%npy
  if (local_settings%npx * local_settings%npy .ne. num_ranks) then
     stop 'pass arguments so that decomposition in x * y = number of ranks'
  end if

  ! For now set simple initial data (rank dependent)
  local_data%Temp = dble(my_rank)

  call adios2_init( adios, MPI_COMM_WORLD, ierr )
!  call adios2_init( adios, 'example.xml', MPI_COMM_WORLD, ierr )
  call adios2_declare_io( ioPut, adios, 'TempWrite', ierr )
     !declare hdf5 engine: try sst at some point? or compare to bp?
  ! note if container not set up properly, cant do this!
#ifdef HDF5
  call adios2_set_engine(ioPut, 'HDF5', ierr)
#endif
  ! also set params - threads etc??

  !adios variables have global shape and local start and count
  icount = (/ local_settings%ndx                       , local_settings%ndy     /)
  istart = (/ local_settings%posx * local_settings%ndx , local_settings%posy * local_settings%ndy  /)
  ishape = (/ local_settings%npx  * local_settings%ndx , local_settings%npy  * local_settings%ndy  /)
  call adios2_define_variable( temperature, ioPut, 'temperatures', &
                               adios2_type_dp, 2, &
                               ishape, istart, icount, adios2_constant_dims, &
                               ierr )

  icount = (/ local_settings%ndx+2                         , local_settings%ndy+2     /)
  istart = (/ local_settings%posx * (local_settings%ndx+2) , local_settings%posy * (local_settings%ndy+2)  /)
  ishape = (/ local_settings%npx  * (local_settings%ndx+2) , local_settings%npy  * (local_settings%ndy+2)  /)

  call adios2_define_variable( temperature_full, ioPut, 'temperatures_full', &
                               adios2_type_dp, 2, &
                               ishape, istart, icount, adios2_constant_dims, &
                               ierr )

  ijkcount = (/ local_settings%ndx, local_settings%ndy, local_settings%ndz /)
  ijkstart = (/ local_settings%posx * local_settings%ndx , local_settings%posy * local_settings%ndy, 0  /)
  ijkshape = (/ local_settings%npx  * local_settings%ndx , local_settings%npy  * local_settings%ndy, &
                local_settings%npz  * local_settings%ndz  /)

  call adios2_define_variable( var_3d, ioPut, 'something_3d', &
                               adios2_type_dp, 3, &
                               ijkshape, ijkstart, ijkcount, adios2_constant_dims, &
                               ierr )

#ifdef HDF5
  if (my_rank .eq. 0) write(0,*) "Writing output in HDF5 format"
  call adios2_open( bpWriter, ioPut, 'initial_dat.hdf5', adios2_mode_write, &
                    ierr )
#else
  if (my_rank .eq. 0) write(0,*) "Writing output in BP format"
  call adios2_open( bpWriter, ioPut, 'initial_dat.bp', adios2_mode_write, &
                    ierr )
#endif
  call adios2_begin_step(bpWriter, ierr)
  call adios2_put( bpWriter, temperature,      local_data%Temp,     ierr )
  call adios2_put( bpWriter, temperature_full, local_data%Temp,     ierr )
  call adios2_put( bpWriter, var_3d,           local_data%field_3d, ierr )
  call adios2_end_step(bpWriter, ierr)

  edge_temp = 4.8 !set edges
  ! do more work and write more data
  do t = 1,  local_settings%steps

      call adios2_begin_step(bpWriter, ierr)

      if (my_rank .eq. 0) write(0,*) "Step ", t
      ! increase this to do more work/comms per write
      do iter = 1, 1! local_settings%iterations

         ! operators
         ! 2D
         call apply_diffusion(local_settings, local_data)

         ! 3D
         call apply_column_addition(local_settings, local_data)

         ! mpi
         call exchange(local_settings, local_data)

         ! BC
         call apply_heat(edge_temp, local_settings, local_data)

      end do !iter

      ! Write something and increment step
      call adios2_put( bpWriter, temperature, local_data%Temp, ierr )
      call adios2_put( bpWriter, temperature_full, local_data%Temp, ierr )
      call adios2_put( bpWriter, var_3d, local_data%field_3d, ierr )
      call adios2_end_step(bpWriter, ierr)
  end do !t

  !delete objects
  call clear_data(local_data)

  call adios2_close( bpWriter, ierr )

  if (my_rank .eq. 0) THEN
    call system_clock(time_end)
    write(0,*) 'Wallclock time (s): ', ( dble(time_end) - dble(time_start) )/dble(system_count_rate)
  end if
  call MPI_FINALIZE ( ierr )

end program
