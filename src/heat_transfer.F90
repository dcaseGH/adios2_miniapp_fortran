Module heat_transfer

  use settings, only : run_settings
  use mpi

IMPLICIT NONE

PRIVATE

PUBLIC data_object, apply_diffusion, initialise, clear_data, exchange, apply_heat, apply_column_addition

TYPE data_object
  double precision, Allocatable, Dimension(:,:) :: Temp, temp_Temp
  double precision, Allocatable, Dimension(:,:,:) :: field_3d
  double precision :: omega

end type data_object

CONTAINS

SUBROUTINE apply_diffusion(local_settings, local_data)

   ! Numerical gradient

   TYPE(run_settings), intent(IN)    :: local_settings
   TYPE(data_object), intent(INOUT)  :: local_data

   integer :: i, j

   local_data%temp_Temp = 0.

   do j = 1, local_settings%ndy
      do i = 1, local_settings%ndx

         local_data%temp_Temp(i,j) = (local_data%omega / 4) * ( local_data%temp(i+1,j) + &
                                                                local_data%temp(i-1,j) + &
                                                                local_data%temp(i,j+1) + &
                                                                local_data%temp(i,j-1) ) &
                                   + (1. - local_data%omega) *  local_data%temp(i,j)
      end do !i
   end do !j

   local_data%temp = local_data%temp_Temp

END SUBROUTINE apply_diffusion

SUBROUTINE clear_data(local_data)

   !clear data arrays
   TYPE(data_object), intent(INOUT)  :: local_data

   DEALLOCATE(local_data%Temp)
   DEALLOCATE(local_data%temp_Temp)

   !3d
   DEALLOCATE(local_data%field_3d)

END SUBROUTINE clear_data

SUBROUTINE initialise(local_settings, local_data)

   ! just invent some data

   TYPE(run_settings), intent(IN)    :: local_settings
   TYPE(data_object), intent(INOUT)  :: local_data

   integer :: i, j
   double precision :: hx, hy, x, y

   local_data%omega = 0.8 !because it needs to be something...

   !2d
   ALLOCATE(local_data%Temp(0:local_settings%ndx+1,0:local_settings%ndy+1))
   ALLOCATE(local_data%temp_Temp(0:local_settings%ndx+1,0:local_settings%ndy+1))

   !3d
   ALLOCATE(local_data%field_3d(local_settings%ndx, local_settings%ndy, local_settings%ndz))

   !2d
   local_data%temp      = 0.
   local_data%temp_Temp = 0.

   hx = 2.0 * 4.0 * atan(1.0) / local_settings%gndx
   hy = 2.0 * 4.0 * atan(1.0) / local_settings%gndy
   !write(0,*) 'Debug', local_settings%posy, shape(local_data%Temp), local_settings%ndx, local_settings%ndy

   do j=0, local_settings%ndy + 1
      y = 0.0 + hy * (j + local_settings%posy * local_settings%ndy - 1)

      do i=0, local_settings%ndx + 1
         x = 0. + hx * (i + local_settings%posx *local_settings%ndy - 1)
         local_data%Temp(i,j) = cos(8 * x) + cos(6 * x) - cos(4 * x) + cos(2 * x) &
                              - cos(x) + sin(8 * y) - sin(6 * y) + sin(4 * y) &
                              - sin(2 * y) + sin(y)
      end do
   end do

   !3D
   local_data%field_3d = 4.2

END SUBROUTINE initialise

SUBROUTINE exchange(local_settings, local_data)

   TYPE(run_settings), intent(IN)    :: local_settings
   TYPE(data_object), intent(INOUT)  :: local_data
   integer status(MPI_STATUS_SIZE)

   !mpi datatypes are integers in fortran?
   integer :: ierr, mpi_vector, tag

   call MPI_Type_vector(local_settings%ndy+2, 1, local_settings%ndx + 2, MPI_DOUBLE_PRECISION, mpi_vector, ierr)
   call MPI_Type_commit(mpi_vector, ierr)

    !Exchange ghost cells, in the order left-right-up-down

    ! send to left + receive from right
    tag = 1
    if (local_settings% rank_left .ge. 0) THEN
       call MPI_SEND( local_data%Temp(1,:), local_settings%ndy+2, MPI_DOUBLE_PRECISION, local_settings%rank_left, &
                      tag, MPI_COMM_WORLD, ierr )
    end if
    if (local_settings%rank_right .ge. 0) THEN
       call MPI_RECV(local_data%Temp(local_settings%ndx+1,:), local_settings%ndy+2 , MPI_DOUBLE_PRECISION, &
                     local_settings%rank_right, &
                     tag, MPI_COMM_WORLD, status, ierr )
    end if

    ! send to right + receive from left
    tag = 2
    if (local_settings%rank_right .ge. 0) THEN
       call MPI_SEND( local_data%Temp(local_settings%ndx,:), local_settings%ndy+2, MPI_DOUBLE_PRECISION, &
                      local_settings%rank_right, &
                      tag, MPI_COMM_WORLD, ierr )
    end if
    if (local_settings%rank_left .ge. 0) THEN
      call MPI_Recv(local_data%Temp(0,:), local_settings%ndy+2, MPI_DOUBLE_PRECISION, &
                     local_settings%rank_left, tag, MPI_COMM_WORLD, status, ierr)
    end if

    !send down + receive from above
    tag = 3
    if (local_settings%rank_down .ge. 0) THEN
       call MPI_SEND( local_data%Temp(:, 1),                   local_settings%ndx+2, MPI_DOUBLE_PRECISION, &
                      local_settings%rank_down, tag, MPI_COMM_WORLD, ierr)
    end if
    if (local_settings%rank_up .ge. 0) THEN
       call MPI_Recv(local_data%Temp(:, local_settings%ndy+1), local_settings%ndx+2, MPI_DOUBLE_PRECISION, &
                     local_settings%rank_up, tag, MPI_COMM_WORLD, status, ierr)
    end if

    !send up + receive from below
    tag = 4
    if (local_settings%rank_up .ge. 0) THEN
       call MPI_SEND( local_data%Temp(:, local_settings%ndy),  local_settings%ndx+2, MPI_DOUBLE_PRECISION, &
                      local_settings%rank_up, tag, MPI_COMM_WORLD, ierr)
    end if
    if (local_settings%rank_down .ge. 0) THEN
       call MPI_Recv(local_data%Temp(:, 0),                    local_settings%ndx+2, MPI_DOUBLE_PRECISION, &
                     local_settings%rank_down, tag, MPI_COMM_WORLD, status, ierr)
    end if

   call MPI_Type_free(mpi_vector, ierr)


END SUBROUTINE exchange

SUBROUTINE apply_column_addition(local_settings, local_data)

   ! Add something to the 3D field
   ! Add more to higher points

   TYPE(run_settings), intent(IN)    :: local_settings
   TYPE(data_object), intent(INOUT)  :: local_data
   Integer :: i

   do i=1, local_settings%ndz
      local_data%field_3d(:,:,i) = local_data%field_3d(:,:,i) + DBLE(i) / 10.00
   end do

END SUBROUTINE apply_column_addition

SUBROUTINE apply_heat(edge_temp, local_settings, local_data)

   ! Apply heat to set edge_temp of the cells outside the interested zone
   ! Want to find points which are globally on edges, i.e. = or num_ranks * ndx + 1

   double precision, intent(in) :: edge_temp
   TYPE(run_settings), intent(IN)    :: local_settings
   TYPE(data_object), intent(INOUT)  :: local_data

   !
   if (local_settings%posx .eq. 0) THEN
      local_data%Temp(0,:) = edge_temp
   end if

   if (local_settings%posx .eq. (local_settings%npx-1)) THEN
      local_data%Temp(local_settings%ndx+1,:) = edge_temp
   end if

   if (local_settings%posy .eq. 0) THEN
      local_data%Temp(:,0) = edge_temp
   end if

   if (local_settings%posy .eq. (local_settings%npy-1)) THEN
      local_data%Temp(:,local_settings%ndy+1) = edge_temp
   end if

END SUBROUTINE apply_heat

End Module heat_transfer

