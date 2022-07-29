Module Settings

IMPLICIT NONE

PRIVATE

PUBLIC run_settings!, define_local_settings


type run_settings
    !npx y number processors
    ! ndx y size local dimension
    integer :: npx, npy, npz, ndx, ndy, ndz, steps, iterations, &
               posx, posy, posz, rank_up, rank_down, rank_left, rank_right, & !mpi stuff
               gndx, gndy, gndz !global x and y dimensions and processor numbers

    ! variable grid
    ! assume variables are 2d and have same grid for now (x by y) - but need to calculate offsets
    ! now add 3d vars - but mpi distribution still only in x and y

    ! file names

    contains
    procedure :: define_local_settings


end type run_settings

contains

Subroutine define_local_settings(self, my_rank)!, num_ranks)

  ! define local and global dimensions, steps, mpi stuff etc
  ! pass in command line args and parse them

  class(run_settings), intent (INOUT) :: self
  character(len=12), dimension(:), allocatable :: args
  integer, intent(IN) :: my_rank!, num_ranks
  integer :: itemp, num_args, i

  write(0,*) "Defining local settings"

  self%ndx                  = 10 !n points in x
  self%ndy                  = 5
  self%ndz                  = 12
  self%steps                = 10 !steps of sim
  self%iterations           = 1 !n times apply operator per step

  num_args = command_argument_count()
  allocate(args(num_args))
  do i = 1, num_args
     call get_command_argument(i, args(i))
  end do

  !Use ridiculous fortran 'internal IO'???
  read( args(1), '(i6)' ) itemp
  self%npx       = itemp
  read( args(2), '(i6)' ) itemp
  self%npy       = itemp
  self%npz       = 1 ! do not split along z


  !calculate global array size and the local offsets in that global space
  self%gndx = self%npx * self%ndx
  self%gndy = self%npy * self%ndy
  self%gndz = self%npz * self%ndz
  self%posx = MODULO(my_rank,  self%npx)
  self%posy = my_rank / self%npx
  self%posz = 1

  ! determine neighbours
  if (self%posy .eq. 0) THEN
      self%rank_down = -1
  else
      self%rank_down = my_rank - 1 * self%npx
  end if

  if (self%posy .eq. (self%npy - 1)) THEN
      self%rank_up = -1
  else
      self%rank_up = my_rank + 1 * self%npx
  end if

  if (self%posx .eq. 0) THEN
      self%rank_left = -1
  else
      self%rank_left = my_rank - 1
  end if

  if (self%posx .eq. (self%npx - 1)) THEN
      self%rank_right = -1
  else
      self%rank_right = my_rank + 1
  end if


End Subroutine define_local_settings


End Module Settings
