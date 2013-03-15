program main
! Daniel Buscombe Aug 2011
  implicit none

  character ( len = 32 ) :: input1
  character ( len = 32 ) :: input2
  character ( len = 32 ) :: input3
  character ( len = 32 ) :: input4
  integer ( kind = 4 ) numgrains
  integer ( kind = 4 ) model
  integer ( kind = 4 ) batch
  integer ( kind = 4 ) it_max
  character ( len = 80 ) :: file_out_name

  call timestamp ( )

  call get_command_argument(1, input1)
  read(input1,*) numgrains

  call get_command_argument(2, input2)
  read(input2,*) model

  call get_command_argument(3, input3)
  read(input3,*) batch

  call get_command_argument(4, input4)
  read(input4,*) it_max

  call get_command_argument(5, file_out_name)

  call gen_points (numgrains, model, batch, it_max, file_out_name)
  call timestamp ( )

  stop
end


subroutine gen_points (numgrains, model, batch, it_max, file_out_name)

!*****************************************************************************
! generates points in 3D.

  ! passed variables
  character*(*) file_out_name
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) batch
  integer ( kind = 4 ) model

  ! internal variables
  integer ( kind = 4 ), parameter :: dim_num = 3 
  real ( kind = 8 ) r(dim_num,numgrains)
  real ( kind = 8 ) energy
  character ( len = 80 ) init_string
  real ( kind = 8 ) it_diff
  integer ( kind = 4 ) it_fixed
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) sample_num
  character ( len = 80 ) sample_string
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_init

  !write ( *, * ) '  Compute points in 3D.'
  !write ( *, *) 'number of grains: '
  !write ( *, *) numgrains
  !write ( *, *) 'model: '
  !if (model == 1) then
  !   write ( *, *) 'pvt'
  !else if (model == 2) then
  !   write ( *, *) 'cvt'
  !else
  !   write ( *, *) 'halton'
  !end if
  !write ( *, *) '            '
  !write ( *, *) 'batch size: '
  !write ( *, *) batch
  !write ( *, *) 'it_max: '
  !write ( *, *) it_max
  !write ( *, *) 'output file '
  !write ( *, *) file_out_name

  init_string = 'uniform'
  it_fixed = 1
  sample_num = 10000
  sample_string = 'uniform'
  seed = 123456789
  seed_init = seed

  model = model-2

  call cvt ( dim_num, numgrains, batch, model, model, sample_num, it_max, it_fixed, &
    seed, r, it_num, it_diff, energy )

  call r8mat_write ( file_out_name, dim_num, numgrains, r )

  return
end

