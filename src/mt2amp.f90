program mt2amp
  use cls_radiation, only: radiation
  implicit none
  type(radiation) :: rad
  double precision :: mrr, mtt, mff, mrt, mrf, mtf
  double precision :: mxx, myy, mzz, mxy, mxz, myz
  double precision :: azi, inc
  double precision, parameter :: d_theta = 2.d0
  double precision, parameter :: d_phi = 2.d0
  integer, parameter :: n_theta = nint(90.d0 / d_theta)
  integer, parameter :: n_phi = nint(360.d0 / d_phi)
  double precision, parameter :: pi = acos(-1.d0)
  double precision :: theta, phi
  double precision :: amp(3)
  double precision :: northing, easting
  integer :: i, j, n_args
  logical :: output_all = .true.

  n_args = command_argument_count()
  
  if (n_args /= 6 .and. n_args /= 8) then
     print *, 'Usage: mt2amp [Mrr] [Mtt] [Mff] [Mrt] [Mrf] [Mtf] ([azi] [inc])'
     stop
  end if



  block
    character(len=256) :: arg
    call get_command_argument(1, arg)
    read(arg, *) mrr
    call get_command_argument(2, arg)
    read(arg, *) mtt
    call get_command_argument(3, arg)
    read(arg, *) mff
    call get_command_argument(4, arg)
    read(arg, *) mrt
    call get_command_argument(5, arg)
    read(arg, *) mrf
    call get_command_argument(6, arg)
    read(arg, *) mtf

    if (n_args == 8) then
       output_all = .false.
       call get_command_argument(7, arg)
       read(arg, *) azi
       call get_command_argument(8, arg)
       read(arg, *) inc

       !if (inc > 90.d0) then
       !   azi = azi + 180.d0
       !else if (inc < 0.0 .or. inc > 180.0) then
       !   print *, 'Error: inc must be in the range [0, 180]'
       !   stop
       !end if
       
    end if
  end block

  mxx = mtt
  myy = mff
  mzz = mrr
  mxy = -mtf
  mxz = mrt
  myz = -mrf

  rad = radiation(mxx, myy, mzz, mxy, mxz, myz)

  if (output_all) then
     do i = 1, n_phi
        phi = (i - 0.5d0) * d_phi / 180.d0 * pi
        do j = 1, n_theta
           theta = (j - 0.5d0) * d_theta / 180.d0 * pi
           call rad%calc_radiation(phi, theta, "P", amp)
           call rad%calc_projection_point(phi, theta, northing, easting)
           
           print *, easting, northing, amp(1)
           
        end do
     end do
  else
     phi = azi / 180.d0 * pi
     theta = inc / 180.d0 * pi
     call rad%calc_radiation(phi, theta, "P", amp)
     print *, amp(1)
  end if
     
   stop
end program mt2amp

