program CMTbyAmp
  use cls_radiation, only: radiation
  use cls_observation, only: observation
  implicit none

  ! Note the coordinate definition. For f-net CMT solution,
  ! the x-axis is pointing to the north, the y-axis is pointing to
  ! the east, and the z-axis is pointing downward.
  double precision :: m_xx = -4.7005d0, m_yy = -0.6250d0, m_zz = 5.3255d0
  double precision :: m_xy = 2.6538d0, m_xz = 2.5805d0, m_yz = -2.0563d0
  double precision :: phi, theta, amp(3)
  double precision :: easting, northing
  double precision, parameter :: pi = acos(-1.0d0)
  integer :: i, j
  type(radiation) :: rad
  type(observation) :: obs
  character(*), parameter :: pol_file = "/home/akuhara/Develop/CMTbyAmp/Data/obs_less.dat"
  character(*), parameter :: sta_file = "/home/akuhara/Develop/CMTbyAmp/Data/sta.dat"


  rad = radiation(m_xx, m_yy, m_zz, m_xy, m_xz, m_yz)

  obs = observation(sta_file = sta_file, pol_file=pol_file)

  do i = 1, 360
     phi = i * pi / 180.d0
     do j = 1, 90
        theta = dble(j) * pi / 180.d0
        call rad%calc_radiation(phi, theta, "P", amp)
        ! In case of F-net definition, the x-axis is pointing to the north
        call rad%calc_projection_point(phi, theta, northing, easting)

        write(111,*) easting, northing, amp(1), amp(2), amp(3)
     end do
  end do
  
end program CMTbyAmp
  
