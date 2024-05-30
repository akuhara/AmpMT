program CMTbyAmp
  use mod_mpi
  use mod_random, only: init_random
  use cls_radiation, only: radiation
  use cls_observation, only: observation
  use cls_moment, only: moment
  use cls_model, only: model
  
  implicit none
  double precision, parameter :: pi = acos(-1.0d0)

  ! Note the coordinate definition. For f-net CMT solution,
  ! the x-axis is pointing to the north, the y-axis is pointing to
  ! the east, and the z-axis is pointing downward.
  double precision :: m_xx = -4.7005d0, m_yy = -0.6250d0, m_zz = 5.3255d0
  double precision :: m_xy = 2.6538d0, m_xz = 2.5805d0, m_yz = -2.0563d0
  !double precision :: u =  pi * 3.d0 / 8.d0, v = -1.d0 / 9.d0, k = 4.d0 * pi / 5.d0
  !double precision :: s = - pi / 2.d0, h = 3.d0 / 4.d0
  
  double precision :: phi, theta, amp(3)
  double precision :: easting, northing

  double precision, allocatable :: azi(:,:), inc(:,:)
  character(1), allocatable :: pol(:,:)
  character(16), allocatable :: stations(:)
  integer :: i, j, n_sta, n_evt
  integer :: ierr, rank, n_procs
  type(radiation) :: rad
  type(observation) :: obs
  type(moment), allocatable :: mt(:)
  type(model) :: u, v, k, s, h
  character(*), parameter :: pol_file = "/home/akuhara/Develop/CMTbyAmp/Data/obs_less.dat"
  character(*), parameter :: sta_file = "/home/akuhara/Develop/CMTbyAmp/Data/station.list"


  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call mpi_comm_size(mpi_comm_world, n_procs, ierr)

  call init_random(14141424, 22341, 100984, 92842433, rank)
  
  obs = observation(sta_file = sta_file, pol_file=pol_file)
  n_sta = obs%get_n_stations()
  n_evt = obs%get_n_events()

  azi = obs%get_azi()
  inc = obs%get_inc()
  pol = obs%get_pol()
  stations = obs%get_stations()

  do i = 1, n_evt
     print *, "Event ", i, "-----------------------------------"
     do j = 1, n_sta
        if (pol(j,i) /= "-") then
           print *, stations(j), azi(j,i), inc(j,i), pol(j,i)
        end if
     end do
  end do

  ! Initialize moment
  allocate(mt(n_evt))
  mt(1) = moment(u=0.d0, v=0.d0, k=0.d0, s=0.d0, h=0.d0) !dummy
  
  ! set model parameters 
  u = model(nx=n_evt)
  v = model(nx=n_evt)
  k = model(nx=n_evt)
  s = model(nx=n_evt)
  h = model(nx=n_evt)
  do i = 1, n_evt
     call u%set_prior(i, mt(1)%get_u_min(), mt(1)%get_u_max(), 1)
     call v%set_prior(i, mt(1)%get_v_min(), mt(1)%get_v_max(), 1)
     call k%set_prior(i, mt(1)%get_k_min(), mt(1)%get_k_max(), 1)
     call s%set_prior(i, mt(1)%get_s_min(), mt(1)%get_s_max(), 1)
     call h%set_prior(i, mt(1)%get_h_min(), mt(1)%get_h_max(), 1)
  end do
  call u%generate_model()
  call v%generate_model()
  call k%generate_model()
  call s%generate_model()
  call h%generate_model()
     
  
  rad = radiation(m_xx, m_yy, m_zz, m_xy, m_xz, m_yz)
  
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



  call mpi_finalize(ierr)
  
end program CMTbyAmp
  
