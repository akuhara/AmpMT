!***********************************************************************

!  CMTbyAmp: A program to estimate CMT solutions using amplitude ratios

!* Note: 1st axis: North, 2nd axis: East, 3rd axis: Down
  !        This is different from Tape and Tape (2012) convention
  !        but consistent with the F-net definition
  !        Azimuth is measured clockwise from the North
  !        Inclination is measured from the vertical downward

  !      Lower hemisphere equal area projection is used (Schmidt net)
program cmtbyamp
  use mod_mpi
  use mod_random
  use cls_radiation, only: radiation
  use cls_observation, only: observation
  use cls_moment, only: moment
  use cls_model, only: model
  use cls_parallel, only: parallel
  use cls_mcmc, only: mcmc
  use cls_forward, only: forward
  use cls_output, only: output
  implicit none
  double precision, parameter :: pi = acos(-1.0d0)


  double precision :: m(3,3)
  double precision :: phi, theta, amp(3)
  double precision :: easting, northing

  double precision, allocatable :: azi(:,:), inc(:,:)
  character(1), allocatable :: pol(:,:)
  character(16), allocatable :: stations(:)
  integer :: i, j, n_sta, n_evt, i_sta, i_evt
  integer :: ierr, rank, n_procs
  integer :: n_chain = 5, n_iter = 200000, n_cool = 1, n_burn = 100000
  integer :: n_interval = 5000
  double precision :: temp_high = 300.d0, temp
  type(radiation) :: rad
  type(observation) :: obs
  type(moment), allocatable :: mt(:)
  type(model) :: u, v, k, s, h, q
  type(model) :: u_perturb, v_perturb, k_perturb, s_perturb, h_perturb, q_perturb
  type(parallel) :: pt
  type(mcmc) :: mc
  type(forward) :: fwd
  type(output), allocatable :: out(:)
  character(*), parameter :: pol_file = "/home/akuhara/Develop/CMTbyAmp/Data/obs_least.dat"
  character(*), parameter :: sta_file = "/home/akuhara/Develop/CMTbyAmp/Data/station.list"
  logical :: verb
  double precision, parameter :: eps = epsilon(1.d0)
  double precision :: log_prior_ratio, log_likelihood
  logical :: prior_ok
  integer :: id
  
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
  call mpi_comm_size(mpi_comm_world, n_procs, ierr)

  if (rank == 0) then
     verb = .true.
  else
     verb = .false.
  end if
  
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

  fwd = forward(obs=obs)

  allocate(out(n_evt))
  do i = 1, n_evt
     out(i) = output()
  end do

  
  ! Initialize moment
  allocate(mt(n_evt))
  mt(1) = moment(u=3.d0*pi/8.d0, v=-1.d0/9.d0, k=4.d0*pi/5.d0, &
       s=-pi/2.d0, h=3.d0/4.d0)
  m = mt(1)%get_moment()


  pt = parallel(n_proc = n_procs, rank = rank, n_chain = n_chain, verb = verb)
  do j = 1, n_chain
     ! set model parameters 
     u = model(nx=n_evt)
     v = model(nx=n_evt)
     k = model(nx=n_evt)
     s = model(nx=n_evt)
     h = model(nx=n_evt)
     q = model(nx=n_evt * n_sta)
     do i = 1, n_evt
        call u%set_prior(i, mt(1)%get_u_min(), mt(1)%get_u_max(), 1)
        call v%set_prior(i, mt(1)%get_v_min(), mt(1)%get_v_max(), 1)
        call k%set_prior(i, mt(1)%get_k_min(), mt(1)%get_k_max(), 1)
        call s%set_prior(i, mt(1)%get_s_min(), mt(1)%get_s_max(), 1)
        call h%set_prior(i, mt(1)%get_h_min(), mt(1)%get_h_max(), 1)
        do i_sta = 1, n_sta
           call q%set_prior((i-1)*n_sta + i_sta, 0.7d0, 1.d0, 1)
        end do
        
        call u%set_perturb(i, (mt(1)%get_u_max() - mt(1)%get_u_min()) / 50.d0)
        call v%set_perturb(i, (mt(1)%get_v_max() - mt(1)%get_v_min()) / 50.d0)
        call k%set_perturb(i, (mt(1)%get_k_max() - mt(1)%get_k_min()) / 50.d0)
        call s%set_perturb(i, (mt(1)%get_s_max() - mt(1)%get_s_min()) / 50.d0)
        call h%set_perturb(i, (mt(1)%get_h_max() - mt(1)%get_h_min()) / 50.d0)
        do i_sta = 1, n_sta
           call q%set_perturb((i-1)*n_sta + i_sta, 1.d0 / 50.d0)
        end do
     end do
     call u%generate_model()
     call v%generate_model()
     call k%generate_model()
     call s%generate_model()
     call h%generate_model()
     call q%generate_model()

     if (verb) then
        print *, "Chain ", j
        print *, "u: ", u%get_all_x()
        print *, "v: ", v%get_all_x()
        print *, "k: ", k%get_all_x()
        print *, "s: ", s%get_all_x()
        print *, "h: ", h%get_all_x()
     end if


     mc = mcmc(&
          u = u, &
          v = v, &
          k = k, &
          s = s, &
          h = h, &
          q = q, &
          n_iter = n_iter)


     
     ! Set temperatures
     if (j <= n_cool) then
        call mc%set_temp(1.d0)
     else
        temp = exp((rand_u() * (1.d0 - eps) + eps) * log(temp_high))

        call mc%set_temp(temp)
     end if
     call pt%set_mc(j, mc)
  end do
  
  ! Main MCMC loop

  do i = 1, n_iter

     if (mod(i, n_interval) == 0) then
        if (rank == 0) print *, "Iter ", i, " / ", n_iter
     end if
     do j = 1, n_chain
        mc = pt%get_mc(j)
        call mc%propose_model(u_perturb, v_perturb, k_perturb, &
             s_perturb, h_perturb, q_perturb, log_prior_ratio, prior_ok, id)


        if (prior_ok) then
           call fwd%forward_calculation_full(u_perturb, v_perturb, k_perturb, &
                s_perturb, h_perturb, q_perturb, log_likelihood)
        end if
        
        call mc%judge_model(u_perturb, v_perturb, k_perturb, &
             s_perturb, h_perturb, q_perturb, log_likelihood, log_prior_ratio, prior_ok)
        call pt%set_mc(j, mc)

        ! recording
        if (mc%get_temp() < 1.d0 + eps .and. mod(i, n_interval) == 1 .and. &
             i > n_burn) then
           !write(222,*) i, j, mc%write_out_q()
           !write(333,*) i, j, mc%write_out_u()
           !write(444,*) i, j, mc%write_out_v()

           u = mc%get_u()
           v = mc%get_v()
           k = mc%get_k()
           s = mc%get_s()
           h = mc%get_h()

           do i_evt = 1, n_evt
              mt(i_evt) = moment(u=u%get_x(i_evt), v=v%get_x(i_evt), &
                   k=k%get_x(i_evt), s=s%get_x(i_evt), h=h%get_x(i_evt))
              call out(i_evt)%count_pol(mt(i_evt))
              call out(i_evt)%count_source_type(mt(i_evt))
              call out(i_evt)%count_fault_type(mt(i_evt))
              call out(i_evt)%count_principal_axes(mt(i_evt))
              call out(i_evt)%save_likelihood(mc%get_log_likelihood())
              call out(i_evt)%save_moment(mt(i_evt))
           end do
           !write(555,*) i, j, (mt(i_evt)%get_strike(), i_evt = 1, n_evt)
           !write(666,*) i, j, (mt(i_evt)%get_dip(), i_evt = 1, n_evt)
           !write(777,*) i, j, (mt(i_evt)%get_rake(), i_evt = 1, n_evt)
           
        end if
     end do
     call pt%swap_temperature(verb=verb)
  end do


  block
    double precision, allocatable :: pol_cum_all(:,:)
    double precision, allocatable :: pol_cum(:,:)
    integer :: n_pol_cum, n_pol_cum_all
    integer :: dim(2)
    character(255) :: out_file
    
    
    if (rank==0) print *, "Writing output"
    
    do i = 1, n_evt

       ! polarization
       n_pol_cum = out(i)%get_n_pol_cum()
       pol_cum = out(i)%get_pol_cum()
       dim = shape(pol_cum)
       allocate(pol_cum_all(dim(1), dim(2)))
       
       call mpi_reduce(pol_cum, pol_cum_all, size(pol_cum), &
            mpi_double_precision, mpi_sum, 0, MPI_COMM_WORLD, ierr)
       call mpi_reduce(n_pol_cum, n_pol_cum_all, 1, &
            mpi_integer, mpi_sum, 0, MPI_COMM_WORLD, ierr)
       
       
       if (rank == 0) then
          call out(i)%set_pol_cum(pol_cum_all)
          call out(i)%set_n_pol_cum(n_pol_cum_all)
          write(out_file, "(A, I5.5, A)") "pol_", i, ".dat"
          call out(i)%write_pol(out_file)
          
       end if

       deallocate(pol_cum_all)

       ! source type

       write(out_file, "(A, I5.5, A)") "source_", i, ".dat"

       if (rank == 0) then
          call out(i)%write_source_type(out_file, "replace")
       end if
       call mpi_barrier(MPI_COMM_WORLD, ierr)
       
       do j = 1, n_procs-1
          if (rank == j) then
             call out(i)%write_source_type(out_file, "append")
          end if
          call mpi_barrier(MPI_COMM_WORLD, ierr)
       end do

       ! fault type
       write(out_file, "(A, I5.5, A)") "fault_", i, ".dat"

       if (rank == 0) then
          call out(i)%write_fault_type(out_file, "replace")
       end if
       call mpi_barrier(MPI_COMM_WORLD, ierr)

       do j = 1, n_procs-1
          if (rank == j) then
             call out(i)%write_fault_type(out_file, "append")
          end if
          call mpi_barrier(MPI_COMM_WORLD, ierr)
       end do

       ! principal axes
       write(out_file, "(A, I5.5, A)") "axes_", i, ".dat"

       if (rank == 0) then
          call out(i)%write_principal_axes(out_file, "replace")
       end if
       call mpi_barrier(MPI_COMM_WORLD, ierr)

       do j = 1, n_procs-1
          if (rank == j) then
             call out(i)%write_principal_axes(out_file, "append")
          end if
          call mpi_barrier(MPI_COMM_WORLD, ierr)
       end do

       !summary

       write(out_file, "(A, I5.5, A)") "summary_", i, ".dat"

       if (rank == 0) then
          call out(i)%write_summary(out_file, "replace")
       end if
       call mpi_barrier(MPI_COMM_WORLD, ierr)

       do j = 1, n_procs-1
          if (rank == j) then
             call out(i)%write_summary(out_file, "append")
          end if
          call mpi_barrier(MPI_COMM_WORLD, ierr)
       end do
    end do
  end block
    
!  rad = radiation(m_xx, m_yy, m_zz, m_xy, m_xz, m_yz)
!  
!  do i = 1, 360
!     phi = i * pi / 180.d0
!     do j = 1, 90
!        theta = dble(j) * pi / 180.d0
!        call rad%calc_radiation(phi, theta, "P", amp)
!        ! In case of F-net definition, the x-axis is pointing to the north
!        call rad%calc_projection_point(phi, theta, northing, easting)
!
!        write(111,*) easting, northing, amp(1), amp(2), amp(3)
!       
!     end do
!  end do



  call mpi_finalize(ierr)
  
end program cmtbyamp
  
