module cls_forward
  use cls_observation, only: observation
  use cls_moment, only: moment
  use cls_radiation, only: radiation
  use cls_model, only: model
  implicit none

  type forward
     private

     type(observation) :: obs
     double precision :: m(3,3)
     integer :: n_stations
     integer :: n_events
     
   contains
     
     procedure :: forward_calculation_full => forward_forward_calculation_full
     
  end type forward
  
  interface forward
     module procedure init_forward
  end interface forward

contains

  !-----------------------------------------------------------------------

  type(forward) function init_forward(obs) result(self)
    type(observation), intent(in) :: obs

    self%obs = obs
    self%n_stations = obs%get_n_stations()
    self%n_events = obs%get_n_events()

  end function init_forward

  !-----------------------------------------------------------------------

  subroutine forward_forward_calculation_full(self, u, v, k, s, h, q, log_likelihood)
    class(forward), intent(inout) :: self
    type(model), intent(in) :: u, v, k, s, h, q
    double precision, intent(out) :: log_likelihood
    type(moment) :: mt
    double precision :: m(3,3)
    type(radiation) :: rad
    integer :: i_sta, i_evt
    double precision :: azi(self%n_stations, self%n_events)
    double precision :: inc(self%n_stations, self%n_events)
    character(1) :: pol(self%n_stations, self%n_events)
    double precision, parameter :: pi = acos(-1.d0)
    double precision :: amp(3)
    double precision :: tmp

    

    pol = self%obs%get_pol()
    azi = self%obs%get_azi()
    inc = self%obs%get_inc()

    log_likelihood = 0.d0
    
    ! Calculate the forward model
    
    do i_evt = 1, self%n_events
       mt = moment(u=u%get_x(i_evt), v=v%get_x(i_evt), k=k%get_x(i_evt), &
            s=s%get_x(i_evt), h=h%get_x(i_evt))
       self%m = mt%get_moment()
       rad = radiation(m_xx=self%m(1,1), m_yy=self%m(2,2), m_zz=self%m(3,3), &
            m_xy=self%m(1,2), m_xz=self%m(1,3), m_yz=self%m(2,3))
       
       
       do i_sta = 1, self%n_stations
          if (pol(i_sta, i_evt) == '-') cycle
          call rad%calc_radiation(phi = azi(i_sta, i_evt) * pi / 180.d0, &
               theta = inc(i_sta, i_evt) * pi / 180.d0, &
               phase = "P", amp=amp)

          tmp = q%get_x((i_evt-1)*self%n_stations + i_sta)
          if (amp(1) > 0.d0 .and. pol(i_sta, i_evt) == "U") then
             log_likelihood = log_likelihood + &
                  & log(tmp)
          else if (amp(1) < 0.d0 .and. pol(i_sta, i_evt) == "D") then
             log_likelihood = log_likelihood + &
                  & log(tmp)
          else
             log_likelihood = log_likelihood + log(1.d0 - tmp)
          end if
       end do
    end do
    
  end subroutine forward_forward_calculation_full
end module cls_forward
