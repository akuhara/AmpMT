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
     logical :: sample_prior
     logical :: use_amp
     double precision :: amp_ignore = 100.d0
     
     
   contains
     
     procedure :: forward_calculation_full => forward_forward_calculation_full
     
  end type forward
  
  interface forward
     module procedure init_forward
  end interface forward

contains

  !-----------------------------------------------------------------------

  type(forward) function init_forward(obs, sample_prior, use_amp) result(self)
    type(observation), intent(in) :: obs
    logical, intent(in) :: sample_prior, use_amp

    self%obs = obs
    self%sample_prior = sample_prior
    self%use_amp = use_amp
    self%n_stations = obs%get_n_stations()
    self%n_events = obs%get_n_events()
    
    
  end function init_forward

  !-----------------------------------------------------------------------

  subroutine forward_forward_calculation_full(self, u, v, k, s, h, q, &
       & log_likelihood, a_mle)
    class(forward), intent(inout) :: self
    type(model), intent(in) :: u, v, k, s, h, q
    double precision, intent(out) :: log_likelihood
    double precision, intent(out) :: a_mle
    type(moment) :: mt
    double precision :: m(3,3)
    type(radiation) :: rad
    integer :: i_sta, i_evt
    double precision :: azi(self%n_stations, self%n_events)
    double precision :: inc(self%n_stations, self%n_events)
    character(1) :: pol(self%n_stations, self%n_events)
    double precision :: amp_obs(self%n_stations, self%n_events)
    double precision, parameter :: pi = acos(-1.d0)
    double precision :: amp(3)
    double precision :: tmp
    
    

    integer :: n_amp
    double precision :: oo, ss, os

    pol = self%obs%get_pol()
    azi = self%obs%get_azi()
    inc = self%obs%get_inc()
    amp_obs = self%obs%get_amp()

    log_likelihood = 0.d0
    if (self%sample_prior) return
    
    ! Calculate the forward model
    
    do i_evt = 1, self%n_events
       mt = moment(u=u%get_x(i_evt), v=v%get_x(i_evt), k=k%get_x(i_evt), &
            s=s%get_x(i_evt), h=h%get_x(i_evt))
       self%m = mt%get_moment()
       rad = radiation(m_xx=self%m(1,1), m_yy=self%m(2,2), m_zz=self%m(3,3), &
            m_xy=self%m(1,2), m_xz=self%m(1,3), m_yz=self%m(2,3))
       
       n_amp = 0
       os = 0.d0
       ss = 0.d0
       oo = 0.d0
       a_mle = 0.d0
       do i_sta = 1, self%n_stations
          if (pol(i_sta, i_evt) == '.') cycle
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

          if (self%use_amp) then
             if (abs(amp_obs(i_sta, i_evt)) > self%amp_ignore) cycle
             ss = ss + amp(1) * amp(1)
             oo = oo + amp_obs(i_sta, i_evt) * amp_obs(i_sta, i_evt)
             os = os + amp(1) * amp_obs(i_sta, i_evt)
             n_amp = n_amp + 1
          end if
       end do
       
       if (self%use_amp) then
          log_likelihood = log_likelihood - &
               & 0.5d0 * n_amp * log(1.d0 - os*os / (ss * oo))
          !log_likelihood = log_likelihood
          a_mle = os / ss
       end if
       
    end do
    
  end subroutine forward_forward_calculation_full
end module cls_forward
