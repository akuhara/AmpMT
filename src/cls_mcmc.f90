module cls_mcmc
  use mod_random
  use cls_model
  
  implicit none 
  
  type mcmc 
     private
     integer :: n_events, n_sta

     ! Model parameters
     type(model) :: u, v, k, s, h, q

     
     integer :: n_proposal_type = 6
     integer :: n_iter
     
     integer :: i_iter
     integer, allocatable :: n_propose(:)
     integer, allocatable :: n_accept(:)
     character(6), allocatable :: label(:)
     integer :: i_proposal_type
     double precision :: log_likelihood
     double precision :: temp = 1.d0
     double precision :: a_mle
     logical :: is_accepted

     logical :: dc_only = .false.

     double precision :: p_u, p_v, p_k, p_s, p_h, p_q
     
     
   contains
     procedure :: propose_model => mcmc_propose_model
     procedure :: judge_model   => mcmc_judge_model
     procedure :: one_step_summary => mcmc_one_step_summary
     procedure :: set_log_likelihood => mcmc_set_log_likelihood
     procedure :: set_temp => mcmc_set_temp
     procedure :: get_temp => mcmc_get_temp
     procedure :: set_a_mle => mcmc_set_a_mle
     procedure :: get_a_mle => mcmc_get_a_mle
     procedure :: get_log_likelihood => mcmc_get_log_likelihood
     procedure :: get_n_propose => mcmc_get_n_propose
     procedure :: get_n_accept => mcmc_get_n_accept
     procedure :: get_n_iter => mcmc_get_n_iter
     procedure :: get_u => mcmc_get_u
     procedure :: get_v => mcmc_get_v
     procedure :: get_k => mcmc_get_k
     procedure :: get_s => mcmc_get_s
     procedure :: get_h => mcmc_get_h
     procedure :: get_q => mcmc_get_q
     procedure :: get_label => mcmc_get_label

     procedure :: write_out_u => mcmc_write_out_u
     procedure :: write_out_v => mcmc_write_out_v
     procedure :: write_out_k => mcmc_write_out_k
     procedure :: write_out_s => mcmc_write_out_s
     procedure :: write_out_h => mcmc_write_out_h
     procedure :: write_out_q => mcmc_write_out_q
     
  end type mcmc
  
  interface mcmc
     module procedure init_mcmc
  end interface mcmc


contains

  !---------------------------------------------------------------------
  
  type(mcmc) function init_mcmc(u, v, k, s, h, q, n_iter, dc_only) result(self)
    type(model), intent(in) :: u, v, k, s, h, q
    integer, intent(in) :: n_iter
    logical, intent(in), optional :: dc_only

    if (present(dc_only)) then
       self%dc_only = dc_only
    end if
    
    self%u = u
    self%v = v
    self%k = k
    self%s = s
    self%h = h
    self%q = q

    self%n_events = u%get_nx()
    self%n_sta = q%get_nx() / self%n_events

    allocate(self%n_propose(self%n_proposal_type))
    allocate(self%n_accept(self%n_proposal_type))
    allocate(self%label(self%n_proposal_type))
    self%label = [ character(6) :: "u", "v", "k", "s", "h", "q" ]
    self%n_propose = 0
    self%n_accept  = 0
    self%n_iter    = n_iter
    self%i_iter    = 0
    self%log_likelihood = -9.d+300

    ! set proposal probability
    if (.not. dc_only) then
       self%p_u = 1.d0 / dble(self%n_sta + 5)
       self%p_v = 1.d0 / dble(self%n_sta + 5)
       self%p_k = 1.d0 / dble(self%n_sta + 5)
       self%p_s = 1.d0 / dble(self%n_sta + 5)
       self%p_h = 1.d0 / dble(self%n_sta + 5)
       self%p_q = dble(self%n_sta) / dble(self%n_sta + 5)
    else
       self%p_u = 0.d0
       self%p_v = 0.d0
       self%p_k = 1.d0 / dble(self%n_sta + 3)
       self%p_s = 1.d0 / dble(self%n_sta + 3)
       self%p_h = 1.d0 / dble(self%n_sta + 3)
       self%p_q = dble(self%n_sta) / dble(self%n_sta + 3)
    end if
    
    print *, self%p_u, self%p_v, self%p_k, self%p_s, self%p_h, self%p_q
    print *, "n_events = ", self%n_events
    print *, "n_sta = ", self%n_sta

    self%a_mle = 0.d0
    
    return 
  end function init_mcmc

  !---------------------------------------------------------------------
  
  subroutine mcmc_propose_model(self, u_proposed, v_proposed, &
       & k_proposed, s_proposed, h_proposed, q_proposed, &
       & log_prior_ratio, prior_ok, evt_id)
    class(mcmc), intent(inout) :: self
    type(model), intent(out) :: u_proposed, v_proposed, k_proposed, &
         s_proposed, h_proposed, q_proposed
    double precision, intent(out) :: log_prior_ratio
    integer, intent(out) :: evt_id
    integer :: itype, icmp, id, sta_id
    logical, intent(out) :: prior_ok
    double precision :: a_select


    u_proposed = self%u
    v_proposed = self%v
    k_proposed = self%k
    s_proposed = self%s
    h_proposed = self%h
    q_proposed = self%q
    
    a_select = rand_u()

    prior_ok = .true. 
    evt_id = -999

    evt_id = int(rand_u() * self%n_events) + 1
    if (a_select < self%p_u) then
       ! Perturb u
       call u_proposed%perturb(evt_id, log_prior_ratio, prior_ok)
       self%i_proposal_type = 1
    else if (a_select < self%p_u + self%p_v) then
       ! Perturb v
       call v_proposed%perturb(evt_id, log_prior_ratio, prior_ok)
       self%i_proposal_type = 2
    else if (a_select < self%p_u + self%p_v + self%p_k) then
       ! Perturb k
       call k_proposed%perturb(evt_id, log_prior_ratio, prior_ok)
       self%i_proposal_type = 3
    else if (a_select < self%p_u + self%p_v + self%p_k + self%p_s) then
       ! Perturb s
       call s_proposed%perturb(evt_id, log_prior_ratio, prior_ok)
       self%i_proposal_type = 4
    else if (a_select < &
         self%p_u + self%p_v + self%p_k + self%p_s + self%p_h) then
       ! Perturb h
       call h_proposed%perturb(evt_id, log_prior_ratio, prior_ok)
       self%i_proposal_type = 5
    else
       ! Perturb q
       sta_id = int(rand_u() * self%n_sta) + 1
       call q_proposed%perturb((evt_id-1) * self%n_sta + sta_id, &
            & log_prior_ratio, prior_ok)
       self%i_proposal_type = 6
    end if
    
    
    return
  end subroutine mcmc_propose_model

  !---------------------------------------------------------------------

  subroutine mcmc_judge_model(self, u, v, k, s, h, q, &
       & log_likelihood, a_mle, log_prior_ratio, prior_ok)
    class(mcmc), intent(inout) :: self
    type(model), intent(in) :: u, v, k, s, h, q
    double precision, intent(in) :: log_likelihood, log_prior_ratio
    double precision, intent(in) :: a_mle
    logical, intent(in) :: prior_ok
    double precision :: ratio
    double precision :: r
    double precision, parameter :: eps = epsilon(1.d0)
    
    if (self%temp < 1.d0 + eps) then
       self%n_propose(self%i_proposal_type) = &
            & self%n_propose(self%i_proposal_type) + 1
    end if
    

    self%is_accepted = .false.
    if (prior_ok) then
       ratio = (log_likelihood - self%log_likelihood) / self%temp
       ratio = ratio + log_prior_ratio 
       
       r = rand_u()
       if (r >= eps) then
          if (log(r) <= ratio) then
             self%is_accepted = .true.
          end if
       end if
    end if
    
    

    if (self%is_accepted) then
       ! Accept model
       self%u = u
       self%v = v
       self%k = k
       self%s = s
       self%h = h
       self%q = q
       self%log_likelihood = log_likelihood
       self%a_mle = a_mle
       if (self%temp < 1.d0 + eps) then
          self%n_accept(self%i_proposal_type) = &
               & self%n_accept(self%i_proposal_type) + 1
       end if
    end if

    ! Adds iteration counter
    self%i_iter = self%i_iter + 1


    return 
  end subroutine mcmc_judge_model
  
  !---------------------------------------------------------------------

  subroutine mcmc_one_step_summary(self)
    class(mcmc), intent(in) :: self
    write(*,*)"Iteration    : ", self%i_iter, "/", self%n_iter
    write(*,*)"Proposal type: ", self%i_proposal_type
    write(*,*)"Accepted     : ", self%is_accepted
    write(*,*)"Likelihood   : ", self%log_likelihood
    write(*,*)
  end subroutine mcmc_one_step_summary

  !---------------------------------------------------------------------

  subroutine mcmc_set_log_likelihood(self, log_likelihood)
    class(mcmc), intent(inout) :: self
    double precision, intent(in) :: log_likelihood

    self%log_likelihood = log_likelihood

    return 
  end subroutine mcmc_set_log_likelihood
  
  !---------------------------------------------------------------------
  
  subroutine mcmc_set_temp(self, temp)
    class(mcmc), intent(inout) :: self
    double precision, intent(in) :: temp

    self%temp = temp

    return 
  end subroutine mcmc_set_temp

  !---------------------------------------------------------------------

  double precision function mcmc_get_temp(self) result(temp)
    class(mcmc), intent(in) :: self

    temp = self%temp
    
    return 
  end function mcmc_get_temp

  !---------------------------------------------------------------------

  subroutine mcmc_set_a_mle(self, a_mle)
    class(mcmc), intent(inout) :: self
    double precision, intent(in) :: a_mle
    
    self%a_mle = a_mle
    
    return 
  end subroutine mcmc_set_a_mle
  
  !---------------------------------------------------------------------

  double precision function mcmc_get_a_mle(self) result(a_mle)
    class(mcmc), intent(in) :: self
    
    a_mle = self%a_mle
    
    return 
  end function mcmc_get_a_mle

  !---------------------------------------------------------------------

  double precision function mcmc_get_log_likelihood(self) result(l)
    class(mcmc), intent(in) :: self
    
    l = self%log_likelihood
    
    return 
  end function mcmc_get_log_likelihood
  
  !---------------------------------------------------------------------

  function mcmc_get_n_accept(self) result(n_accept)
    class(mcmc), intent(in) :: self
    integer :: n_accept(self%n_proposal_type)

    n_accept = self%n_accept

    return 
  end function mcmc_get_n_accept
  
  !---------------------------------------------------------------------

  function mcmc_get_n_propose(self) result(n_propose)
    class(mcmc), intent(in) :: self
    integer :: n_propose(self%n_proposal_type)

    n_propose = self%n_propose

    return 
  end function mcmc_get_n_propose
  
  !---------------------------------------------------------------------

  function mcmc_get_n_iter(self) result(n_iter)
    class(mcmc), intent(in) :: self
    integer :: n_iter

    n_iter = self%n_iter

    return 
  end function mcmc_get_n_iter
  
  !---------------------------------------------------------------------
  
  function mcmc_get_u(self) result(u)
    class(mcmc), intent(in) :: self
    type(model) :: u
    
    u = self%u
    
    return 
  end function mcmc_get_u

  !---------------------------------------------------------------------

  function mcmc_get_v(self) result(v)
    class(mcmc), intent(in) :: self
    type(model) :: v
    
    v = self%v
    
    return 
  end function mcmc_get_v

  !---------------------------------------------------------------------

  function mcmc_get_k(self) result(k)
    class(mcmc), intent(in) :: self
    type(model) :: k
    
    k = self%k
    
    return 
  end function mcmc_get_k

  !---------------------------------------------------------------------

  function mcmc_get_s(self) result(s)
    class(mcmc), intent(in) :: self
    type(model) :: s
    
    s = self%s
    
    return 
  end function mcmc_get_s

  !---------------------------------------------------------------------

  function mcmc_get_h(self) result(h)
    class(mcmc), intent(in) :: self
    type(model) :: h
    
    h = self%h
    
    return 
  end function mcmc_get_h

  !---------------------------------------------------------------------

  function mcmc_get_q(self) result(q)
    class(mcmc), intent(in) :: self
    type(model) :: q
    
    q = self%q
    
    return 
  end function mcmc_get_q
  
  !---------------------------------------------------------------------
    
  function mcmc_get_label(self) result(label)
    class(mcmc), intent(in) :: self
    character(6), dimension(self%n_proposal_type) :: label
    
    label = self%label
    
    return
  end function mcmc_get_label
    
  !---------------------------------------------------------------------

  function mcmc_write_out_u(self) result(out)
    class(mcmc), intent(inout) :: self
    double precision :: out(self%n_events)
    
    out = self%u%get_all_x()
    
    return 
  end function mcmc_write_out_u

  !---------------------------------------------------------------------

  function mcmc_write_out_v(self) result(out)
    class(mcmc), intent(inout) :: self
    double precision :: out(self%n_events)
    
    out = self%v%get_all_x()
    
    return 
  end function mcmc_write_out_v

  !---------------------------------------------------------------------

  function mcmc_write_out_k(self) result(out)
    class(mcmc), intent(inout) :: self
    double precision :: out(self%n_events)
    
    out = self%k%get_all_x()
    
    return 
  end function mcmc_write_out_k

  !---------------------------------------------------------------------

  function mcmc_write_out_s(self) result(out)
    class(mcmc), intent(inout) :: self
    double precision :: out(self%n_events)
    
    out = self%s%get_all_x()
    
    return 
  end function mcmc_write_out_s

  !---------------------------------------------------------------------

  function mcmc_write_out_h(self) result(out)
    class(mcmc), intent(inout) :: self
    double precision :: out(self%n_events)
    
    out = self%h%get_all_x()
    
    return 
  end function mcmc_write_out_h

  !---------------------------------------------------------------------

  function mcmc_write_out_q(self) result(out)
    class(mcmc), intent(inout) :: self
    double precision :: out(self%n_sta * self%n_events)
    
    out = self%q%get_all_x()
    
    return 
  end function mcmc_write_out_q
  
end module cls_mcmc
