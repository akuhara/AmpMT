module cls_moment
  implicit none
  double precision, private, parameter :: pi = acos(-1.d0)

  private :: v2gamma, h2theta, u2beta, f, f_prime, make_lambda
  private :: make_rota
  
  type moment
     private
     
     double precision :: u, v, k, s, h
     double precision :: gamma, theta, beta
     double precision :: lambda(3,3)
     double precision :: rota(3,3)
     double precision :: moment(3,3)
     double precision :: u_min = 0.d0, u_max = 0.75d0 * pi
     double precision :: v_min = -1.0 / 3.0, v_max = 1.0 / 3.0
     double precision :: k_min = 0.d0, k_max = 2.d0 * pi
     double precision :: s_min = -0.5d0 * pi, s_max = 0.5d0 * pi
     double precision :: h_min = 0.d0 * pi, h_max = 1.d0
  
   contains
     procedure :: construct_moment => moment_construct_moment
     procedure :: get_moment => moment_get_moment
     procedure :: get_u => moment_get_u
     procedure :: get_v => moment_get_v
     procedure :: get_k => moment_get_k
     procedure :: get_s => moment_get_s
     procedure :: get_h => moment_get_h
     procedure :: set_u => moment_set_u
     procedure :: set_v => moment_set_v
     procedure :: set_k => moment_set_k
     procedure :: set_s => moment_set_s
     procedure :: set_h => moment_set_h
     procedure :: get_u_min => moment_get_u_min
     procedure :: get_u_max => moment_get_u_max
     procedure :: get_v_min => moment_get_v_min
     procedure :: get_v_max => moment_get_v_max
     procedure :: get_k_min => moment_get_k_min
     procedure :: get_k_max => moment_get_k_max
     procedure :: get_s_min => moment_get_s_min
     procedure :: get_s_max => moment_get_s_max
     procedure :: get_h_min => moment_get_h_min
     procedure :: get_h_max => moment_get_h_max
     procedure :: get_strike => moment_get_strike
     procedure :: get_dip => moment_get_dip
     procedure :: get_rake => moment_get_rake
     procedure :: get_eigenvalues => moment_get_eigenvalues
     procedure :: get_principal_axes => moment_get_principal_axes
     procedure :: check_eigen => moment_check_eigen

  end type moment

  interface moment
     module procedure init_moment
  end interface moment
  
contains

  !---------------------------------------------------------------------

  function make_rota(k, s, theta) result(rota)
    double precision, intent(in) :: k, s, theta
    double precision :: rota(3,3)
    double precision :: z_k(3,3), x_theta(3,3), z_s(3,3), y_p4(3,3)
    double precision :: cos_k, sin_k, cos_theta, sin_theta, cos_s, sin_s
    double precision, parameter :: cos_p4 = cos(-pi / 4.d0), sin_p4 = sin(-pi / 4.d0)

    cos_k = cos(-k)
    sin_k = sin(-k)
    cos_theta = cos(theta)
    sin_theta = sin(theta)
    cos_s = cos(s)
    sin_s = sin(s)

    z_k = 0.d0
    z_k(1,1) = cos_k
    z_k(1,2) = -sin_k
    z_k(2,1) = sin_k
    z_k(2,2) = cos_k
    z_k(3,3) = 1.d0

    x_theta = 0.d0
    x_theta(1,1) = 1.d0
    x_theta(2,2) = cos_theta
    x_theta(2,3) = -sin_theta
    x_theta(3,2) = sin_theta
    x_theta(3,3) = cos_theta

    z_s = 0.d0
    z_s(1,1) = cos_s
    z_s(1,2) = -sin_s
    z_s(2,1) = sin_s
    z_s(2,2) = cos_s
    z_s(3,3) = 1.d0
    
    y_p4 = 0.d0
    y_p4(1,1) = cos_p4
    y_p4(1,3) = sin_p4
    y_p4(2,2) = 1.d0
    y_p4(3,1) = -sin_p4
    y_p4(3,3) = cos_p4

    rota = matmul(z_k, matmul(x_theta, matmul(z_s, y_p4)))
    
  end function make_rota
  
  !---------------------------------------------------------------------

  
  function make_lambda(beta, gamma) result(lambda)
    double precision, intent(in) :: beta, gamma
    double precision :: lambda(3,3)
    double precision :: a1, a2, a3
    double precision, parameter :: s2 = sqrt(2.d0)
    double precision, parameter :: s3 = sqrt(3.d0)

    a1 = sin(beta) * cos(gamma)
    a2 = sin(beta) * sin(gamma)
    a3 = cos(beta)


    lambda = 0.d0
    lambda(1,1) =  s3 * a1      -a2 + s2 * a3
    lambda(2,2) =         2.d0 * a2 + s2 * a3
    lambda(3,3) = -s3 * a1      -a2 + s2 * a3

    lambda = lambda / sqrt(6.d0)
    
  end function make_lambda
  
  !---------------------------------------------------------------------
  
  double precision function f(x, a)
    double precision, intent(in) :: x, a
    
    f = 0.75d0 * x - 0.5d0 * sin(2.d0 * x) + sin(4.d0 * x) / 16.d0 - a
    
  end function f

  !---------------------------------------------------------------------

  double precision function f_prime(x)
    double precision, intent(in) :: x
    
    f_prime = 0.75d0 - cos(2.d0 * x) + cos(4.d0 * x) / 4.d0
    
  end function f_prime
    
  !---------------------------------------------------------------------

  double precision function u2beta(u) result(beta)
    double precision, intent(in) :: u
    double precision :: beta1, beta2, beta_m
    ! Newton's method
    !beta = u * 4.d0 / 3.d0
    !
    !do while (abs(f(beta, u)) > 1.d-12)
    !   print *, " -- beta= ", beta, " f(beta) = ", f(beta, u), " f'(beta) = ", f_prime(beta)
    !  beta = beta - f(beta, u) / f_prime(beta)
    !end do

    ! Bisection method
    beta1 = 0.d0
    beta2 = pi

    
    do while (abs(beta1 - beta2) > 1.d-12)
       beta_m = (beta1 + beta2) / 2.d0
       if (f(beta_m, u) * f(beta1, u) > 0.d0) then
          beta1 = beta_m
       else if (f(beta_m, u) * f(beta2, u) > 0.d0) then
          beta2 = beta_m
       else
          exit
       end if
    end do

    beta = beta_m
    
  end function u2beta
  
  !---------------------------------------------------------------------

  double precision function v2gamma(v) result(gamma)
    double precision, intent(in) :: v
    
    gamma = asin(3.d0 * v) / 3.d0
    
  end function v2gamma
  
  !---------------------------------------------------------------------

  double precision function h2theta(h) result(theta)
    double precision, intent(in) :: h
    
    theta = acos(h)
    
  end function h2theta
  
  !---------------------------------------------------------------------

  type(moment) function init_moment(u, v, k, s, h) result(self)
    double precision, intent(in) :: u, v, k, s, h
    
    self%u = u
    self%v = v
    self%k = k
    self%s = s
    self%h = h
    
    call self%construct_moment()
    
  end function init_moment
  
  !---------------------------------------------------------------------
  
  subroutine moment_construct_moment(self)
    class(moment), intent(inout) :: self

    !print *, "u = ", self%u
    !print *, "v = ", self%v
    !print *, "k = ", self%k
    !print *, "s = ", self%s
    !print *, "h = ", self%h
    
    self%gamma = v2gamma(self%v)
    !print *, "gamma = ", self%gamma
    self%theta = h2theta(self%h)
    !print *, "theta = ", self%theta
    self%beta = u2beta(self%u)
    !print *, "beta = ", self%beta
    self%lambda = make_lambda(self%beta, self%gamma)
    !print *, "lambda = ", self%lambda
    self%rota = make_rota(self%k, self%s, self%theta)
    !print *, "rota = ", self%rota    
    self%moment = matmul(self%rota,matmul(self%lambda, transpose(self%rota)))
    !print *, "moment = ", self%moment
    
    
  end subroutine moment_construct_moment
  
  !---------------------------------------------------------------------

  function moment_get_moment(self) result(m)
    class(moment), intent(in) :: self
    double precision :: m(3,3)
    
    m = self%moment
    
  end function moment_get_moment

  !---------------------------------------------------------------------

  double precision function moment_get_u(self) result(u)
    class(moment), intent(in) :: self
    
    u = self%u
    
  end function moment_get_u

  !---------------------------------------------------------------------

  double precision function moment_get_v(self) result(v)
    class(moment), intent(in) :: self
    
    v = self%v
    
  end function moment_get_v

  !---------------------------------------------------------------------

  double precision function moment_get_k(self) result(k)
    class(moment), intent(in) :: self
    
    k = self%k
    
  end function moment_get_k

  !---------------------------------------------------------------------

  double precision function moment_get_s(self) result(s)
    class(moment), intent(in) :: self
    
    s = self%s
    
  end function moment_get_s

  !---------------------------------------------------------------------

  double precision function moment_get_h(self) result(h)
    class(moment), intent(in) :: self
    
    h = self%h
    
  end function moment_get_h

  !---------------------------------------------------------------------

  double precision function moment_get_strike(self) result(strike)
    class(moment), intent(in) :: self
    
    strike = self%k * 180.d0 / pi
    
  end function moment_get_strike
  
  !---------------------------------------------------------------------

  double precision function moment_get_dip(self) result(dip)
    class(moment), intent(in) :: self
    
    dip = self%s * 180.d0 / pi
    
  end function moment_get_dip

  !---------------------------------------------------------------------

  double precision function moment_get_rake(self) result(rake)
    class(moment), intent(in) :: self
    
    rake = self%theta * 180.d0 / pi
    
  end function moment_get_rake

  !---------------------------------------------------------------------
  
  subroutine moment_set_u(self, u)
    class(moment), intent(inout) :: self
    double precision, intent(in) :: u
    
    self%u = u
    
  end subroutine moment_set_u

  !---------------------------------------------------------------------

  subroutine moment_set_v(self, v)
    class(moment), intent(inout) :: self
    double precision, intent(in) :: v
    
    self%v = v
    
  end subroutine moment_set_v

  !---------------------------------------------------------------------

  subroutine moment_set_k(self, k)
    class(moment), intent(inout) :: self
    double precision, intent(in) :: k
    
    self%k = k
    
  end subroutine moment_set_k

  !---------------------------------------------------------------------

  subroutine moment_set_s(self, s)
    class(moment), intent(inout) :: self
    double precision, intent(in) :: s
    
    self%s = s
    
  end subroutine moment_set_s

  !---------------------------------------------------------------------

  subroutine moment_set_h(self, h)
    class(moment), intent(inout) :: self
    double precision, intent(in) :: h
    
    self%h = h
    
  end subroutine moment_set_h

  !---------------------------------------------------------------------

  double precision function moment_get_u_min(self) result(u_min)
    class(moment), intent(in) :: self
    
    u_min = self%u_min
    
  end function moment_get_u_min

  !---------------------------------------------------------------------

  double precision function moment_get_u_max(self) result(u_max)
    class(moment), intent(in) :: self
    
    u_max = self%u_max
    
  end function moment_get_u_max

  !---------------------------------------------------------------------

  double precision function moment_get_v_min(self) result(v_min)
    class(moment), intent(in) :: self
    
    v_min = self%v_min
    
  end function moment_get_v_min

  !---------------------------------------------------------------------

  double precision function moment_get_v_max(self) result(v_max)
    class(moment), intent(in) :: self
    
    v_max = self%v_max
    
  end function moment_get_v_max

  !---------------------------------------------------------------------

  double precision function moment_get_k_min(self) result(k_min)
    class(moment), intent(in) :: self
    
    k_min = self%k_min
    
  end function moment_get_k_min

  !---------------------------------------------------------------------

  double precision function moment_get_k_max(self) result(k_max)
    class(moment), intent(in) :: self
    
    k_max = self%k_max
    
  end function moment_get_k_max

  !---------------------------------------------------------------------

  double precision function moment_get_s_min(self) result(s_min)
    class(moment), intent(in) :: self
    
    s_min = self%s_min
    
  end function moment_get_s_min

  !---------------------------------------------------------------------

  double precision function moment_get_s_max(self) result(s_max)
    class(moment), intent(in) :: self
    
    s_max = self%s_max
    
  end function moment_get_s_max

  !---------------------------------------------------------------------

  double precision function moment_get_h_min(self) result(h_min)
    class(moment), intent(in) :: self
    
    h_min = self%h_min
    
  end function moment_get_h_min

  !---------------------------------------------------------------------

  double precision function moment_get_h_max(self) result(h_max)
    class(moment), intent(in) :: self
    
    h_max = self%h_max
    
  end function moment_get_h_max

  !---------------------------------------------------------------------

  function moment_get_eigenvalues(self) result(eigen_v)
    class(moment), intent(in) :: self
    double precision :: eigen_v(3)

    eigen_v(1:3) = [self%lambda(1,1), self%lambda(2,2), self%lambda(3,3)]
    
  end function moment_get_eigenvalues

  !---------------------------------------------------------------------

  function moment_get_principal_axes(self) result(principal_axes)
    class(moment), intent(in) :: self
    double precision :: principal_axes(3,3)

    principal_axes = self%rota

    !Note that the principal axes are the columns of the rotation matrix
    ! but in a reverse order, i.e. the last column is the first principal axis
    ! See Section 3.4 in Tape and Tape (2015).

    
  end function moment_get_principal_axes

  !---------------------------------------------------------------------

  subroutine moment_check_eigen(self)
    class(moment), intent(in) :: self
    double precision :: eigen_v(3)
    double precision :: principal_axes(3,3)
    double precision :: moment(3,3)
    

    eigen_v = self%get_eigenvalues()
    principal_axes = self%get_principal_axes()
    moment = self%get_moment()

    print *, "Eigenvalues: ", eigen_v
    print *, "Principal axes: ", principal_axes
    print *, "Moment: ", moment

    print *, "Check eigenvalues: ", matmul(transpose(principal_axes), matmul(moment, principal_axes))
    
  end subroutine moment_check_eigen
  
end module cls_moment
  

