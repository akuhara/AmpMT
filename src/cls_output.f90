module cls_output
  use cls_radiation, only: radiation
  use cls_moment, only: moment
  implicit none

  type output
     private
     
     double precision :: d_theta = 2.0d0
     double precision :: d_phi = 2.0d0
     integer :: n_theta, n_phi
     integer :: n_mod_max = 0
     
     double precision, allocatable :: pol_cum(:,:)
     integer :: n_pol_cum = 0
     
     double precision, allocatable :: source_a(:), source_b(:)
     double precision, allocatable :: fault_h(:), fault_v(:)
     integer :: n_source_type = 0
     integer :: n_fault_type = 0
     
   contains 
     procedure :: count_pol => output_count_pol
     procedure :: write_pol => output_write_pol
     procedure :: count_source_type => output_count_source_type
     procedure :: get_pol_cum => output_get_pol_cum
     procedure :: count_fault_type => output_count_fault_type
     
     procedure :: get_n_pol_cum => output_get_n_pol_cum
     procedure :: set_pol_cum => output_set_pol_cum
     procedure :: set_n_pol_cum => output_set_n_pol_cum
     procedure :: get_n_theta => output_get_n_theta
     procedure :: get_n_phi => output_get_n_phi
     procedure :: add_source_type => output_add_source_type
     procedure :: write_source_type => output_write_source_type
     procedure :: write_append_source_type => output_write_append_source_type
     procedure :: write_fault_type => output_write_fault_type
     procedure :: write_append_fault_type => output_write_append_fault_type
     
  end type output
    
  
  interface output
       module procedure :: init_output
  end interface output

contains

  !---------------------------------------------------------------------
  
  type(output) function init_output(n_iter, n_chain, n_cool, &
       n_burn, n_interval) result(self)

    integer, intent(in) :: n_iter, n_chain, n_cool, n_burn, n_interval

    self%n_mod_max = n_chain * n_cool * (n_iter - n_burn) / n_interval
    
    self%n_phi = 360.0d0 / self%d_phi
    self%n_theta = 90.0d0 / self%d_theta
    
    
    allocate(self%pol_cum(self%n_phi, self%n_theta))
    self%pol_cum = 0.0d0

    
    allocate(self%source_a(self%n_mod_max))
    allocate(self%source_b(self%n_mod_max))
    allocate(self%fault_h(self%n_mod_max))
    allocate(self%fault_v(self%n_mod_max))
    
  end function init_output
  
!---------------------------------------------------------------------

  subroutine output_count_pol(self, mt)
    class(output), intent(inout) :: self
    type(moment), intent(in) :: mt
    type(radiation) :: rad
    double precision :: m(3,3)
    integer :: i, j
    double precision :: phi, theta
    double precision :: r
    double precision, parameter :: pi = acos(-1.0d0)
    double precision :: amp(3)

    m = mt%get_moment()

    rad = radiation(m(1,1), m(2,2), m(3,3), m(1,2), m(1,3), m(2,3))
    
    do i = 1, self%n_theta
       theta = (i-0.5d0) * self%d_theta * pi / 180.d0
       do j = 1, self%n_phi
          phi = (j-0.5d0) * self%d_phi * pi / 180.0d0
          call rad%calc_radiation(phi, theta, "P", amp)
          if (amp(1) > 0.d0) then
             self%pol_cum(j, i) = self%pol_cum(j, i) + 1.d0
          else
             self%pol_cum(j, i) = self%pol_cum(j, i) - 1.d0
          end if
       end do
    end do
    
    self%n_pol_cum = self%n_pol_cum + 1
    
  end subroutine output_count_pol

  !---------------------------------------------------------------------

  subroutine output_write_pol(self, filename)
    class(output), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer :: io, ierr
    integer :: i, j
    double precision :: phi, theta
    double precision :: point(2)
    double precision :: easting, northing
    double precision, parameter :: pi = acos(-1.0d0)
    

    open(newunit=io, file=filename, status='replace', iostat=ierr)
    if (ierr /= 0) then
       print *, 'Error opening file'
       error stop
    end if

    print *, 'Writing to file ', filename
    do i = 1, self%n_theta
       theta = (i-0.5d0) * self%d_theta * pi / 180.d0
       
       do j = 1, self%n_phi
          phi = (j-0.5d0) * self%d_phi * pi / 180.0d0
          point = projection_point(phi, theta)
          northing = point(1)
          easting = point(2)
          write(io, *) easting, northing, &
               self%pol_cum(j, i) / (self%n_pol_cum)
       end do
    end do
    close(io)
  end subroutine output_write_pol
    
    
  
  !---------------------------------------------------------------------

  subroutine output_count_source_type(self, mt)
    class(output), intent(inout) :: self
    type(moment), intent(in) :: mt
    double precision :: l(3)
    double precision :: a, b, l1, l2, l3, l_sum, f
    double precision, parameter :: pi = acos(-1.0d0)
    double precision, parameter :: s3 = sqrt(3.d0)

    ! Source type plot from Aso et al. (2016): see their type j
    
    l = mt%get_eigenvalues()
    l1 = l(1)
    l2 = l(2)
    l3 = l(3)
    l_sum = l1 + l2 + l3
    f = sqrt(1.d0 - abs(l_sum) / s3)
    

    a = 6.d0 / pi * atan((l1 - 2.d0 * l2 + l3) / (s3 * (l1 - l3))) * f
    b = (l_sum / s3) / (1.d0 + f)

    self%n_source_type = self%n_source_type + 1
    
    self%source_a(self%n_source_type) = a
    self%source_b(self%n_source_type) = b
    
    
  end subroutine output_count_source_type
  
  !---------------------------------------------------------------------

  subroutine output_count_fault_type(self, mt)
    class(output), intent(inout) :: self
    type(moment), intent(in) :: mt
    double precision :: axes(3,3)
    double precision :: p_dip, b_dip, t_dip
    double precision, parameter :: pi = acos(-1.0d0)
    double precision, parameter :: rad45 = 45.d0 * pi / 180.d0
    double precision, parameter :: rad35_264 = 35.264d0 * pi / 180.d0
    double precision :: psi, h, v, denom
    

    axes = mt%get_principal_axes()

    
    p_dip = atan(abs(axes(1,3)) / sqrt(axes(1,1)**2 + axes(1,2)**2))
    b_dip = atan(abs(axes(2,3)) / sqrt(axes(2,1)**2 + axes(2,2)**2))
    t_dip = atan(abs(axes(3,3)) / sqrt(axes(3,1)**2 + axes(3,2)**2))


    psi = atan(sin(t_dip) / sin(p_dip)) - rad45

    denom = sin(rad35_264) * sin(b_dip) + cos(rad35_264) * cos(b_dip) * cos(psi)
    denom = denom * 3.d0

    h = sqrt(2.d0) * cos(b_dip) * sin(psi) / denom
    v = sqrt(2.d0) * (cos(rad35_264) * sin(b_dip) - sin(rad35_264) * cos(b_dip) * cos(psi)) / &
         denom

    !print *, "h=", h, "v=", v
    self%n_fault_type = self%n_fault_type + 1
    self%fault_h(self%n_fault_type) = h
    self%fault_v(self%n_fault_type) = v

    !print *, "p_dip=", p_dip * 180.d0 / pi, "b_dip=", b_dip * 180.d0 / pi, &
    !        "t_dip=", t_dip * 180.d0 / pi, "psi=", psi * 180.d0 / pi, "h=", h, "v=", v, "denom=", denom
    !
    !print *, "validate ", sin(b_dip)**2 + sin(p_dip)**2 + sin(t_dip)**2     
  end subroutine output_count_fault_type
    
  
  !---------------------------------------------------------------------
  

  function projection_point(phi, theta) result (point)
    double precision, intent(in) :: phi, theta
    double precision :: point(2)
    double precision :: x_on_sphere, y_on_sphere, z_on_sphere

    x_on_sphere = sin(theta)*cos(phi)
    y_on_sphere = sin(theta)*sin(phi)
    z_on_sphere = -cos(theta)

    ! stereographic projection
    if (z_on_sphere == 1.d0) then
       point(1) = 0.d0
       point(2) = 0.d0
    else
       point(1) = x_on_sphere * sqrt(2.d0 / (1.d0 - z_on_sphere))
       point(2) = y_on_sphere * sqrt(2.d0 / (1.d0 - z_on_sphere))
    end if

    point(1:2) = point(1:2) / sqrt(2.d0)
  end function projection_point
    
  !---------------------------------------------------------------------

  function output_get_pol_cum(self) result(pol_cum)
    class(output), intent(in) :: self
    double precision :: pol_cum(self%n_phi, self%n_theta)
    
    pol_cum = self%pol_cum
  end function output_get_pol_cum

  !---------------------------------------------------------------------

  function output_get_n_pol_cum(self) result(n_pol_cum)
    class(output), intent(in) :: self
    integer :: n_pol_cum
    
    n_pol_cum = self%n_pol_cum
  end function output_get_n_pol_cum

  !---------------------------------------------------------------------

  subroutine output_set_pol_cum(self, pol_cum)
    class(output), intent(inout) :: self
    double precision, intent(in) :: pol_cum(self%n_phi, self%n_theta)
    
    self%pol_cum = pol_cum
  end subroutine output_set_pol_cum


  !---------------------------------------------------------------------

  subroutine output_set_n_pol_cum(self, n_pol_cum)
    class(output), intent(inout) :: self
    integer, intent(in) :: n_pol_cum
    
    self%n_pol_cum = n_pol_cum
  end subroutine output_set_n_pol_cum

  !---------------------------------------------------------------------

  function output_get_n_theta(self) result(n_theta)
    class(output), intent(in) :: self
    integer :: n_theta
    
    n_theta = self%n_theta
  end function output_get_n_theta

  !---------------------------------------------------------------------

  function output_get_n_phi(self) result(n_phi)
    class(output), intent(in) :: self
    integer :: n_phi
    
    n_phi = self%n_phi
  end function output_get_n_phi

  !---------------------------------------------------------------------

  subroutine output_add_source_type(self, a, b)
    class(output), intent(inout) :: self
    double precision, intent(in) :: a(:), b(:)

    if (size(a) /= size(b)) then
       print *, 'Error: source type a and b have different sizes'
       error stop
    end if
    
    self%source_a = [self%source_a, a]
    self%source_b = [self%source_b, b]
  end subroutine output_add_source_type

  !---------------------------------------------------------------------

  subroutine output_write_source_type(self, filename)
    class(output), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer :: io, ierr
    integer :: i
    
    open(newunit=io, file=filename, status='replace', iostat=ierr)
    if (ierr /= 0) then
       print *, 'Error opening file'
       error stop
    end if

    print *, 'Writing to file ', filename
    do i = 1, self%n_source_type
       write(io, *) self%source_a(i), self%source_b(i)
    end do
    close(io)
  end subroutine output_write_source_type

  !---------------------------------------------------------------------

  subroutine output_write_append_source_type(self, filename)
    class(output), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer :: io, ierr
    integer :: i
    
    open(newunit=io, file=filename, status='old', iostat=ierr, position='append')
    if (ierr /= 0) then
       print *, 'Error opening file'
       error stop
    end if

    print *, 'Writing to file ', filename
    write(io, *)"TEST"
    do i = 1, self%n_source_type
       write(io, *) self%source_a(i), self%source_b(i)
    end do
    close(io)
  end subroutine output_write_append_source_type

  !---------------------------------------------------------------------

  subroutine output_write_fault_type(self, filename)
    class(output), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer :: io, ierr
    integer :: i
    
    open(newunit=io, file=filename, status='replace', iostat=ierr)
    if (ierr /= 0) then
       print *, 'Error opening file'
       error stop
    end if

    print *, 'Writing to file ', filename
    do i = 1, self%n_fault_type
       write(io, *) self%fault_h(i), self%fault_v(i)
    end do
    close(io)
  end subroutine output_write_fault_type

  !---------------------------------------------------------------------

  subroutine output_write_append_fault_type(self, filename)
    class(output), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer :: io, ierr
    integer :: i
    
    open(newunit=io, file=filename, status='old', iostat=ierr, position='append')
    if (ierr /= 0) then
       print *, 'Error opening file'
       error stop
    end if

    print *, 'Writing to file ', filename
    do i = 1, self%n_fault_type
       write(io, *) self%fault_h(i), self%fault_v(i)
    end do
    close(io)
  end subroutine output_write_append_fault_type

end module cls_output
