module cls_output
  use cls_radiation, only: radiation
  use cls_moment, only: moment
  implicit none
  
  type output
     private
     
     double precision :: d_theta = 2.0d0
     double precision :: d_phi = 2.0d0
     integer :: n_theta, n_phi
     
     double precision, allocatable :: pol_cum(:,:)
     integer :: n_pol_cum = 0
     integer :: n_stations = 0
     
     double precision, allocatable :: source_a(:), source_b(:)
     double precision, allocatable :: fault_h(:), fault_v(:)
     double precision, allocatable :: p_dip(:), b_dip(:), t_dip(:)
     double precision, allocatable :: p_strike(:), b_strike(:), t_strike(:)
     double precision, allocatable :: likelihood(:), a_mle(:)
     double precision, allocatable :: mxx(:), myy(:), mzz(:), mxy(:), mxz(:), myz(:)
     double precision, allocatable :: q(:)
     integer :: n_source_type = 0
     integer :: n_fault_type = 0
     
     
   contains 
     procedure :: count_pol => output_count_pol
     procedure :: write_pol => output_write_pol
     procedure :: count_source_type => output_count_source_type
     procedure :: get_pol_cum => output_get_pol_cum
     procedure :: count_fault_type => output_count_fault_type
     procedure :: count_principal_axes => output_count_principal_axes
     procedure :: save_likelihood => output_save_likelihood
     procedure :: save_a_mle => output_save_a_mle
     procedure :: save_moment => output_save_moment
     procedure :: save_q => output_save_q
     procedure :: get_n_pol_cum => output_get_n_pol_cum
     procedure :: set_pol_cum => output_set_pol_cum
     procedure :: set_n_pol_cum => output_set_n_pol_cum
     procedure :: get_n_theta => output_get_n_theta
     procedure :: get_n_phi => output_get_n_phi
     procedure :: add_source_type => output_add_source_type
     procedure :: write_source_type => output_write_source_type
     procedure :: write_fault_type => output_write_fault_type
     procedure :: write_principal_axes => output_write_principal_axes
     procedure :: write_summary => output_write_summary
     procedure :: write_q => output_write_q
     
  end type output
  
  
  interface output
     module procedure :: init_output
  end interface output
  
contains
  
  !---------------------------------------------------------------------
  
  type(output) function init_output() result(self)
    
    
    self%n_phi = 360.0d0 / self%d_phi
    self%n_theta = 90.0d0 / self%d_theta
    
    
    allocate(self%pol_cum(self%n_phi, self%n_theta))
    
    self%pol_cum = 0.0d0

    
    self%source_a = [ double precision :: ]
    self%source_b = [ double precision :: ]
    self%fault_h = [ double precision :: ]
    self%fault_v = [ double precision :: ]
    self%p_dip = [ double precision :: ]
    self%b_dip = [ double precision :: ]
    self%t_dip = [ double precision :: ]
    self%p_strike = [ double precision :: ]
    self%b_strike = [ double precision :: ]
    self%t_strike = [ double precision :: ]
    self%likelihood = [ double precision :: ]
    self%a_mle = [ double precision :: ]
    self%mxx = [ double precision :: ]
    self%myy = [ double precision :: ]
    self%mzz = [ double precision :: ]
    self%mxy = [ double precision :: ]
    self%mxz = [ double precision :: ]
    self%myz = [ double precision :: ]
    self%q = [ double precision :: ]
    
  end function init_output
  
  !---------------------------------------------------------------------

  subroutine output_save_likelihood(self, likelihood)
    class(output), intent(inout) :: self
    double precision, intent(in) :: likelihood
    
    self%likelihood = [self%likelihood, likelihood]
    
  end subroutine output_save_likelihood

  !---------------------------------------------------------------------

  subroutine output_save_a_mle(self, a_mle)
    class(output), intent(inout) :: self
    double precision, intent(in) :: a_mle
    
    self%a_mle = [self%a_mle, a_mle]
    
  end subroutine output_save_a_mle
  
  !---------------------------------------------------------------------

  subroutine output_save_moment(self, mt)
    class(output), intent(inout) :: self
    type(moment), intent(in) :: mt
    double precision :: m(3,3)
    
    m = mt%get_moment()
    
    self%mxx = [self%mxx, m(1,1)]
    self%myy = [self%myy, m(2,2)]
    self%mzz = [self%mzz, m(3,3)]
    self%mxy = [self%mxy, m(1,2)]
    self%mxz = [self%mxz, m(1,3)]
    self%myz = [self%myz, m(2,3)]
    
  end subroutine output_save_moment

  !---------------------------------------------------------------------

  subroutine output_save_q(self, q)
    class(output), intent(inout) :: self
    double precision, intent(in) :: q(:)

    self%n_stations = size(q)
    self%q = [self%q, q(1:self%n_stations)]
    
    
  end subroutine output_save_q
  
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
    double precision :: a, b, l1, l2, l3
    ! double precision :: f, l_sum
    double precision, parameter :: pi = acos(-1.0d0)
    double precision, parameter :: s3 = sqrt(3.d0)
    double precision :: m_iso, m_dc, m_clvd

    ! Source type plot from Aso et al. (2016): see their type j
    
    l = mt%get_eigenvalues()
    l1 = l(1)
    l2 = l(2)
    l3 = l(3)
    !l_sum = l1 + l2 + l3
    !f = sqrt(1.d0 - abs(l_sum) / s3)
    

    
    !a = 6.d0 / pi * atan((l1 - 2.d0 * l2 + l3) / (s3 * (l1 - l3))) * f
    !b = (l_sum / s3) / (1.d0 + f)

    m_dc = 0.5d0 * (l1 - l3 - abs(l1 + l3 - 2.d0 * l2))
    m_clvd = 2.d0 * (l1 + l3 - 2.d0 * l2) / 3.d0 
    m_iso = (l1 + l2 + l3) / 3.d0
    
    a = m_clvd / (m_dc + abs(m_clvd) + abs(m_iso))
    b = m_iso / (m_dc + abs(m_clvd) + abs(m_iso))

    self%n_source_type = self%n_source_type + 1

    self%source_a = [self%source_a, a]
    self%source_b = [self%source_b, b]
    !self%source_a(self%n_source_type) = a
    !self%source_b(self%n_source_type) = b
    
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
    double precision :: l, n, zt, zp, zb
    

    axes = mt%get_principal_axes()

    if (axes(3,1) < 0.d0) then
       axes(1,1) = -axes(1,1)
       axes(2,1) = -axes(2,1)
       axes(3,1) = -axes(3,1)
    end if
    if (axes(3,2) < 0.d0) then
      axes(1,2) = -axes(1,2)
       axes(2,2) = -axes(2,2)
       axes(3,2) = -axes(3,2)
    end if
   if (axes(3,3) < 0.d0) then
       axes(1,3) = -axes(1,3)
       axes(2,3) = -axes(2,3)
       axes(3,3) = -axes(3,3)
    end if
    t_dip = asin((axes(3,1)) / sqrt(axes(1,1)**2 + axes(2,1)**2 + axes(3,1)**2))
    b_dip = asin((axes(3,2)) / sqrt(axes(1,2)**2 + axes(2,2)**2 + axes(3,2)**2))
    p_dip = asin((axes(3,3)) / sqrt(axes(1,3)**2 + axes(2,3)**2 + axes(3,3)**2))
    
    !psi = atan(sin(t_dip) / sin(p_dip)) - rad45
    !
    !denom = sin(rad35_264) * sin(b_dip) + cos(rad35_264) * cos(b_dip) * cos(psi)
    !denom = denom * 3.d0
    !
    !h = sqrt(2.d0) * cos(b_dip) * sin(psi) / denom
    !v = sqrt(2.d0) * (cos(rad35_264) * sin(b_dip) - sin(rad35_264) * cos(b_dip) * cos(psi)) / &
    !     denom


    zp = sin(p_dip)
    zb = sin(b_dip)
    zt = sin(t_dip)
    l = 2.d0 * sin(0.5d0 * acos((zt + zp + zb) / sqrt(3.d0)))
    n = sqrt(2.d0 * ((zb-zp)**2 + (zb-zt)**2 + (zt-zp)**2))
    h = sqrt(3.d0) * l / n * (zt - zp)
    v = l / n * (2.d0 * zb - zp - zt)


    self%n_fault_type = self%n_fault_type + 1
    
    self%fault_h = [self%fault_h, h]
    self%fault_v = [self%fault_v, v]
    
    
  end subroutine output_count_fault_type
  
  
  !---------------------------------------------------------------------
  
  subroutine output_count_principal_axes(self, mt)
    class(output), intent(inout) :: self
    type(moment), intent(in) :: mt
    double precision :: axes(3,3)
    double precision :: p_dip, b_dip, t_dip, p_strike, b_strike, t_strike
    double precision, parameter :: pi = acos(-1.0d0)
    double precision, parameter :: rad2deg = 180.d0 / pi
    axes = mt%get_principal_axes()

    if (axes(3,1) < 0.d0) then
       axes(1,1) = -axes(1,1)
       axes(2,1) = -axes(2,1)
       axes(3,1) = -axes(3,1)
    end if
    if (axes(3,2) < 0.d0) then
      axes(1,2) = -axes(1,2)
       axes(2,2) = -axes(2,2)
       axes(3,2) = -axes(3,2)
    end if
   if (axes(3,3) < 0.d0) then
       axes(1,3) = -axes(1,3)
       axes(2,3) = -axes(2,3)
       axes(3,3) = -axes(3,3)
    end if
    t_dip = asin((axes(3,1)) / sqrt(axes(1,1)**2 + axes(2,1)**2 + axes(3,1)**2)) &
         * rad2deg
    b_dip = asin((axes(3,2)) / sqrt(axes(1,2)**2 + axes(2,2)**2 + axes(3,2)**2)) &
         * rad2deg
    p_dip = asin((axes(3,3)) / sqrt(axes(1,3)**2 + axes(2,3)**2 + axes(3,3)**2)) &
         * rad2deg
    

    t_strike = atan2(axes(2,1), axes(1,1)) * rad2deg
    b_strike = atan2(axes(2,2), axes(1,2)) * rad2deg
    p_strike = atan2(axes(2,3), axes(1,3)) * rad2deg

    

    self%p_dip = [self%p_dip, p_dip]
    self%b_dip = [self%b_dip, b_dip]
    self%t_dip = [self%t_dip, t_dip]
    self%p_strike = [self%p_strike, p_strike]
    self%b_strike = [self%b_strike, b_strike]
    self%t_strike = [self%t_strike, t_strike]
    
    
  end subroutine output_count_principal_axes
  
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

  subroutine output_write_source_type(self, filename, mode)
    class(output), intent(inout) :: self
    character(len=*), intent(in) :: filename
    character(len=*), intent(in), optional :: mode
    integer :: io, ierr
    integer :: i

    if (present(mode)) then
       if (mode == 'replace') then
          open(newunit=io, file=filename, status='replace', iostat=ierr)
       else if (mode == 'append') then
          open(newunit=io, file=filename, status='old', iostat=ierr, &
               position='append')
       else
          print *, 'Error: mode must be either "replace" or "append"'
          error stop
       end if
    else
       open(newunit=io, file=filename, status='replace', iostat=ierr)
    end if
    
    if (ierr /= 0) then
       print *, 'Error opening file'
       error stop
    end if

    print *, 'Writing to file ', filename
    do i = 1, size(self%source_a)
       write(io, *) self%source_a(i), self%source_b(i)
    end do
    close(io)
  end subroutine output_write_source_type

  !---------------------------------------------------------------------
  
  subroutine output_write_fault_type(self, filename, mode)
    class(output), intent(inout) :: self
    character(len=*), intent(in) :: filename
    character(len=*), intent(in), optional :: mode
    integer :: io, ierr
    integer :: i

    if (present(mode)) then
       if (mode == 'replace') then
          open(newunit=io, file=filename, status='replace', iostat=ierr)
       else if (mode == 'append') then
          open(newunit=io, file=filename, status='old', iostat=ierr, &
               position='append')
       else
          print *, 'Error: mode must be either "replace" or "append"'
          error stop
       end if
    else
       open(newunit=io, file=filename, status='replace', iostat=ierr)
    end if
    
    if (ierr /= 0) then
       print *, 'Error opening file'
       error stop
    end if

    print *, 'Writing to file ', filename
    do i = 1, size(self%fault_h)
       write(io, *) self%fault_h(i), self%fault_v(i)
    end do
    close(io)
  end subroutine output_write_fault_type

  !---------------------------------------------------------------------

  subroutine output_write_principal_axes(self, filename, mode)
    class(output), intent(inout) :: self
    character(len=*), intent(in) :: filename
    character(len=*), intent(in), optional :: mode
    integer :: io, ierr
    integer :: i

    if (present(mode)) then
       if (mode == 'replace') then
          open(newunit=io, file=filename, status='replace', iostat=ierr)
       else if (mode == 'append') then
          open(newunit=io, file=filename, status='old', iostat=ierr, &
               position='append')
       else
          print *, 'Error: mode must be either "replace" or "append"'
          error stop
       end if
    else
       open(newunit=io, file=filename, status='replace', iostat=ierr)
    end if

    if (ierr /= 0) then
       print *, 'Error opening file'
       error stop
    end if

    print *, 'Writing to file ', filename
    do i = 1, size(self%p_dip)
       write(io, *) self%t_strike(i), self%t_dip(i), self%b_strike(i), &
            self%b_dip(i), self%p_strike(i), self%p_dip(i)
       
       
    end do
    close(io)

  end subroutine output_write_principal_axes
  
  !--------------------------------------------------------------
  
  subroutine output_write_summary(self, filename, mode)
    class(output), intent(inout) :: self
    character(len=*), intent(in) :: filename
    character(len=*), intent(in), optional :: mode
    integer :: io, ierr
    integer :: i
    
    if (present(mode)) then
       if (mode == 'replace') then
          open(newunit=io, file=filename, status='replace', iostat=ierr)
       else if (mode == 'append') then
          open(newunit=io, file=filename, status='old', iostat=ierr, &
               position='append')
       else
          print *, 'Error: mode must be either "replace" or "append"'
          error stop
       end if
    else
       open(newunit=io, file=filename, status='replace', iostat=ierr)
    end if

    if (ierr /= 0) then
       print *, 'Error opening file'
       error stop
    end if

    print *, 'Writing to file ', filename

    do i = 1, size(self%likelihood)
       ! mrr, mtt, mff, mrt, mrf, mtf
       write(io, '(E16.3,6F8.3,6F9.2,5F12.7)') self%likelihood(i), self%mzz(i), self%mxx(i), self%myy(i), &
            self%mxz(i), -self%myz(i), -self%mxy(i), &
            self%t_strike(i), self%t_dip(i), &
            self%b_strike(i), self%b_dip(i), &
            self%p_strike(i), self%p_dip(i), &
            self%fault_h(i), self%fault_v(i), &
            self%source_a(i), self%source_b(i), &
            self%a_mle(i)
    end do
    close(io)
    
  end subroutine output_write_summary

  !--------------------------------------------------------------

  subroutine output_write_q(self, filename, mode)
    class(output), intent(inout) :: self
    character(len=*), intent(in) :: filename
    character(len=*), intent(in), optional :: mode
    character(len=20) :: fmt
    integer :: io, ierr
    integer :: i, j
    
    if (present(mode)) then
       if (mode == 'replace') then
          open(newunit=io, file=filename, status='replace', iostat=ierr)
       else if (mode == 'append') then
          open(newunit=io, file=filename, status='old', iostat=ierr, &
               position='append')
       else
          print *, 'Error: mode must be either "replace" or "append"'
          error stop
       end if
    else
       open(newunit=io, file=filename, status='replace', iostat=ierr)
    end if

    if (ierr /= 0) then
       print *, 'Error opening file'
       error stop
    end if

    print *, 'Writing to file ', trim(filename)
    write(fmt, '(A, I0, A)') '(', self%n_stations, 'F8.3)'
    if (size(self%q) == 0) then
       print *, 'No data to write'
       return
    end if
    do i = 1, size(self%q) / self%n_stations
       write(io, fmt) (self%q((i-1)*self%n_stations + j), j= 1, self%n_stations)
    end do
    close(io)
    
  end subroutine output_write_q
  
end module cls_output
