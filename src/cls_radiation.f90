module cls_radiation
  implicit none
  
  type radiation
     private

     ! Moment tensor components (for visualization purposes, not used in the code)
     double precision :: m_xx, m_yy, m_zz, m_xy, m_xz, m_yz

     ! Moment tensor components (for the code)
     double precision :: mt(3,3)
     

   contains
     procedure :: get_mt => radiation_get_mt
     procedure :: calc_radiation => radiation_calc_radiation
     procedure :: calc_projection_point => radiation_calc_projection_point

  end type radiation

  interface radiation
     module procedure init_radiation
  end interface radiation

contains

  !-----------------------------------------------------------------------

  type(radiation) function init_radiation(m_xx, m_yy, m_zz, m_xy, m_xz, m_yz) result(self)
    double precision, intent(in) :: m_xx, m_yy, m_zz, m_xy, m_xz, m_yz

    
    ! Set moment tensor components
    self%m_xx = m_xx
    self%m_yy = m_yy
    self%m_zz = m_zz
    self%m_xy = m_xy
    self%m_xz = m_xz
    self%m_yz = m_yz

    self%mt(1,1) = m_xx
    self%mt(2,2) = m_yy
    self%mt(3,3) = m_zz
    self%mt(1,2) = m_xy
    self%mt(1,3) = m_xz
    self%mt(2,3) = m_yz
    self%mt(2,1) = m_xy
    self%mt(3,1) = m_xz
    self%mt(3,2) = m_yz

    

  end function init_radiation

  !-----------------------------------------------------------------------

  subroutine radiation_calc_radiation(self, phi, theta, phase, amp)
    class(radiation), intent(inout) :: self
    double precision, intent(in) :: phi, theta
    character(len=*), intent(in) :: phase

    integer :: i, j
    double precision, intent(out) :: amp(3)
    double precision :: gamma2(3,3), amp_xyz(3)
    double precision :: e_r(3), e_theta(3), e_phi(3)
    
    
    ! Check phase validity
    if (phase /= "P" .and. phase /= "S") then
       print *, "Error: phase must be P or S"
       error stop
    end if

    ! Convert spherical coordinates to cartesian coordinates
    e_r(1) = sin(theta)*cos(phi)
    e_r(2) = sin(theta)*sin(phi)
    e_r(3) = cos(theta)

    e_theta(1) = cos(theta)*cos(phi)
    e_theta(2) = cos(theta)*sin(phi)
    e_theta(3) = -sin(theta)

    e_phi(1) = -sin(phi)
    e_phi(2) = cos(phi)
    e_phi(3) = 0.d0
    


    ! concurrent loop for i=1, 3 and j=1, 3
    do concurrent(i = 1:3, j = 1:3)
       gamma2(i,j) = e_r(i) * e_r(j)
    end do

    ! Calculate radiation pattern
    amp = 0.d0
    if (phase == "P") then
       !print *, "P-wave radiation pattern"

       do i = 1, 3
          amp_xyz(i) = e_r(i) * sum(gamma2(:,:)*self%mt(:,:))
       end do
       
    else
       print *, "S-wave radiation pattern"
    end if
    

    ! Convert to spherical coordinates
    amp(1) =  dot_product(amp_xyz, e_r) ! P-wave
    amp(2) =  dot_product(amp_xyz, e_theta) ! SV-wave
    amp(3) =  dot_product(amp_xyz, e_phi) ! SH-wave

    
  end subroutine radiation_calc_radiation

  !-----------------------------------------------------------------------

  subroutine radiation_calc_projection_point(self, phi, theta, x, y)
    class(radiation), intent(inout) :: self
    double precision, intent(in) :: phi, theta
    double precision, intent(out) :: x, y
    double precision :: x_on_sphere, y_on_sphere, z_on_sphere

    x_on_sphere = sin(theta)*cos(phi)
    y_on_sphere = sin(theta)*sin(phi)
    z_on_sphere = -cos(theta)

    ! stereographic projection
    if (z_on_sphere == 1.d0) then
       x = 0.d0
       y = 0.d0
    else 
       x = x_on_sphere * sqrt(2.d0 / (1.d0 - z_on_sphere))
       y = y_on_sphere * sqrt(2.d0 / (1.d0 - z_on_sphere))
    end if
    
    
    
  end subroutine radiation_calc_projection_point

  !-----------------------------------------------------------------------
  function radiation_get_mt(self) result(mt)
    class(radiation), intent(in) :: self
    double precision :: mt(6)

    mt(1) = self%m_xx
    mt(2) = self%m_yy
    mt(3) = self%m_zz
    mt(4) = self%m_xy
    mt(5) = self%m_xz
    mt(6) = self%m_yz

  end function radiation_get_mt
    
     
     
end module cls_radiation
