module cls_observation
  implicit none

  type observation
     private

     character(len=256) :: pol_file, sta_file
     
     integer :: n_events
     integer :: n_stations
     integer, allocatable :: n_obs(:)
 
     double precision, allocatable :: azi(:,:), inc(:,:)
     character(len=1), allocatable :: pol(:,:)
     character(len=16), allocatable :: sta(:)
   contains
     procedure :: read_polarity_data => observation_read_polarity_data
     procedure :: read_sta_file => observation_read_sta_file
     procedure :: get_n_events => observation_get_n_events
     procedure :: get_n_stations => observation_get_n_stations
     procedure :: get_n_obs_single => observation_get_n_obs_single
     procedure :: get_n_obs_all => observation_get_n_obs_all
     procedure :: get_azi_single => observation_get_azi_single
     procedure :: get_azi_all => observation_get_azi_all
     procedure :: get_inc_single => observation_get_inc_single
     procedure :: get_inc_all => observation_get_inc_all
     procedure :: get_pol_single => observation_get_pol_single
     procedure :: get_pol_all => observation_get_pol_all
     procedure :: get_sta_id => observation_get_sta_id
     procedure :: get_stations => observation_get_stations
     
     
     generic :: get_n_obs => get_n_obs_single, get_n_obs_all
     generic :: get_azi => get_azi_single, get_azi_all
     generic :: get_inc => get_inc_single, get_inc_all
     generic :: get_pol => get_pol_single, get_pol_all
     
     
   end type observation
  
  interface observation
     module procedure init_observation
  end interface observation
     

contains


  !-----------------------------------------------------------------------
  
  type(observation) function init_observation(sta_file, pol_file) &
       result(self)
    character(*), intent(in) :: sta_file, pol_file
    
    self%sta_file = sta_file
    self%pol_file = pol_file


    call self%read_sta_file()
    call self%read_polarity_data()


    

    

  end function init_observation

  !-----------------------------------------------------------------------

  subroutine observation_read_polarity_data(self)
    class(observation), intent(inout) :: self
    integer :: io, ierr

    open(newunit=io, file=self%pol_file, status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
       write(*,*) 'Error opening file ', trim(self%pol_file)
       stop
       
    end if

    block
      integer :: id, n_obs, i, n_obs_max, i_evt, i_obs, sta_id
      double precision :: lat, lon, dep, azi, inc
      character(len=256) :: time, dummy
      character(len=16) :: sta
      character(len=1) :: pol
      logical :: debug = .false.

      ! Fist, cound the number of events and each number of observations
      self%n_events = 0
      n_obs_max = 0
      self%n_obs = [integer :: ]
      do
         read(io, *, iostat=ierr) id, time, lat, lon, dep
         if (ierr /= 0) exit
         self%n_events = self%n_events + 1
         read(io, *, iostat=ierr) dummy, n_obs
         if (n_obs > n_obs_max) then
            n_obs_max = n_obs
         end if
         self%n_obs = [self%n_obs, n_obs]
         do i = 1, n_obs
            read(io,*)
         end do
      end do

      ! Next read the azimuth and inclination for each observation
      allocate(self%azi(self%n_stations, self%n_events))
      allocate(self%inc(self%n_stations, self%n_events))
      allocate(self%pol(self%n_stations, self%n_events))
      self%azi = 0.0
      self%inc = 0.0
      self%pol = '-'
      
      rewind(io)

      
      do i_evt = 1, self%n_events
          read(io, *) sta_id, time, lat, lon, dep
          read(io, *) dummy, n_obs
          do i_obs = 1, n_obs
             read(io, *) sta, azi, inc, pol
             sta_id = self%get_sta_id(sta)
             self%azi(sta_id, i_evt) = azi
             self%inc(sta_id, i_evt) = inc
             self%pol(sta_id, i_evt) = pol
             
          end do
      end do

      close(io)

      ! Check
      if (debug) then
         do i_evt = 1, self%n_events
            do i_obs = 1, self%n_obs(i_evt)
               print *, i_evt, self%azi(i_obs, i_evt), &
                    & self%inc(i_obs, i_evt), self%pol(i_obs, i_evt)
            end do
         end do
      end if
      
    end block
    
  end subroutine observation_read_polarity_data


  !-----------------------------------------------------------------------

  subroutine observation_read_sta_file(self)
    class(observation), intent(inout) :: self
    integer :: io, ierr, i
    character(len=16) :: sta
    double precision :: lat, lon, dep
    
    open(newunit=io, file=self%sta_file, status='old', iostat=ierr)

    if (ierr /= 0) then
       write(*,*) 'Error opening file ', trim(self%sta_file)
       error stop
    end if

    self%n_stations = 0
    
    do
       read(io, *, iostat=ierr)
       if (ierr /= 0) exit
       self%n_stations = self%n_stations + 1
    end do

    allocate(self%sta(self%n_stations))

    rewind(io)

    do i = 1, self%n_stations
       read(io, *) self%sta(i), lat, lon, dep
    end do
    
    close(io)
    
  end subroutine observation_read_sta_file

  !-----------------------------------------------------------------------
  
  integer function observation_get_n_events(self) result(n_events)
    class(observation), intent(in) :: self

    n_events = self%n_events
    
  end function observation_get_n_events

  !-----------------------------------------------------------------------

  integer function observation_get_n_stations(self) result(n_stations)
    class(observation), intent(in) :: self

    n_stations = self%n_stations
    
  end function observation_get_n_stations
  
  !-----------------------------------------------------------------------

  integer function observation_get_n_obs_single(self, i_evt) result(n_obs)
    class(observation), intent(in) :: self
    integer, intent(in) :: i_evt

    n_obs = self%n_obs(i_evt)
    
  end function observation_get_n_obs_single

  !-----------------------------------------------------------------------

  function observation_get_n_obs_all(self) result(n_obs)
    class(observation), intent(in) :: self
    integer :: n_obs(self%n_events)
    
    n_obs = self%n_obs
    
  end function observation_get_n_obs_all
  
  !-----------------------------------------------------------------------

  function observation_get_azi_single(self, i_evt) result(azi)
    class(observation), intent(in) :: self
    integer, intent(in) :: i_evt
    double precision :: azi(self%n_stations)

    if (i_evt > self%n_events .or. i_evt < 1) then
       print *, 'Error: event number ', i_evt, ' is out of range'
       error stop
    end if
    azi(1:self%n_stations) = self%azi(1:self%n_stations, i_evt)

  end function observation_get_azi_single

  !-----------------------------------------------------------------------

  function observation_get_azi_all(self) result(azi)
    class(observation), intent(in) :: self
    double precision :: azi(self%n_stations, self%n_events)
    integer :: i_evt

    azi(1:self%n_stations, 1:self%n_events) &
         = self%azi(1:self%n_stations, 1:self%n_events)

  end function observation_get_azi_all
  

  !-----------------------------------------------------------------------

  function observation_get_inc_single(self, i_evt) result(inc)
    class(observation), intent(in) :: self
    integer, intent(in) :: i_evt
    double precision :: inc(self%n_stations)

    if (i_evt > self%n_events .or. i_evt < 1) then
       print *, 'Error: event number ', i_evt, ' is out of range'
       error stop
    end if
    inc(1:self%n_stations) = self%inc(1:self%n_stations, i_evt)

  end function observation_get_inc_single

  !-----------------------------------------------------------------------

  function observation_get_inc_all(self) result(inc)
    class(observation), intent(in) :: self
    double precision :: inc(self%n_stations, self%n_events)

    inc(1:self%n_stations, 1:self%n_events) &
         = self%inc(1:self%n_stations, 1:self%n_events)

  end function observation_get_inc_all

  !-----------------------------------------------------------------------

  function observation_get_pol_single(self, i_evt) result(pol)
    class(observation), intent(in) :: self
    integer, intent(in) :: i_evt
    character(len=1) :: pol(self%n_stations)

    if (i_evt > self%n_events .or. i_evt < 1) then
       print *, 'Error: event number ', i_evt, ' is out of range'
       error stop
    end if
    pol(1:self%n_stations) = self%pol(1:self%n_stations, i_evt)

  end function observation_get_pol_single

  !-----------------------------------------------------------------------

  function observation_get_pol_all(self) result(pol)
    class(observation), intent(in) :: self
    character(len=1) :: pol(self%n_stations, self%n_events)

    pol(1:self%n_stations, 1:self%n_events) &
         = self%pol(1:self%n_stations, 1:self%n_events)

  end function observation_get_pol_all
  
  !-----------------------------------------------------------------------

  integer function observation_get_sta_id(self, sta) result(id)
    class(observation), intent(in) :: self
    character(len=16), intent(in) :: sta
    integer :: i
    
    id = -1
    do i = 1, self%n_stations
       if (self%sta(i) == sta) then
          id = i
          exit
       end if
    end do
    
  end function observation_get_sta_id

  !-----------------------------------------------------------------------

  function observation_get_stations(self) result(sta)
    class(observation), intent(in) :: self
    character(len=16) :: sta(self%n_stations)
    integer :: i
    
    sta(1:self%n_stations) = self%sta(1:self%n_stations)
    
  end function observation_get_stations
  
end module cls_observation

    
  
