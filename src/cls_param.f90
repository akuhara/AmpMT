!=======================================================================
!   CMTbyAmp
!   Copyright (C) 2024 Takeshi Akuhara
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!
!   Contact information
!
!   Email  : akuhara @ eri. u-tokyo. ac. jp 
!   Address: Earthquake Research Institute, The Univesity of Tokyo
!           1-1-1, Yayoi, Bunkyo-ku, Tokyo 113-0032, Japan
!
!=======================================================================
module cls_param
  use cls_line_text
  implicit none 
  
  type param
     private
     character(line_max) :: param_file
     
     ! # of processes for MPI
     integer :: n_procs
     
     ! Input files
     character(line_max) :: station_file
     character(line_max) :: polarity_file
     
     ! MCMC
     integer :: n_iter, n_burn, n_interval
     integer :: n_chains
     integer :: n_cool
     double precision :: temp_high

     ! Mode
     logical :: dc_only
     logical :: sample_prior
     logical :: use_amp
          
     ! Verbosity
     logical :: verb
     
   contains
     procedure :: read_param_file => param_read_param_file
     procedure :: get_n_procs => param_get_n_procs
     procedure :: get_n_iter => param_get_n_iter
     procedure :: get_n_burn => param_get_n_burn
     procedure :: get_n_interval => param_get_n_interval
     procedure :: get_n_chains => param_get_n_chains
     procedure :: get_n_cool => param_get_n_cool
     procedure :: get_temp_high => param_get_temp_high
     procedure :: get_dc_only => param_get_dc_only
     procedure :: get_sample_prior => param_get_sample_prior
     procedure :: get_use_amp => param_get_use_amp
     procedure :: get_station_file => param_get_station_file
     procedure :: get_polarity_file => param_get_polarity_file
     procedure :: set_value => param_set_value
     
  end type param
  
  interface param
     module procedure init_param
  end interface param

contains
  
  !---------------------------------------------------------------------
  
  type(param) function init_param(param_file, verb) result(self)
    character(len=*), intent(in) :: param_file
    logical, intent(in)      :: verb
    
    self%verb = verb
    
    
    ! Read parmeter file
    
    if (self%verb) then
       write(*,'(3A)')"<< Reading parameters from ", &
            & trim(param_file), " >>"
    end if    
    self%param_file = param_file
    call self%read_param_file()
    if (self%verb) then
       write(*,*)
    end if
    
    return 
  end function init_param

  !---------------------------------------------------------------------
    
  subroutine param_read_param_file(self)
    class(param), intent(inout) :: self
    character(len=line_max) :: line, name, val
    integer :: ierr, io
    type(line_text) :: lt
    logical :: is_ok

    open(newunit = io, file = self%param_file, &
         & status = 'old', iostat = ierr)
    if (ierr /= 0) then
       if (self%verb) then
          write(0,*) "ERROR: cannot open ", trim(self%param_file)
       end if
       stop
    end if
    
    do 
       read(io, '(a)', iostat=ierr) line
       if (ierr /= 0) then
          exit
       end if
       lt = init_line_text(line)
       call lt%read_value(name, val, is_ok)
       if (is_ok) then
          call self%set_value(name, val)
       end if
    end do
    close(io)
    

    return 
  end subroutine param_read_param_file
  
  !---------------------------------------------------------------------  
  
  subroutine param_set_value(self, name, val)
    class(param), intent(inout) :: self
    character(len=*), intent(in) :: name, val
    
    if (self%verb) then
       write(*,*)trim(name), " <- ", trim(val)
    end if

    if (name == "station_file") then
       self%station_file = val
    else if (name == "polarity_file") then
       self%polarity_file = val
    else if (name == "dc_only") then
       read(val,*) self%dc_only
    else if (name == "sample_prior") then
       read(val,*) self%sample_prior
    else if (name == "use_amp") then
        read(val,*) self%use_amp
    else if (name == "n_procs") then
       read(val,*) self%n_procs
    else if (name == "n_iter") then
       read(val,*) self%n_iter
    else if (name == "n_burn") then
       read(val,*) self%n_burn
    else if (name == "n_interval") then
       read(val,*) self%n_interval
    else if (name == "n_chains") then
       read(val,*) self%n_chains
    else if (name == "n_cool") then
       read(val,*) self%n_cool
    else if (name == "temp_high") then
       read(val,*) self%temp_high
    else
       if (self%verb) then
          write(0,*)"ERROR: Invalid parameter name"
          write(0,*)"        : ", name, "  (?)"
       end if
       stop
    end if
    
    
    return 
  end subroutine param_set_value

  !---------------------------------------------------------------------

  integer function param_get_n_procs(self) result(n_procs)
    class(param), intent(in) :: self
    
    n_procs = self%n_procs

    return 
  end function param_get_n_procs

  !---------------------------------------------------------------------

  integer function param_get_n_iter(self) result(n_iter)
    class(param), intent(in) :: self
    
    n_iter = self%n_iter

    return 
  end function param_get_n_iter

  !---------------------------------------------------------------------

  integer function param_get_n_interval(self) result(n_interval)
    class(param), intent(in) :: self
    
    n_interval = self%n_interval

    return 
  end function param_get_n_interval

  !---------------------------------------------------------------------

  integer function param_get_n_burn(self) result(n_burn)
    class(param), intent(in) :: self
    
    n_burn = self%n_burn

    return 
  end function param_get_n_burn

  !---------------------------------------------------------------------

  integer function param_get_n_chains(self) result(n_chains)
    class(param), intent(in) :: self
    
    n_chains = self%n_chains

    return 
  end function param_get_n_chains

  !---------------------------------------------------------------------
  
  integer function param_get_n_cool(self) result(n_cool)
    class(param), intent(in) :: self
    
    n_cool = self%n_cool

    return 
  end function param_get_n_cool

  !---------------------------------------------------------------------

  logical function param_get_dc_only(self) result(dc_only)
    class(param), intent(in) :: self
    
    dc_only = self%dc_only

    return 
  end function param_get_dc_only

  !---------------------------------------------------------------------

  logical function param_get_sample_prior(self) result(sample_prior)
    class(param), intent(in) :: self
    
    sample_prior = self%sample_prior

    return 
  end function param_get_sample_prior

  !---------------------------------------------------------------------

  logical function param_get_use_amp(self) result(use_amp)
    class(param), intent(in) :: self
    
    use_amp = self%use_amp

    return 
  end function param_get_use_amp
  
  !---------------------------------------------------------------------

  double precision function param_get_temp_high(self) result(temp_high)
    class(param), intent(in) :: self
    
    temp_high = self%temp_high

    return 
  end function param_get_temp_high
  
  !---------------------------------------------------------------------

  character(len=line_max) function param_get_station_file(self) result(station_file)
    class(param), intent(in) :: self
    
    station_file = self%station_file

    return 
  end function param_get_station_file

  !---------------------------------------------------------------------

  character(len=line_max) function param_get_polarity_file(self) result(polarity_file)
    class(param), intent(in) :: self
    
    polarity_file = self%polarity_file

    return 
  end function param_get_polarity_file

  !---------------------------------------------------------------------
  
end module cls_param
