module cls_forward
  implicit none

  type forward
     private

     character(len=10) :: data_type = "P_pol"

   contains

  end type forward

  interface forward
     module procedure init_forward
  end interface forward

contains

  !-----------------------------------------------------------------------

  type(forward) function init_forward(data_type) result(self)
    character(*), intent(in) :: data_type

    self%data_type = data_type

  end function init_forward
    
end module cls_forward
