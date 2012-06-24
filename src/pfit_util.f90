#include "util.h"

module mod_util
  ! Some general-purpose utility functions
  implicit none

  interface to_str
     module procedure itos
     module procedure r8tos
     module procedure iarraytos
  end interface

contains

  subroutine assert_msg(val,str,file,line,msg)
    implicit none
    logical :: val
    character(len=*) :: str, file, msg
    integer :: line

    if (.not. val) then
       write (0,'(3a,i0,2a)') 'Error: assertion failed at ', file, ', line ', &
            line, ': ', str
       if (len(msg) > 0) then
          write (0,'(2a)') 'Message: ', msg
       end if
       stop
    end if
  end subroutine assert_msg


  function itos(n)
    ! Convert integer to a string
    implicit none
    integer, intent(in) :: n
    ! Is there any way to return a string of exactly the right length?
    character(len=32)   :: itos

    write (itos,'(i0)') n
  end function itos


  function r8tos(r,fmt)
    ! Convert double-precision to a string
    implicit none
    real*8, intent(in) :: r
    character(len=*), optional, intent(in) :: fmt
    character(len=32)   :: r8tos

    if (present(fmt)) then
       write(r8tos,'(' // trim(adjustl(fmt)) // ')') r
    else
       write (r8tos, '(es25.15)') r
    end if
  end function r8tos

  function iarraytos(v) result(str)
    ! Convert array of integers to a string
    implicit none
    integer, intent(in) :: v(:)

    character(len=2+size(v)*8) :: str
    integer :: i

    str = ' '
    if (size(v) == 0) then
       write (str,'(a)') '[]'
       return
    end if

    write (str,'(a,i0)') '[', v(1)
    do i=2,size(v)
       write (str(index(str,' '):),'(a,i0)') ',', v(i)
    end do
    write (str(index(str,' '):),'(a)') ']'
  end function iarraytos


end module mod_util
