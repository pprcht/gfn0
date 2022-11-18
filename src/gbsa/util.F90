module solvation_solv_util

  implicit none

  character,parameter :: space = ' '
#ifdef _WIN32
  character,parameter :: pathdel = ';'
  character,parameter :: pathsep = '\'
#else
  character,parameter :: pathdel = ':'
  character,parameter :: pathsep = '/'
#endif

contains
!========================================================================================!

  subroutine rdpath(path,arg,fname,ex)
    implicit none
    character(len=*),intent(in)  :: arg
    character(len=*),intent(in)  :: path
    character(len=:),allocatable,intent(out) :: fname
    logical,intent(out),optional :: ex

!* temporary variables
    character(len=:),allocatable :: scratch1
    character(len=:),allocatable :: scratch2
    character(len=:),allocatable :: fpath
    logical :: exist
    integer :: i

    scratch1 = path
    do
      i = index(scratch1,pathdel)
      if (i .eq. 0) then
        scratch2 = scratch1
      else
        scratch2 = scratch1(:i - 1)
        scratch1 = scratch1(i + 1:)
      end if
      fpath = scratch2//pathsep//arg
      inquire (file=fpath,exist=exist)
      if (exist) exit
!     print*,fpath,space,scratch1,exist
      if (i .eq. 0) exit
    end do

!  print*,fpath,exist

    if (exist) fname = fpath
    if (present(ex)) ex = exist

  end subroutine rdpath

!========================================================================================!

  function lowercase(s)
    implicit none
    character(len=*),intent(in) :: s
    character(len=:),allocatable :: sout
    character(len=:),allocatable :: lowerCase
    integer :: ic,i
    character(26),Parameter :: high = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26),Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
    sout = s
    do i = 1,LEN_TRIM(s)
      ic = INDEX(high,s(i:i))
      if (ic > 0) sout(i:i) = low(ic:ic)
    end do
    call move_alloc(sout,lowerCase)
  end function lowercase

!========================================================================================!

end module solvation_solv_util
