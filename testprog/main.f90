program gfn0_main
    use iso_fortran_env, only: wp=>real64
    use testmol
    implicit none
    
    integer :: nat
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)

    integer :: i,j,k,l

!========================================================================================!

    nat = testnat
    allocate(at(nat),xyz(3,nat))
    
    at = testat
    xyz = testxyz


    write(*,*) nat
    write(*,*)
    do i=1,nat
      write(*,'(2x,i3,3x,3f16.6)') at(i),xyz(1:3,i)
    enddo
   

    deallocate(xyz,at)


!=======================================================================================!
end program gfn0_main
