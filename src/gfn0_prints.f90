!================================================================================!
! This file is part of gfn0.
!
! Copyright (C) 2022-2023 Philipp Pracht
!
! gfn0 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! gfn0 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with gfn0.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!
module gfn0_prints
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use gfn0_types
  use gfn0_parameter
  use gfn0_basisset
  use wfn_module
  implicit none
  public

  !>  convert Hartree to eV and back
  real(wp),parameter,private :: autoev = 27.21138505_wp
  real(wp),parameter,private :: evtoau = 1.0_wp/autoev


contains
!========================================================================================!
  subroutine pr_wfn_param(iunit,xtbData,basis,wfn)
    integer :: iunit
    type(TxTBData_mod)  :: xtbData
    type(Tbasisset)     :: basis
    type(TWavefunction) :: wfn
    real(wp) :: intcut
    character(len=*),parameter :: scifmt = &
                                  '(10x,":",2x,a,e22.7,1x,a,1x,":")'
    character(len=*),parameter :: dblfmt = &
                                  '(10x,":",2x,a,f18.7,5x,a,1x,":")'
    character(len=*),parameter :: intfmt = &
                                  '(10x,":",2x,a,i18,      10x,":")'
    character(len=*),parameter :: chrfmt = &
                                  '(10x,":",2x,a,a18,      10x,":")'

    intcut = 25.0_wp - 10.0_wp * log10(xtbData%acc)
    intcut = max(20.0_wp,intcut)

    write (iunit,'(/,10x,51("."))')
    write (iunit,'(10x,":",22x,a,22x,":")') "SETUP"
    write (iunit,'(10x,":",49("."),":")')
    write (iunit,intfmt) "# basis functions  ",basis%nbf
    write (iunit,intfmt) "# atomic orbitals  ",basis%nao
    write (iunit,intfmt) "# shells           ",basis%nshell
    write (iunit,intfmt) "# electrons        ",wfn%nel
    write (iunit,dblfmt) "electronic temp.   ",xtbdata%etemp,"K   "
    write (iunit,dblfmt) "accuracy           ",xtbData%acc,"    "
    write (iunit,scifmt) "-> integral cutoff ",intcut,"    "
    write (iunit,'(10x,51("."))')
    write (iunit,*)
    return
  end subroutine pr_wfn_param
!========================================================================================!

subroutine pr_orbital_eigenvalues(iunit,wfn,range)
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: range
   type(TWavefunction),intent(in) :: wfn
   character(len=*),parameter :: hlfmt = '(    a24,f21.7,1x,"Eh",f18.4,1x,"eV")'
   integer :: maxorb,minorb,iorb
   real(wp) :: gap

   minorb = max(wfn%ihomoa - (range+1), 1)
   maxorb = min(wfn%ihomoa +  range, wfn%nao)
   gap = wfn%emo(wfn%ihomoa+1) - wfn%emo(wfn%ihomoa)

   write(iunit,'(a)')
   write(iunit,'(a10,a14,a21,a21)') "#","Occupation","Energy/Eh","Energy/eV"
   write(iunit,'(6x,61("-"))')
   if (minorb .gt. 1) then
      call write_line(1,wfn%focc,wfn%emo,wfn%ihomo)
      if (minorb .gt. 2) &
         write(iunit,'(a10,a14,a21,a21)') "...","...","...","..."
   endif
   do iorb = minorb,maxorb
      call write_line(iorb,wfn%focc,wfn%emo,wfn%ihomo)
   enddo
   if (maxorb .lt. wfn%nao) then
      if (maxorb .lt. wfn%nao-1) then
         if (wfn%focc(maxorb) > 1.0e-7_wp) then
            write(iunit,'(a10,a14,a21,a21)') "...","...","...","..."
         else
            write(iunit,'(a10,a14,a21,a21)') "...",   "","...","..."
         endif
      endif
      call write_line(wfn%nao,wfn%focc,wfn%emo,wfn%ihomo)
   endif
   write(iunit,'(6x,61("-"))')
   write(iunit,hlfmt) "HL-Gap",gap*evtoau,gap
   !write(iunit,hlfmt) "Fermi-level",(wfn%efa+wfn%efb)/2*evtoau,(wfn%efa+wfn%efb)/2
contains
subroutine write_line(iorb,focc,emo,ihomo)
   integer, intent(in) :: iorb
   integer, intent(in) :: ihomo
   real(wp),intent(in) :: focc(:)
   real(wp),intent(in) :: emo (:)
   character(len=*),parameter :: mofmt = '(i10,f14.4,f21.7,f21.4)'
   character(len=*),parameter :: vofmt = '(i10,14x,  f21.7,f21.4)'
   if (focc(iorb) < 1.0e-7_wp) then
      write(iunit,vofmt,advance='no') iorb,             emo(iorb)*evtoau, emo(iorb)
   else
      write(iunit,mofmt,advance='no') iorb, focc(iorb), emo(iorb)*evtoau, emo(iorb)
   endif
   if (iorb == ihomo) then
      write(iunit,'(1x,"(HOMO)")')
   elseif (iorb == ihomo+1) then
      write(iunit,'(1x,"(LUMO)")')
   else
      write(iunit,'(a)')
   endif
end subroutine write_line
end subroutine pr_orbital_eigenvalues

!========================================================================================!

end module gfn0_prints
