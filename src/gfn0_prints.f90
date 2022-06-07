

module gfn0_prints
   use iso_fortran_env, only: wp=>real64, stdout=>output_unit
   use gfn0_types
   use gfn0_parameter
   use gfn0_basisset
   use wavefunction_module
   implicit none
   public

contains
!========================================================================================!
subroutine pr_wfn_param(xtbData,basis,wfn)
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
  
   intcut=25.0_wp-10.0_wp*log10(xtbData%acc)
   intcut=max(20.0_wp,intcut)
   
   write(stdout,*)
      write(stdout,'(/,10x,51("."))')
      write(stdout,'(10x,":",22x,a,22x,":")') "SETUP"
      write(stdout,'(10x,":",49("."),":")')
      write(stdout,intfmt) "# basis functions  ",basis%nbf
      write(stdout,intfmt) "# atomic orbitals  ",basis%nao
      write(stdout,intfmt) "# shells           ",basis%nshell
      write(stdout,intfmt) "# electrons        ",wfn%nel
      write(stdout,dblfmt) "electronic temp.   ",xtbdata%etemp, "K   "
      write(stdout,dblfmt) "accuracy           ",xtbData%acc, "    "
      write(stdout,scifmt) "-> integral cutoff ",intcut,  "    "
      write(stdout,'(10x,51("."))')
    write(stdout,*)
    return
end subroutine pr_wfn_param
!========================================================================================!

end module gfn0_prints
