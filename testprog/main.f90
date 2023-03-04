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
program gfn0_main_tester
    use iso_fortran_env, only: wp=>real64,stdout=>output_unit
    use testmol
    use gfn0_module
    use gfn0_interface
    implicit none
    
    integer :: nat
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    integer :: chrg
    integer :: uhf

    integer :: i,j,k,l

!========================================================================================!
   type(TBasisset)    :: basis   !< basis set info
   type(TxTBData_mod) :: xtbData !< main data/parameter frame 
   type(Twavefunction) :: wfn    !< wavefunction data

   type(gfn0_data) :: gdat   !< gfn0 wrapper of types
   type(gfn0_results) :: res   

   real(wp) :: energy
   real(wp),allocatable :: gradient(:,:)

   real(wp),allocatable :: cn(:)
   real(wp),allocatable :: dcndr(:,:,:)

   real(wp),allocatable :: qat(:)
   real(wp),allocatable :: dqdr(:,:,:)

   real(wp) :: ies,edisp,erep,esrb,eel,gnorm

   logical :: fail

   real(wp),allocatable :: occ(:)
   integer :: nao,nel
   integer,allocatable :: active(:)


!========================================================================================!
    fail = .false.

    nat = testnat
    allocate(at(nat),xyz(3,nat))
    
    at = testat
    xyz = testxyz
    uhf = 0
    chrg = 0

    energy = 0.0_wp
    ies    = 0.0_wp
    edisp  = 0.0_wp
    erep   = 0.0_wp
    esrb   = 0.0_wp
    eel    = 0.0_wp
    gnorm  = 0.0_wp
    allocate(gradient(3,nat),source=0.0_wp)


    write(*,*) nat
    write(*,*)
    do i=1,nat
      write(*,'(2x,i3,3x,3f16.6)') at(i),xyz(1:3,i)
    enddo
    call writetestcoord()

    call initGFN0Params(nat, at, xyz, chrg, uhf, basis, xtbData)


    call xtbData%writeinfo(stdout)

!=======================================================================================!
    write(*,*)

!>--- CN and CN gradient
    write(*,*) 'Calculating CN'
    allocate(cn(nat),dcndr(3,nat,nat), source=0.0_wp)
    call gfn0_getCN(nat,at,xyz,cn,dcndr)

!>--- EEQ charges, gradients and IES energy
    write(*,*) 'Calculating EEQ charges'
    allocate(qat(nat),dqdr(3,nat,nat))
    call gfn0_electrostatics(nat, at, xyz, chrg, xtbData, &
    & cn, dcndr, ies, gradient, qat, dqdr)
    write(*,'(3x,a5,1x,f16.8)') 'IES',ies
 
    write(*,*)
    write(*,'(a5,a8,a8)') 'at','CN','q'
    do i=1,nat
    write(*,'(i5,f8.3,f8.3)') at(i), cn(i), qat(i)
    enddo

!>--- D4 two-body dispersion energy
    write(*,*) 
    write(*,*) 'Calculating Dispersion'
    call gfn0_dispersion(nat, at, xyz, chrg, xtbData, &
    &  cn, dcndr, edisp, gradient)
    write(*,'(3x,a5,1x,f16.8)') 'Edisp',edisp

!>--- Repulsion energy
    write(*,*)
    write(*,*) 'Calculating Repulsion'
    call gfn0_repulsion(nat, at, xyz, xtbData%repulsion, erep, gradient)
    write(*,'(3x,a5,1x,f16.8)') 'Erep',erep

!>--- SRB energy
    write(*,*)
    write(*,*) 'Calculating SRB correction'
    call gfn0_shortranged(nat, at, xyz, xtbData%srb, cn, dcndr, &
    &                     esrb, gradient)
    write(*,'(3x,a5,1x,f16.8)') 'Esrb',esrb

!>--- QM part
    write(*,*)
    write(*,*) 'Calculating EHT energy'
    call wfnsetup(xtbData, basis, nat, at, uhf, chrg, wfn)
    call pr_wfn_param(stdout,xtbdata,basis,wfn)
    call gfn0_eht(nat, at, xyz, xtbData, basis, cn, dcndr, qat, dqdr, &
   &                wfn, eel, gradient, fail)
    write(*,'(3x,a5,1x,f16.8)') 'Eel',eel

!>--- Summaryund total energy

   energy = ies+edisp+erep+esrb+eel
   gnorm  = sqrt(sum( gradient**2 ))
   write(*,*)
   write(*,*) "Summary"
   write(*,'(3x,a5,1x,f16.8)') 'IES',ies
   write(*,'(3x,a5,1x,f16.8)') 'Edisp',edisp
   write(*,'(3x,a5,1x,f16.8)') 'Erep',erep
   write(*,'(3x,a5,1x,f16.8)') 'Esrb',esrb
   write(*,'(3x,a5,1x,f16.8)') 'Eel',eel
   write(*,*)'--------------------------'
   write(*,'(3x,a5,1x,f16.8)') 'Etot',energy
   write(*,'(3x,a5,1x,f16.8)') 'gnorm',gnorm


   deallocate(dqdr,qat,dcndr,cn)

!=======================================================================================!
   write(*,*)
   write(*,*) '!===========================================================!'
   write(*,*)
   write(*,*) 'API call'

   energy   = 0.0_wp
   gradient = 0.0_wp 
   
   call gfn0_init(nat,at,xyz,chrg,uhf,gdat)  
   call gfn0_print_summary(stdout,gdat)
   call gfn0_singlepoint(nat,at,xyz,chrg,uhf,gdat,energy,gradient,fail,res)
  
   call res%print(stdout)   
   write(*,'(3x,a5,15x,l)') 'fail?',fail


!=======================================================================================!
   write(*,*)
   write(*,*) '!===========================================================!'
   write(*,*)
   write(*,*) 'API call with GBSA(H2O)'
   energy   = 0.0_wp
   gradient = 0.0_wp

   call gfn0_init(nat,at,xyz,chrg,uhf,gdat,solv='h2o',alpb=.false.)
   call gfn0_singlepoint(nat,at,xyz,chrg,uhf,gdat,energy,gradient,fail,res)
   call gfn0_print_summary(stdout,gdat,res)
   write(*,'(3x,a5,15x,l)') 'fail?',fail


!=======================================================================================!
   write(*,*)
   write(*,*) '!===========================================================!'
   write(*,*)
   write(*,*) 'API call with ALPB(H2O)'
   energy   = 0.0_wp
   gradient = 0.0_wp

   call gfn0_init(nat,at,xyz,chrg,uhf,gdat,solv='h2o',alpb=.true.)
   call gfn0_singlepoint(nat,at,xyz,chrg,uhf,gdat,energy,gradient,fail,res)
   call gfn0_print_summary(stdout,gdat,res)
   write(*,'(3x,a5,15x,l)') 'fail?',fail

!=======================================================================================!
   deallocate(gradient)
   deallocate(xyz,at)
!=======================================================================================!
end program gfn0_main_tester
