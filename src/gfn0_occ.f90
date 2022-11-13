!================================================================================!
! This file is part of gfn0.
!
! Copyright (C) 2022 Philipp Pracht
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
module gfn0_module
  use iso_fortran_env,only:wp => real64,stdout => output_unit

  use gfn0_types
  use gfn0_parameter
  use gfn0_basisset
  use wfn_module
  use gfn0_qm
  implicit none
  private

  !> from this module
  public :: gfn0_eht_occ

  !> privat params
  real(wp),private,parameter :: autoev = 27.21138505_wp
  real(wp),private,parameter :: evtoau = 1.0_wp / autoev

contains

!========================================================================================!

  subroutine gfn0_eht_occ(nat,at,xyz,nlev,occ,xtbData,basis,cn,dcndr,qat,dqdr, &
     &                wfn,eels,gradients,fail)
!> subroutine gfn0_eht_occ
!> Set up the H0 hamiltonian and overlap,
!> and from that get the "QM" part of the energy
!> Takes an occupation matrix as input and can calculate
!> the energy for multiple (nlev) occupations
!> No Fermi smearing is conducted!
    use wfn_module
    use math_wrapper,only:contract
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer, intent(in)  :: nlev
    real(wp),intent(in) :: occ(basis%nao,nlev)
    real(wp),intent(in) :: cn(:)
    real(wp),intent(in) :: dcndr(:,:,:)
    real(wp),intent(in) :: qat(:)
    real(wp),intent(in) :: dqdr(:,:,:)
    type(TBasisset),intent(in)   :: basis
    type(TxTBData_mod),intent(in) :: xtbData
    type(Twavefunction),intent(inout) :: wfn
    !> OUTPUT
    real(wp),intent(out) :: eels(nlev)
    real(wp),intent(inout) :: gradients(3,nat,nlev)
    logical,intent(out)  :: fail

    !> LOCAL
    integer :: i,j,k
    real(wp) :: intcut,scfconv
    real(wp) :: et,efa,efb,nfoda,nfodb,ga,gb
    integer :: maxshell,nao,naop
    real(wp),allocatable :: selfEnergy(:,:)
    real(wp),allocatable :: dSEdcn(:,:)
    real(wp),allocatable :: dSEdq(:,:)
    real(wp),allocatable :: S(:,:)
    real(wp),allocatable :: H0(:)
    real(wp),allocatable :: dHdcn(:)
    real(wp),allocatable :: dHdq(:)
    real(wp),allocatable :: Pew(:,:)
    real(wp),allocatable :: tmp(:)
    real(wp),allocatable :: tmpgrd(:,:)

    eels = 0.0_wp
    fail = .false.

    !>---  Numerical stuff and cutoffs
    intcut = 25.0_wp - 10.0_wp * log10(xtbData%acc)
    intcut = max(20.0_wp,intcut)
    scfconv = 1.e-6_wp * xtbData%acc

    !>--- allocation
    nao = basis%nao
    naop = basis%nao * (basis%nao + 1) / 2
    allocate (H0(naop),source=0.0_wp)
    allocate (S(nao,nao),source=0.0_wp)
    maxshell = maxval(xtbData%nshell)
    allocate (selfEnergy(maxshell,nat),source=0.0_wp)
    allocate (dSEdcn(maxshell,nat),source=0.0_wp)
    allocate (dSEdq(maxshell,nat),source=0.0_wp)

    et = xtbData%etemp
    wfn%q = qat

    !>--- Setup, AO integrals of S and H0
    call getSelfEnergy(xtbData%hamiltonian,xtbData%nShell,at,cn,wfn%q, &
       & selfEnergy,dSEdcn,dSEdq)
    call build_SH0(xtbData%nShell,xtbData%hamiltonian,selfEnergy, &
          & nat,at,basis,nao,xyz,intcut,S,H0)

    !>--- Diagonalization
    call solve(nao,H0,S,wfn%C,wfn%P,wfn%emo,fail)
    if (fail) return

!==========! GRADIENT(S) !=========!
    allocate (dHdcn(nat),dHdq(nat),pew(nao,nao),tmp(nao), &
    &      source=0.0_wp)
    allocate(tmpgrd(3,nat), source= 0.0_wp)

    do i=1,nlev
      !>--- setup 
      dHdcn = 0.0_wp
      dHdq  = 0.0_wp
      tmp(:) = occ(:,i)    
      tmpgrd(:,:) = gradients(:,:,i)    
      
      !>--- No Fermi Smearing this time!
      call dmat(nao,tmp,wfn%C,wfn%P)
      eels(i) = sum(tmp * wfn%emo) * evtoau
  
      
      !> setup energy weighted density matrix = Pew for gradient calculation
      tmp = wfn%focc * wfn%emo * evtoau
      call dmat(nao,tmp,wfn%C,Pew)
      call build_dSH0(xtbData%nShell,xtbData%hamiltonian,selfEnergy, &
          & dSEdcn,dSEdq,nat,basis,intcut,nao,at,xyz, &
          & wfn%P,Pew,tmpgrd,dHdcn,dHdq)
  
      !>--- Gradient Contraction
      call contract(dcndr,dhdcn,tmpgrd,beta=1.0_wp)
      call contract(dqdr,dhdq,tmpgrd,beta=1.0_wp)
  
      !>--- back to output
      gradients(:,:,i) = tmpgrd(:,:)
    enddo

    !>--- deallocation
    deallocate(tmpgrd)
    deallocate (tmp,Pew,dHdq,dHdcn)
    deallocate (dSEdq,dSEdcn,selfEnergy)
    deallocate (S,H0)

    return
  end subroutine gfn0_eht_occ
!========================================================================================!



!========================================================================================!
end module gfn0_module
