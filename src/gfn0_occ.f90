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
module gfn0_occ
  use iso_fortran_env,only:wp => real64,stdout => output_unit,stderr => error_unit

  use gfn0_types
  use gfn0_parameter
  use gfn0_basisset
  use wfn_module
  use gfn0_qm
  use gfn0_math_wrapper,only:contract
  implicit none
  private

  !> from this module
  public :: gfn0_eht_occ
  public :: generate_config

  !> privat params
  real(wp),private,parameter :: autoev = 27.21138505_wp
  real(wp),private,parameter :: evtoau = 1.0_wp / autoev

contains

!========================================================================================!

  subroutine gfn0_eht_occ(nat,at,xyz,occ,xtbData,basis,cn,dcndr,qat,dqdr, &
     &                wfn,eel,gradient,fail)
!> subroutine gfn0_eht_occ
!> Set up the H0 hamiltonian and overlap,
!> and from that get the "QM" part of the energy
!> Takes an occupation matrix as input and can calculate
!> the energy for given occupations
!> No Fermi smearing is conducted!
    !> INPUT
    type(TBasisset),intent(in)   :: basis
    type(TxTBData_mod),intent(in) :: xtbData
    type(Twavefunction),intent(inout) :: wfn
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(in) :: occ(basis%nao)
    real(wp),intent(in) :: cn(:)
    real(wp),intent(in) :: dcndr(:,:,:)
    real(wp),intent(in) :: qat(:)
    real(wp),intent(in) :: dqdr(:,:,:)
    !> OUTPUT
    real(wp),intent(out) :: eel
    real(wp),intent(inout) :: gradient(3,nat)
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

    eel = 0.0_wp
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
    if(allocated(wfn%S)) wfn%S = S

    !>--- Diagonalization
    call solve(nao,H0,S,wfn%C,wfn%P,wfn%emo,fail)
    if (fail) return

!==========! GRADIENT(S) !=========!
    allocate (dHdcn(nat),dHdq(nat),pew(nao,nao),tmp(nao), source=0.0_wp)

      !>--- setup
      dHdcn = 0.0_wp
      dHdq = 0.0_wp
      tmp(:) = occ(:)

      !>--- No Fermi Smearing this time!
      call dmat(nao,tmp,wfn%C,wfn%P)
      eel = sum(tmp * wfn%emo) * evtoau

      !> setup energy weighted density matrix = Pew for gradient calculation
      tmp(:) = occ(:) * wfn%emo * evtoau
      Pew = 0.0_wp
      call dmat(nao,tmp,wfn%C,Pew)
      call build_dSH0(xtbData%nShell,xtbData%hamiltonian,selfEnergy, &
          & dSEdcn,dSEdq,nat,basis,intcut,nao,at,xyz, &
          & wfn%P,Pew,gradient,dHdcn,dHdq)

      !>--- Gradient Contraction
      call contract(dcndr,dhdcn,gradient,beta=1.0_wp)
      call contract(dqdr,dhdq,gradient,beta=1.0_wp)

    !>--- deallocation
    deallocate (tmp,Pew,dHdcn,dHdq)
    deallocate (dSEdq,dSEdcn,selfEnergy)
    deallocate (S,H0)

    return
  end subroutine gfn0_eht_occ
!========================================================================================!

  subroutine generate_config(nel,nao,occ,active)
    implicit none
    !> INPUT
    integer,intent(in) :: nel
    integer,intent(in) :: nao
    integer,intent(in) :: active(:)
    !> OUTPUT
    real(wp),intent(out) :: occ(nao)
    !> LOCAL
    integer :: i,l,dl
    integer :: nel_active,dnel

    occ = 0.0_wp

    l = size(active,1)
    nel_active = sum(active)

    if ((nel_active > nel) .or. (l > nao)) then
      error stop 'error in generate_config()'
    end if
 
    do i=1,l
      if((active(i) > 2) .or. (active(i)<0))then
       write (stderr,*) 'error in generate_config()'
       write (stderr,*) 'unphysical occupation in ',i
       error stop
      endif
    enddo


    dnel = nel - nel_active
    if (MOD(dnel,2) .ne. 0) then
      write (stderr,*) 'error in generate_config()'
      write (stderr,*) 'mismatching number of electrons'
      error stop
    end if

    dl = dnel / 2
    if ((dl + l) > nao) then
      error stop 'error in generate_config()'
    end if
    !>--- fill doubly occupied levels first
    do i = 1,dl
      occ(i) = 2.0_wp
    end do
    !>--- add the active (user-set) configuration
    do i = dl + 1,dl + l
      occ(i) = float(active(i - dl))
    end do

  end subroutine generate_config

!========================================================================================!
end module gfn0_occ
