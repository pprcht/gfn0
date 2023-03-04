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
!--------------------------------------------------------------------------------!
!> The original (unmodified) source code can be found under the GNU LGPL 3.0 license
!> Copyright (C) 2019-2020 Sebastian Ehlert
!> at https://github.com/grimme-lab/xtb
!================================================================================!
module wfn_module
!> A module that defines the wavefunction type "wfn"
   use iso_fortran_env, only : wp=>real64
   implicit none

   public :: TWavefunction

   private

   type :: TWavefunction
      integer :: n = 0
      integer :: nel = 0
      integer :: nopen = 0
      integer :: nao = 0
      integer :: nshell = 0
      real(wp),allocatable :: P(:,:)    ! density matrix
      real(wp),allocatable :: q(:)      ! partial charges
      real(wp),allocatable :: qsh(:)    ! shell charges
      real(wp),allocatable :: dipm(:,:) ! dipole moments
      real(wp),allocatable :: qp(:,:)   ! quadrupole moments
      real(wp),allocatable :: wbo(:,:)  ! wiberg bond orders
      integer :: ihomo = 0,ihomoa = 0,ihomob = 0 ! HOMO position
      real(wp) :: efa = 0.0_wp, efb = 0.0_wp
      real(wp),allocatable :: focc(:)   ! fractional occupation
      real(wp),allocatable :: focca(:)  ! for alpha space
      real(wp),allocatable :: foccb(:)  ! for beta space
      real(wp),allocatable :: emo(:)    ! orbital energies
      real(wp),allocatable :: C(:,:)    ! molecular orbitals
      real(wp),allocatable :: S(:,:)    ! overlap (not always allocated)
   contains
   procedure :: allocate => allocate_wavefunction
   procedure :: deallocate => deallocate_wavefunction
   procedure :: refresh_occu
   end type TWavefunction

contains

subroutine allocate_wavefunction(self,n,nshell,nao)
   class(TWavefunction),intent(inout) :: self
   integer,intent(in) :: n,nshell,nao
   self%n = n
   self%nshell = nshell
   self%nao = nao
   self%ihomo = 0
   self%ihomoa = 0
   self%ihomob = 0
   call self%deallocate
   allocate( self%P(nao,nao),  source = 0.0_wp )
   allocate( self%q(n),        source = 0.0_wp )
   allocate( self%qsh(nshell), source = 0.0_wp )
   allocate( self%dipm(3,n),   source = 0.0_wp )
   allocate( self%qp(6,n),     source = 0.0_wp )
   allocate( self%wbo(n,n),    source = 0.0_wp )
   allocate( self%focc(nao),   source = 0.0_wp )
   allocate( self%focca(nao),  source = 0.0_wp )
   allocate( self%foccb(nao),  source = 0.0_wp )
   allocate( self%emo(nao),    source = 0.0_wp )
   allocate( self%C(nao,nao),  source = 0.0_wp )
   !> overlap not automatically allocated!
end subroutine allocate_wavefunction

subroutine deallocate_wavefunction(self)
   class(TWavefunction),intent(inout) :: self
   if(allocated(self%P))    deallocate(self%P)
   if(allocated(self%q))    deallocate(self%q)
   if(allocated(self%qsh))  deallocate(self%qsh)
   if(allocated(self%dipm)) deallocate(self%dipm)
   if(allocated(self%qp))   deallocate(self%qp)
   if(allocated(self%wbo))  deallocate(self%wbo)
   if(allocated(self%focc)) deallocate(self%focc)
   if(allocated(self%focca))deallocate(self%focca)
   if(allocated(self%foccb))deallocate(self%foccb)
   if(allocated(self%emo))  deallocate(self%emo)
   if(allocated(self%C))    deallocate(self%C)
   if(allocated(self%S))    deallocate(self%S)
end subroutine deallocate_wavefunction

subroutine refresh_occu(self,nel,uhf)
!(ndim,   nel,    nopen,    ihomoa,    ihomob,    focca,    foccb)
!(wfn%nao,wfn%nel,wfn%nopen,wfn%ihomoa,wfn%ihomob,wfn%focca,wfn%foccb)
   class(TWavefunction),intent(inout) :: self
   integer,intent(in)  :: nel  !> total # of electrons
   integer,intent(in)  :: uhf  !> uhf input parameter
    integer  :: nopen
    integer  :: ndim
    integer  :: ihomoa
    integer  :: ihomob
    integer,allocatable  :: focc(:)
    integer  :: i,na,nb,ihomo

    self%nel = nel
    self%nopen = uhf
    if (self%nopen == 0 .and. mod(self%nel,2) /= 0) self%nopen = 1
        
    nopen = self%nopen
    ndim = self%nao  
    allocate(focc(ndim), source=0)
    focc = 0
    self%focca = 0.0_wp
    self%foccb = 0.0_wp
!>--- even nel
    if (mod(nel,2) .eq. 0) then
      ihomo = nel / 2
      do i = 1,ihomo
        focc(i) = 2
      end do
      if (2 * ihomo .ne. nel) then
        ihomo = ihomo + 1
        focc(ihomo) = 1
        if (nopen .eq. 0) nopen = 1
      end if
      if (nopen .gt. 1) then
        do i = 1,nopen / 2
          focc(ihomo - i + 1) = focc(ihomo - i + 1) - 1
          focc(ihomo + i) = focc(ihomo + i) + 1
        end do
      end if
!>--- odd nel
    else
      na = nel / 2 + (nopen - 1) / 2 + 1
      nb = nel / 2 - (nopen - 1) / 2
      do i = 1,na
        focc(i) = focc(i) + 1
      end do
      do i = 1,nb
        focc(i) = focc(i) + 1
      end do
    end if

    do i = 1,ndim
      if (focc(i) .eq. 2) then
        self%focca(i) = 1.0d0
        self%foccb(i) = 1.0d0
      end if
      if (focc(i) .eq. 1) self%focca(i) = 1.0d0
    end do

    ihomoa = 0
    ihomob = 0
    do i = 1,ndim
      if (self%focca(i) .gt. 0.99) ihomoa = i
      if (self%foccb(i) .gt. 0.99) ihomob = i
    end do

    self%ihomoa = ihomoa
    self%ihomo  = ihomoa
    self%ihomob = ihomob
    self%focc   = self%focca + self%foccb
    deallocate(focc)

end subroutine refresh_occu


end module wfn_module
