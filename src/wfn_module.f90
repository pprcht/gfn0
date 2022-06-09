
!> The original source code can be found under the GNU LGPL 3.0 license
!> at https://github.com/grimme-lab/xtb

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
   contains
   procedure :: allocate => allocate_wavefunction
   procedure :: deallocate => deallocate_wavefunction
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
end subroutine deallocate_wavefunction

end module wfn_module
