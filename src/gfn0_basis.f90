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
!> The original source code can be found under the GNU LGPL 3.0 license
!> at https://github.com/grimme-lab/xtb
!================================================================================!
module gfn0_basisset
   use iso_fortran_env, only : wp=>real64

   use gfn0_types
   use slater_module

   implicit none

   public :: TBasisset
   public :: newBasisset

   private

   real(wp),private,parameter :: pi=3.1415926535897932385_wp

   type :: TBasisset
      integer  :: maxao = 0
      integer  :: n = 0
      integer  :: nbf = 0
      integer  :: nao = 0
      integer  :: nshell = 0
      ! ------------------------------------------------------------------------
      !  This variables describe the ATOMS
      integer, allocatable :: caoshell(:,:)  ! atom number -> basis function
      integer, allocatable :: saoshell(:,:)  ! atom number -> AO
      integer, allocatable :: fila(:,:)      ! atom number -> basis function
      integer, allocatable :: fila2(:,:)     ! atom number -> AO
      integer, allocatable :: shells(:,:)    ! atom number -> shell
      ! ------------------------------------------------------------------------
      !  This variables describe the indivdual SHELLS
      integer, allocatable :: lsh(:)   ! shell -> azimudal quantum number
      integer, allocatable :: ash(:)   ! shell -> atom number
      real(wp),allocatable :: zeta(:)  ! shell -> exponent
      real(wp),allocatable :: level(:) ! shell -> level energy
      real(wp),allocatable :: minalp(:)! shell -> most diffuse exponent
      integer, allocatable :: sh2bf(:,:) ! shell -> BF
      integer, allocatable :: sh2ao(:,:) ! shell -> AO
      integer, allocatable :: valsh(:) ! shell -> shell type
      ! ------------------------------------------------------------------------
      !  This variables describe the SPHERICAL atomic orbitals
      real(wp),allocatable :: aoexp(:)    ! AO -> exponent
      integer, allocatable :: ao2sh(:)    ! AO -> shell
      integer, allocatable :: lao2(:)     ! AO -> azimudal quantum number
      integer, allocatable :: aoat2(:)    ! AO -> atom number
      real(wp),allocatable :: hdiag2(:)   ! AO -> level energy
      integer, allocatable :: valao2(:)   ! AO -> AO type
      ! ------------------------------------------------------------------------
      !  This variables describe the CARTESIAN basis functions
      integer, allocatable :: primcount(:) ! BF -> primitive count
      real(wp),allocatable :: hdiag(:)     ! BF -> level energy
      integer, allocatable :: nprim(:)     ! BF -> primitive number
      integer, allocatable :: aoat(:)      ! BF -> atom
      integer, allocatable :: valao(:)     ! BF -> BF type
      integer, allocatable :: lao(:)       ! BF -> azimudal quantum number
      ! ------------------------------------------------------------------------
      !  This variables decribe the PRIMITIV basis functions
      real(wp),allocatable :: alp(:)  ! primitive -> primitive exponent
      real(wp),allocatable :: cont(:) ! primitive -> contraction coeffient
   contains
   procedure :: allocate => allocate_basisset
   procedure :: deallocate => deallocate_basisset
   end type TBasisset

contains

!========================================================================================!
!========================================================================================!

subroutine allocate_basisset(self,n,nbf,nao,nshell)
   implicit none
   class(TBasisset),intent(inout) :: self
   integer,intent(in) :: n,nbf,nao,nshell
   self%n=n
   self%nbf=nbf
   self%nao=nao
   self%nshell=nshell
   call self%deallocate
   allocate( self%shells(2,n),    source = 0 )
   allocate( self%sh2ao(2,nshell),source = 0 )
   allocate( self%sh2bf(2,nshell),source = 0 )
   allocate( self%minalp(nshell), source = 0.0_wp )
   allocate( self%level(nshell),  source = 0.0_wp )
   allocate( self%zeta(nshell),   source = 0.0_wp )
   allocate( self%valsh(nshell),  source = 0 )
   allocate( self%hdiag(nbf),     source = 0.0_wp )
   allocate( self%alp(9*nbf),     source = 0.0_wp )
   allocate( self%cont(9*nbf),    source = 0.0_wp )
   allocate( self%hdiag2(nao),    source = 0.0_wp )
   allocate( self%aoexp(nao),     source = 0.0_wp )
   allocate( self%ash(nao),       source = 0 )
   allocate( self%lsh(nao),       source = 0 )
   allocate( self%ao2sh(nao),     source = 0 )
   allocate( self%nprim(nbf),     source = 0 )
   allocate( self%primcount(nbf), source = 0 )
   allocate( self%caoshell(5,n),  source = 0 )
   allocate( self%saoshell(5,n),  source = 0 )
   allocate( self%fila(2,n),      source = 0 )
   allocate( self%fila2(2,n),     source = 0 )
   allocate( self%lao(nbf),       source = 0 )
   allocate( self%aoat(nbf),      source = 0 )
   allocate( self%valao(nbf),     source = 0 )
   allocate( self%lao2(nao),      source = 0 )
   allocate( self%aoat2(nao),     source = 0 )
   allocate( self%valao2(nbf),    source = 0 )
end subroutine allocate_basisset

subroutine deallocate_basisset(self)
   implicit none
   class(TBasisset),intent(inout) :: self
   if(allocated(self%shells))    deallocate(self%shells)
   if(allocated(self%sh2ao))     deallocate(self%sh2ao)
   if(allocated(self%sh2bf))     deallocate(self%sh2bf)
   if(allocated(self%minalp))    deallocate(self%minalp)
   if(allocated(self%valsh))     deallocate(self%valsh)
   if(allocated(self%hdiag))     deallocate(self%hdiag)
   if(allocated(self%alp))       deallocate(self%alp)
   if(allocated(self%cont))      deallocate(self%cont)
   if(allocated(self%hdiag2))    deallocate(self%hdiag2)
   if(allocated(self%aoexp))     deallocate(self%aoexp)
   if(allocated(self%ash))       deallocate(self%ash)
   if(allocated(self%lsh))       deallocate(self%lsh)
   if(allocated(self%ao2sh))     deallocate(self%ao2sh)
   if(allocated(self%nprim))     deallocate(self%nprim)
   if(allocated(self%primcount)) deallocate(self%primcount)
   if(allocated(self%caoshell))  deallocate(self%caoshell)
   if(allocated(self%saoshell))  deallocate(self%saoshell)
   if(allocated(self%fila))      deallocate(self%fila)
   if(allocated(self%fila2))     deallocate(self%fila2)
   if(allocated(self%lao))       deallocate(self%lao)
   if(allocated(self%lao2))      deallocate(self%lao2)
   if(allocated(self%aoat))      deallocate(self%aoat2)
   if(allocated(self%valao))     deallocate(self%valao)
   if(allocated(self%valao2))    deallocate(self%valao2)
   if(allocated(self%zeta))      deallocate(self%zeta)
   if(allocated(self%level))     deallocate(self%level)
end subroutine deallocate_basisset

!========================================================================================!
subroutine newBasisset(xtbData,n,at,basis,ok)
   type(TxTBData_mod), intent(in) :: xtbData
   type(TBasisset),intent(inout) :: basis
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   logical, intent(out) :: ok

   integer  :: elem,valao
   integer  :: i,j,m,l,iat,ati,ish,ibf,iao,ipr,p,nprim,thisprimR,idum,npq,npqR,pqn
   real(wp) :: a(10),c(10),zeta,k1,k2,split1,pp,zqfR,zcnfR,qi,level
   real(wp) :: aR(10),cR(10),ss
   real(wp) :: as(10),cs(10)
   integer :: info

   call xbasis0(xtbData,n,at,basis)

   basis%hdiag(1:basis%nbf)=1.d+42

   ibf=0
   iao=0
   ipr=0
   ish=0
   ok=.true.

   atoms: do iat=1,n
      ati = at(iat)
      basis%shells(1,iat)=ish+1
      basis%fila  (1,iat)=ibf+1
      basis%fila2 (1,iat)=iao+1
      shells: do m=1,xtbData%nShell(ati)
         ish = ish+1
         ! principle QN
         npq=xtbData%hamiltonian%principalQuantumNumber(m,ati)
         l=xtbData%hamiltonian%angShell(m,ati)

         level = xtbData%hamiltonian%selfEnergy(m,ati)
         zeta  = xtbData%hamiltonian%slaterExponent(m,ati)
         valao = xtbData%hamiltonian%valenceShell(m,ati)
         if (valao /= 0) then
            nprim = xtbData%hamiltonian%numberOfPrimitives(m,ati)
         else
            thisprimR = xtbData%hamiltonian%numberOfPrimitives(m,ati)
         end if

         basis%lsh(ish) = l
         basis%ash(ish) = iat
         basis%sh2bf(1,ish) = ibf
         basis%sh2ao(1,ish) = iao
         basis%caoshell(m,iat)=ibf
         basis%saoshell(m,iat)=iao

         ! add new shellwise information, for easier reference
         basis%level(ish) = level
         basis%zeta (ish) = zeta
         basis%valsh(ish) = valao

         ! H-He
         if(l.eq.0.and.ati.le.2.and.valao/=0)then
            ! s
            call slaterToGauss(nprim, npq, l, zeta, a, c, .true., info)
            basis%minalp(ish) = minval(a(:nprim))

            ibf =ibf+1
            basis%primcount(ibf) = ipr
            basis%valao    (ibf) = valao
            basis%aoat     (ibf) = iat
            basis%lao      (ibf) = 1
            basis%nprim    (ibf) = nprim
            basis%hdiag    (ibf) = level

            do p=1,nprim
               ipr=ipr+1
               basis%alp (ipr)=a(p)
               basis%cont(ipr)=c(p)
            enddo

            iao = iao+1
            basis%valao2(iao) = valao
            basis%aoat2 (iao) = iat
            basis%lao2  (iao) = 1
            basis%hdiag2(iao) = level
            basis%aoexp (iao) = zeta
            basis%ao2sh (iao) = ish
         endif

         if(l.eq.0.and.ati.le.2.and.valao==0)then
            ! diff s
            call slaterToGauss(thisprimR, npq, l, zeta, aR, cR, .true., info)
            call atovlp(0,nprim,thisprimR,a,aR,c,cR,ss)
            basis%minalp(ish) = min(minval(a(:nprim)),minval(aR(:thisprimR)))

            ibf =ibf+1
            basis%primcount(ibf) = ipr
            basis%valao    (ibf) = valao
            basis%aoat     (ibf) = iat
            basis%lao      (ibf) = 1
            basis%nprim    (ibf) = thisprimR+nprim
            basis%hdiag    (ibf) = level

            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               basis%alp (ipr)=aR(p)
               basis%cont(ipr)=cR(p)
            enddo
            do p=1,nprim
               ipr=ipr+1
               basis%alp (ipr)=a(p)
               basis%cont(ipr)=-ss*c(p)
            enddo
            call atovlp(0,basis%nprim(ibf),basis%nprim(ibf), &
               &        basis%alp(idum),basis%alp(idum), &
               &        basis%cont(idum),basis%cont(idum),ss)
            do p=1,basis%nprim(ibf)
               basis%cont(idum-1+p)=basis%cont(idum-1+p)/sqrt(ss)
            enddo

            iao = iao+1
            basis%valao2(iao) = valao
            basis%aoat2 (iao) = iat
            basis%lao2  (iao) = 1
            basis%hdiag2(iao) = level
            basis%aoexp (iao) = zeta
            basis%ao2sh (iao) = ish
         endif

         ! p polarization
         if(l.eq.1.and.ati.le.2)then
            call slaterToGauss(nprim, npq, l, zeta, a, c, .true., info)
            basis%minalp(ish) = minval(a(:nprim))
            do j=2,4
               ibf=ibf+1
               basis%primcount(ibf) = ipr
               basis%aoat     (ibf) = iat
               basis%lao      (ibf) = j
               basis%valao    (ibf) = -valao
               basis%nprim    (ibf) = nprim
               basis%hdiag    (ibf) = level

               do p=1,nprim
                  ipr=ipr+1
                  basis%alp (ipr)=a(p)
                  basis%cont(ipr)=c(p)
               enddo

               iao = iao+1
               basis%valao2(iao) = -valao
               basis%aoat2 (iao) = iat
               basis%lao2  (iao) = j
               basis%hdiag2(iao) = level
               basis%aoexp (iao) = zeta
               basis%ao2sh (iao) = ish
            enddo
         endif

         ! general sp
         if(l.eq.0.and.ati.gt.2 .and. valao/=0)then
            ! s
            call slaterToGauss(nprim, npq, l, zeta, as, cs, .true., info)
            basis%minalp(ish) = minval(as(:nprim))

            ibf=ibf+1
            basis%primcount(ibf) = ipr
            basis%valao    (ibf) = valao
            basis%aoat     (ibf) = iat
            basis%lao      (ibf) = 1
            basis%nprim    (ibf) = nprim
            basis%hdiag    (ibf) = level

            do p=1,nprim
               ipr=ipr+1
               basis%alp (ipr)=as(p)
               basis%cont(ipr)=cs(p)
            enddo

            iao = iao+1
            basis%valao2(iao) = valao
            basis%aoat2 (iao) = iat
            basis%lao2  (iao) = 1
            basis%hdiag2(iao) = level
            basis%aoexp (iao) = zeta
            basis%ao2sh (iao) = ish
         endif
         ! p
         if(l.eq.1.and.ati.gt.2)then
            call slaterToGauss(nprim, npq, l, zeta, a, c, .true., info)
            basis%minalp(ish) = minval(a(:nprim))
            do j=2,4
               ibf=ibf+1
               basis%primcount(ibf) = ipr
               basis%valao    (ibf) = valao
               basis%aoat     (ibf) = iat
               basis%lao      (ibf) = j
               basis%nprim    (ibf) = nprim
               basis%hdiag    (ibf) = level

               do p=1,nprim
                  ipr=ipr+1
                  basis%alp (ipr)=a(p)
                  basis%cont(ipr)=c(p)
               enddo

               iao = iao+1
               basis%valao2(iao) = valao
               basis%aoat2 (iao) = iat
               basis%lao2  (iao) = j
               basis%hdiag2(iao) = level
               basis%aoexp (iao) = zeta
               basis%ao2sh (iao) = ish
            enddo
         endif

         ! DZ s
         if(l.eq.0 .and. ati > 2 .and. valao==0)then
            call slaterToGauss(thisprimR, npq, l, zeta, aR, cR, .true., info)
            call atovlp(0,nprim,thisprimR,as,aR,cs,cR,ss)
            basis%minalp(ish) = min(minval(as(:nprim)),minval(aR(:thisprimR)))

            ibf=ibf+1
            basis%primcount(ibf) = ipr
            basis%valao    (ibf) = valao
            basis%aoat     (ibf) = iat
            basis%lao      (ibf) = 1
            basis%nprim    (ibf) = thisprimR+nprim
            basis%hdiag    (ibf) = level

            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               basis%alp (ipr)=aR(p)
               basis%cont(ipr)=cR(p)
            enddo
            do p=1,nprim
               ipr=ipr+1
               basis%alp (ipr)=as(p)
               basis%cont(ipr)=-ss*cs(p)
            enddo
            call atovlp(0,basis%nprim(ibf),basis%nprim(ibf), &
               &        basis%alp(idum),basis%alp(idum), &
               &        basis%cont(idum),basis%cont(idum),ss)
            do p=1,basis%nprim(ibf)
               basis%cont(idum-1+p)=basis%cont(idum-1+p)/sqrt(ss)
            enddo

            iao = iao+1
            basis%valao2(iao) = valao
            basis%aoat2 (iao) = iat
            basis%lao2  (iao) = 1
            basis%hdiag2(iao) = level
            basis%aoexp (iao) = zeta
            basis%ao2sh (iao) = ish
         endif

         ! d
         if(l.eq.2)then
            call set_d_function(basis,iat,ish,iao,ibf,ipr, &
               &                npq,l,nprim,zeta,level,valao)
         endif

         ! f
         if(l.eq.3)then
            call set_f_function(basis,iat,ish,iao,ibf,ipr, &
               &                npq,l,nprim,zeta,level,1)
         endif

         basis%sh2bf(2,ish) = ibf-basis%sh2bf(1,ish)
         basis%sh2ao(2,ish) = iao-basis%sh2ao(1,ish)
      enddo shells
      basis%shells(2,iat)=ish
      basis%fila  (2,iat)=ibf
      basis%fila2 (2,iat)=iao
   enddo atoms

   ok = all(basis%alp(:ipr) > 0.0_wp) .and. basis%nbf == ibf .and. basis%nao == iao

end subroutine newBasisset


! ========================================================================
!> determine basisset limits
subroutine xbasis0(xtbData,n,at,basis)
   type(TxTBData_mod), intent(in) :: xtbData
   type(TBasisset),intent(inout) :: basis
   integer,intent(in)  :: n
   integer,intent(in)  :: at(n)
   integer :: nbf
   integer :: nao
   integer :: nshell

   integer i,j,k,l

   call dim_basis(xtbData,n,at,nshell,nao,nbf)

   call basis%allocate(n,nbf,nao,nshell)

end subroutine xbasis0

subroutine dim_basis(xtbData,n,at,nshell,nao,nbf)
   type(TxTBData_mod), intent(in) :: xtbData
   integer,intent(in)  :: n
   integer,intent(in)  :: at(n)
   integer,intent(out) :: nshell
   integer,intent(out) :: nao
   integer,intent(out) :: nbf

   integer i,j,k,l

   nao=0
   nbf=0
   nshell=0

   do i=1,n
      k=0
      do j=1,xtbData%nShell(at(i))
         l = xtbData%hamiltonian%angShell(j,at(i))
         k = k + 1
         nshell=nshell+1
         select case(l)
         case(0) ! s
            nbf = nbf+1
            nao = nao+1
         case(1) ! p
            nbf = nbf+3
            nao = nao+3
         case(2) ! d
            nbf = nbf+6
            nao = nao+5
         case(3) ! f
            nbf = nbf+10
            nao = nao+7
         case(4) ! g
            nbf = nbf+15
            nao = nao+9
         end select
      enddo
      if(k.eq.0) then
         write(*,*) 'no basis found for atom', i,' Z= ',at(i)
         error stop 
      endif
   enddo

end subroutine dim_basis


! ------------------------------------------------------------------------
!  Helper functions

subroutine atovlp(l,npri,nprj,alpa,alpb,conta,contb,ss)
   integer l,npri,nprj
   real(wp) alpa(*),alpb(*)
   real(wp) conta(*),contb(*)
   real(wp) ss

   integer ii,jj
   real(wp) ab,s00,sss,ab05

   SS=0.0_wp
   do ii=1,npri
      do jj=1,nprj
         ab =1./(alpa(ii)+alpb(jj))
         s00=(pi*ab)**1.50_wp
         if(l.eq.0)then
            sss=s00
         endif
         if(l.eq.1)then
            ab05=ab*0.5_wp
            sss=s00*ab05
         endif
         SS=SS+SSS*conta(ii)*contb(jj)
      enddo
   enddo

end subroutine atovlp

subroutine set_s_function(basis,iat,ish,iao,ibf,ipr,npq,l,nprim,zeta,level,valao)
   type(TBasisset), intent(inout) :: basis
   integer, intent(in)    :: iat
   integer, intent(in)    :: ish
   integer, intent(inout) :: ibf
   integer, intent(inout) :: iao
   integer, intent(inout) :: ipr
   integer, intent(in)    :: npq
   integer, intent(in)    :: l
   integer, intent(in)    :: nprim
   integer, intent(in)    :: valao
   real(wp),intent(in)    :: zeta
   real(wp),intent(in)    :: level
   integer  :: p
   real(wp) :: alp(10),cont(10)
   integer :: info

   call slaterToGauss(nprim, npq, l, zeta, alp, cont, .true., info)
   basis%minalp(ish) = minval(alp(:nprim))

   ibf = ibf+1
   basis%primcount(ibf) = ipr
   basis%valao    (ibf) = valao
   basis%aoat     (ibf) = iat
   basis%lao      (ibf) = 1
   basis%nprim    (ibf) = nprim
   basis%hdiag    (ibf) = level

   do p=1,nprim
      ipr = ipr+1
      basis%alp (ipr)=alp (p)
      basis%cont(ipr)=cont(p)
   enddo

   iao = iao+1
   basis%valao2(iao) = valao
   basis%aoat2 (iao) = iat
   basis%lao2  (iao) = 1
   basis%hdiag2(iao) = level
   basis%aoexp (iao) = zeta
   basis%ao2sh (iao) = ish
end subroutine set_s_function

subroutine set_p_function(basis,iat,ish,iao,ibf,ipr,npq,l,nprim,zeta,level,valao)
   type(TBasisset), intent(inout) :: basis
   integer, intent(in)    :: iat
   integer, intent(in)    :: ish
   integer, intent(inout) :: ibf
   integer, intent(inout) :: iao
   integer, intent(inout) :: ipr
   integer, intent(in)    :: npq
   integer, intent(in)    :: l
   integer, intent(in)    :: nprim
   integer, intent(in)    :: valao
   real(wp),intent(in)    :: zeta
   real(wp),intent(in)    :: level
   integer  :: j,p
   real(wp) :: alp(10),cont(10)
   integer :: info

   call slaterToGauss(nprim, npq, l, zeta, alp, cont, .true., info)
   basis%minalp(ish) = minval(alp(:nprim))

   do j = 2, 4

      ibf = ibf+1
      basis%primcount(ibf) = ipr
      basis%valao    (ibf) = valao
      basis%aoat     (ibf) = iat
      basis%lao      (ibf) = j
      basis%nprim    (ibf) = nprim
      basis%hdiag    (ibf) = level

      do p=1,nprim
         ipr = ipr+1
         basis%alp (ipr)=alp (p)
         basis%cont(ipr)=cont(p)
      enddo

      iao = iao+1
      basis%valao2(iao) = valao
      basis%aoat2 (iao) = iat
      basis%lao2  (iao) = j
      basis%hdiag2(iao) = level
      basis%aoexp (iao) = zeta
      basis%ao2sh (iao) = ish

   enddo
end subroutine set_p_function

subroutine set_d_function(basis,iat,ish,iao,ibf,ipr,npq,l,nprim,zeta,level,valao)
   type(TBasisset), intent(inout) :: basis
   integer, intent(in)    :: iat
   integer, intent(in)    :: ish
   integer, intent(inout) :: ibf
   integer, intent(inout) :: iao
   integer, intent(inout) :: ipr
   integer, intent(in)    :: npq
   integer, intent(in)    :: l
   integer, intent(in)    :: nprim
   integer, intent(in)    :: valao
   real(wp),intent(in)    :: zeta
   real(wp),intent(in)    :: level
   integer  :: j,p
   real(wp) :: alp(10),cont(10)
   real(wp) :: trafo(5:10) = &
      & [1.0_wp, 1.0_wp, 1.0_wp, sqrt(3.0_wp), sqrt(3.0_wp), sqrt(3.0_wp)]
   integer :: info

   call slaterToGauss(nprim, npq, l, zeta, alp, cont, .true., info)
   basis%minalp(ish) = minval(alp(:nprim))

   do j = 5, 10

      ibf = ibf+1
      basis%primcount(ibf) = ipr
      basis%valao    (ibf) = valao
      basis%aoat     (ibf) = iat
      basis%lao      (ibf) = j
      basis%nprim    (ibf) = nprim
      basis%hdiag    (ibf) = level

      do p=1,nprim
         ipr = ipr+1
         basis%alp (ipr)=alp (p)
         basis%cont(ipr)=cont(p)*trafo(j)
      enddo

      if (j .eq. 5) cycle

      iao = iao+1
      basis%valao2(iao) = valao
      basis%aoat2 (iao) = iat
      basis%lao2  (iao) = j-1
      basis%hdiag2(iao) = level
      basis%aoexp (iao) = zeta
      basis%ao2sh (iao) = ish

   enddo
end subroutine set_d_function

subroutine set_f_function(basis,iat,ish,iao,ibf,ipr,npq,l,nprim,zeta,level,valao)
   type(TBasisset), intent(inout) :: basis
   integer, intent(in)    :: iat
   integer, intent(in)    :: ish
   integer, intent(inout) :: ibf
   integer, intent(inout) :: iao
   integer, intent(inout) :: ipr
   integer, intent(in)    :: npq
   integer, intent(in)    :: l
   integer, intent(in)    :: nprim
   integer, intent(in)    :: valao
   real(wp),intent(in)    :: zeta
   real(wp),intent(in)    :: level
   integer  :: j,p
   real(wp) :: alp(10),cont(10)
   real(wp) :: trafo(11:20) = &
      & [1.0_wp, 1.0_wp, 1.0_wp, sqrt(5.0_wp), sqrt(5.0_wp), &
      &  sqrt(5.0_wp), sqrt(5.0_wp), sqrt(5.0_wp), sqrt(5.0_wp), sqrt(15.0_wp)]
   integer :: info

   call slaterToGauss(nprim, npq, l, zeta, alp, cont, .true., info)
   basis%minalp(ish) = minval(alp(:nprim))

   do j = 11, 20

      ibf = ibf+1
      basis%primcount(ibf) = ipr
      basis%valao    (ibf) = valao
      basis%aoat     (ibf) = iat
      basis%lao      (ibf) = j
      basis%nprim    (ibf) = nprim
      basis%hdiag    (ibf) = level

      do p=1,nprim
         ipr = ipr+1
         basis%alp (ipr)=alp (p)
         basis%cont(ipr)=cont(p)*trafo(j)
      enddo

      if (j.ge.11 .and. j.le.13) cycle

      iao = iao+1
      basis%valao2(iao) = valao
      basis%aoat2 (iao) = iat
      basis%lao2  (iao) = j-3
      basis%hdiag2(iao) = level
      basis%aoexp (iao) = zeta
      basis%ao2sh (iao) = ish

   enddo
end subroutine set_f_function

!========================================================================================!
!========================================================================================!
end module gfn0_basisset
