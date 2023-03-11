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
module gfn0_qm
  use iso_fortran_env,only:wp => real64,stdout => output_unit

  use gfn0_types
  use gfn0_parameter
  use gfn0_basisset
  use wfn_module

  implicit none
  public :: wfnsetup
  public :: getSelfEnergy
  public :: build_SH0
  public :: fermismear
  interface fermismear
    module procedure :: fermismear_original
    module procedure :: fermismear_nmax
  end interface fermismear
  public :: occ_nmax

  private :: lin
  private :: ncore

  real(wp),private,parameter :: pi = 3.1415926535897932385_wp
  real(wp),private,parameter :: autoev = 27.21138505_wp
  real(wp),private,parameter :: evtoau = 1.0_wp / autoev

contains
!========================================================================================!
  subroutine wfnsetup(xtbData,basis,nat,at,uhf,chrg,wfn)
    type(TxtbData_mod),intent(in)  :: xtbData
    type(TBasisset),intent(in) :: basis
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    integer,intent(in) :: uhf
    integer,intent(in) :: chrg
    type(TWavefunction),intent(inout) :: wfn
    integer,allocatable :: z(:)

    call wfn%allocate(nat,basis%nshell,basis%nao)
    allocate (z(nat),source=0)
    call getZcore(nat,at,z)
    wfn%nel = sum(z) - chrg
    wfn%nel = max(0,wfn%nel)
    deallocate (z)

    wfn%nopen = uhf
    if (wfn%nopen == 0 .and. mod(wfn%nel,2) /= 0) wfn%nopen = 1

    if (wfn%nel .ne. 0) then
      call occu(wfn%nao,wfn%nel,wfn%nopen,wfn%ihomoa,wfn%ihomob,wfn%focca,wfn%foccb)
      wfn%focc = wfn%focca + wfn%foccb
      wfn%ihomo = wfn%ihomoa
!      call wfn%refresh_occu(wfn%nel, uhf)
    else
      wfn%focc = 0.0_wp
      wfn%ihomo = 0
      wfn%ihomoa = 0
      wfn%nopen = 0
    end if

    return
  end subroutine wfnsetup

!========================================================================================!
  subroutine getSelfEnergy(hData,nShell,at,cn,qat,selfEnergy,dSEdcn,dSEdq)
    type(THamiltonianData),intent(in) :: hData
    integer,intent(in) :: nShell(:)
    integer,intent(in) :: at(:)
    real(wp),intent(in),optional :: cn(:)
    real(wp),intent(in),optional :: qat(:)
    real(wp),intent(out) :: selfEnergy(:,:)
    real(wp),intent(out),optional :: dSEdcn(:,:)
    real(wp),intent(out),optional :: dSEdq(:,:)

    integer :: iAt,iZp,iSh

    selfEnergy(:,:) = 0.0_wp
    if (present(dSEdcn)) dSEdcn(:,:) = 0.0_wp
    if (present(dSEdq)) dSEdq(:,:) = 0.0_wp

    !>--- parametrized self-energy
    do iAt = 1,size(cn)
      iZp = at(iAt)
      do iSh = 1,nShell(iZp)
        selfEnergy(iSh,iAt) = hData%selfEnergy(iSh,iZp)
      end do
    end do
    !>--- CN dependent contribution
    if (present(dSEdcn) .and. present(cn)) then
      do iAt = 1,size(cn)
        iZp = at(iAt)
        do iSh = 1,nShell(iZp)
          selfEnergy(iSh,iAt) = selfEnergy(iSh,iAt) &
             & - hData%kCN(iSh,iZp) * cn(iAt)
          dSEdcn(iSh,iAt) = -hData%kCN(iSh,iZp)
        end do
      end do
    end if
    !>--- atomic charge (EEQ) dependent contribution
    if (present(dSEdq) .and. present(qat)) then
      do iAt = 1,size(cn)
        iZp = at(iAt)
        do iSh = 1,nShell(iZp)
          selfEnergy(iSh,iAt) = selfEnergy(iSh,iAt) &
             & - hData%kQShell(iSh,iZp) * qat(iAt) - hData%kQAtom(iZp) * qat(iAt)**2
          dSEdq(iSh,iAt) = -hData%kQShell(iSh,iZp) &
             & - hData%kQAtom(iZp) * 2 * qat(iAt)
        end do
      end do
    end if

  end subroutine getSelfEnergy
!========================================================================================!

  subroutine build_SH0(nShell,hData,selfEnergy,nat,at,basis,nao, &
        & xyz,intcut,sint,h0)
    use overlap_module
    implicit none

    integer,intent(in) :: nShell(:)
    type(THamiltonianData),intent(in) :: hData
    !type(tb_wsc), intent(in) :: wsc
    type(TBasisset),intent(in) :: basis
    real(wp),intent(in) :: selfEnergy(:,:)

    integer,intent(in)  :: nat
    integer,intent(in)  :: nao
    integer,intent(in)  :: at(nat)
    real(wp),intent(in)  :: xyz(3,nat)
    !real(wp),intent(in)  :: lattice(3,3)
    real(wp),intent(in)  :: intcut

    !> overlap matrix S
    real(wp),intent(out) :: sint(nao,nao)
    !> H0 matrix
    real(wp),intent(out) :: h0(nao * (1 + nao) / 2)

    integer,parameter :: llao(0:3) = (/1,3,6,10/)
    integer,parameter :: llao2(0:3) = (/1,3,5,7/)

    integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,ij
    real(wp) tmp1,tmp2,tmp3,tmp4,step,step2
    real(wp) dx,dy,dz,s00r,s00l,s00,alpj
    real(wp) skj,r1,r2,tt,t1,t2,t3,t4,f,ci,cc,cj,alpi,rab2,ab,est

    real(wp) :: den,den2,den4,enpoly
    real(wp) :: zi,zj,zetaij
    real(wp) :: hii,hjj,hav
    real(wp) :: shpoly,km
    logical  :: valaoi,valaoj

    real(wp) ri(3),rj(3),t(3),f1,f2
    real(wp),parameter ::point(3) = 0.0_wp
    real(wp) dtmp(3),qtmp(6),ss(6,6),dd(3,6,6),qq(6,6,6),tmp(6,6)
    real(wp) sstmp(6,6) !TESTST
    integer ip,jp,iat,jat,ati,atj,ish,jsh,icao,jcao,iao,jao,jshmax,il,jl
    integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
    integer itt(0:3)

    integer kat

    parameter(itt=(/0,1,4,10/))
    real(wp) :: saw(10)

    sint = 0.0_wp
    h0 = 0.0_wp

    !$omp parallel do default(none) schedule(dynamic) &
    !$omp private(iat,jat,ij,ati,cc,ci,rab2,atj,ish,ishtyp,valaoi,valaoj, &
    !$omp&        ri,rj,icao,naoi,iptyp,jsh,jshmax,jshtyp,jcao,naoj,jptyp, &
    !$omp&        ss,saw,est,alpi,alpj,ab,iprim,jprim,ip,jp,km,shpoly, &
    !$omp&        mli,mlj,tmp,zi,zj,zetaij,enpoly,iao,jao, &
    !$omp&        ii,jj,k,den,den2,den4,i,j,il,jl,hii,hjj,hav,t) &
    !$omp shared(sint,h0) &
    !$omp shared(basis,at,nShell,hData,xyz,intcut,nat,selfEnergy)
    do iat = 1,nat
      ri = xyz(:,iat)
      ati = at(iat)
      do jat = 1,iat - 1
        atj = at(jat)

        den = hData%electronegativity(ati) - hData%electronegativity(atj)
        den2 = den**2
        den4 = den2**2

        ishells: do ish = 1,nShell(ati)
          ishtyp = hData%angShell(ish,ati)
          icao = basis%caoshell(ish,iat)
          naoi = llao(ishtyp)
          iptyp = itt(ishtyp)
          jshmax = nShell(atj)
          if (iat == jat) jshmax = ish

          jshells: do jsh = 1,jshmax
            jshtyp = hData%angShell(jsh,atj)
            jcao = basis%caoshell(jsh,jat)
            naoj = llao(jshtyp)
            jptyp = itt(jshtyp)

            il = ishtyp + 1
            jl = jshtyp + 1
            !>--- diagonals are the same for all H0 elements
            hii = selfEnergy(ish,iat)
            hjj = selfEnergy(jsh,jat)

            !>--- evaluate the EN polynom for this shells
            enpoly = (1.0_wp + hData%enScale(jshtyp,ishtyp) * den2 &
            &  + hData%enScale4 * hData%enScale(jshtyp,ishtyp) * den4)

            !>--- scale the two shells depending on their exponent
            zi = hData%slaterExponent(ish,ati)
            zj = hData%slaterExponent(jsh,atj)
            zetaij = 2 * sqrt(zi * zj) / (zi + zj)

            !>--- now do the real magic (called EHT enhancement factor)
            km = hData%kScale(jshtyp,ishtyp) * hData%pairParam(ati,atj) &
                &    * zetaij * enpoly

            !>--- check for valence orbitals
            valaoi = hData%valenceShell(ish,ati) .eq. 0
            valaoj = hData%valenceShell(jsh,atj) .eq. 0
            !>--- and scale appropiately
            if (valaoi) then
              if (valaoj) then
                km = 0.0_wp
              else
                km = km * hData%kDiff
              end if
            else
              if (valaoj) km = km * hData%kDiff
            end if

            !>--- averaged H0 element (without overlap contribution!)
            hav = 0.5_wp * km * (hii + hjj)

            !wscAt: do kat = 1, wsc%itbl(jat,iat)
            ss = 0.0_wp
            rj = xyz(:,jat)

            !>--- distance dependent polynomial
            shpoly = shellPoly(hData%shellPoly(il,ati),hData%shellPoly(jl,atj),&
               &      hData%atomicRad(ati),hData%atomicRad(atj),ri,rj)

            !>--- get overlap integral
            call get_overlap(icao,jcao,naoi,naoj,ishtyp,jshtyp,ri,rj,point, &
               &             intcut,basis%nprim,basis%primcount, &
               &             basis%alp,basis%cont,ss)

            !>--- transform from CAO to SAO
            call dtrf2(ss,ishtyp,jshtyp)

            do ii = 1,llao2(ishtyp)
              iao = ii + basis%saoshell(ish,iat)
              do jj = 1,llao2(jshtyp)
                jao = jj + basis%saoshell(jsh,jat)
                if (jao > iao) cycle
                ij = lin(iao,jao)
                sint(jao,iao) = sint(jao,iao) + ss(jj,ii)
                !> Hamiltonian construction
                H0(ij) = H0(ij) + hav * shpoly * ss(jj,ii)
              end do
            end do

            !enddo wscAt
          end do jshells
        end do ishells
      end do
    end do
    !$omp parallel do default(none) shared(nao, sint) private(iao, jao)
    do iao = 1,nao
      do jao = 1,iao - 1
        sint(iao,jao) = sint(jao,iao)
      end do
    end do

    !> diagonal elements
    do iat = 1,nat
      ati = at(iat)
      do ish = 1,nShell(ati)
        iao = 1 + basis%saoshell(ish,iat)
        ishtyp = hData%angShell(ish,ati)
        il = ishtyp + 1
        do iao = 1,llao2(ishtyp)
          i = iao + basis%saoshell(ish,iat)
          sint(i,i) = 1.0_wp + sint(i,i)
          !> H0 is packed, note i*(i-1)/2+i = i*(1+i)/2
          ii = i * (1 + i) / 2

          !> calculate environment dependent shift
          hii = selfEnergy(ish,iat)
          H0(ii) = hii
        end do
      end do
    end do

  end subroutine build_SH0
!> DERIVATIVES of SH0
  subroutine build_dSH0(nShell,hData,selfEnergy,dSEdcn,dSEdq,nat,basis, &
        & thr,nao,at,xyz,P,Pew,g,dHdcn,dHdq)
    use overlap_module
    implicit none

    !> INPUT
    integer,intent(in) :: nShell(:)
    type(THamiltonianData),intent(in) :: hData
    real(wp),intent(in) :: selfEnergy(:,:)
    real(wp),intent(in) :: dSEdcn(:,:)
    real(wp),intent(in) :: dSEdq(:,:)
    integer,intent(in)      :: nat
    !type(tb_wsc), intent(in) :: wsc
    type(TBasisset),intent(in) :: basis
    real(wp),intent(in)      :: thr
    integer,intent(in)      :: nao
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    !real(wp),intent(in) :: lattice(3,3)
    real(wp),intent(in) :: P(nao,nao)
    real(wp),intent(in) :: Pew(nao,nao)
    !> OUTPUT
    real(wp),intent(inout) :: dHdq(nat)
    real(wp),intent(inout) :: dHdcn(nat)
    real(wp),intent(inout) :: g(3,nat)
    !real(wp),intent(inout) :: sigma(3,3)

    integer,parameter :: llao(0:3) = (/1,3,6,10/)
    integer,parameter :: llao2(0:3) = (/1,3,5,7/)

    integer itt(0:3)
    parameter(itt=(/0,1,4,10/))
    real(wp) tmp1,tmp2,tmp3,step,step2,step3,s00r,s00l,s00,alpj
    real(wp) skj,r1,r2,tt,t1,t2,t3,t4
    real(wp) thr2,f,ci,cc,cj,alpi,rab2,ab,est
    real(wp) f1,f2,tmp(6,6),rij(3),ri(3),rj(3),rij2,t(3)
    real(wp) stmp,ral(3,3),rar(3,3),rbl(3,3),pre
    real(wp) dtmp,qtmp,rbr(3,3),r2l(3),r2r(3),qqa(6,6,6,3)
    real(wp) ss(6,6,3),ddc(3,6,6,3),qqc(6,6,6,3),dda(3,6,6,3)
    integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,ij,lin,jshmax
    integer ip,jp,iat,jat,kat,ati,atj,ish,jsh,icao,jcao,iao,jao,ixyz,jxyz
    integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
    real(wp) :: dumdum(3),dum,sdq(6,6),sdqg(3,6,6)
    real(wp),parameter :: point(3) = 0.0_wp

    integer  :: il,jl
    real(wp) :: Hij,Pii,Pij,HPij,H0sr
    real(wp) :: den,den2,den4,enpoly
    real(wp) :: zi,zj,zetaij
    real(wp) :: hii,hjj,hav
    real(wp) :: shpoly,dshpoly(3),km
    logical  :: valaoi,valaoj
    real(wp) :: g_xyz(3)

    !$omp parallel default(none) &
    !$omp private(iat,jat,ixyz,ati,cc,ci,rab2,atj,ish,ishtyp,valaoi,valaoj,g_xyz, &
    !$omp&        ri,rj,icao,naoi,iptyp,jsh,jshmax,jshtyp,jcao,naoj,jptyp,rij2, &
    !$omp&        sdq,sdqg,est,alpi,alpj,ab,iprim,jprim,ip,jp,km,shpoly,dshpoly, &
    !$omp&        mli,mlj,dum,dumdum,tmp,stmp,rij,zi,zj,zetaij,enpoly,iao,jao,t, &
    !$omp&        ii,jj,k,den,den2,den4,i,j,il,jl,hii,hjj,hav,Pij,Hij,HPij,H0sr) &
    !$omp reduction (+:g,dhdcn,dhdq) &
    !$omp shared(basis,at,nShell,hData,xyz,thr,nat,selfEnergy, &
    !$omp        dSEdcn,dSEdq,P,Pew)
    !$omp do schedule(runtime)
    do iat = 1,nat
      ri = xyz(:,iat)
      ati = at(iat)
      do jat = 1,iat - 1
        atj = at(jat)

        den = hData%electronegativity(ati) - hData%electronegativity(atj)
        den2 = den**2
        den4 = den2**2

        ishells: do ish = 1,nShell(ati)
          ishtyp = hData%angShell(ish,ati)
          icao = basis%caoshell(ish,iat)
          naoi = llao(ishtyp)
          iptyp = itt(ishtyp)
          jshmax = nShell(atj)
          if (iat == jat) jshmax = ish

          jshells: do jsh = 1,jshmax
            jshtyp = hData%angShell(jsh,atj)
            jcao = basis%caoshell(jsh,jat)
            naoj = llao(jshtyp)
            jptyp = itt(jshtyp)

            il = ishtyp + 1
            jl = jshtyp + 1
            !>--- diagonals are the same for all H0 elements
            hii = selfEnergy(ish,iat)
            hjj = selfEnergy(jsh,jat)

            !>--- evaluate the EN polynom for this shells
            enpoly = (1.0_wp + hData%enScale(jl - 1,il - 1) * den2 &
               & + hData%enScale4 * hData%enScale(jl - 1,il - 1) * den4)

            !>--- we scale the two shells depending on their exponent
            zi = hData%slaterExponent(ish,ati)
            zj = hData%slaterExponent(jsh,atj)
            zetaij = 2 * sqrt(zi * zj) / (zi + zj)

            !>--- now do the real magic (called EHT enhancement factor)
            km = hData%kScale(jl - 1,il - 1) * hData%pairParam(ati,atj) * zetaij * enpoly

            !>--- check for valence orbitals
            valaoi = hData%valenceShell(ish,ati) .eq. 0
            valaoj = hData%valenceShell(jsh,atj) .eq. 0
            !>--- and scale appropiately
            if (valaoi) then
              if (valaoj) then
                km = 0.0_wp
              else
                km = km * hData%kDiff
              end if
            else
              if (valaoj) km = km * hData%kDiff
            end if

            !>--- averaged H0 element (without overlap contribution!)
            hav = 0.5_wp * km * (hii + hjj)

            !  wscAt: do kat = 1,wsc%itbl(jat,iat)
            !     t = wsc%lattr(:,kat,jat,iat)
            !     rj = xyz(:,jat) + matmul(lattice,t)
            rj = xyz(:,jat)
            rij = ri - rj
            rij2 = sum(rij**2)

            ! distance dependent polynomial
            call dshellPoly(hData%shellPoly(il,ati),hData%shellPoly(jl,atj),&
               & hData%atomicRad(ati),hData%atomicRad(atj),rij2,ri,rj,&
               & shpoly,dshpoly)

            call get_grad_overlap(icao,jcao,naoi,naoj,ishtyp,jshtyp,ri,rj, &
               &                  point,thr,basis%nprim,basis%primcount, &
               &                  basis%alp,basis%cont,sdq,sdqg)

            call dtrf2(sdq,ishtyp,jshtyp)
            do ixyz = 1,3
              !>--- transform from CAO to SAO, only transform the gradient
              tmp = sdqg(ixyz,:,:)
              call dtrf2(tmp,ishtyp,jshtyp)
              sdqg(ixyz,:,:) = tmp
            end do

            do ii = 1,llao2(ishtyp)
              iao = ii + basis%saoshell(ish,iat)
              do jj = 1,llao2(jshtyp)
                jao = jj + basis%saoshell(jsh,jat)

                Pij = P(jao,iao) * evtoau

                !>--- Hamiltonian element without overlap
                Hij = hav * shpoly
                HPij = Hij * Pij

                !>--- overlap contribution
                stmp = 2 * (HPij - Pew(jao,iao))

                !>--- distance dependent polynomial contribution
                H0sr = 2 * HPij * sdq(jj,ii)
                g_xyz = stmp * sdqg(:,jj,ii) + H0sr * dshpoly / shpoly

                !>--- distribute to atom i
                g(:,iat) = g(:,iat) + g_xyz

                !>--- distribute to atom j
                g(:,jat) = g(:,jat) - g_xyz

                        !!>--- reduce to strain
                !sigma(:,1) = sigma(:,1) + g_xyz*rij(1)
                !sigma(:,2) = sigma(:,2) + g_xyz*rij(2)
                !sigma(:,3) = sigma(:,3) + g_xyz*rij(3)

                !>--- Hamiltonian without Hav
                HPij = km * shpoly * Pij * sdq(jj,ii)

                !>--- save dE/dCN for CNi
                dhdcn(iat) = dhdcn(iat) + HPij * dSEdcn(ish,iat)
                !>--- save dE/dCN for CNj
                dhdcn(jat) = dhdcn(jat) + HPij * dSEdcn(jsh,jat)

                !>--- save dE/dq for qi
                dhdq(iat) = dhdq(iat) + HPij * dSEdq(ish,iat)
                !>--- save dE/dq for qj
                dhdq(jat) = dhdq(jat) + HPij * dSEdq(jsh,jat)

              end do
            end do

            !enddo wscAt

          end do jshells
        end do ishells

      end do !> jat
    end do  !> iat
    !$omp end do
    !$omp end parallel
    !>-- diagonal contributions
    do iat = 1,nat
      ati = at(iat)
      do ish = 1,nShell(ati)
        iao = 1 + basis%saoshell(ish,iat)
        ishtyp = hData%angShell(ish,ati)
        il = ishtyp + 1
        do iao = 1,llao2(ishtyp)
          i = iao + basis%saoshell(ish,iat)
          ii = i * (1 + i) / 2
          !>--- H0 is packed, note i*(i-1)/2+i = i*(1+i)/2

          Pii = P(i,i) * evtoau

          !>--- save dE/dCN for CNi
          dhdcn(iat) = dhdcn(iat) + Pii * dSEdcn(ish,iat)
          !>--- save dE/dq for qi
          dhdq(iat) = dhdq(iat) + Pii * dSEdq(ish,iat)
        end do
      end do
    end do

  end subroutine build_dSH0

!=========================================================================================!

!! ========================================================================
!  eigenvalue solver
!! ========================================================================
  subroutine solve(ndim,H,S,X,P,e,fail)
    use gfn0_math_wrapper,only:lapack_sygvd
    !> INPUT
    integer,intent(in)   :: ndim
    real(wp),intent(in)  :: H(ndim * (ndim + 1) / 2)
    real(wp),intent(in)   :: S(ndim,ndim)
    !> OUTPUT
    real(wp),intent(out)  :: X(ndim,ndim)
    real(wp),intent(out)  :: P(ndim,ndim)
    real(wp),intent(out)  :: e(ndim)
    logical,intent(out)  :: fail
    !> LOCAL
    integer :: i,j,k,info

    !>--- SETUP
    fail = .false.

    do i = 1,ndim
      do j = 1,i
        k = j + i * (i - 1) / 2
        X(j,i) = H(k)
        X(i,j) = X(j,i)
      end do
    end do
    P = S
    !>--- DIAG IN NON-ORTHORGONAL BASIS
    call lapack_sygvd(1,'v','u',ndim,X,ndim,P,ndim,e,info)
    if (info .ne. 0) then
      fail = .true.
      return
    end if

    return
  end subroutine solve

!=========================================================================================!

  subroutine setzshell(xtbData,n,at,nshell,z,zsh,e)
    type(TxTBData_mod),intent(in) :: xtbData
    integer,intent(in)  :: n
    integer,intent(in)  :: at(n)
    integer,intent(in)  :: nshell
    real(wp),intent(in)  :: z(n)
    real(wp),intent(out) :: zsh(nshell)
    real(wp),intent(out) :: e

    real(wp) ntot,fracz
    integer i,j,k,m,l,ll(0:3),iat,lll,iver
    data ll/1,3,5,7/

    k = 0
    e = 0.0_wp
    do i = 1,n
      iat = at(i)
      ntot = -1.d-6
      do m = 1,xtbData%nShell(iat)
        l = xtbData%hamiltonian%angShell(m,iat)
        k = k + 1
        zsh(k) = xtbData%hamiltonian%referenceOcc(m,iat)
        ntot = ntot + zsh(k)
        if (ntot .gt. z(i)) zsh(k) = 0
        e = e + xtbData%hamiltonian%selfEnergy(m,iat) * zsh(k)
      end do
    end do

  end subroutine setzshell

!! ========================================================================
!  S(R) enhancement factor
!! ========================================================================
  pure function shellPoly(iPoly,jPoly,iRad,jRad,xyz1,xyz2)
    real(wp),intent(in) :: iPoly,jPoly
    real(wp),intent(in) :: iRad,jRad
    real(wp),intent(in) :: xyz1(3),xyz2(3)
    real(wp) :: shellPoly
    real(wp) :: rab,k1,rr,r,rf1,rf2,dx,dy,dz,a

    a = 0.5_wp       ! R^a dependence 0.5 in GFN1

    dx = xyz1(1) - xyz2(1)
    dy = xyz1(2) - xyz2(2)
    dz = xyz1(3) - xyz2(3)

    rab = (dx * dx + dy * dy + dz * dz)
    rab = sqrt(rab)

    ! this sloppy conv. factor has been used in development, keep it
    rr = jRad + iRad

    r = rab / rr

    rf1 = 1.0_wp + (0.01_wp * iPoly) * (r**a)
    rf2 = 1.0_wp + (0.01_wp * jPoly) * (r**a)

    shellPoly = rf1 * rf2

  end function shellPoly

!! ========================================================================
!  derivative of S(R) enhancement factor
!! ========================================================================
  pure subroutine dshellPoly(iPoly,jPoly,iRad,jRad,rab2,xyz1,xyz2,rf,dxyz)
    real(wp),intent(in)  :: iPoly,jPoly
    real(wp),intent(in)  :: iRad,jRad
    real(wp),intent(in)  :: rab2
    real(wp),intent(out) :: dxyz(3),rf
    real(wp),intent(in)  :: xyz1(3),xyz2(3)
    real(wp) :: rab,k1,k2,rr,r,a,dum,rf1,rf2,dx,dy,dz
    real(wp) :: t14,t15,t17,t22,t20,t23,t10,t11,t35,t13

    a = 0.5_wp   ! R^a dependence 0.5 in GFN1

    dx = xyz1(1) - xyz2(1)
    dy = xyz1(2) - xyz2(2)
    dz = xyz1(3) - xyz2(3)

    rab = sqrt(rab2)

    ! this sloppy conv. factor has been used in development, keep it
    r = iRad + jRad

    rr = rab / r

    k1 = iPoly * 0.01_wp
    k2 = jPoly * 0.01_wp

    t14 = rr**a
    t15 = k1 * t14
    t17 = 1.0_wp / rab2
    t22 = rr**a
    t23 = k2 * t22
    rf = (1.0_wp + t15) * (1.0_wp + k2 * t22)
    dxyz(:) = (t15 * (1.0_wp + t23) + (1.0_wp + t15) * k2 * t22) * a * t17*[dx,dy,dz]

  end subroutine dshellPoly

!=========================================================================================!

  subroutine occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
    integer  :: nel
    integer  :: nopen
    integer  :: ndim
    integer  :: ihomoa
    integer  :: ihomob
    real(wp) :: focca(ndim)
    real(wp) :: foccb(ndim)
    integer  :: focc(ndim)
    integer  :: i,na,nb,ihomo

    focc = 0
    focca = 0
    foccb = 0
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
        focca(i) = 1.0d0
        foccb(i) = 1.0d0
      end if
      if (focc(i) .eq. 1) focca(i) = 1.0d0
    end do

    ihomoa = 0
    ihomob = 0
    do i = 1,ndim
      if (focca(i) .gt. 0.99) ihomoa = i
      if (foccb(i) .gt. 0.99) ihomob = i
    end do

  end subroutine occu

!=========================================================================================!
  subroutine occ_nmax(nao,occ,nmaxa,nmaxb,ihomoa,ihomob)
    implicit none
    !> INPUT
    integer,intent(in) :: nao
    real(wp),intent(in) :: occ(nao)
    !> OUTPUT
    integer,allocatable,intent(out) :: nmaxa(:)
    integer,allocatable,intent(out) :: nmaxb(:)
    integer,intent(out) :: ihomoa,ihomob
    !> LOCAL
    integer,allocatable :: iocc(:),iocca(:),ioccb(:)
    integer :: i,na,nb

    if (.not. allocated(nmaxa)) allocate (nmaxa(nao))
    if (.not. allocated(nmaxb)) allocate (nmaxb(nao))

    allocate (iocc(nao),source=0)
    allocate (iocca(nao),ioccb(nao),source=0)
    iocc(:) = nint(occ(:))

    ihomoa = 0
    ihomob = 0
    na = 0
    nb = 0
    nmaxa(:) = 1
    nmaxb(:) = 1
    do i = 1,nao
      if (iocc(i) == 2) then
        iocca(i) = 1
        na = na + 1
        ioccb(i) = 1
        nb = nb + 1
      elseif (iocc(i) == 1) then
        iocca(i) = 1
        na = na + 1
        ioccb(i) = 0
      else
        iocca(i) = 0
        ioccb(i) = 0
      end if
    end do
    do i = 1,nao
      if (iocca(i) == 1) ihomoa = i
      if (ioccb(i) == 1) ihomob = i
    end do
    do i = 1,nao
      if (i < ihomoa) then
        if (iocca(i) == 0) nmaxa(i) = 0
      end if
      if (i < ihomob) then
        if (ioccb(i) == 0) nmaxb(i) = 0
      end if
    end do

    deallocate (ioccb,iocca)
    deallocate (iocc)
  end subroutine occ_nmax

!=========================================================================================!

  subroutine fermismear_original(prt,norbs,nel,t,eig,occ,fod,e_fermi,s)
    integer,intent(in)  :: norbs
    integer,intent(in)  :: nel
    real(wp),intent(in)  :: eig(norbs)
    real(wp),intent(out) :: occ(norbs)
    real(wp),intent(in)  :: t
    real(wp),intent(out) :: fod
    real(wp),intent(out) :: e_fermi
    logical,intent(in)  :: prt

    real(wp) :: bkt,occt,total_number
    real(wp) :: total_dfermi,dfermifunct,fermifunct,s,change_fermi

    real(wp),parameter :: autoev = 27.21138505_wp
    real(wp),parameter :: kB = 3.166808578545117e-06_wp
    real(wp),parameter :: boltz = kB * autoev
    real(wp),parameter :: thr = 1.d-9
    integer :: ncycle,i,j,m,k,i1,i2

    bkt = boltz * t

    e_fermi = 0.5 * (eig(nel) + eig(nel + 1))
    occt = nel

    do ncycle = 1,200
      total_number = 0.0
      total_dfermi = 0.0
      do i = 1,norbs
        fermifunct = 0.0
        if ((eig(i) - e_fermi) / bkt .lt. 50) then
          fermifunct = 1.0 / (exp((eig(i) - e_fermi) / bkt) + 1.0)
          dfermifunct = exp((eig(i) - e_fermi) / bkt) / &
          &       (bkt * (exp((eig(i) - e_fermi) / bkt) + 1.0)**2)
        else
          dfermifunct = 0.0
        end if
        occ(i) = fermifunct
        total_number = total_number + fermifunct
        total_dfermi = total_dfermi + dfermifunct
      end do
      change_fermi = (occt - total_number) / total_dfermi
      e_fermi = e_fermi + change_fermi
      if (abs(occt - total_number) .le. thr) exit
    end do

    fod = 0
    s = 0
    do i = 1,norbs
      if (occ(i) .gt. thr .and. 1.0d00 - occ(i) .gt. thr) &
      &   s = s + occ(i) * log(occ(i)) + (1.0d0 - occ(i)) * log(1.0d00 - occ(i))
      if (eig(i) .lt. e_fermi) then
        fod = fod + 1.0d0 - occ(i)
      else
        fod = fod + occ(i)
      end if
    end do
    s = s * kB * t

    if (prt) then
      write (*,'('' t,e(fermi),nfod : '',2f10.3,f10.6)') t,e_fermi,fod
    end if

  end subroutine fermismear_original

!=========================================================================================!

  subroutine fermismear_nmax(prt,norbs,ihomo,nmax,t,eig,occ,ef,TS)
    !> INPUT
    logical,intent(in)  :: prt    !> print statement
    integer,intent(in)  :: norbs  !> number of orbitals
    integer,intent(in)  :: ihomo  !> position of the H(S)OMO
    integer,intent(in)  :: nmax(norbs) !> maximum allowed occupation
    real(wp),intent(in)  :: t     !> (electronic) temperature
    real(wp),intent(in)  :: eig(norbs) !> orbital energies
    !> OUTPUT
    real(wp),intent(out) :: occ(norbs)  !> orbital occupations
    real(wp),intent(out),optional :: ef !> fermi temperature
    real(wp),intent(out),optional :: TS !> electronic free energy
    !> LOCAL
    real(wp) :: e_fermi,s
    real(wp) :: bkt,occt,total_number,fod
    real(wp) :: total_dfermi,dfermifunct,fermifunct,change_fermi

    real(wp),parameter :: autoev = 27.21138505_wp
    real(wp),parameter :: kB = 3.166808578545117e-06_wp
    real(wp),parameter :: boltz = kB * autoev
    real(wp),parameter :: thr = 1.d-9
    integer :: ncycle,i,j,m,k,i1,i2,nholes
    real(wp) :: deigkt
    !integer :: nfermi,ncycle2
    !real(wp),allocatable :: efermis(:),occsum(:)

    bkt = boltz * t

    e_fermi = 0.5 * (eig(ihomo) + eig(ihomo + 1))
    nholes = count((nmax(:) < 1),1)
    occt = ihomo - nholes

    do ncycle = 1,200
        total_number = 0.0
        total_dfermi = 0.0
        do i = 1,norbs
          fermifunct = 0.0
          dfermifunct = 0.0

          deigkt = (eig(i) - e_fermi) / bkt

          if ((deigkt .lt. 50)) then
            if ((nmax(i) > 0)) then !> smear electrons
              fermifunct = 1.0 / (exp(deigkt) + 1.0)
              dfermifunct = exp(deigkt) / &
              &       (bkt * (exp(deigkt) + 1.0)**2)
            else !> smear holes
              fermifunct = 1.0 - (1.0 / (exp(deigkt) + 1.0))
              dfermifunct = -(exp(deigkt) / &
              &       (bkt * (exp(deigkt) + 1.0)**2))
            end if
          end if
          occ(i) = fermifunct
          total_number = total_number + fermifunct
          total_dfermi = total_dfermi + dfermifunct
        end do
        change_fermi = (occt - total_number) / total_dfermi
        e_fermi = e_fermi + change_fermi
        if (abs(occt - total_number) .le. thr) exit
      end do

    fod = 0
    s = 0
    do i = 1,norbs
      if (occ(i) .gt. thr .and. 1.0d00 - occ(i) .gt. thr) &
      &   s = s + occ(i) * log(occ(i)) + (1.0d0 - occ(i)) * log(1.0d00 - occ(i))
      if (eig(i) .lt. e_fermi) then
        fod = fod + 1.0d0 - occ(i)
      else
        fod = fod + occ(i)
      end if
    end do
    s = s * kB * t

    if (prt) then
      write (*,'('' t,e(fermi),nfod : '',2f10.3,f10.6)') t,e_fermi,fod
    end if

    if (present(ef)) ef = e_fermi
    if (present(TS)) TS = s

  end subroutine fermismear_nmax
!=========================================================================================!
  subroutine dmat(ndim,focc,C,P)
    !> density matrix
    !> C: MO coefficient
    !> P  dmat
    use gfn0_math_wrapper,only:gemm
    integer,intent(in)  :: ndim
    real(wp),intent(in)  :: focc(:)
    real(wp),intent(in)  :: C(:,:)
    real(wp),intent(out) :: P(:,:)
    integer :: i,m
    real(wp),allocatable :: Ptmp(:,:)

    allocate (Ptmp(ndim,ndim))
    Ptmp = 0.0_wp

    do m = 1,ndim
      do i = 1,ndim
        Ptmp(i,m) = C(i,m) * focc(m)
      end do
    end do
    call gemm(C,Ptmp,P,transb='t')

    deallocate (Ptmp)
  end subroutine dmat

!=========================================================================================!

  pure elemental integer function lin(i1,i2)
    integer,intent(in) :: i1,i2
    integer :: idum1,idum2
    idum1 = max(i1,i2)
    idum2 = min(i1,i2)
    lin = idum2 + idum1 * (idum1 - 1) / 2
    return
  end function lin
!=========================================================================================!
  pure elemental integer function ncore(at)
    integer,intent(in) :: at
    if (at .le. 2) then
      ncore = 0
    elseif (at .le. 10) then
      ncore = 2
    elseif (at .le. 18) then
      ncore = 10
    elseif (at .le. 29) then   !zn
      ncore = 18
    elseif (at .le. 36) then
      ncore = 28
    elseif (at .le. 47) then
      ncore = 36
    elseif (at .le. 54) then
      ncore = 46
    elseif (at .le. 71) then
      ncore = 54
    elseif (at .le. 79) then
      ncore = 68
    elseif (at .le. 86) then
      ncore = 78
    end if
  end function ncore
!=========================================================================================!
  pure subroutine getZcore(nat,at,z)
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    integer,intent(out) :: z(nat)
    integer :: i
    z = 0
    do i = 1,nat
      z(i) = at(i) - ncore(at(i))
    end do
  end subroutine getZcore

!=========================================================================================!
end module gfn0_qm
