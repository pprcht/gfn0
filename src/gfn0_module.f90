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
module gfn0_module
  use iso_fortran_env,only:wp => real64,stdout => output_unit

  use gfn0_types
  use gfn0_parameter
  use gfn0_basisset
  use wfn_module
  use gfn0_qm
  use gfn0_prints
  use gfn0_gbsa
  use gfn0_occ
  implicit none
  private

  !> interface to other modules
  public :: TxTBData_mod
  public :: TBasisset
  public :: Twavefunction,wfnsetup
  public :: gfn0_partials
  public :: pr_wfn_param
  public :: pr_orbital_eigenvalues
  public :: gfn0_eht_occ
  public :: generate_config
  public :: gfn0_gbsa_init
  public :: TBorn
  public :: gfn0_solvation

  !> from this module
  public :: initGFN0Params
  public :: gfn0_getCN
  public :: gfn0_electrostatics
  public :: gfn0_dispersion
  public :: gfn0_repulsion
  public :: gfn0_shortranged
  public :: gfn0_eht
  public :: gfn0_eht_var

  !> privat params
  real(wp),private,parameter :: autoev = 27.21138505_wp
  real(wp),private,parameter :: evtoau = 1.0_wp / autoev

contains

!========================================================================================!

  subroutine initGFN0Params(nat,at,xyz,chrg,uhf,basis,xtbData)
!> subroutine initGFN0Params
!> Initialize parametrization and basisset for calculation
!> Note, this information is independent of
!> the molecular geometry, charge and uhf configuration
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in) :: chrg
    integer,intent(in) :: uhf

    type(TBasisset),intent(out)   :: basis   !< basis set info
    type(TxTBData_mod),intent(out) :: xtbData !< main data/parameter frame

    logical :: ok

    !> Get the parametrization
    call initGFN0(xtbData)

    !> Get the basisset
    call newBasisset(xtbData,nat,at,basis,ok)

    return
  end subroutine initGFN0Params

!========================================================================================!

  subroutine gfn0_getCN(nat,at,xyz,cn,dcndr)
!> subroutine gfn0_getCN
!> Calculate coordination number and its Cartesian derivative
    use cn_module
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    !> Coordination number
    real(wp),intent(out) :: cn(:)
    !> Derivative of the coordination number w.r.t. Cartesian coordinates
    real(wp),intent(out) :: dcndr(:,:,:)

    call getCoordinationNumber(nat,at,xyz,40.0_wp,cnType%erf,cn,dcndr)
    call cutCoordinationNumber(nat,cn,dcndr,maxCN=8.0_wp)
    return
  end subroutine gfn0_getCN

!========================================================================================!

  subroutine gfn0_electrostatics(nat,at,xyz,chrg,xtbData, &
     & cn,dcndr,energy,gradient,qat,dqdr,gbsa)
!> subroutine gfn0_electrostatics
!> Calculate atomic charges using the charge equilibirum model (EEQ).
!> Also returns the Cartesian derivative dqdr, as well as the
!> respective isotropic electrostatic energies (IES) and its gradient.
    use eeq_module
    !> INPUTS
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in) :: chrg
    type(TxTBData_mod),intent(in) :: xtbData
    !> Coordination number
    real(wp),intent(in) :: cn(:)
    !> Derivative of the coordination number w.r.t. Cartesian coordinates
    real(wp),intent(in) :: dcndr(:,:,:)
    !> gbsa object
    type(TBorn),intent(inout),optional :: gbsa
    !> OUTPUTS
    !> Electrostatic energy
    real(wp),intent(inout),optional :: energy
    !> Molecular gradient
    real(wp),intent(inout),optional :: gradient(:,:)
    !> Atomic partial charges
    real(wp),intent(out),optional :: qat(:)
    !> Derivative of the partial charges w.r.t. Cartesian coordinates
    real(wp),intent(out),optional :: dqdr(:,:,:)

    if (.not. present(gbsa)) then
      call chargeEquilibration(nat,at,xyz,chrg,cn,dcndr, &
         & xtbData%coulomb%electronegativity,  & !chi
         & xtbData%coulomb%kcn,                & !kcn
         & xtbData%coulomb%chemicalHardness,   & !gam
         & xtbData%coulomb%chargeWidth,        & !rad
         & energy,gradient,qat,dqdr)
    else
      call chargeEquilibration(nat,at,xyz,chrg,cn,dcndr, &
         & xtbData%coulomb%electronegativity,  & !chi
         & xtbData%coulomb%kcn,                & !kcn
         & xtbData%coulomb%chemicalHardness,   & !gam
         & xtbData%coulomb%chargeWidth,        & !rad
         & energy,gradient,qat,dqdr,gbsa=gbsa)
    end if
    return
  end subroutine gfn0_electrostatics
!========================================================================================!

  subroutine gfn0_dispersion(nat,at,xyz,chrg,xtbData, &
     &  cn,dcndr,energy,gradient)
!> subroutine gfn0_dispersion
!> Calculate two-body D4 dispersion interaction energy and the
!> respective Cartesian gradient.
!> The D4 calculation in GFN0 uses a SEPARATE set of
!> EEQ charges that are set up with the erf-type CN used
!> in the rest of GFN0 terms, AND a set of CNs only
!> for the D4 calculation.
    use eeq_module
    use dftd4param,only:D4gam,D4chi,D4kcn,D4alp
    use gfn0_dftd4,only:d4_gradient
    use cn_module
    !> INPUTS
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in) :: chrg
    type(TxTBData_mod),intent(in) :: xtbData !< main data/parameter frame
    !> Coordination number
    real(wp),intent(in) :: cn(:)
    !> Derivative of the coordination number w.r.t. Cartesian coordinates
    real(wp),intent(in) :: dcndr(:,:,:)

    !> OUTPUTS
    !> Electrostatic energy
    real(wp),intent(inout) :: energy
    !> Molecular gradient
    real(wp),intent(inout) :: gradient(:,:)

    !> Atomic partial charges
    real(wp),allocatable :: qat(:)
    !> Derivative of the partial charges w.r.t. Cartesian coordinates
    real(wp),allocatable :: dqdr(:,:,:)
    !> Coordination number (only for D4 calc.)
    real(wp),allocatable :: ccn(:)
    !> Derivative of the coordination number w.r.t. Cartesian coordinates (only for D4 calc.)
    real(wp),allocatable :: dccndr(:,:,:)

    !> dummy argument
    real(wp) :: edum
    real(wp),allocatable :: gdum(:,:)
    integer :: i

    allocate (qat(nat),dqdr(3,nat,nat),source=0.0_wp)
    allocate (ccn(nat),dccndr(3,nat,nat),source=0.0_wp)
    allocate (gdum(3,nat),source=0.0_wp)

    !> get EEQ charges for D4 parametrization, but GFN0 erf-type CN!
    call chargeEquilibration(nat,at,xyz,chrg,cn,dcndr, &
       & D4chi,D4kcn,D4gam,D4alp,&
       & edum,gdum,qat,dqdr)

    !> get D4 covalent coordination number (ccn & dccndr)
    call getCoordinationNumber(nat,at,xyz,40.0_wp,cnType%cov,ccn,dccndr)

    !> calculate two-body D4 contribution
    call d4_gradient(nat,at,xyz,xtbData%dispersion,60.0_wp, &
       &  ccn,dccndr,qat,dqdr,energy,gradient)

    deallocate (gdum)
    deallocate (dccndr,ccn,dqdr,qat)

  end subroutine gfn0_dispersion

!========================================================================================!

  subroutine gfn0_repulsion(nat,at,xyz,repData,erep,gradient)
!> subroutine gfn0_repulsion
!> calculate the repulsion energy and respective gradients
!> This is a classical term only dependent on the
!> molecular geometry and fully parametrized
    implicit none

    !> number of atoms of the molecule
    integer,intent(in) :: nat
    !> atom types by atomic number
    integer,intent(in) :: at(nat)
    !> Cartesian coordinates
    real(wp),intent(in) :: xyz(3,nat)
    !> Repulsion potential parameters
    type(TRepulsionData),intent(in) :: repData

    !> Cartesian gradient
    real(wp),intent(inout) :: gradient(:,:)
    real(wp),intent(out) :: erep

    integer :: i,j,k,lin
    integer :: iat,jat,ati,atj
    real(wp) :: ri(3),dr(3)
    real(wp) :: rij(3)
    real(wp) :: r2,r,r2top34,r2top34top2
    real(wp) :: den2,den4
    real(wp) :: alpha
    real(wp) :: repab
    real(wp) :: expterm
    real(wp) :: dtmp
    real(wp),parameter :: rthr = 1600.0_wp
    real(wp) :: w,t(3)
    integer  :: latrep(3),tx,ty,tz,itr

    w = 1.0_wp
    !> initialize
    erep = 0.0_wp
    do i = 1,nat
      ati = at(i)
      do j = 1,i
        rij = xyz(:,i) - xyz(:,j)
        r2 = sum(rij**2)
        if (r2 .gt. rthr .or. r2 .lt. 1.0e-6_wp) cycle
        atj = at(j)
        r = sqrt(r2)
        den2 = (repData%electronegativity(ati) - repData%electronegativity(atj))**2
        den4 = den2**2
        alpha = sqrt(repData%alpha(ati) * repData%alpha(atj)) &
                * (1.0_wp + (0.01_wp * den2 + 0.01_wp * den4) * repData%enScale)
        repab = repData%zeff(ati) * repData%zeff(atj)
        r2top34 = r2**0.75_wp
        r2top34top2 = r2top34**2
        expterm = exp(-alpha * r2top34) * repab
        !> save repulsion energy
        erep = erep + expterm / r * w
        !> save repulsion gradient
        dtmp = expterm * (1.5_wp * alpha * r2top34 + 1) / r2top34top2 * w
        gradient(:,i) = gradient(:,i) - dtmp * rij
        gradient(:,j) = gradient(:,j) + dtmp * rij
      end do !> j atom
    end do !> i atom
    return
  end subroutine gfn0_repulsion

!========================================================================================!

  subroutine gfn0_shortranged(nat,at,xyz,srb, &
     &       cn,dcndr,esrb,gradient)
!> subroutine gfn0_shortranranged
!> short-ranged bond (SRB) correction
!> The SRB term is fully classical and depends on the
!> atomic CNs.
    use gfn0_srb
    use gfn0_math_wrapper,only:contract
    implicit none

    !> number of atoms of the molecule
    integer,intent(in) :: nat
    !> atom types by atomic number
    integer,intent(in) :: at(nat)
    !> Cartesian coordinates
    real(wp),intent(in) :: xyz(3,nat)
    !> Coordination number and gradient
    real(wp),intent(in) :: cn(:)
    real(wp),intent(in) :: dcndr(:,:,:)
    !> SRB param object
    type(TShortRangeData),intent(in) :: srb
    !> SRB energy and gradient
    real(wp),intent(out) :: esrb
    real(wp),intent(inout) :: gradient(:,:)

    integer :: i,j,k,wsAt
    integer :: iat,jat,ati,atj
    integer :: nsrb,itr
    real(wp) :: den
    real(wp) :: expterm
    ! distances
    real(wp) :: rij(3)
    real(wp) :: r2,r,dr,rab
    real(wp) :: dtmp,pre
    real(wp) :: w
    ! allocatables
    real(wp),allocatable :: dEdr0(:)
    real(wp),allocatable :: rab0(:)
    real(wp),allocatable :: drab0dr(:,:,:)
    integer,allocatable :: srblist(:,:)

    !> initialize
    esrb = 0.0_wp

    w = 1.0_wp

    call build_srblist(nat,at,xyz,nsrb,srblist)
    if (nsrb .eq. 0) return

    !> get memory
    allocate (dEdr0(nsrb),source=0.0_wp)
    allocate (drab0dr(3,nat,nsrb),source=0.0_wp)
    allocate (rab0(nsrb),source=0.0_wp)
    !> get approximated distances rab and gradients
    call approx_rab(nat,at,xyz,cn,dcndr,nsrb,srblist, &
    &               srb%shift,rab0,drab0dr)
    do i = 1,nsrb
      iat = srblist(1,i)
      jat = srblist(2,i)
      ati = at(iat)
      atj = at(jat)
      den = paulingEN(ati) - paulingEN(atj)
      pre = srb%steepness * (1.0_wp + srb%enScale * den**2)
      rij = xyz(:,iat) - xyz(:,jat)
      rab = norm2(rij)
      dr = rab - rab0(i)
      expterm = srb%prefactor * exp(-pre * dr**2)
      !> save SRB energy
      esrb = esrb + expterm * w
      dtmp = 2.0_wp * pre * dr * expterm * w
      gradient(:,iat) = gradient(:,iat) - dtmp * rij / rab
      gradient(:,jat) = gradient(:,jat) + dtmp * rij / rab
      !> three body gradient
      dEdr0(i) = dEdr0(i) + dtmp
    end do !> i

    !> construct gradient
    call contract(drab0dr,dEdr0,gradient,beta=1.0_wp)

    return
  end subroutine gfn0_shortranged

!========================================================================================!

  subroutine gfn0_eht(nat,at,xyz,xtbData,basis,cn,dcndr,qat,dqdr, &
     &                wfn,eel,gradient,fail)
!> subroutine gfn0_eht
!> Set up the H0 hamiltonian and overlap,
!> and from that get the "QM" part of the energy
    use wfn_module
    use gfn0_math_wrapper,only:contract
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(in) :: cn(:)
    real(wp),intent(in) :: dcndr(:,:,:)
    real(wp),intent(in) :: qat(:)
    real(wp),intent(in) :: dqdr(:,:,:)
    type(TBasisset),intent(in)   :: basis
    type(TxTBData_mod),intent(in) :: xtbData
    type(Twavefunction),intent(inout) :: wfn
    !> OUTPUT
    real(wp),intent(out) :: eel
    real(wp),intent(inout) :: gradient(:,:)
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
    if (allocated(wfn%S)) wfn%S = S

    !>--- Diagonalization
    call solve(nao,H0,S,wfn%C,wfn%P,wfn%emo,fail)
    if (fail) return

    !>--- Fermi Smearing and electronic energy
    if (et .gt. 0.1_wp) then
      if (wfn%ihomoa + 1 .le. nao) then
        call fermismear(.false.,nao,wfn%ihomoa,et,wfn%emo,wfn%focca,nfoda,efa,ga)
      end if
      if (wfn%ihomob + 1 .le. nao) then
        call fermismear(.false.,nao,wfn%ihomob,et,wfn%emo,wfn%foccb,nfodb,efb,gb)
      end if
      wfn%focc = wfn%focca + wfn%foccb
    end if
    call dmat(nao,wfn%focc,wfn%C,wfn%P)
    eel = sum(wfn%focc * wfn%emo) * evtoau + ga + gb

    !>--- Gradient Setup
    allocate (dHdcn(nat),dHdq(nat),pew(nao,nao),tmp(nao), &
       &      source=0.0_wp)
    tmp = wfn%focc * wfn%emo * evtoau
    !> setup energy weighted density matrix = Pew for gradient calculation
    call dmat(nao,tmp,wfn%C,Pew)
    call build_dSH0(xtbData%nShell,xtbData%hamiltonian,selfEnergy, &
        & dSEdcn,dSEdq,nat,basis,intcut,nao,at,xyz, &
        & wfn%P,Pew,gradient,dHdcn,dHdq)

    !>--- Gradient Contraction
    call contract(dcndr,dhdcn,gradient,beta=1.0_wp)
    call contract(dqdr,dhdq,gradient,beta=1.0_wp)

    !>--- deallocation
    deallocate (tmp,Pew,dHdq,dHdcn)
    deallocate (dSEdq,dSEdcn,selfEnergy)
    deallocate (S,H0)

    return
  end subroutine gfn0_eht
!========================================================================================!
  subroutine gfn0_eht_var(nat,at,xyz,xtbData,basis,wfn,part, &
     &                    eel,gradient,fail,update,occ)
!> subroutine gfn0_eht_var
!> Set up the H0 hamiltonian and overlap,
!> and from that get the "QM" part of the energy
    use wfn_module
    use gfn0_math_wrapper,only:contract
    !> INPUT
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(TBasisset),intent(in)   :: basis
    type(TxTBData_mod),intent(in) :: xtbData
    type(Twavefunction),intent(inout) :: wfn
    type(gfn0_partials),intent(inout) :: part
    logical,intent(in),optional :: update
    real(wp),intent(in),optional :: occ(basis%nao)
    !> OUTPUT
    real(wp),intent(out) :: eel
    real(wp),intent(inout) :: gradient(:,:)
    logical,intent(out)  :: fail

    !> LOCAL
    integer :: i,j,k
    real(wp) :: intcut,scfconv
    real(wp) :: et,efa,efb,nfoda,nfodb,ga,gb
    integer :: maxshell,nao,naop
    real(wp),allocatable :: S(:,:)
    real(wp),allocatable :: H0(:)
    real(wp),allocatable :: dHdcn(:)
    real(wp),allocatable :: dHdq(:)
    real(wp),allocatable :: Pew(:,:)
    real(wp),allocatable :: tmp(:)
    integer :: ihomoa,ihomob
    integer,allocatable  :: nmaxa(:),nmaxb(:)
    logical :: refresh

    eel = 0.0_wp
    fail = .false.
    if(present(update))then
      refresh = update
    else
      refresh = .true.
    endif 

    !>---  Numerical stuff and cutoffs
    intcut = 25.0_wp - 10.0_wp * log10(xtbData%acc)
    intcut = max(20.0_wp,intcut)
    scfconv = 1.e-6_wp * xtbData%acc

    !>--- allocation
    nao = basis%nao
    naop = basis%nao * (basis%nao + 1) / 2
    if (refresh) then
      allocate (H0(naop),source=0.0_wp)
      allocate (S(nao,nao),source=0.0_wp)
      maxshell = maxval(xtbData%nshell)

      if (.not. allocated(part%selfEnergy)) &
      &  allocate (part%selfEnergy(maxshell,nat),source=0.0_wp)
      if (.not. allocated(part%dSEdcn)) &
      &  allocate (part%dSEdcn(maxshell,nat),source=0.0_wp)
      if (.not. allocated(part%dSEdq)) &
      &  allocate (part%dSEdq(maxshell,nat),source=0.0_wp)

      wfn%q = part%qat

      !>--- Setup, AO integrals of S and H0
      call getSelfEnergy(xtbData%hamiltonian,xtbData%nShell,at,part%cn,wfn%q, &
         & part%selfEnergy,part%dSEdcn,part%dSEdq)
      call build_SH0(xtbData%nShell,xtbData%hamiltonian,part%selfEnergy, &
            & nat,at,basis,nao,xyz,intcut,S,H0)
      if (allocated(wfn%S)) wfn%S = S

      !>--- Diagonalization
      call solve(nao,H0,S,wfn%C,wfn%P,wfn%emo,fail)
      deallocate (S,H0)
      if (fail) return
    end if

    !>--- Electronic energy
    et = xtbData%etemp
    eel = 0.0_wp
    allocate (dHdcn(nat),dHdq(nat),Pew(nao,nao),tmp(nao),source=0.0_wp)
    if (present(occ)) then
      if(et .gt. 0.1_wp)then
        call occ_nmax(nao,occ,nmaxa,nmaxb,ihomoa,ihomob)
        if(ihomoa + 1 .le. nao )then
          call fermismear(.false.,nao,ihomoa,nmaxa,et,wfn%emo,wfn%focca,TS=ga)
        endif
        if(ihomob + 1 .le. nao )then
          call fermismear(.false.,nao,ihomob,nmaxb,et,wfn%emo,wfn%foccb,TS=gb)
        endif
        wfn%focc = wfn%focca + wfn%foccb
        !eel = ga + gb
      else  
        wfn%focc = occ
      endif
      !write(*,*) sum(occ),ihomoa,ihomob
    else
      if (et .gt. 0.1_wp) then
        if (wfn%ihomoa + 1 .le. nao) then
          call fermismear(.false.,nao,wfn%ihomoa,et,wfn%emo,wfn%focca,nfoda,efa,ga)
        end if
        if (wfn%ihomob + 1 .le. nao) then
          call fermismear(.false.,nao,wfn%ihomob,et,wfn%emo,wfn%foccb,nfodb,efb,gb)
        end if
        wfn%focc = wfn%focca + wfn%foccb
      end if
      eel = ga + gb
    end if
    !write(*,*) wfn%focc
    tmp = wfn%focc
    call dmat(nao,tmp,wfn%C,wfn%P)
    eel = eel + sum(tmp * wfn%emo) * evtoau
    tmp(:) = wfn%focc(:) * wfn%emo * evtoau


    !>--- Gradient Setup
    !> setup energy weighted density matrix = Pew for gradient calculation
    call dmat(nao,tmp,wfn%C,Pew)
    call build_dSH0(xtbData%nShell,xtbData%hamiltonian, &
        & part%selfEnergy,part%dSEdcn,part%dSEdq, &
        & nat,basis,intcut,nao,at,xyz, &
        & wfn%P,Pew,gradient,dHdcn,dHdq)

    !>--- Gradient Contraction
    call contract(part%dcndr,dhdcn,gradient,beta=1.0_wp)
    call contract(part%dqdr,dhdq,gradient,beta=1.0_wp)

    !>--- deallocation
    if(allocated(nmaxa)) deallocate(nmaxa)
    if(allocated(nmaxb)) deallocate(nmaxb)
    deallocate (tmp,Pew,dHdq,dHdcn)
    !deallocate (dSEdq,dSEdcn,selfEnergy)

    return
  end subroutine gfn0_eht_var
!========================================================================================!
end module gfn0_module
