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

!====================================================!
! module gfn0_interface
! An interface to GFN0 calculations
!====================================================!
module gfn0_interface
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use gfn0_module
  implicit none
  private

  !> module types
  public :: gfn0_data
  type :: gfn0_data
    type(TBasisset),allocatable     :: basis   !< basis set info
    type(TxTBData_mod),allocatable  :: xtbData !< main data/parameter frame
    type(Twavefunction),allocatable :: wfn     !< wavefunction data
    type(TBorn),allocatable         :: gbsa    !< gbsa/alpb model
    type(gfn0_partials),allocatable :: part    !< partial derivatives
    logical :: update = .true.
    real(wp),allocatable :: refxyz(:,:)
    real(wp) :: eclassic
    real(wp),allocatable :: gclassic(:,:)
  end type gfn0_data

  public :: gfn0_results
  type :: gfn0_results
    real(wp) :: etot = 0.0_wp
    real(wp) :: ies = 0.0_wp
    real(wp) :: edisp = 0.0_wp
    real(wp) :: erep = 0.0_wp
    real(wp) :: esrb = 0.0_wp
    real(wp) :: esolv = 0.0_wp
    real(wp) :: eel = 0.0_wp
    real(wp) :: gnorm = 0.0_wp
  contains
    procedure :: print => print_gfn0_res
  end type gfn0_results

  !> module routines
  public :: gfn0_init
  public :: gfn0_singlepoint
  interface gfn0_singlepoint
    module procedure :: gfn0_singlepoint_full
    module procedure :: gfn0_singlepoint_wrap
    module procedure :: gfn0_singlepoint_short
  end interface gfn0_singlepoint
  public :: gfn0_occ_singlepoint
  interface gfn0_occ_singlepoint
    module procedure :: gfn0_occ_singlepoint_wrap
  end interface gfn0_occ_singlepoint
 
  public :: gfn0_print_summary

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine print_gfn0_res(self,iunit)
    class(gfn0_results) :: self
    integer,intent(in) :: iunit
    write (iunit,*)
    write (iunit,*) "GFN0 calculation summary"
    write (iunit,'(3x,a5,1x,f16.8)') 'IES',self%ies
    write (iunit,'(3x,a5,1x,f16.8)') 'Edisp',self%edisp
    write (iunit,'(3x,a5,1x,f16.8)') 'Erep',self%erep
    write (iunit,'(3x,a5,1x,f16.8)') 'Esrb',self%esrb
    write (iunit,'(3x,a5,1x,f16.8)') 'Esolv',self%esolv
    write (iunit,'(3x,a5,1x,f16.8)') 'Eel',self%eel
    write (iunit,*) '--------------------------'
    write (iunit,'(3x,a5,1x,f16.8)') 'Etot',self%etot
    write (iunit,'(3x,a5,1x,f16.8)') 'gnorm',self%gnorm
    write (iunit,*)
  end subroutine print_gfn0_res

!========================================================================================!

  subroutine gfn0_print_summary(iunit,gdat,res)
    implicit none
    integer, intent(in) :: iunit
    type(gfn0_data),intent(in)    :: gdat
    type(gfn0_results),intent(in),optional :: res
    
    call pr_wfn_param(iunit,gdat%xtbData,gdat%basis,gdat%wfn) 
    if(allocated(gdat%gbsa))then
    call gdat%gbsa%info(iunit)
    endif
    if(allocated(gdat%wfn))then
     call pr_orbital_eigenvalues(iunit,gdat%wfn,9)
    endif
    if(present(res))then
    call res%print(iunit)
    endif

  end subroutine gfn0_print_summary

!========================================================================================!

  subroutine gfn0_init(nat,at,xyz,chrg,uhf,gdat,solv,alpb)
    implicit none
    !> INPUT
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in)  :: uhf
    integer,intent(in)  :: chrg
    character(len=*),intent(in),optional :: solv
    logical,intent(in),optional :: alpb
    !> OUTPUT
    type(gfn0_data),intent(inout) :: gdat

    if (allocated(gdat%basis)) deallocate (gdat%basis)
    if (allocated(gdat%xtbData)) deallocate (gdat%xtbData)
    if (allocated(gdat%wfn)) deallocate (gdat%wfn)
    if (allocated(gdat%gbsa)) deallocate (gdat%gbsa)
    if (allocated(gdat%part)) deallocate (gdat%part)

    allocate (gdat%basis)
    allocate (gdat%xtbData)
    allocate (gdat%wfn)
    call initGFN0Params(nat,at,xyz,chrg,uhf,gdat%basis,gdat%xtbData)
    call wfnsetup(gdat%xtbData,gdat%basis,nat,at,uhf,chrg,gdat%wfn)

    if (present(solv)) then
      allocate(gdat%gbsa)
      if (present(alpb)) then
        call gfn0_gbsa_init(nat,at,alpb,solv,gdat%gbsa)
      else
        call gfn0_gbsa_init(nat,at,.false.,solv,gdat%gbsa)
      end if
    end if

    allocate(gdat%part)
 
    if(allocated(gdat%refxyz))deallocate(gdat%refxyz)
    allocate(gdat%refxyz(3,nat), source=0.0_wp)
    !gdat%refxyz = xyz   
    gdat%update = .true.
    if(allocated(gdat%gclassic))deallocate(gdat%gclassic)
    allocate(gdat%gclassic(3,nat), source=0.0_wp)

    return
  end subroutine gfn0_init
!========================================================================================!

  subroutine gfn0_singlepoint_classical(nat,at,xyz,chrg,uhf,basis,xtbData,wfn,part,gbsa, &
  &          energy,gradient,fail,res)
    implicit none
    !> INPUT
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in)  :: uhf
    integer,intent(in)  :: chrg
    type(TBasisset),intent(inout)     :: basis   !< basis set info
    type(TxTBData_mod),intent(inout)  :: xtbData !< main data/parameter frame
    type(Twavefunction),intent(inout) :: wfn     !< wavefunction data
    type(gfn0_partials),intent(inout) :: part
    type(TBorn),allocatable,intent(inout) :: gbsa
    type(gfn0_results),intent(inout),optional :: res
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,nat)
    logical,intent(out)  :: fail
    !> LOCAL
    real(wp),allocatable :: cn(:)
    real(wp),allocatable :: dcndr(:,:,:)

    real(wp),allocatable :: qat(:)
    real(wp),allocatable :: dqdr(:,:,:)

    real(wp) :: ies,edisp,erep,esrb,eel,esolv,gnorm
    energy = 0.0_wp
    gradient = 0.0_wp
    fail = .false.
    ies = 0.0_wp
    edisp = 0.0_wp
    erep = 0.0_wp
    esrb = 0.0_wp
    esolv = 0.0_wp
    eel = 0.0_wp
    gnorm = 0.0_wp
!>--- CN and CN gradient
    if(.not.allocated(part%cn)) allocate (part%cn(nat),source=0.0_wp)
    if(.not.allocated(part%dcndr)) allocate (part%dcndr(3,nat,nat),source=0.0_wp)
    call gfn0_getCN(nat,at,xyz,part%cn,part%dcndr)

    if(.not.allocated(part%qat)) allocate (part%qat(nat),source=0.0_wp)
    if(.not.allocated(part%dqdr)) allocate (part%dqdr(3,nat,nat),source=0.0_wp)
!>--- Solvation+EEQ (if required)
    if (allocated(gbsa)) then
      call gfn0_solvation(nat,at,xyz,gbsa,esolv,gradient)
      call gfn0_electrostatics(nat,at,xyz,chrg,xtbData, &
      & part%cn,part%dcndr,ies,gradient,part%qat,part%dqdr,gbsa)
    else
!>--- EEQ charges, gradients and IES energy (gasphase)
      call gfn0_electrostatics(nat,at,xyz,chrg,xtbData, &
      & part%cn,part%dcndr,ies,gradient,part%qat,part%dqdr)
    end if

!>--- D4 two-body dispersion energy
    call gfn0_dispersion(nat,at,xyz,chrg,xtbData, &
    &  part%cn,part%dcndr,edisp,gradient)

!>--- Repulsion energy
    call gfn0_repulsion(nat,at,xyz,xtbData%repulsion,erep,gradient)

!>--- SRB energy
    call gfn0_shortranged(nat,at,xyz,xtbData%srb,part%cn,part%dcndr, &
    &                     esrb,gradient)

!>--- Add up total energy
    energy = ies + edisp + erep + esrb + esolv
    gnorm = sqrt(sum(gradient**2))

    if (present(res)) then
      res%etot = energy
      res%ies = ies
      res%edisp = edisp
      res%erep = erep
      res%esrb = esrb
      res%esolv = esolv
      res%eel = eel
      res%gnorm = gnorm
    end if

    !deallocate (dqdr,qat,dcndr,cn)
  end subroutine gfn0_singlepoint_classical

!========================================================================================!

  subroutine gfn0_singlepoint_full(nat,at,xyz,chrg,uhf,basis,xtbData,wfn,part,gbsa, &
  &          energy,gradient,fail,res)
    implicit none
    !> INPUT
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in)  :: uhf
    integer,intent(in)  :: chrg
    type(TBasisset),intent(inout)     :: basis   !< basis set info
    type(TxTBData_mod),intent(inout)  :: xtbData !< main data/parameter frame
    type(Twavefunction),intent(inout) :: wfn     !< wavefunction data
    type(gfn0_partials),intent(inout) :: part
    type(TBorn),allocatable,intent(inout) :: gbsa
    type(gfn0_results),intent(inout),optional :: res
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,nat)
    logical,intent(out)  :: fail
    !> LOCAL
    real(wp) :: ies,edisp,erep,esrb,eel,esolv,gnorm

    energy = 0.0_wp
    gradient = 0.0_wp
    fail = .false.
    eel = 0.0_wp
    gnorm = 0.0_wp

!>--- Classical part
    call gfn0_singlepoint_classical(nat,at,xyz,chrg,uhf,basis,xtbData,wfn,part,gbsa, &
  &          energy,gradient,fail,res=res)

!>--- QM part
    !call wfnsetup(xtbData,basis,nat,at,uhf,chrg,wfn)
    call gfn0_eht(nat,at,xyz,xtbData,basis,part%cn,part%dcndr,part%qat,part%dqdr, &
   &                wfn,eel,gradient,fail)

!>--- Add up total energy
    energy = energy + eel
    gnorm = sqrt(sum(gradient**2))

    if (present(res)) then
      res%etot = energy
      res%eel = eel
      res%gnorm = gnorm
    end if

    !deallocate (dqdr,qat,dcndr,cn)
  end subroutine gfn0_singlepoint_full

!========================================================================================!

  subroutine gfn0_singlepoint_wrap(nat,at,xyz,chrg,uhf,gdat, &
&          energy,gradient,fail,res)
    implicit none
    !> INPUT
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in)  :: uhf
    integer,intent(in)  :: chrg
    type(gfn0_data),intent(inout) :: gdat
    type(gfn0_results),intent(inout),optional :: res
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,nat)
    logical,intent(out)  :: fail
    !> LOCAL
    real(wp) :: eel,gnorm
    
!>--- Check if geometry has changed
    call check_geo(nat,xyz,gdat%refxyz,gdat%update)

    energy = 0.0_wp
    gradient = 0.0_wp
    fail = .false.
    eel = 0.0_wp
    gnorm = 0.0_wp

!>--- Classical part
    if(gdat%update)then
    call gfn0_singlepoint_classical(nat,at,xyz,chrg,uhf, &
  &  gdat%basis, gdat%xtbData, gdat%wfn, gdat%part, gdat%gbsa, &
  &          gdat%eclassic, gdat%gclassic, fail, res=res)
    endif
    energy = gdat%eclassic
    gradient = gdat%gclassic

!>--- QM part
!    call gfn0_eht_occ(nat,at,xyz,occ, gdat%xtbData, gdat%basis, &
!    &   gdat%part%cn, gdat%part%dcndr, gdat%part%qat, gdat%part%dqdr, &
!    &   gdat%wfn, eel,gradient,fail)

    call gfn0_eht_var(nat,at,xyz, gdat%xtbData, gdat%basis, gdat%wfn, gdat%part, &
     &                    eel,gradient,fail, update=gdat%update)


!>--- Add  EHT energy to total
      energy = energy + eel
      gnorm = sqrt(sum(gradient**2))

    if (present(res)) then
      res%etot = energy
      res%eel = eel
      res%gnorm = gnorm
    end if

    gdat%update = .false.
  end subroutine gfn0_singlepoint_wrap

!========================================================================================!
!> shortcut version of the singlepoint call that always sets up the entire parametrization
  subroutine gfn0_singlepoint_short(nat,at,xyz,chrg,uhf,energy,gradient,fail,solv,alpb)
    implicit none
    !> INPUT
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in)  :: uhf
    integer,intent(in)  :: chrg
    character(len=*),intent(in),optional :: solv
    logical,intent(in),optional :: alpb
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,nat)
    logical,intent(out)  :: fail
    !> LOCAL
    type(gfn0_data) :: gdat

    if (present(solv)) then
      if (present(alpb)) then
        call gfn0_init(nat,at,xyz,chrg,uhf,gdat,solv,alpb)
      else
        call gfn0_init(nat,at,xyz,chrg,uhf,gdat,solv)
      end if
    else
     call gfn0_init(nat,at,xyz,chrg,uhf,gdat) 
    end if
    call gfn0_singlepoint_wrap(nat,at,xyz,chrg,uhf,gdat,energy,gradient,fail)

  end subroutine gfn0_singlepoint_short

!========================================================================================!

  subroutine gfn0_occ_singlepoint_wrap(nat,at,xyz,chrg,uhf,occ,gdat, &
&          energy,gradient,fail,res)
    implicit none
    !> INPUT
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in)  :: uhf
    integer,intent(in)  :: chrg
    type(gfn0_data),intent(inout) :: gdat
    real(wp),intent(in) :: occ(gdat%basis%nao)
    type(gfn0_results),intent(inout),optional :: res
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,nat)
    logical,intent(out)  :: fail
    !> LOCAL
    real(wp) :: ies,edisp,erep,esrb,eel,esolv,gnorm
    integer :: i

!>--- Check if geometry has changed
    call check_geo(nat,xyz,gdat%refxyz,gdat%update)

    energy = 0.0_wp
    gradient = 0.0_wp
    fail = .false.
    eel = 0.0_wp
    gnorm = 0.0_wp

!>--- Classical part
    if(gdat%update)then
    call gfn0_singlepoint_classical(nat,at,xyz,chrg,uhf, &
  &  gdat%basis, gdat%xtbData, gdat%wfn, gdat%part, gdat%gbsa, &
  &          gdat%eclassic, gdat%gclassic, fail, res=res)
    endif
    energy = gdat%eclassic
    gradient = gdat%gclassic

!>--- QM part
!    call gfn0_eht_occ(nat,at,xyz,occ, gdat%xtbData, gdat%basis, &
!    &   gdat%part%cn, gdat%part%dcndr, gdat%part%qat, gdat%part%dqdr, &
!    &   gdat%wfn, eel,gradient,fail)

    call gfn0_eht_var(nat,at,xyz, gdat%xtbData, gdat%basis, gdat%wfn, gdat%part, &
     &                    eel,gradient,fail, update=gdat%update, occ=occ)
 

!>--- Add  EHT energy to total
      energy = energy + eel
      gnorm = sqrt(sum(gradient**2))

    if (present(res)) then
      res%etot = energy
      res%eel = eel
      res%gnorm = gnorm
    end if


    gdat%update = .false. 
  end subroutine gfn0_occ_singlepoint_wrap

!========================================================================================!

  subroutine check_geo(nat,xyz,refxyz,update)
    implicit none
    integer,intent(in) :: nat
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(inout) :: refxyz(3,nat)
    logical,intent(out)  :: update
    real(wp),allocatable :: diff(:,:)
    real(wp),parameter :: thr = 1.0d-8
    integer :: i,j
    update = .false. 
    !allocate(diff(3,nat),source=0.0_wp)
    !diff(:,:) = abs( xyz(:,:) - refxyz(:,:)) 
    !if(any(diff.gt.thr).or.nat==1) then
    !  update = .true.
    !  refxyz = xyz
    !endif
    !deallocate(diff)
    if(nat==1)then
       update = .true. 
       return 
    endif 
    ILOOP : do i=1,nat
      JLOOP : do j=1,3
         if(abs( xyz(j,i) - refxyz(j,i)).gt.thr)then
            update = .true.
            refxyz = xyz
            exit ILOOP
         endif
      enddo JLOOP
    enddo ILOOP
  end subroutine check_geo

!========================================================================================!
end module gfn0_interface
