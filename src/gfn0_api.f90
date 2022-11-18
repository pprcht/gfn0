!====================================================!
! module tblite_api
! An interface to GFN0 calculations
!====================================================!

module gfn0_api
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
  public :: gfn0_setup
  public :: gfn0_singlepoint
  interface gfn0_singlepoint
    module procedure :: gfn0_singlepoint_full
    module procedure :: gfn0_singlepoint_wrap
    module procedure :: gfn0_singlepoint_short
  end interface gfn0_singlepoint
  public :: gfn0_occ_singlepoint
  interface gfn0_occ_singlepoint
    module procedure :: gfn0_occ_singlepoint_full
    module procedure :: gfn0_occ_singlepoint_wrap
  end interface gfn0_occ_singlepoint

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine print_gfn0_res(self,iunit)
    class(gfn0_results) :: self
    integer,intent(in) :: iunit
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
  end subroutine print_gfn0_res

!========================================================================================!

  subroutine gfn0_setup(nat,at,xyz,chrg,uhf,gdat,solv,alpb)
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

    return
  end subroutine gfn0_setup
!========================================================================================!

  subroutine gfn0_singlepoint_full(nat,at,xyz,chrg,uhf,basis,xtbData,wfn,gbsa, &
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
    allocate (cn(nat),dcndr(3,nat,nat),source=0.0_wp)
    call gfn0_getCN(nat,at,xyz,cn,dcndr)

    allocate (qat(nat),dqdr(3,nat,nat))
!>--- Solvation+EEQ (if required)
    if (allocated(gbsa)) then
      call gfn0_solvation(nat,at,xyz,gbsa,esolv,gradient)
      call gfn0_electrostatics(nat,at,xyz,chrg,xtbData, &
      & cn,dcndr,ies,gradient,qat,dqdr,gbsa)
    else
!>--- EEQ charges, gradients and IES energy (gasphase)
      call gfn0_electrostatics(nat,at,xyz,chrg,xtbData, &
      & cn,dcndr,ies,gradient,qat,dqdr)
    end if

!>--- D4 two-body dispersion energy
    call gfn0_dispersion(nat,at,xyz,chrg,xtbData, &
    &  cn,dcndr,edisp,gradient)

!>--- Repulsion energy
    call gfn0_repulsion(nat,at,xyz,xtbData%repulsion,erep,gradient)

!>--- SRB energy
    call gfn0_shortranged(nat,at,xyz,xtbData%srb,cn,dcndr, &
    &                     esrb,gradient)

!>--- QM part
    !call wfnsetup(xtbData,basis,nat,at,uhf,chrg,wfn)
    call gfn0_eht(nat,at,xyz,xtbData,basis,cn,dcndr,qat,dqdr, &
   &                wfn,eel,gradient,fail)

!>--- Add up total energy
    energy = ies + edisp + erep + esrb + esolv + eel
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

    deallocate (dqdr,qat,dcndr,cn)
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

    if (present(res)) then
      call gfn0_singlepoint_full(nat,at,xyz,chrg,uhf,&
  &        gdat%basis,gdat%xtbData,gdat%wfn,gdat%gbsa, &
  &        energy,gradient,fail,res)
    else
      call gfn0_singlepoint_full(nat,at,xyz,chrg,uhf, &
  &          gdat%basis,gdat%xtbData,gdat%wfn,gdat%gbsa, &
  &          energy,gradient,fail)
    end if

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
        call gfn0_setup(nat,at,xyz,chrg,uhf,gdat,solv,alpb)
      else
        call gfn0_setup(nat,at,xyz,chrg,uhf,gdat,solv)
      end if
    else
     call gfn0_setup(nat,at,xyz,chrg,uhf,gdat) 
    end if
    call gfn0_singlepoint_wrap(nat,at,xyz,chrg,uhf,gdat,energy,gradient,fail)

  end subroutine gfn0_singlepoint_short

!========================================================================================!

  subroutine gfn0_occ_singlepoint_full(nat,at,xyz,chrg,uhf,nlev,occ, &
  &          basis,xtbData,wfn,gbsa, &
  &          energies,gradients,fail,res)
    implicit none
    !> INPUT
    type(TBasisset),intent(inout)     :: basis   !< basis set info
    type(TxTBData_mod),intent(inout)  :: xtbData !< main data/parameter frame
    type(Twavefunction),intent(inout) :: wfn     !< wavefunction data
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in)  :: nlev
    real(wp),intent(in) :: occ(basis%nao,nlev)
    integer,intent(in)  :: uhf
    integer,intent(in)  :: chrg
    type(gfn0_results),intent(inout),optional :: res
    type(TBorn),allocatable,intent(inout) :: gbsa
    !> OUTPUT
    real(wp),intent(out) :: energies(nlev)
    real(wp),intent(out) :: gradients(3,nat,nlev)
    logical,intent(out)  :: fail
    !> LOCAL
    real(wp),allocatable :: cn(:)
    real(wp),allocatable :: dcndr(:,:,:)

    real(wp),allocatable :: qat(:)
    real(wp),allocatable :: dqdr(:,:,:)

    real(wp),allocatable :: eels(:)
    real(wp),allocatable :: gradient(:,:)

    real(wp) :: ies,edisp,erep,esrb,eel,esolv,energy,gnorm
    integer :: i

    energies = 0.0_wp
    gradients = 0.0_wp
    fail = .false.
    ies = 0.0_wp
    edisp = 0.0_wp
    erep = 0.0_wp
    esrb = 0.0_wp
    esolv = 0.0_wp
    eel = 0.0_wp
    gnorm = 0.0_wp

    allocate (gradient(3,nat),source=0.0_wp)

!>--- CN and CN gradient
    allocate (cn(nat),dcndr(3,nat,nat),source=0.0_wp)
    call gfn0_getCN(nat,at,xyz,cn,dcndr)

    allocate (qat(nat),dqdr(3,nat,nat))
!>--- Solvation+EEQ (if required)
    if (allocated(gbsa)) then
      call gfn0_solvation(nat,at,xyz,gbsa,esolv,gradient)
      call gfn0_electrostatics(nat,at,xyz,chrg,xtbData, &
      & cn,dcndr,ies,gradient,qat,dqdr,gbsa)
    else
!>--- EEQ charges, gradients and IES energy (gasphase)
      call gfn0_electrostatics(nat,at,xyz,chrg,xtbData, &
      & cn,dcndr,ies,gradient,qat,dqdr)
    end if

!>--- D4 two-body dispersion energy
    call gfn0_dispersion(nat,at,xyz,chrg,xtbData, &
    &  cn,dcndr,edisp,gradient)

!>--- Repulsion energy
    call gfn0_repulsion(nat,at,xyz,xtbData%repulsion,erep,gradient)

!>--- SRB energy
    call gfn0_shortranged(nat,at,xyz,xtbData%srb,cn,dcndr, &
    &                     esrb,gradient)

!>--- QM part
    !call wfnsetup(xtbData,basis,nat,at,uhf,chrg,wfn)
    allocate (eels(nlev),source=0.0_wp)
    !    call gfn0_eht(nat,at,xyz,xtbData,basis,cn,dcndr,qat,dqdr, &
    !   &                wfn,eel,gradient,fail)
    call gfn0_eht_occ(nat,at,xyz,nlev,occ,xtbData,basis,cn,dcndr,qat,dqdr, &
   &                wfn,eels,gradients,fail)

!>--- Add up total energies
    do i = 1,nlev
      energies(i) = ies + edisp + erep + esrb + esolv + eels(i)
      gradients(:,:,i) = gradients(:,:,i) + gradient(:,:)
      energy = energy + energies(i) / float(nlev)
      gnorm = gnorm + sqrt(sum(gradient**2)) / float(nlev)
    end do

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

    deallocate (eels,dqdr,qat,dcndr,cn,gradient)
  end subroutine gfn0_occ_singlepoint_full

!========================================================================================!

  subroutine gfn0_occ_singlepoint_wrap(nat,at,xyz,chrg,uhf,nlev,occ,gdat, &
&          energies,gradients,fail,res)
    implicit none
    !> INPUT
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    integer,intent(in)  :: uhf
    integer,intent(in)  :: chrg
    type(gfn0_data),intent(inout) :: gdat
    integer,intent(in)  :: nlev
    real(wp),intent(in) :: occ(gdat%basis%nao,nlev)
    type(gfn0_results),intent(inout),optional :: res
    !> OUTPUT
    real(wp),intent(out) :: energies(nlev)
    real(wp),intent(out) :: gradients(3,nat,nlev)
    logical,intent(out)  :: fail

    if (present(res)) then
      call gfn0_occ_singlepoint_full(nat,at,xyz,chrg,uhf,nlev,occ, &
  &   gdat%basis,gdat%xtbData,gdat%wfn,gdat%gbsa,energies,gradients,fail,res)
    else
      call gfn0_occ_singlepoint_full(nat,at,xyz,chrg,uhf,nlev,occ, &
  &   gdat%basis,gdat%xtbData,gdat%wfn,gdat%gbsa,energies,gradients,fail)
    end if

  end subroutine gfn0_occ_singlepoint_wrap

!========================================================================================!

end module gfn0_api
