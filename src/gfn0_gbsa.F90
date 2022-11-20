!====================================================!
! module gfn0_gbsa
! An interface to GBSA and ALPB calculations for GFN0
!
! NOTE: GBSA or ALPB were never parametrized for GFN0
!       Hence, GBSA parameters for GFN2-xTB and
!       ALPB parameters for GFN-FF are used instead.
!====================================================!

module gfn0_gbsa
  use iso_fortran_env,only:wp => real64,stdout => output_unit
#ifdef WITH_GBSA
  use solvation_solv_input, only : TSolvInput
  use solvation_solv_state, only : solutionState
  use solvation_solv_kernel, only : gbKernel
  use solvation_solv_model, only: TSolvModel, init, newBornModel, info
  use solvation_solv_gbsa, only: TBorn
#endif
  implicit none
  private

#ifndef WITH_GBSA
  !> TBorn placeholder if compiled without GBSA support
  type :: TBorn
    integer :: dummy
  end type TBorn
#endif

  public :: TBorn, gfn0_gbsa_init, gfn0_solvation

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!


subroutine gfn0_gbsa_init(nat,at,alpb,solv,gbsa)
!> set up the GBSA (or ALPB) parametrization
   implicit none
   !> INPUT
   integer,intent(in)          :: nat      !> number of atoms
   integer,intent(in)          :: at(nat)  !> atom types
   character(len=*),intent(in) :: solv     !> which solvent to use
   logical,intent(in)          :: alpb     !> use ALPB instead of GBSA?
   !> OUTPUT
   type(TBorn),intent(out)     :: gbsa     !> gbsa type returned 
   !> LOCAL
#ifdef WITH_GBSA
   type(TSolvInput) :: input
   type(TSolvModel) :: model

   if(alpb)then
   input = TSolvInput(solvent=trim(solv), alpb=.true., kernel=gbKernel%p16)
   else
   input = TSolvInput(solvent=trim(solv), alpb=.false., kernel=gbKernel%still)
   endif
   call init(model, input, 0)
   !call info(model,stdout)
   call newBornModel(model, gbsa, at)
#else
   gbsa%dummy = 0
#endif

   return
end subroutine gfn0_gbsa_init

subroutine gfn0_gbsa_print(gbsa,iunit)
   implicit none
   integer :: iunit
   type(TBorn),intent(in)     :: gbsa     !> gbsa type
#ifdef WITH_GBSA
   call gbsa%info(iunit)
#endif
end subroutine gfn0_gbsa_print

!========================================================================================!

subroutine gfn0_solvation(nat,at,xyz,gbsa,esolv,gradient)
!> get energy and gradient contribution from SASA term
   implicit none
   !> INPUT
   integer,intent(in)          :: nat        !> number of atoms
   integer,intent(in)          :: at(nat)    !> atom types 
   real(wp),intent(in)         :: xyz(3,nat) !> coordinats (bohr)
   type(TBorn),intent(inout)   :: gbsa       !> gbsa type returned
   !> OUTPUT
   real(wp),intent(out)        :: esolv      !> solv contribution to E
   real(wp),intent(inout)      :: gradient(3,nat) !> atomic gradient (added to)

#ifdef WITH_GBSA
   call gbsa%update(at, xyz)
   esolv = gbsa%gsasa + gbsa%gshift
   gradient = gradient + gbsa%dsdr
#else
   esolv = 0.0_wp
#endif
end subroutine gfn0_solvation
!========================================================================================!
end module gfn0_gbsa
