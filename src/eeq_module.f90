


!> The original (unmodified) source code can be found under the GNU LGPL 3.0 license
!> Copyright (C) 2019-2020 Sebastian Ehlert
!> at https://github.com/grimme-lab/xtb

!> Implementation of the electronegativity equilibration model
module eeq_module
   !use xtb_mctc_blas, only : mctc_dot, mctc_symv
   !use xtb_mctc_lapack, only : lapack_sytrf, lapack_sytri
   !use xtb_mctc_la, only : contract
   use iso_fortran_env, only: wp=>real64
   use math_wrapper
   implicit none
   private

   public :: chargeEquilibration
   !public :: TENEquilibration, init

   real(wp),private,parameter :: pi=3.1415926535897932385_wp
   real(wp),private,parameter :: sqrtpi  = sqrt(pi)


   !!> Electronegativity equilibration model
   !type :: TENEquilibration

   !   !> Electronegativity of each species
   !   real(wp), allocatable :: chi(:, :)

   !   !> Coordination number dependence of each species
   !   real(wp), allocatable :: kcn(:, :)

   !   !> Chemical hardness of each species
   !   real(wp), allocatable :: gam(:, :)

   !   !!> Solver for the linear system
   !   !procedure(solve), nopass, pointer :: solve_lineq => null()

   !   !!> Inversion of the linear system (for response)
   !   !procedure(solve), nopass, pointer :: invert_lineq => null()

   !contains

   !   !> Evaluate charge equilibration model
   !   procedure :: chargeEquilibration

   !end type TENEquilibration


contains


subroutine chargeEquilibration(nat, at, xyz, chrg, cn, dcndr, &
      & chi, kcn, gam, rad,&
      & energy, gradient, qat, dqdr)

   !> Source for error generation
   character(len=*), parameter :: source = 'eeq_chargeEquilibration'

   !!> Instance of the equilibration model
   !class(TENEquilibration), intent(inout) :: self
   !!> Calculation environment
   !type(TEnvironment), intent(inout) :: env
   !!> Molecular structure data
   !type(TMolecule), intent(in) :: mol
   !!> Instance of the Coulomb evaluator
   !class(TCoulomb), intent(inout) :: coulomb

   !> Number of atoms
   integer, intent(in) :: nat
   !> Atom types via the atomic number
   integer, intent(in) :: at(nat)
   !> Cartesian coordinates (a.u.)
   real(wp),intent(in) :: xyz(3,nat)
   !> Molecular charge
   integer, intent(in) :: chrg

   !> Electronegativity of each species 
   real(wp), intent(in) :: chi(:)
   !> Coordination number dependence of each species
   real(wp), intent(in) :: kcn(:)
   !> Chemical hardness of each species
   real(wp), intent(in) :: gam(:)
   !> Shell/Atomic hardnesses for each species (referred as 'alpha' in older versions)
   real(wp), intent(in) :: rad(:)
   !> Coordination number
   real(wp), intent(in) :: cn(:)
   !> Derivative of the coordination number w.r.t. Cartesian coordinates
   real(wp), intent(in) :: dcndr(:, :, :)

   !> Electrostatic energy
   real(wp), intent(inout), optional :: energy
   !> Molecular gradient
   real(wp), intent(inout), optional :: gradient(:, :)
   !> Atomic partial charges
   real(wp), intent(out), optional :: qat(:)
   !> Derivative of the partial charges w.r.t. Cartesian coordinates
   real(wp), intent(out), optional :: dqdr(:, :, :)

   real(wp), allocatable :: jmat(:, :), inv(:, :)
   real(wp), allocatable :: djdr(:, :, :), djdtr(:, :)
   real(wp), allocatable :: xvec(:), qvec(:), dxdcn(:), shift(:)
   real(wp), allocatable :: dxdr(:, :, :)

   logical :: deriv, response, exitRun
   integer :: nsh, ndim
   integer :: iat, ish, ii, iid
   real(wp) :: tmp

   deriv = present(gradient)
   response = present(dqdr)
   !nsh = sum(coulomb%itbl(2, :))
   nsh = nat
   ndim = nsh + 1
   allocate(jmat(ndim, ndim))
   allocate(xvec(ndim))
   allocate(qvec(ndim))
   allocate(dxdcn(nsh))

   !call getCoulombMatrix(nat, xyz, jmat)
   call getCoulombMatrixSmeared(nat, at, xyz, rad, jmat)

   !do iat = 1, mol%n
   !   ii = coulomb%itbl(1, iat)
   !   iid = mol%id(iat)
   !   do ish = 1, coulomb%itbl(2, iat)
   !      tmp = self%kcn(ish, iid) / (sqrt(cn(iat)) + 1.0e-14_wp)
   !      xvec(ii+ish) = -self%chi(ish, iid) + tmp*cn(iat)
   !      dxdcn(ii+ish) = 0.5_wp*tmp
   !      jmat(ii+ish, ii+ish) = jmat(ii+ish, ii+ish) + self%gam(ish, iid)
   !      jmat(ii+ish, ndim) = 1.0_wp
   !      jmat(ndim, ii+ish) = 1.0_wp
   !   end do
   !end do
   do iat = 1, nat
      iid = at(iat) !> atom type
      tmp = kcn(iid) / (sqrt(cn(iat)) + 1.0e-14_wp)
      xvec(iat) = -chi(iid) + tmp*cn(iat)
      dxdcn(iat) = 0.5_wp*tmp
      jmat(iat, iat) = jmat(iat, iat) + gam(iid)
      jmat(iat, ndim) = 1.0_wp
      jmat(ndim, iat) = 1.0_wp
   end do
   xvec(ndim) = float(chrg)
   jmat(ndim, ndim) = 0.0_wp

   inv = jmat
   if (response) then
      !call self%invert_lineq(env, inv, xvec, qvec)
      call solve_sytri(inv, xvec, qvec)
   else
      !call self%solve_lineq(env, inv, xvec, qvec)
      call solve_sysv(inv, xvec, qvec)
   end if

   !call env%check(exitRun)
   !if (exitRun) then
   !   call env%error("Solving linear equations failed", source)
   !   return
   !end if

   if (present(qat)) then
      qat(:) = 0.0_wp
      !do iat = 1, mol%n
      !   ii = coulomb%itbl(1, iat)
      !   do ish = 1, coulomb%itbl(2, iat)
      !      qat(iat) = qat(iat) + qvec(ii+ish)
      !   end do
      !end do
      do iat = 1, nat
        qat(iat) = qat(iat) + qvec(iat)
      end do
   end if

   if (present(energy)) then
      shift = xvec
      call mctc_symv(jmat, qvec, shift, alpha=0.5_wp, beta=-1.0_wp)
      energy = energy + mctc_dot(qvec, shift)
   end if

   if (deriv .or. response) then
      allocate(djdr(3, nat, ndim))
      allocate(djdtr(3, ndim))
      allocate(dxdr(3, nat, ndim))
      !call getCoulombDerivs(nat, xyz, qvec, djdr, djdtr)
      call getCoulombDerivsSmeared(nat, at, xyz, rad, qvec, djdr, djdtr)
      do iat = 1, nat
        dxdr(:, :, iat) = -dxdcn(iat)*dcndr(:, :, iat)
      end do
      dxdr(:, :, ndim) = 0.0_wp

      if (present(gradient)) then
         call contract(djdr, qvec, gradient, beta=1.0_wp) !>contract312
         call contract(dxdr, qvec, gradient, beta=1.0_wp) !>contract312
      end if

      if (response) then
         do iat = 1, nat
           djdr(:, iat, iat) = djdr(:, iat, iat) + djdtr(:, iat)
         end do

         if (present(dqdr)) then
            dqdr(:, :, :) = 0.0_wp
            call contract(djdr, inv(:, :nsh), dqdr, alpha=-1.0_wp, beta=0.0_wp) !>contract323
            call contract(dxdr, inv(:, :nsh), dqdr, alpha=-1.0_wp, beta=1.0_wp) !>contract323
         end if

      end if
   end if

end subroutine chargeEquilibration


!========================================================================================!
subroutine solve_sysv(mat, rhs, vec)

   !> Source for error generation
   character(len=*), parameter :: source = 'xtb_eeq_sysv'

   real(wp), intent(inout) :: mat(:, :)

   real(wp), intent(in) :: rhs(:)

   real(wp), intent(out), target :: vec(:)

   integer :: lwork, ndim, info
   real(wp), pointer :: ptr(:, :)
   real(wp) :: test(1)
   real(wp), allocatable :: work(:)
   integer, allocatable :: ipiv(:)

   ndim = size(mat, dim=1)
   if (size(mat, dim=2) < ndim .or. size(vec) < ndim .or. size(rhs) < ndim) then
      write(0,'(a,1x,a)')"A not so carefully crafted algorithm did mess up", source
      return
   end if

   vec(:ndim) = rhs(:ndim)
   allocate(ipiv(ndim))
   ptr(1:ndim, 1:1) => vec(1:ndim)

   ! assume work space query, set best value to test after first dsysv call
   call dsysv('l', ndim, 1, mat, ndim, ipiv, ptr, ndim, test, -1, info)
   lwork = int(test(1))
   allocate(work(lwork))

   call dsysv('l', ndim, 1, mat, ndim, ipiv, ptr, ndim, work, lwork, info)

   if (info > 0) then
      write(0,'(a,1x,a)')"LAPACK linear equation solver failed", source
      return
   end if

end subroutine solve_sysv


!========================================================================================!
subroutine solve_sytri(mat, rhs, vec)

   !> Source for error generation
   character(len=*), parameter :: source = 'xtb_eeq_sytri'

   !!> Calculation environment
   !type(TEnvironment), intent(inout) :: env

   real(wp), intent(inout) :: mat(:, :)

   real(wp), intent(in) :: rhs(:)

   real(wp), intent(out), target :: vec(:)

   integer :: ndim
   logical :: exitRun

   ndim = size(mat, dim=1)
   if (size(mat, dim=2) < ndim .or. size(vec) < ndim .or. size(rhs) < ndim) then
      write(0,'(a,1x,a)')"A not so carefully crafted algorithm did mess up", source 
      return
   end if

   call invert_sytri(mat)

   call mctc_symv(mat, rhs, vec)

end subroutine solve_sytri

!========================================================================================!
subroutine invert_sytri(mat)

   !> Source for error generation
   character(len=*), parameter :: source = 'eeq_sytri'

   real(wp), intent(inout) :: mat(:, :)

   integer :: ii, jj
   integer :: lwork, ndim, info
   real(wp) :: test(1)
   real(wp), allocatable :: work(:)
   integer, allocatable :: ipiv(:)

   ndim = size(mat, dim=1)
   if (size(mat, dim=2) < ndim) then
      write(0,'(a,1x,a)')"A not so carefully crafted algorithm did mess up", source
      return
   end if

   allocate(ipiv(ndim))

   ! assume work space query, set best value to test after first dsytrf call
   call lapack_sytrf('L', ndim, mat, ndim, ipiv, test, -1, info)
   lwork = int(test(1))
   allocate(work(lwork))

   ! Bunch-Kaufman factorization A = L*D*L**T
   call lapack_sytrf('L', ndim, mat, ndim, ipiv, work, lwork, info)
   if(info > 0)then
      write(0,'(a,1x,a)') 'Bunch-Kaufman factorization failed', source
      return
   endif

   ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in matrix
   ! matrix is overwritten with lower triangular part of A⁻¹
   call lapack_sytri('L', ndim, mat, ndim, ipiv, work, info)
   if (info > 0) then
      write(0,'(a,1x,a)') 'Bunch-Kaufman inversion failed', source
      return
   endif

   ! symmetrizes A⁻¹ matrix from lower triangular part of inverse matrix
   do ii = 1, ndim
      do jj = ii+1, ndim
         mat(ii,jj) = mat(jj,ii)
      enddo
   enddo

end subroutine invert_sytri

!========================================================================================!
!> subroutine getCoulombMatrix
!> Compute all elements |r_ij|⁻¹
subroutine getCoulombMatrix(nat, xyz, jmat)
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp),intent(in) :: xyz(3,nat)
   !> Coulomb matrix
   real(wp), intent(out) :: jmat(:, :)

   integer :: iat, jat, ish, jsh, ii, jj
   real(wp) :: r1, rterm
   jmat(:, :) = 0.0_wp

   !$omp parallel do default(none) shared(nat, xyz, jmat) &
   !$omp private(iat, jat, r1, rterm)
   do iat = 1, nat
      do jat = 1, iat-1
         r1 = sqrt(sum((xyz(:, jat) - xyz(:, iat))**2))
         rterm = 1.0_wp/r1
         jmat(jat, iat) = rterm
         jmat(iat, jat) = rterm
      end do
   end do
   !$omp end parallel do
end subroutine getCoulombMatrix
!> subroutine getCoulombMatrixSmeard
!> Compute all elements of the Coulomb matrix with smeard Gaussian charges
subroutine getCoulombMatrixSmeared(nat, at, xyz, rad, jmat)

   !> Number of atoms
   integer, intent(in) :: nat
   !> Atomic types by atom number
   integer, intent(in) :: at(nat)
   !> Cartesian coordinates
   real(wp),intent(in) :: xyz(3,nat)

   !> Shell/Atomic hardnesses for each species
   real(wp), intent(in) :: rad(:)

   !> Coulomb matrix
   real(wp), intent(out) :: jmat(:, :)

   integer :: iat, jat, ish, jsh, ii, jj, iid, jid
   real(wp) :: r1, rterm, gij

   jmat(:, :) = 0.0_wp

   !$omp parallel do default(none) shared(nat, xyz, rad, jmat) &
   !$omp private(iat, jat, iid, jid, r1, rterm, gij)
   do iat = 1, nat
      iid = at(iat)
      do jat = 1, iat-1
         jid = at(jat)
         r1 = sqrt(sum((xyz(:, jat) - xyz(:, iat))**2))
         gij = 1.0_wp/sqrt(rad(iid)**2 + rad(jid)**2)
         rterm = erf(gij*r1)/r1
         jmat(jat, iat) = rterm
         jmat(iat, jat) = rterm
      end do
      gij = sqrt(0.5_wp)/rad(iid)
      jmat(iat, iat) = 2.0_wp/sqrtpi*gij
   end do
   !$omp end parallel do

end subroutine getCoulombMatrixSmeared


!========================================================================================!
!> subroutine getCoulombDerivs
!> Compute all derivaties ∂(q_j|r_ij|⁻¹)/∂r_i 
subroutine getCoulombDerivs(nat, xyz, qvec, djdr, djdtr)

   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp),intent(in) :: xyz(3,nat)

   !> Charges
   real(wp), intent(in) :: qvec(:)

   !> Derivative of Coulomb matrix w.r.t. Cartesian coordinates
   real(wp), intent(out) :: djdr(:, :, :)

   !> Trace derivative of Coulomb matrix
   real(wp), intent(out) :: djdtr(:, :)

   integer :: iat, jat, ish, jsh, ii, jj
   real(wp) :: r1, vec(3), dG(3), dS(3, 3)

   djdr(:, :, :) = 0.0_wp
   djdtr(:, :) = 0.0_wp

   !$omp parallel do default(none) reduction(+:djdr, djdtr) &
   !$omp shared(nat, xyz, qvec) &
   !$omp private(iat, jat, r1, vec, dG, dS)
   do iat = 1, nat
      do jat = 1, iat-1
         vec(:) = xyz(:, jat) - xyz(:, iat)
         r1 = norm2(vec)
         dG(:) = -vec/r1**3
         dS(:, :) = 0.5_wp * spread(dG, 1, 3) * spread(vec, 2, 3)
         djdr(:, iat, jat) = djdr(:, iat, jat) - dG*qvec(iat)
         djdr(:, jat, iat) = djdr(:, jat, iat) + dG*qvec(jat)
         djdtr(:, jj+jsh) = djdtr(:, jat) + dG*qvec(iat)
         djdtr(:, ii+ish) = djdtr(:, iat) - dG*qvec(jat)
      end do
   end do
   !$omp end parallel do

end subroutine getCoulombDerivs
!> subroutine getCoulombDerivsSmeared
!> Compute all Coulomb matrix derivatives with smeared Gaussian charges
subroutine getCoulombDerivsSmeared(nat, at, xyz, rad, qvec, djdr, djdtr)

   !> Number of atoms
   integer, intent(in) :: nat
   !> Atomic types by atom number
   integer, intent(in) :: at(nat)
   !> Cartesian coordinates
   real(wp),intent(in) :: xyz(3,nat)

   !> Shell/Atomic hardnesses for each species
   real(wp), intent(in) :: rad(:)

   !> Charges
   real(wp), intent(in) :: qvec(:)

   !> Derivative of Coulomb matrix w.r.t. Cartesian coordinates
   real(wp), intent(out) :: djdr(:, :, :)

   !> Trace derivative of Coulomb matrix
   real(wp), intent(out) :: djdtr(:, :)

   integer :: iat, jat, ish, jsh, ii, jj, iid, jid
   real(wp) :: r2, g1, gij, vec(3), dG(3), dS(3, 3)

   djdr(:, :, :) = 0.0_wp
   djdtr(:, :) = 0.0_wp

   !$omp parallel do default(none) reduction(+:djdr, djdtr) &
   !$omp shared(nat, at, xyz, qvec, rad) &
   !$omp private(iat, jat, iid, jid, r2, g1, gij, vec, dG, dS)
   do iat = 1, nat
      iid = at(iat)
      do jat = 1, iat-1
         jid = at(jat)
         vec(:) = xyz(:, jat) - xyz(:, iat)
         r2 = sum(vec**2)
         gij = 1.0_wp/(rad(iid)**2 + rad(jid)**2)
         g1 = erf(sqrt(gij*r2))/sqrt(r2)
         dG(:) = (2*sqrt(gij)*exp(-gij*r2)/sqrtpi - g1) * vec/r2
         dS(:, :) = 0.5_wp * spread(dG, 1, 3) * spread(vec, 2, 3)
         djdr(:, iat, jat) = djdr(:, iat, jat) - dG*qvec(iat)
         djdr(:, jat, iat) = djdr(:, jat, iat) + dG*qvec(jat)
         djdtr(:, jat) = djdtr(:, jat) + dG*qvec(iat)
         djdtr(:, iat) = djdtr(:, iat) - dG*qvec(jat)
      end do
   end do
   !$omp end parallel do

end subroutine getCoulombDerivsSmeared


!========================================================================================!
end module eeq_module
