module math_wrapper
    use iso_fortran_env, only: wp=>real64
    implicit none
    public

    interface contract
       module procedure contract312
       module procedure contract323
    end interface contract

contains
!========================================================================================!
subroutine contract312(amat, bvec, cvec, alpha, beta)

   real(wp), intent(in), contiguous, target :: amat(:, :, :)
   real(wp), intent(in), contiguous :: bvec(:)
   real(wp), intent(inout), contiguous, target :: cvec(:, :)
   real(wp), intent(in), optional :: alpha
   real(wp), intent(in), optional :: beta
   real(wp), pointer :: aptr(:, :)
   real(wp), pointer :: cptr(:)

   integer :: nn, mm
   real(wp) :: aa, bb
   !> BLAS
   external :: dgemv

   aa = 1.0_wp
   if (present(alpha)) aa = alpha
   bb = 0.0_wp
   if (present(beta)) bb = beta

   nn = size(amat, dim=1)*size(amat, dim=2)
   mm = size(amat, dim=3)
   aptr(1:size(amat, dim=1) * size(amat, dim=2), 1:size(amat, dim=3)) => amat
   cptr(1:size(cvec, dim=1) * size(cvec, dim=2)) => cvec

   call dgemv('n', nn, mm, aa, aptr, nn, bvec, 1, bb, cptr, 1)

end subroutine contract312
!========================================================================================!
subroutine contract323(amat, bmat, cmat, alpha, beta)

   real(wp), intent(in), contiguous, target :: amat(:, :, :)
   real(wp), intent(in), contiguous :: bmat(:, :)
   real(wp), intent(inout), contiguous, target :: cmat(:, :, :)
   real(wp), intent(in), optional :: alpha
   real(wp), intent(in), optional :: beta
   real(wp), pointer :: aptr(:, :)
   real(wp), pointer :: cptr(:, :)
   
   integer :: nn, mm, kk
   real(wp) :: aa, bb
   !> BLAS
   external :: dgemm

   aa = 1.0_wp
   if (present(alpha)) aa = alpha
   bb = 0.0_wp
   if (present(beta)) bb = beta

   nn = size(amat, dim=1)*size(amat, dim=2)
   mm = size(amat, dim=3)
   kk = size(cmat, dim=3)
   aptr(1:size(amat, dim=1) * size(amat, dim=2), 1:size(amat, dim=3)) => amat
   cptr(1:size(cmat, dim=1) * size(cmat, dim=2), 1:size(cmat, dim=3)) => cmat

   call dgemm('n', 'n', nn, kk, mm, aa, aptr, nn, bmat, mm, bb, cptr, nn)

end subroutine contract323
!========================================================================================!





end module math_wrapper
