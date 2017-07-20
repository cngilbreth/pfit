! pfit.f90: Module for polynomial least-squares fitting
! http://infty.net/pfit/pfit.html
! v0.8.1
!
! Copyright (c) 2010-2013 Christopher N. Gilbreth
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module mod_pfit

  integer, parameter, private :: rk = selected_real_kind(p=15)

contains


  subroutine pfit(x,y,sig,a,cov,coeff)
    ! Fit data to a polynomial a_0 + a_1 x + ... + a_d x**d
    ! Inputs:
    !   x(1:npt)         - abscissas
    !   y(1:npt)         - data values
    !   sig(1:npt)       - data errors
    ! Outputs:
    !   a(1:d+1)         - max. likelihood parameters
    !   coeff(1:npt,1:d+1) - coefficients giving the max. likelihood parameters
    !                    in terms of the data:
    !                      a(i) = \Sum_{j} coeff(j,i) * y(j)
    !   cov(1:n,1:d+1)   - Covariance matrix, cov(i,j) = Cov(a(i),a(j))
    !                    The estimated error in a(i) is sqrt(Cov(a(i),a(i))).
    ! Notes:
    !   This routine uses a QR decomposition method, which should be more
    !   numerically stable than solving the normal equations.
    implicit none
    real(rk), intent(in)  :: x(:), y(:), sig(:)
    real(rk), intent(out) :: a(:)
    real(rk), intent(out), optional :: cov(:,:), coeff(:,:)

    real(rk), allocatable :: work(:), C(:,:), Q(:,:), R(:,:), b(:)
    integer :: ipiv(size(a)), lwork
    integer :: i,j,k,d,npt,ifail
    real(rk) :: coeff1(size(x),size(a))

    npt = size(x) ! Number of data points
    d = size(a)-1 ! Max degree of polynomial

    if (size(y) .ne. npt) stop "Error 1 in pfit"
    if (size(sig) .ne. npt) stop "Error 2 in pfit"
    if (d+1 .gt. npt) stop "Error 4 in pfit"
    if (present(coeff)) then
       if (size(coeff,1) .ne. npt) stop "Error 6 in pfit"
       if (size(coeff,2) .ne. d+1) stop "Error 5 in pfit"
    end if
    if (present(cov)) then
       if (size(cov,1) .ne. d+1) stop "Error 7 in pfit"
       if (size(cov,2) .ne. d+1) stop "Error 8 in pfit"
    end if

    allocate(C(npt,d+1), Q(npt,d+1), R(d+1,d+1), b(d+1))

    ! Vandermonde matrix
    do j=1,d+1
       do i=1,npt
          C(i,j) = x(i)**(j-1)/sig(i)
       end do
    end do

    ! QR decomposition
    call DQRF(C,Q,R,work)

    ! Inversion of R factor
    ! TODO: Should instead be able to efficiently just multiply R^-1 by (Q^T y)
    ! since R is right triangular at Q^T y is a single column vector.
    R = dinverse(R)

    ! Compute max-likelihood parameters
    ! a = R^-1 Q^T y/σ
    b = 0.d0
    do j=1,d+1
       do k=1,npt
          b(j) = b(j) + Q(k,j) * y(k) / sig(k)
       end do
    end do

    a = 0.d0
    do i=1,d+1
       do j=1,d+1
          a(i) = a(i) + R(i,j) * b(j)
       end do
    end do

    ! Compute coefficient matrix coeff such that a(i) = Σ_j coeff(j,i) y(j)
    ! Here a(i) = R^{-1}(i,j) Q(k,j) y(k)/σ(k)
    ! So coeff(k,i) = R^{-1}(i,j) Q(k,j) / σ(k)
    if (present(coeff) .or. present(cov)) then
       coeff1 = 0.d0
       do i=1,d+1
          do j=1,d+1
             do k=1,npt
                coeff1(k,i) = coeff1(k,i) + R(i,j) * Q(k,j) / sig(k)
             end do
          end do
       end do
       if (present(coeff)) coeff = coeff1
    end if

    ! Compute covariance matrix Cov(a(i),a(j)) = Σ_k C(k,i) C(k,j) σ(k)^2
    if (present(cov)) then
       cov = 0.d0
       do j=1,d+1
          do i=1,d+1
             do k=1,npt
                cov(i,j) = cov(i,j) + coeff1(k,i) * coeff1(k,j) * sig(k)**2
             end do
          end do
       end do
    end if
  end subroutine pfit


  subroutine DQRF(A,Q,R,work)
    ! Compute the QR factorization of a general real matrix A:
    !   A = Q R
    ! where Q is unitary and R is upper triangular, using the LAPACK routine
    ! zgeqrf.
    ! Inputs:
    !   A:     Matrix to be factorized, m x n
    ! Ouputs:
    !   Q:     Unitary matrix, m x m
    !   R:     Upper triangular, n x n
    ! Input/output:
    !   work:  real(8) allocatable workspace array. If unallocated, this
    !          routine will allocate it to an appropriate size. If allocated,
    !          it is assumed to be the correct size for this problem.
    implicit none
    real(8), intent(in)  :: A(:,:)
    real(8), intent(out) :: Q(:,:), R(:,:)
    real(8), allocatable :: work(:)

    integer :: m, n, lwork, ierr, i, j
    real(8) :: tau(size(A,2)), qwork(1)
    real(8) :: A1(size(A,1),size(A,2))

    m = size(A,1)
    n = size(A,2)
    if (m .lt. n) stop "Error in DQRF: m < n"
    if (size(Q,1) .ne. m) stop "Error in DQRF (2)"
    if (size(Q,2) .ne. n) stop "Error in DQRF (3)"
    if (size(R,1) .ne. n) stop "Error in DQRF (4)"
    if (size(R,2) .ne. n) stop "Error in DQRF (5)"

    A1 = A
    if (.not. allocated(work)) then
       ! Compute size of workspace
       lwork = -1
       call DGEQRF(m, n, A1, m, TAU, qwork, LWORK, ierr)
       if (ierr .ne. 0) stop "Error calling DGEQRF (1)"
       lwork = qwork(1)
       allocate(work(lwork))
    end if

    lwork = size(work)
    call dgeqrf(m,n,A1,m,tau,work,lwork,ierr)
    if (ierr .ne. 0) stop "Error calling DGEQRF (2)"
    R = 0.d0
    do j=1,n
       do i=1,j
          R(i,j) = A1(i,j)
       end do
    end do
    Q(:,1:n) = A1
    call dorgqr(m,n,n,Q,m,tau,work,lwork,ierr)
    if (ierr .ne. 0) stop "Error calling DORGQR"
  end subroutine DQRF


  function dinverse(A)
    ! Invert a square matrix
    implicit none
    real(8), intent(in)  :: A(:,:)
    real(8) :: dinverse(size(A,1),size(A,1))

    integer :: ipiv(size(A,1)), ierr, lwork
    real*8, allocatable :: work(:)
    real*8 :: work1(1)

    dinverse = A
    call dgetrf(size(A,1), size(A,1), dinverse, size(A,1), ipiv, ierr)
    if (ierr .ne. 0) stop "Error computing LU decomposition for matrix inverse."

    lwork = -1
    call dgetri(size(A,1), dinverse, size(A,1), ipiv, work1, lwork, ierr)
    if (ierr.ne.0) stop "Error allocating space for dgetri"
    lwork = int(work1(1),kind(lwork))
    allocate(work(max(1,lwork)))
    call dgetri(size(A,1), dinverse, size(A,1), ipiv, work, lwork, ierr)
    if (ierr .ne. 0) stop "Error calling zgetri."
  end function dinverse


end module mod_pfit
