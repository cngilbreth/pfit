! Simple command-line polynomial fitting program
! C. N. Gilbreth, 2010 cngilbreth@gmail.com
#include "util.h"

program prog_pfit
  use mod_util
  use mod_options
  implicit none

  type(options_t) :: opts
  character(len=256) :: filename
  real*8, allocatable :: x_raw(:), y_raw(:), errs_raw(:), a(:), cov(:,:)
  real*8, allocatable :: x(:), y(:), errs(:)
  logical, allocatable :: mask(:)
  real*8  :: ub, lb
  integer :: degree, i, ierr, j, index
  logical :: print_corr

  call define_option(opts,"degree",1,&
       abbrev='d', &
       description="Degree of polynomial to fit (required)", &
       required=.true., &
       min=0)
  call define_option(opts,"lb",0.d0,&
       description="Lower bound of fit range")
  call define_option(opts,"ub",1.d0,&
       description="Upper bound of fit range")
  call define_option(opts,"index",0,&
       abbrev="i",&
       description="Index of dataset in input, starting from 0.&
       & (Datasets are separated by >= two blank lines.&
       & Default: use all data.)", &
       min=0)
  call define_flag(opts,"corr",&
       description="Print correlation matrix")

  if (command_argument_count() == 0) then
     call print_help(opts)
     stop
  end if

  call process_command_line(opts,ierr)
  if (ierr .ne. 0) then
     call print_help(opts)
     stop
  end if


  ! First argument is the data file
  call get_arg(opts,1,filename,ierr)
  if (ierr .ne. 0) then
     write (0,*) "Error: not enough arguments."
     call print_help(opts)
     stop
  end if
  call get_option(opts,"degree", degree)
  if (opt_found(opts,"index")) then
     call get_option(opts,"index",index)
  else
     index = -1
  end if
  call read_data(filename,index,x_raw,y_raw,errs_raw)

  if (size(x_raw) == 0) then
     write (*,'(a)') "Error: no data found"
     stop
  else if (size(x_raw) < degree) then
     write (*,'(a,i0,a)') "Error: not enough data for a degree ", degree, " fit."
  end if

  allocate(a(degree+1), cov(degree+1,degree+1))

  ! Process upper and lower bounds

  allocate(mask(size(x_raw)))
  mask = .true.
  if (opt_found(opts,"lb")) then
     call get_option(opts,"lb",lb)
     do i=1,size(x_raw)
        if (x_raw(i) < lb) then
           mask(i) = .false.
        end if
     end do
  end if
  if (opt_found(opts,"ub")) then
     call get_option(opts,"ub",ub)
     do i=1,size(x_raw)
        if (x_raw(i) > ub) then
           mask(i) = .false.
        end if
     end do
  end if

  if (size(x_raw) == 0) then
     write (*,'(a)') "Error: no data found"
     stop
  else if (size(x_raw) < degree) then
     write (*,'(a,i0,a)') "Error: not enough data for a degree ", degree, " fit."
  end if

  if (count(mask) < degree) then
     write (*,'(a)') "Error: not enough data within bounds for this fit."
     stop
  end if

  allocate(x(count(mask)), y(count(mask)), errs(count(mask)))
  x = pack(x_raw,mask)
  y = pack(y_raw,mask)
  errs = pack(errs_raw,mask)

  ! Fit

  call pfit(x,y,errs,a,cov)

  ! Print results

  write (*,'(a)',advance='no') "# Result of fitting data to f(x) = "
  do i=0,degree
     if (i > 0) then
        write (*,'(a)', advance='no') " + "
     end if
     write (*,'(a,i0)',advance='no') "a", i
     if (i > 0) then
        write (*,'(a,i0)',advance='no') " x^", i
     end if
  end do
  write (*,'(a)') ':'
  write (*,'(a,t8,a20,a20)') "#Coeff", "estimate", "error"
  do i=0,degree
     write (*,'(a,i0,t8,es20.10,es20.10)') "a", i, a(i+1), sqrt(cov(i+1,i+1))
  end do

  call get_option(opts,"corr",print_corr)
  if (print_corr) then
     write (*,'(a)') ""
     write (*,'(a)') "# Correlation matrix:"
     write (*,'(a,t4)',advance='no') '#'
     do j=1,size(a)
        write (*,'(a20,tr2)',advance='no') "a"//trim(to_str(j-1))
     end do
     write (*,'(a)') ''
     do i=1,size(a)
        write (*,'(a,t4)',advance='no') "a"//trim(to_str(i-1))
        do j=1,size(a)
           write (*,'(es20.10,tr2)',advance='no') cov(i,j)/sqrt(cov(i,i) * cov(j,j))
        end do
        write (*,'(a)') ''
     end do
  end if

contains

  subroutine print_help(opts)
    implicit none
    type(options_t) :: opts

    write (*,'(a)') "pfit: fit a univariate polynomial f(x) to some data"
    write (*,'(a)') "Usage: pfit [options] <filename>"
    write (*,'(a)') ''
    call print_info(opts)
  end subroutine print_help



  subroutine read_data(filename,index,x,y,errs)
    implicit none
    character(len=*) :: filename
    integer, intent(in) :: index
    real*8, allocatable :: x(:), y(:), errs(:)

    integer, parameter :: MAX_COLS = 16
    integer, parameter :: MAX_LINE_LEN = 1024

    integer :: line
    character(len=MAX_LINE_LEN) :: buf
    real*8  :: vals(MAX_COLS)
    integer :: nvals, char_col, val_col, current_index, n_prev_blank_lines


    open(unit=1,file=filename)
    ! Two passes: (1) count lines, (2) load data
    line = 0
    nvals = 0
    current_index = 0
    n_prev_blank_lines = 0
    do
       line = line + 1
       read (1,'(a'//trim(itos(MAX_LINE_LEN))//')',end=90) buf
       if (is_comment(buf)) then
          cycle
       else if (is_blank(buf)) then
          n_prev_blank_lines = n_prev_blank_lines + 1
          cycle
       else
          if (n_prev_blank_lines >= 2) then
             current_index = current_index + 1
             n_prev_blank_lines = 0
          end if
       end if

       if (current_index == index .or. index < 0) then
          nvals = nvals + 1
       end if
    end do
90  continue
    close(1)

    allocate(x(nvals), y(nvals), errs(nvals))
    x = 0.d0; y = 0.d0; errs = 1.d0
    open(unit=1,file=filename)

    line = 0
    nvals = 0
    current_index = 0
    n_prev_blank_lines = 0
    ! For each line in the file ...
    do
       line = line + 1
       read (1,'(a'//trim(itos(MAX_LINE_LEN))//')',end=91) buf
       if (is_comment(buf)) then
          cycle
       else if (is_blank(buf)) then
          n_prev_blank_lines = n_prev_blank_lines + 1
          cycle
       else
          if (n_prev_blank_lines >= 2) then
             current_index = current_index + 1
             n_prev_blank_lines = 0
          end if
       end if

       if (current_index == index .or. index < 0) then
          ! Read the numerical values ...
          char_col=1; val_col=0
          do while (val_col < MAX_COLS)
             ! FIXME: How can we handle tabs here?
             do while (buf(char_col:char_col) == ' ')
                if (char_col == len(buf)) goto 10
                char_col = char_col + 1
             end do
             val_col = val_col + 1
             ASSERT(val_col <= MAX_COLS, "Too many columns in input; increase MAX_COLS")
             read (buf(char_col:),*) vals(val_col)
             do while (buf(char_col:char_col) .ne. ' ')
                if (char_col == len(buf)) goto 10
                char_col = char_col + 1
             end do
          end do
10        continue

          ! And pick out the ones we'll use for the fit
          ! (Can be generalized later)
          ASSERT(val_col >= 2, "Missing data on line "//itos(line))
          if (val_col >= 2) then
             nvals = nvals + 1
             x(nvals) = vals(1)
             y(nvals) = vals(2)
          end if
          if (val_col >= 3) then
             errs(nvals) = vals(3)
          end if
       end if
    end do
91  continue
  end subroutine read_data


  logical function is_comment(str)
    implicit none
    character(len=*) :: str
    character(len=1) :: first_nonblank_char

    if (len(str) == 0) then
       is_comment = .false.
       return
    end if

    first_nonblank_char = adjustl(str)

    if (first_nonblank_char == '#' .or. first_nonblank_char == '!') then
       is_comment = .true.
    else
       is_comment = .false.
    end if
  end function is_comment


  logical function is_blank(str)
    implicit none
    character(len=*) :: str
    character(len=1) :: first_nonblank_char

    if (len(str) == 0) then
       is_blank = .true.
       return
    end if

    first_nonblank_char = adjustl(str)

    if (first_nonblank_char == ' ') then
       is_blank = .true.
    else
       is_blank = .false.
    end if
  end function is_blank


  subroutine pfit(x,y,sig,a,cov,coeff)
    ! Fit data to a polynomial a_0 + a_1 x + ... + a_(m-1) x**(m-1)
    ! Inputs:
    !   x(1:n)         - abscissas
    !   y(1:n)         - data values
    !   sig(1:n)       - data errors
    ! Outputs:
    !   a(1:m)         - max. likelihood parameters
    !   cov(1:m,1:m)   - covariance matrix
    !   coeff(1:n,1:m) - coefficients giving the max. likelihood parameters
    !                    in terms of the data:
    !                      a(i) = \Sum_{j} coeff(j,i) * y(j)
    ! Requirements:
    !   Each data point y(i) must be pulled from a normal distribution with std.
    !   dev. sig(i). Otherwise the covariance matrix will not give useful error
    !   estimates.
    implicit none
    real*8, intent(in)  :: x(:), y(:), sig(:)
    real*8, intent(out) :: a(:), cov(:,:)
    real*8, intent(out), optional :: coeff(:,:)
    real*8 :: work(size(a)*10)

    integer :: ipiv(size(a)), lwork
    integer :: i,j,k,n,m,ifail

    if (size(x) .ne. size(y)) stop "Error 1 in pfit"
    if (size(sig) .ne. size(y)) stop "Error 2 in pfit"
    if (size(a) .ne. size(cov,dim=1)) stop "Error 3 in pfit"
    if (size(a) .gt. size(x)) stop "Error 4 in pfit"
    if (present(coeff)) then
       if (size(a) .ne. size(coeff,2)) stop "Error 5 in pfit"
       if (size(x) .ne. size(coeff,1)) stop "Error 6 in pfit"
    end if

    n = size(x) ! Number of data points
    m = size(a) ! Number of parameters
    cov = 0.d0; a = 0.d0
    do j=1,m
       do k=1,m
          do i=1,n
             cov(j,k) = cov(j,k) + (x(i)**(j-1) * x(i)**(k-1))/sig(i)**2
          end do
       end do
       do i=1,n
          a(j) = a(j) + y(i) * x(i)**(j-1) / sig(i)**2
       end do
    end do

    ! Invert the matrix
    call dgetrf(m,m,cov,m,ipiv,ifail)
    if (ifail .ne. 0) stop "LU decomposition failed in pfit"
    lwork = size(work)
    call dgetri(m,cov,m,ipiv,work,lwork,ifail)
    if (ifail .ne. 0) stop "Inversion failed in pfit"

    a = matmul(cov,a)

    if (present(coeff)) then
       coeff = 0.d0
       do j=1,m
          do i=1,n
             do k=1,m
                coeff(i,j) = coeff(i,j) + cov(j,k) * x(i)**(k-1) / sig(i)**2
             end do
          end do
       end do
    end if
  end subroutine pfit


end program prog_pfit
