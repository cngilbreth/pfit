! Simple command-line polynomial fitting program
! C. N. Gilbreth, 2010 cngilbreth@gmail.com
#include "util.h"

program prog_pfit
  use util
  use mod_options
  use math, only: pfit
  implicit none

  type(options) :: opts
  integer :: nargs
  character(len=256) :: filename
  real*8, allocatable :: x_raw(:), y_raw(:), errs_raw(:), a(:), cov(:,:)
  real*8, allocatable :: x(:), y(:), errs(:)
  logical, allocatable :: mask(:)
  real*8  :: ub, lb
  integer :: degree, i, ierr, j
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

  nargs = get_num_args(opts)
  if (nargs .ne. 1) then
     call print_help(opts)
     stop
  end if

  ! First argument is the data file
  call get_arg(opts,1,filename,"")
  call read_data(filename,x_raw,y_raw,errs_raw)

  call get_option(opts,"degree", degree)
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
    type(options) :: opts

    write (*,'(a)') "pfit: fit a univariate polynomial f(x) to some data"
    write (*,'(a)') "Usage: pfit [options] <filename>"
    write (*,'(a)') ''
    call print_info(opts)
  end subroutine print_help



  subroutine read_data(filename,x,y,errs)
    implicit none
    character(len=*) :: filename
    real*8, allocatable :: x(:), y(:), errs(:)

    integer, parameter :: MAX_COLS = 16
    integer, parameter :: MAX_LINE_LEN = 1024

    integer :: line
    character(len=MAX_LINE_LEN) :: buf
    real*8  :: vals(MAX_COLS)
    integer :: nvals, char_col, val_col


    open(unit=1,file=filename)
    ! Two passes: (1) count lines, (2) load data
    line = 0
    nvals = 0
    do
       line = line + 1
       read (1,'(a'//trim(itos(MAX_LINE_LEN))//')',end=90) buf
       if (is_comment_or_blank(buf)) cycle
       nvals = nvals + 1
    end do
90  continue
    close(1)

    allocate(x(nvals), y(nvals), errs(nvals))
    x = 0.d0; y = 0.d0; errs = 1.d0
    open(unit=1,file=filename)

    line = 0
    nvals = 0
    ! For each line in the file ...
    do
       line = line + 1
       read (1,'(a'//trim(itos(MAX_LINE_LEN))//')',end=91) buf
       if (is_comment_or_blank(buf)) cycle

       ! Read all the numerical values ...
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
10     continue

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
    end do
91  continue
  end subroutine read_data


end program prog_pfit
