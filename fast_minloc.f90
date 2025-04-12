! by Ivan Pribec, lightly edited by Beliavsky
! https://fortran-lang.discourse.group/t/speed-of-array-intrinsics/3251/16
! Compares different methods for computing the index of the minimum 
! value in a real array, including the Fortran built-in `minloc` function and 
! user-defined parallel versions using OpenMP: a threaded approach, a custom 
! reduction, and SIMD-based reduction.
module fast_minloc
implicit none
private

public :: minloc_thread
public :: minloc_reduce
public :: minloc_simd

type :: minloc_type
    real :: val = huge(1.0)
    integer :: idx = 0
end type
!$omp declare reduction(min:minloc_type:omp_out=minloc_min(omp_out,omp_in))

contains

    pure function minloc_min(a,b) result(c)
        type(minloc_type), intent(in) :: a, b
        type(minloc_type) :: c
        c = merge(a,b,a%val <= b%val)
    end function

    pure function minloc_reduce(a) result(idx)
        real, intent(in) :: a(:)
        integer :: idx

        type(minloc_type) :: t
        integer :: i

        t = minloc_type(huge(a),0)
        !$omp parallel do reduction(min: t)
        do i = 1, size(a)
            t = minloc_min(t,minloc_type(a(i),i))
        end do

        idx = t%idx

    end function

    pure function minloc_simd(a) result(idx)
        real, intent(in) :: a(:)
        integer :: idx

        ! Type needed for custom reduction
        type :: ridx
           real :: val = huge(1.0)
           integer :: idx = 0
        end type

!$omp declare reduction(min: ridx: &
!$omp& omp_out = merge(omp_out,omp_in,omp_out%val < omp_in%val))

        integer :: i
        type(ridx) :: r

        if (size(a) == 0) then
            idx = 0
            return
        else
            r = ridx(a(1),1)            
            !$omp simd reduction(min: r)
            do i = 2, size(a)
                if (a(i) < r%val) r = ridx(a(i),i)
            end do
            idx = r%idx
        end if

    end function

    pure function minloc_thread(a) result(idx)
        !$ use omp_lib
        real, intent(in) :: a(:)

        integer :: n, nt, it, nc, lb, ub, mt
        integer :: idx

        real, allocatable :: b(:)
        integer, allocatable :: t(:)

        n = size(a)

        if (n == 0) then
            idx = 0
            return
        else

        mt = 1
        !$ mt = omp_get_max_threads()
        allocate(b(mt),source=huge(1.0))
        allocate(t(mt),source=0)

        !$omp parallel default(private) shared(a, b, t, n, idx)

            nt = 1
            it = 0
        !$  nt = omp_get_num_threads()
        !$  it = omp_get_thread_num()

            ! Divide array into chunks
            nc = n / nt

            lb = it*nc + 1
            ub = min((it+1)*nc, n)

            t(it+1) = minloc(a(lb:ub),dim=1) + lb - 1
            b(it+1) = a(t(it+1))
            !$omp barrier

            !$omp single
            idx = t(minloc(b(1:nt),dim=1))
            !$omp end single nowait

        !$omp end parallel

        end if

    end function

end module

program main
use, intrinsic :: iso_fortran_env, only: &
   compiler_version, compiler_options
use, intrinsic :: iso_fortran_env, only: int64
use fast_minloc, only: minloc_thread, minloc_reduce, minloc_simd
implicit none

real, allocatable :: a(:)
integer :: idx, i, n
integer(kind=int64) :: t1, t2, rate
character(len=64) :: str
character (len=*), parameter :: fmt_out = "(i10, f14.8, a30, f10.6)"
print "('compiler version: ', a)", compiler_version()
print "('compiler options: ', a)", compiler_options()

call get_command_argument(1,str)
read(str,*) n
print "('n = ',i0, /)", n
allocate(a(n))
call random_number(a)

print "(a10, a14, a30, a10)", "minloc", "minval", "method", "time"
call system_clock(t1,count_rate=rate)
idx = minloc(a,dim=1)
call system_clock(t2)
print fmt_out, idx, a(idx), "built-in", real(t2 - t1)/rate

call system_clock(t1)
idx = minloc_thread(a)
call system_clock(t2)
print fmt_out, idx, a(idx), "threaded", real(t2 - t1)/rate

call system_clock(t1)
idx = minloc_reduce(a)
call system_clock(t2)
print fmt_out, idx, a(idx), "user-defined reduction", real(t2-t1)/rate

call system_clock(t1)
idx = minloc_simd(a)
call system_clock(t2)
print fmt_out, idx, a(idx), "simd user-defined reduction", real(t2 - t1)/rate

end program
