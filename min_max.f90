module min_max_mod
implicit none
integer, parameter :: wp = kind(1.0d0)
contains
pure function min_max(x) result(y)
real(kind=wp), intent(in) :: x(:)
real(kind=wp)             :: y(2)
integer                   :: i,n
n = size(x)
if (n < 1) return
y(1) = x(1)
y(2) = y(1)
do i=2,n
   if (x(i) < y(1)) y(1) = x(i)
   if (x(i) > y(2)) y(2) = x(i)
end do
end function min_max
!
pure function min_max_mean(x) result(y)
real(kind=wp), intent(in) :: x(:)
real(kind=wp)             :: y(3)
integer                   :: i,n
n = size(x)
if (n < 1) return
y(1)   = x(1)
y(2:3) = y(1)
do i=2,n
   if (x(i) < y(1)) y(1) = x(i)
   if (x(i) > y(2)) y(2) = x(i)
   y(3) = y(3) + x(i)
end do
if (n > 0) y(3) = y(3)/n
end function min_max_mean
end module min_max_mod
!
program xmin_max
use min_max_mod
implicit none
integer :: i,n
integer, parameter :: nt = 5, nlen = 30
real(kind=wp), allocatable :: x(:)
real(kind=wp)              :: extremes(2,2),t(nt),xmin_max_mean(3,2)
character (len=nlen) :: labels(nt-1) = [character(len=nlen) :: &
   "minval(x), maxval(x)","min_max(x)", &
   "minval(x),maxval(x),sum(x)/n","min_max_mean(x)"]
character (len=*), parameter :: fmt_cr = "(a20,*(1x,f16.12))"
n = 10**8
allocate (x(n))
call random_number(x)
call cpu_time(t(1))
extremes(:,1) = [minval(x),maxval(x)]
call cpu_time(t(2))
extremes(:,2) = min_max(x)
call cpu_time(t(3))
xmin_max_mean(:,2) = [minval(x),maxval(x),sum(x)/n]
call cpu_time(t(4))
xmin_max_mean(:,1) = min_max_mean(x)
call cpu_time(t(5))
print fmt_cr,"min, max:",extremes(:,1)
print fmt_cr,"min, max:",extremes(:,2)
print fmt_cr,"min, max, mean:",xmin_max_mean(:,1)
print fmt_cr,"min, max, mean:",xmin_max_mean(:,2)
print "(/,a8,5x,a30)","time","task"
print "(f8.4,5x,a30)",(t(i+1)-t(i),trim(labels(i)),i=1,nt-1)
end program xmin_max
