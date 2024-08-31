module binomial_mod
use kind_mod, only: i64, dp
implicit none
private
public :: falling_factorial, factorial, binomial_coeff, negative_binomial, nfail
contains
elemental function falling_factorial(n, k) result(j)
! returns the product n*(n-1)*(n-2)*...*(n-k+1)
integer(kind=i64), intent(in) :: n, k
integer(kind=i64) :: j
integer(kind=i64) :: i
if (n > 0) then
   if (k <= 0) then
      j = 1
   else if (k > n) then
      j = 0
   else
      j = n
      do i=n-1, n-k+1, -1
         j = j*i
      end do
   end if
else if (n == 0) then
   j = merge(1, 0, k < 1)
end if
end function falling_factorial
!
elemental function factorial(n) result(fac)
integer(kind=i64), intent(in) :: n
integer(kind=i64)             :: fac
integer(kind=i64)             :: i
if (n < 2) then
   fac = 1
   return
end if
fac = 2
do i=3, n
   fac = fac*i
end do
end function factorial
!
elemental function binomial_coeff(n, k) result(j)
integer(kind=i64), intent(in) :: n, k
integer(kind=i64) :: j
j = falling_factorial(n, k) / factorial(k)
end function binomial_coeff
!
elemental function negative_binomial(k, r, p) result(y)
! negative binomial probability density at k, given r and p
integer(kind=i64), intent(in)  :: k ! # of failures
integer(kind=i64), intent(in)  :: r ! # of successes
real(kind=dp)    , intent(in)  :: p ! probability of success
real(kind=dp)                  :: y
y = binomial_coeff(k+r-1, k) * ((1-p)**k) * (p**r)
end function negative_binomial
!
function nfail(rmax, p) result(nf)
! simulate the # of failures before 1, 2, ..., rmax successes
integer(kind=i64), intent(in) :: rmax
real(kind=dp)    , intent(in) :: p
integer(kind=i64) :: nf(rmax)
integer(kind=i64) :: i, nsuccess
real(kind=dp) :: x
if (rmax < 1) return
i = 0
nsuccess = 0        
do
   i = i+1
   call random_number(x)
   if (x < p) then
      nsuccess = nsuccess + 1
      nf(nsuccess) = i - nsuccess
      if (nsuccess == rmax) return
   end if   
end do
end function nfail
!
end module binomial_mod
