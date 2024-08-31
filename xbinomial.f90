program xbinomial
use kind_mod, only: i64
use binomial_mod
implicit none
integer, parameter :: nmax=9, kmin=1
integer(kind=i64) :: n, k
!print "(a13, *(1x,i12))", "n", (k, k=kmin,nmax+1)
print "('falling factorials', /, a7, *(1x,i6))", "n", (k, k=kmin,nmax)
do n=0,nmax
   print "(*(1x,i6))",n, (falling_factorial(n, k), k=kmin,n)
end do
print "(/,'factorials:',*(1x,i0))", (factorial(n), n=0,10)
print "(/,'binomial coefficients', /, a7, *(1x,i6))", "n", (k, k=kmin,nmax)
do n=1,nmax
   print "(*(1x,i6))",n, (binomial_coeff(n, k), k=kmin,n)
end do
end program xbinomial
