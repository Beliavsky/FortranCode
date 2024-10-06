program xxrandom
! test random_normal
use kind_mod, only: dp
use random, only: random_normal
implicit none
integer :: i, iter
integer, parameter :: n = 10**7, niter=5
real(kind=dp) :: x(n)
call random_seed()
do iter=1, niter
   x = random_normal(n)
   print "(a15, *(f10.6))", "moments 1 to 4:", [(sum(x**i), i=1, 4)]/n
end do
end program xxrandom
