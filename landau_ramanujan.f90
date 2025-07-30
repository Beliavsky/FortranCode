! https://www.johndcook.com/blog/2025/07/29/counting-sums-of-squares/
! "Counting sums of squares", by John D. Cook
! compute:
!   heuristic = c / sqrt(log(n_current))
!   exact     = (# of integers < n_current that equal i^2 + j^2) / n_current
! print heuristic, exact, and their difference.

program landau_ramanujan
   implicit none
   integer, parameter :: dp    = kind(1.0d0)
   real(dp), parameter :: c    = 0.76422_dp
   integer, parameter :: kmin  = 1
   integer, parameter :: kmax  = 8
   integer, parameter :: maxn  = 10 ** kmax

   integer, allocatable :: sumsq(:)
   integer :: k, n_current
   integer :: i, j, i2, jmax, idx
   integer :: total
   real(dp) :: heuristic, exact, diff

   allocate(sumsq(0:maxn-1))

   print *, " k          n       heuristic      exact     difference"
   print *, '-----------------------------------------------------------'

   do k = kmin, kmax
      n_current = 10 ** k
      sumsq(0:n_current-1) = 0          ! clear only the needed slice

      heuristic = c / sqrt(log(real(n_current, dp)))

      do i = 0, n_current - 1
         i2 = i * i
         if (i2 >= n_current) exit
         jmax = int(sqrt(real(n_current - i2 - 1, dp)))
         do j = 0, jmax
            idx = i2 + j * j
            sumsq(idx) = 1
         end do
      end do

      total = sum(sumsq(0:n_current-1))
      exact = real(total, dp) / real(n_current, dp)
      diff  = exact - heuristic

      print "(i3,1x,i10,*(1x,f13.6))", k, n_current, heuristic, exact, diff
   end do

   deallocate(sumsq)
end program landau_ramanujan
! output:
!   k          n       heuristic      exact     difference
!  -----------------------------------------------------------
!   1         10      0.503629      0.700000      0.196371
!   2        100      0.356119      0.430000      0.073881
!   3       1000      0.290770      0.330000      0.039230
!   4      10000      0.251814      0.274900      0.023086
!   5     100000      0.225230      0.240280      0.015050
!   6    1000000      0.205606      0.216341      0.010735
!   7   10000000      0.190354      0.198546      0.008192
!   8  100000000      0.178060      0.184578      0.006519

