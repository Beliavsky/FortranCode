program xtrig
use kind_mod, only: dp
use trig_mod, only: pi, sin_approx, cos_approx
implicit none
integer, parameter :: n = 12
integer :: i
real(kind=dp) :: x
print "(*(a12))", "x", "sin_approx", "sin", "error", &
                       "cos_approx", "cos", "error"
do i=0, n
   x = i * 0.5_dp * pi / n
   print "(*(f12.6))", x, sin_approx(x), sin(x), sin_approx(x) - sin(x), &
      cos_approx(x), cos(x), cos_approx(x) - cos(x)
end do
end program xtrig
