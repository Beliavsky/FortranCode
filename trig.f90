module trig_mod
! https://www.johndcook.com/blog/2024/08/31/sine-approx/
! Ancient accurate approximation for sine, by John D. Cook
use kind_mod, only: dp
implicit none
private
public :: pi, sin_approx, cos_approx
real(kind=dp), parameter :: pi = 4*atan(1.0_dp), pisq = pi**2
contains
elemental function sin_approx(x) result(y)
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
y = 16*x*(pi - x) / (5*pisq - 4*x*(pi - x))
end function sin_approx
!
elemental function cos_approx(x) result(y)
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
y = (pisq - 4*x**2) / (pisq + x**2)
end function cos_approx
end module trig_mod
