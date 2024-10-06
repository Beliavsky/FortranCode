module random
! A module for random number generation from the following distributions:
!
!     Distribution                    function/subroutine name
!
!     Normal (Gaussian)               random_normal
!     Gamma                           random_gamma
!     Chi-squared                     random_chisq
!     Exponential                     random_exponential
!     Weibull                         random_Weibull
!     Beta                            random_beta
!     t                               random_t
!     Multivariate normal             random_mvnorm
!     Generalized inverse Gaussian    random_inv_gauss
!     Poisson                         random_Poisson
!     Binomial                        random_binomial1   *
!                                     random_binomial2   *
!     Negative binomial               random_neg_binomial
!     von Mises                       random_von_Mises
!     Cauchy                          random_Cauchy
!
!  Generate a random ordering of the integers 1 .. N
!                                     random_order
!     Initialize (seed) the uniform random number generator for ANY compiler
!                                     seed_random_number

!     Lognormal - see note below.

!  ** Two functions are provided for the binomial distribution.
!  if the parameter values remain constant, it is recommended that the
!  first function is used (random_binomial1).   if one or both of the
!  parameters change, use the second function (random_binomial2).

! The compilers own random number generator, subroutine RANdoM_NUMBER(r),
! is used to provide a source of uniformly distributed random numbers.

! N.B. At this stage, only one random number is generated at each call to
!      one of the functions above.

! The module uses the following functions which are included here:
! bin_prob to calculate a single binomial probability
! lngamma  to calculate the logarithm to base e of the gamma function

! Some of the code is adapted from Dagpunar's book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
!
! In most of Dagpunar's routines, there is a test to see whether the value
! of one or two floating-point parameters has changed since the last call.
! These tests have been replaced by using a logical variable FIRST.
! This should be set to .TRUE. on the first call using new values of the
! parameters, and .FALSE. if the parameter values are the same as for the
! previous call.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lognormal distribution
! if X has a lognormal distribution, then log(X) is normally distributed.
! Here the logarithm is the natural logarithm, that is to base e, sometimes
! denoted as ln.  To generate random variates from this distribution, generate
! a random deviate from the normal distribution with mean and variance equal
! to the mean and variance of the logarithms of X, then take its exponential.

! Relationship between the mean & variance of log(X) and the mean & variance
! of X, when X has a lognormal distribution.
! Let m = mean of log(X), and s^2 = variance of log(X)
! Then
! mean of X     = exp(m + 0.5s^2)
! variance of X = (mean(X))^2.[exp(s^2) - 1]

! In the reverse direction (rarely used)
! variance of log(X) = log[1 + var(X)/(mean(X))^2]
! mean of log(X)     = log(mean(X) - 0.5var(log(X))

! N.B. The above formulae relate to population parameters; they will only be
!      approximate if applied to sample values.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Version 1.13, 2 October 2000
! Changes from version 1.01
! 1. The random_order, random_Poisson & random_binomial routines have been
!    replaced with more efficient routines.
! 2. A routine, seed_random_number, has been added to seed the uniform random
!    number generator.   This requires input of the required number of seeds
!    for the particular compiler from a specified I/O unit such as a keyboard.
! 3. Made compatible with Lahey's ELF90.
! 4. Marsaglia & Tsang algorithm used for random_gamma when shape parameter > 1.
! 5. intent for array f corrected in random_mvnorm.

!     Author: Alan Miller
!     e-mail: amiller @ bigpond.net.au

implicit none
public :: dp,random_normal,random_exponential,random_laplace,lngamma,pi, &
          seed_random_number,random_order,random_cauchy,random_binomial1, &
          random_binomial2,random_neg_binomial,random_von_mises,random_poisson, &
          random_inv_gauss,random_weibull,random_beta,random_mvnorm,random_chisq, &
          random_t,random_normal_vec_inline,random_beta_vec,random_shuffle, &
          random_seed_init
real, private      :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0,   &
                      vsmall = tiny(1.0), vlarge = huge(1.0)
private            :: integral
integer, parameter :: dp = selected_real_kind(12, 60)
real(kind=dp), parameter :: pi = 3.141592653589793D0
interface random_laplace
   module procedure random_laplace_scalar,random_laplace_vec
end interface
interface random_normal
   module procedure random_normal_scalar,random_normal_vec,random_normal_mat
end interface
contains
!
function random_laplace_scalar() result(yy)
real(kind=dp) :: yy
real(kind=dp) :: xran
call random_number(xran)
yy = random_exponential()*merge(-1,1,xran>0.5_dp)
end function random_laplace_scalar
!
function random_laplace_vec(n,standard) result(yy)
integer, intent(in) :: n
logical, intent(in), optional :: standard
real(kind=dp) :: yy(n)
real(kind=dp) :: xran
integer :: i
do i=1,n
   call random_number(xran)
   yy(i) = random_exponential()*merge(-1,1,xran>0.5_dp)
end do
if (present(standard)) then
   if (standard) yy = yy / sqrt(2.0_dp)
end if
end function random_laplace_vec
!
function random_normal_mat(n1, n2) result(mat)
! return an n1-by-n2 matrix of random normal variates
integer, intent(in) :: n1, n2
real(kind=dp)       :: mat(n1, n2)
integer             :: i1, i2
do i2=1,n2
   do i1=1,n1
      mat(i1,i2) = random_normal_scalar()
   end do
end do
end function random_normal_mat
!
function random_normal_vec(n) result(vec)
! return n random normal variates
integer, intent(in) :: n
real(kind=dp)       :: vec(n)
integer             :: i
do i=1,n
   vec(i) = random_normal_scalar()
end do
end function random_normal_vec
!
function random_normal_scalar() result(fn_val)

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

real :: fn_val

!     Local variables
real, parameter :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,    &
            r1 = 0.27597, r2 = 0.27846
real :: u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

do
  call random_number(u)
  call random_number(v)
  v = 1.7156 * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = abs(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  if (q < r1) exit
!     Reject P if outside outer ellipse
  if (q > r2) cycle
!     Reject P if outside acceptance region
  if (v**2 < -4.0*log(u)*u**2) exit
end do

!     Return ratio of P's coordinates as the normal deviate
fn_val = v/u
end function random_normal_scalar
!
function random_normal_vec_inline(n) result(fn_val)

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.
integer, intent(in) :: n
real :: fn_val(n)

!     Local variables
real, parameter :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,    &
            r1 = 0.27597, r2 = 0.27846
real :: u, v, x, y, q
integer :: i
!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

do i=1,n
do
  call random_number(u)
  call random_number(v)
  v = 1.7156 * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = abs(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  if (q < r1) exit
!     Reject P if outside outer ellipse
  if (q > r2) cycle
!     Reject P if outside acceptance region
  if (v**2 < -4.0*log(u)*u**2) exit
end do
!     Return ratio of P's coordinates as the normal deviate
fn_val(i) = v/u
end do
end function random_normal_vec_inline
!
function random_gamma(s, first) result(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

!     function GENERATES A RANdoM GAMMA VARIATE.
!     callS EITHER random_gamma1 (S > 1.0)
!     OR random_exponential (S = 1.0)
!     OR random_gamma2 (S < 1.0).

!     S = SHAPE parameter OF DISTRIBUTION (0 < real).

real, intent(IN)    :: s
logical, intent(IN) :: first
real                :: fn_val

if (s <= zero) then
  write(*, *) 'SHAPE parameter VALUE MUST BE POSITIVE'
  stop
end if

if (s > one) then
  fn_val = random_gamma1(s, first)
ELSE if (s < one) then
  fn_val = random_gamma2(s, first)
ELSE
  fn_val = random_exponential()
end if

return
end function random_gamma



function random_gamma1(s, first) result(fn_val)

! Uses the algorithm in
! Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
! gamma variables', Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.

! Generates a random gamma deviate for shape parameter s >= 1.

real, intent(IN)    :: s
logical, intent(IN) :: first
real                :: fn_val

! Local variables
real, save  :: c, d
real        :: u, v, x

if (first) then
  d = s - one/3.
  c = one/SQRT(9.0*d)
end if

! Start of main loop
do

! Generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.

  do
    x = random_normal()
    v = (one + c*x)**3
    if (v > zero) exit
  end do

! Generate uniform variable U

  call random_number(u)
  if (u < one - 0.0331*x**4) then
    fn_val = d*v
    exit
  ELSE if (log(u) < half*x**2 + d*(one - v + log(v))) then
    fn_val = d*v
    exit
  end if
end do

return
end function random_gamma1



function random_gamma2(s, first) result(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! function GENERATES A RANdoM VARIATE IN [0,INFINITY) FROM
! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO
! GAMMA2**(S-1) * exp(-GAMMA2),
! USING A SWITCHING METHOD.

!    S = SHAPE parameter OF DISTRIBUTION
!          (real < 1.0)

real, intent(IN)    :: s
logical, intent(IN) :: first
real                :: fn_val

!     Local variables
real       :: r, x, w
real, save :: a, p, c, uf, vr, d

if (s <= zero .OR. s >= one) then
  write(*, *) 'SHAPE parameter VALUE OUTSIDE PERMITTED RANGE'
  stop
end if

if (first) then                        ! Initialization, if necessary
  a = one - s
  p = a/(a + s*exp(-a))
  if (s < vsmall) then
    write(*, *) 'SHAPE parameter VALUE TOO SMALL'
    stop
  end if
  c = one/s
  uf = p*(vsmall/a)**s
  vr = one - vsmall
  d = a*log(a)
end if

do
  call random_number(r)
  if (r >= vr) then
    cycle
  ELSE if (r > p) then
    x = a - log((one - r)/(one - p))
    w = a*log(x)-d
  ELSE if (r > uf) then
    x = a*(r/p)**c
    w = x
  ELSE
    fn_val = zero
    return
  end if

  call random_number(r)
  if (one-r <= w .AND. r > zero) then
    if (r*(w + one) >= one) cycle
    if (-log(r) <= w) cycle
  end if
  exit
end do

fn_val = x
return

end function random_gamma2



function random_chisq(ndf, first) result(fn_val)

!     Generates a random variate from the chi-squared distribution with
!     ndf degrees of freedom

integer, intent(IN) :: ndf
logical, intent(IN) :: first
real                :: fn_val

fn_val = two * random_gamma(half*ndf, first)
return

end function random_chisq

function random_exponential() result(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! function GENERATES A RANdoM VARIATE IN [0,INFINITY) FROM
! A NEGATIVE expONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
! TO exp(-random_exponential), USING INVERSION.
real  :: fn_val
!     Local variable
real  :: r
do
  call random_number(r)
  if (r > zero) exit
end do
fn_val = -log(r)
return
end function random_exponential



function random_Weibull(a) result(fn_val)

!     Generates a random variate from the Weibull distribution with
!     probability density:
!                      a
!               a-1  -x
!     f(x) = a.x    e

real, intent(IN) :: a
real             :: fn_val

!     For speed, there is no checking that a is not zero or very small.

fn_val = random_exponential() ** (one/a)
return
end function random_Weibull
!
function random_beta_vec(n,aa,bb) result(fn_val)
integer      , intent(in)    :: n
real(kind=dp), intent(IN)    :: aa, bb
real(kind=dp)                :: fn_val(n)
integer                      :: i
if (n < 1) return
fn_val(1) = random_beta(aa,bb,first=.true.)
do i=2,n
   fn_val(i) = random_beta(aa,bb,first=.false.)
end do
end function random_beta_vec
!
function random_beta(aa, bb, first) result(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! function GENERATES A RANdoM VARIATE IN [0,1]
! FROM A BETA DISTRIBUTION WITH DENSITY
! PROPORTIONAL TO BETA**(AA-1) * (1-BETA)**(BB-1).
! USING CHENG'S log logISTIC METHOD.

!     AA = SHAPE parameter FROM DISTRIBUTION (0 < real)
!     BB = SHAPE parameter FROM DISTRIBUTION (0 < real)

real(kind=dp), intent(IN)    :: aa, bb
logical, intent(IN) :: first
real(kind=dp)                :: fn_val

!     Local variables
real(kind=dp), parameter  :: aln4 = 1.3862944_dp
real(kind=dp)             :: a, b, g, r, s, x, y, z
real(kind=dp), save       :: d, f, h, t, c
logical, save    :: swap

if (aa <= zero .OR. bb <= zero) then
  write(*, *) 'IMPERMISSIBLE SHAPE parameter VALUE(S)'
  stop
end if

if (first) then                        ! Initialization, if necessary
  a = aa
  b = bb
  swap = b > a
  if (swap) then
    g = b
    b = a
    a = g
  end if
  d = a/b
  f = a+b
  if (b > one) then
    h = SQRT((two*a*b - f)/(f - two))
    t = one
  ELSE
    h = b
    t = one/(one + (a/(vlarge*b))**b)
  end if
  c = a+h
end if

do
  call random_number(r)
  call random_number(x)
  s = r*r*x
  if (r < vsmall .OR. s <= zero) cycle
  if (r < t) then
    x = log(r/(one - r))/h
    y = d*exp(x)
    z = c*x + f*log((one + d)/(one + y)) - aln4
    if (s - one > z) then
      if (s - s*z > one) cycle
      if (log(s) > z) cycle
    end if
    fn_val = y/(one + y)
  ELSE
    if (4.0_dp*s > (one + one/d)**f) cycle
    fn_val = one
  end if
  exit
end do

if (swap) fn_val = one - fn_val
end function random_beta



function random_t(m) result(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! function GENERATES A RANdoM VARIATE FROM A
! T DISTRIBUTION USING kindERMAN AND MONAHAN'S RATIO METHOD.

!     M = DEGREES OF FREEdoM OF DISTRIBUTION
!           (1 <= 1NTEGER)

integer, intent(IN) :: m
real                :: fn_val

!     Local variables
real, save      :: s, c, a, f, g
real            :: r, x, v

real, parameter :: three = 3.0, four = 4.0, quart = 0.25,   &
                   five = 5.0, sixteen = 16.0
integer         :: mm = 0

if (m < 1) then
  write(*, *) 'IMPERMISSIBLE DEGREES OF FREEdoM'
  stop
end if

if (m /= mm) then                    ! Initialization, if necessary
  s = m
  c = -quart*(s + one)
  a = four/(one + one/s)**c
  f = sixteen/a
  if (m > 1) then
    g = s - one
    g = ((s + one)/g)**c*SQRT((s+s)/g)
  ELSE
    g = one
  end if
  mm = m
end if

do
  call random_number(r)
  if (r <= zero) cycle
  call random_number(v)
  x = (two*v - one)*g/r
  v = x*x
  if (v > five - a*r) then
    if (m >= 1 .AND. r*(v + three) > f) cycle
    if (r > (one + v/s)**c) cycle
  end if
  exit
end do

fn_val = x
return
end function random_t



subroutine random_mvnorm(n, h, d, f, first, x, ier)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! N.B. An extra argument, ier, has been added to Dagpunar's routine

!     subroutine GENERATES AN N VARIATE RANdoM NORMAL
!     VECTOR USING A CHOLESKY DECOMPOSITION.

! ARGUMENTS:
!        N = NUMBER OF VARIATES IN VECTOR
!           (INput,integer >= 1)
!     H(J) = J'TH ELEMENT OF VECTOR OF MEANS
!           (INput,real)
!     X(J) = J'TH ELEMENT OF DELIVERED VECTOR
!           (OUTput,real)
!
!    D(J*(J-1)/2+I) = (I,J)'TH ELEMENT OF VARIANCE MATRIX (J> = I)
!            (INput,real)
!    F((J-1)*(2*N-J)/2+I) = (I,J)'TH ELEMENT OF LOWER TRIANGULAR
!           DECOMPOSITION OF VARIANCE MATRIX (J <= I)
!            (OUTput,real)

!    FIRST = .TRUE. if THIS IS THE FIRST call OF THE ROUTINE
!    OR if THE DISTRIBUTION HAS CHANGED SINCE THE LAST call OF THE ROUTINE.
!    OTHERWISE SET TO .FALSE.
!            (INput,logical)

!    ier = 1 if the input covariance matrix is not +ve definite
!        = 0 otherwise

integer, intent(IN)   :: n
real, intent(IN)      :: h(:), d(:)   ! d(n*(n+1)/2)
real, intent(IN OUT)  :: f(:)         ! f(n*(n+1)/2)
real, intent(OUT)     :: x(:)
logical, intent(IN)   :: first
integer, intent(OUT)  :: ier

!     Local variables
integer       :: j, i, m
real          :: y, v
integer, save :: n2

if (n < 1) then
  write(*, *) 'SIZE OF VECTOR IS NON POSITIVE'
  stop
end if

ier = 0
if (first) then                        ! Initialization, if necessary
  n2 = 2*n
  if (d(1) < zero) then
    ier = 1
    return
  end if

  f(1) = SQRT(d(1))
  y = one/f(1)
  do j = 2,n
    f(j) = d(1+j*(j-1)/2) * y
  end do

  do i = 2,n
    v = d(i*(i-1)/2+i)
    do m = 1,i-1
      v = v - f((m-1)*(n2-m)/2+i)**2
    end do

    if (v < zero) then
      ier = 1
      return
    end if

    v = SQRT(v)
    y = one/v
    f((i-1)*(n2-i)/2+i) = v
    do j = i+1,n
      v = d(j*(j-1)/2+i)
      do m = 1,i-1
        v = v - f((m-1)*(n2-m)/2+i)*f((m-1)*(n2-m)/2 + j)
      end do ! m = 1,i-1
      f((i-1)*(n2-i)/2 + j) = v*y
    end do ! j = i+1,n
  end do ! i = 2,n
end if

x(1:n) = h(1:n)
do j = 1,n
  y = random_normal()
  do i = j,n
    x(i) = x(i) + f((j-1)*(n2-j)/2 + i) * y
  end do ! i = j,n
end do ! j = 1,n

return
end subroutine random_mvnorm



function random_inv_gauss(h, b, first) result(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! function GENERATES A RANdoM VARIATE IN [0,INFINITY] FROM
! A REparameterISED GENERALISED INVERSE GAUSSIAN (GIG) DISTRIBUTION
! WITH DENSITY PROPORTIONAL TO  GIG**(H-1) * exp(-0.5*B*(GIG+1/GIG))
! USING A RATIO METHOD.

!     H = parameter OF DISTRIBUTION (0 <= real)
!     B = parameter OF DISTRIBUTION (0 < real)

real, intent(IN)    :: h, b
logical, intent(IN) :: first
real                :: fn_val

!     Local variables
real            :: ym, xm, r, w, r1, r2, x
real, save      :: a, c, d, e
real, parameter :: quart = 0.25

if (h < zero .OR. b <= zero) then
  write(*, *) 'IMPERMISSIBLE DISTRIBUTION parameter VALUES'
  stop
end if

if (first) then                        ! Initialization, if necessary
  if (h > quart*b*SQRT(vlarge)) then
    write(*, *) 'THE RATIO H:B IS TOO SMALL'
    stop
  end if
  e = b*b
  d = h + one
  ym = (-d + SQRT(d*d + e))/b
  if (ym < vsmall) then
    write(*, *) 'THE VALUE OF B IS TOO SMALL'
    stop
  end if

  d = h - one
  xm = (d + SQRT(d*d + e))/b
  d = half*d
  e = -quart*b
  r = xm + one/xm
  w = xm*ym
  a = w**(-half*h) * SQRT(xm/ym) * exp(-e*(r - ym - one/ym))
  if (a < vsmall) then
    write(*, *) 'THE VALUE OF H IS TOO LARGE'
    stop
  end if
  c = -d*log(xm) - e*r
end if

do
  call random_number(r1)
  if (r1 <= zero) cycle
  call random_number(r2)
  x = a*r2/r1
  if (x <= zero) cycle
  if (log(r1) < d*log(x) + e*(x + one/x) + c) exit
end do

fn_val = x

return
end function random_inv_gauss



function random_Poisson(mu, first) result(ival)
!**********************************************************************
!     Translated to Fortran 90 by Alan Miller from:
!                           RANLIB
!
!     Library of Fortran Routines for Random Number Generation
!
!                    Compiled and Written by:
!
!                         Barry W. Brown
!                          James Lovato
!
!             Department of Biomathematics, Box 237
!             The University of Texas, M.D. Anderson Cancer Center
!             1515 Holcombe Boulevard
!             Houston, TX      77030
!
! This work was supported by grant CA-16672 from the National Cancer Institute.

!                    GENerate POIsson random deviate

!                            function

! Generates a single random deviate from a Poisson distribution with mean mu.

!                            Arguments

!     mu --> The mean of the Poisson distribution from which
!            a random deviate is to be generated.
!                              real mu

!                              Method

!     For details see:

!               Ahrens, J.H. and Dieter, U.
!               Computer Generation of Poisson Deviates
!               From Modified Normal Distributions.
!               ACM Trans. Math. Software, 8, 2
!               (June 1982),163-179

!     TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
!     COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL

!     SEPARATION OF CASES A AND B

!     .. Scalar Arguments ..
real, intent(IN)    :: mu
logical, intent(IN) :: first
integer             :: ival
!     ..
!     .. Local Scalars ..
real          :: b1, b2, c, c0, c1, c2, c3, del, difmuk, e, fk, fx, fy, g,  &
                 omega, px, py, t, u, v, x, xx
real, save    :: s, d, p, q, p0
integer       :: j, k, kflag
logical, save :: full_init
integer, save :: l, m
!     ..
!     .. Local Arrays ..
real, save    :: pp(35)
!     ..
!     .. Data statements ..
real, parameter :: a0 = -.5, a1 = .3333333, a2 = -.2500068, a3 = .2000118,  &
                   a4 = -.1661269, a5 = .1421878, a6 = -.1384794,   &
                   a7 = .1250060

real, parameter :: fact(10) = (/ 1., 1., 2., 6., 24., 120., 720., 5040.,  &
                                 40320., 362880. /)

!     ..
!     .. Executable Statements ..
if (mu > 10.0) then
!     C A S E  A. (RECALCULATION OF S, D, L if MU HAS CHANGED)

  if (first) then
    s = SQRT(mu)
    d = 6.0*mu*mu

!             THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
!             PROBABILITIES FK WHENEVER K >= M(MU). L=ifIX(MU-1.1484)
!             IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .

    l = mu - 1.1484
    full_init = .false.
  end if


!     STEP N. NORMAL SAMPLE - random_normal() FOR STANDARD NORMAL DEVIATE

  g = mu + s*random_normal()
  if (g > 0.0) then
    ival = g

!     STEP I. IMMEDIATE ACCEPTANCE if ival IS LARGE ENOUGH

    if (ival>=l) return

!     STEP S. SQUEEZE ACCEPTANCE - SAMPLE U

    fk = ival
    difmuk = mu - fk
    call random_number(u)
    if (d*u >= difmuk*difmuk*difmuk) return
  end if

!     STEP P. PREPARATIONS FOR STEPS Q AND H.
!             (RECALCULATIONS OF parameterS if NECESSARY)
!             .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
!             THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
!             APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
!             C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-function.

  if (.NOT. full_init) then
    omega = .3989423/s
    b1 = .4166667E-1/mu
    b2 = .3*b1*b1
    c3 = .1428571*b1*b2
    c2 = b2 - 15.*c3
    c1 = b1 - 6.*b2 + 45.*c3
    c0 = 1. - b1 + 3.*b2 - 15.*c3
    c = .1069/mu
    full_init = .true.
  end if

  if (g < 0.0) GO TO 50

!             'subroutine' F IS callED (KFLAG=0 FOR CORRECT return)

  kflag = 0
  GO TO 70

!     STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)

  40 if (fy-u*fy <= py*exp(px-fx)) return

!     STEP E. expONENTIAL SAMPLE - random_exponential() FOR STANDARD expONENTIAL
!             DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
!             (if T <= -.6744 then PK < FK FOR ALL MU >= 10.)

  50 e = random_exponential()
  call random_number(u)
  u = u + u - one
  t = 1.8 + SIGN(e, u)
  if (t <= (-.6744)) GO TO 50
  ival = mu + s*t
  fk = ival
  difmuk = mu - fk

!             'subroutine' F IS callED (KFLAG=1 FOR CORRECT return)

  kflag = 1
  GO TO 70

!     STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)

  60 if (c*abs(u) > py*exp(px+e) - fy*exp(fx+e)) GO TO 50
  return

!     STEP F. 'subroutine' F. CALCULATION OF PX, PY, FX, FY.
!             CASE ival < 10 USES FACTORIALS FROM TABLE FACT

  70 if (ival>=10) GO TO 80
  px = -mu
  py = mu**ival/fact(ival+1)
  GO TO 110

!             CASE ival >= 10 USES POLYNOMIAL APPROXIMATION
!             A0-A7 FOR ACCURACY WHEN ADVISABLE
!             .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)

  80 del = .8333333E-1/fk
  del = del - 4.8*del*del*del
  v = difmuk/fk
  if (abs(v)>0.25) then
    px = fk*log(one + v) - difmuk - del
  ELSE
    px = fk*v*v* (((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0) - del
  end if
  py = .3989423/SQRT(fk)
  110 x = (half - difmuk)/s
  xx = x*x
  fx = -half*xx
  fy = omega* (((c3*xx + c2)*xx + c1)*xx + c0)
  if (kflag <= 0) GO TO 40
  GO TO 60

!---------------------------------------------------------------------------
!     C A S E  B.    mu < 10
!     START NEW TABLE AND CALCULATE P0 if NECESSARY

ELSE
  if (first) then
    m = MAX(1, INT(mu))
    l = 0
    p = exp(-mu)
    q = p
    p0 = p
  end if

!     STEP U. UNifORM SAMPLE FOR INVERSION METHOD

  do
    call random_number(u)
    ival = 0
    if (u <= p0) return

!     STEP T. TABLE COMPARISON UNTIL THE end PP(L) OF THE
!             PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
!             (0.458=PP(9) FOR MU=10)

    if (l == 0) GO TO 150
    j = 1
    if (u > 0.458) j = MIN(l, m)
    do k = j, l
      if (u <= pp(k)) GO TO 180
    end do
    if (l == 35) cycle

!     STEP C. CREATION OF NEW POISSON PROBABILITIES P
!             AND THEIR CUMULATIVES Q=PP(K)

    150 l = l + 1
    do k = l, 35
      p = p*mu / k
      q = q + p
      pp(k) = q
      if (u <= q) GO TO 170
    end do
    l = 35
  end do

  170 l = k
  180 ival = k
  return
end if

return
end function random_Poisson



function random_binomial1(n, p, first) result(ival)

! function GENERATES A RANdoM BINOMIAL VARIATE USING C.D.Kemp's method.
! This algorithm is suitable when many random variates are required
! with the SAME parameter values for n & p.

!    P = BERNOULLI SUCCESS PROBABILITY
!           (0 <= real <= 1)
!    N = NUMBER OF BERNOULLI TRIALS
!           (1 <= integer)
!    FIRST = .TRUE. for the first call using the current parameter values
!          = .FALSE. if the values of (n,p) are unchanged from last call

! Reference: Kemp, C.D. (1986). `A modal method for generating binomial
!            variables', Commun. Statist. - Theor. Meth. 15(3), 805-813.

integer, intent(IN) :: n
real, intent(IN)    :: p
logical, intent(IN) :: first
integer             :: ival

!     Local variables

integer         :: ru, rd
integer, save   :: r0
real            :: u, pd, pu
real, save      :: odds_ratio, p_r
real, parameter :: zero = 0.0, one = 1.0

if (first) then
  r0 = (n+1)*p
  p_r = bin_prob(n, p, r0)
  odds_ratio = p / (one - p)
end if

call random_number(u)
u = u - p_r
if (u < zero) then
  ival = r0
  return
end if

pu = p_r
ru = r0
pd = p_r
rd = r0
do
  rd = rd - 1
  if (rd >= 0) then
    pd = pd * (rd+1) / (odds_ratio * (n-rd))
    u = u - pd
    if (u < zero) then
      ival = rd
      return
    end if
  end if

  ru = ru + 1
  if (ru <= n) then
    pu = pu * (n-ru+1) * odds_ratio / ru
    u = u - pu
    if (u < zero) then
      ival = ru
      return
    end if
  end if
end do
!     This point should not be reached, but just in case:
ival = r0
end function random_binomial1

function bin_prob(n, p, r) result(fn_val)
!     Calculate a binomial probability

integer, intent(IN) :: n, r
real, intent(IN)    :: p
real                :: fn_val

!     Local variable
real                :: one = 1.0

fn_val = exp( lngamma(DBLE(n+1)) - lngamma(DBLE(r+1)) - lngamma(DBLE(n-r+1)) &
              + r*log(p) + (n-r)*log(one - p) )
end function bin_prob

function lngamma(x) result(fn_val)

! Logarithm to base e of the gamma function.
!
! Accurate to about 1.e-14.
! Programmer: Alan Miller

! Latest revision of Fortran 77 version - 28 February 1988

real (dp), intent(IN) :: x
real (dp)             :: fn_val

!       Local variables

real (dp) :: a1 = -4.166666666554424D-02, a2 = 2.430554511376954D-03,  &
             a3 = -7.685928044064347D-04, a4 = 5.660478426014386D-04,  &
             temp, arg, product, lnrt2pi = 9.189385332046727D-1
logical   :: reflect

!       lngamma is not defined if x = 0 or a negative integer.

if (x > 0.d0) GO TO 10
if (x /= INT(x)) GO TO 10
fn_val = 0.d0
return

!       if x < 0, use the reflection formula:
!               gamma(x) * gamma(1-x) = pi * cosec(pi.x)

10 reflect = (x < 0.d0)
if (reflect) then
  arg = 1.d0 - x
ELSE
  arg = x
end if

!       Increase the argument, if necessary, to make it > 10.

product = 1.d0
20 if (arg <= 10.d0) then
  product = product * arg
  arg = arg + 1.d0
  GO TO 20
end if

!  Use a polynomial approximation to Stirling's formula.
!  N.B. The real Stirling's formula is used here, not the simpler, but less
!       accurate formula given by De Moivre in a letter to Stirling, which
!       is the one usually quoted.

arg = arg - 0.5D0
temp = 1.d0/arg**2
fn_val = lnrt2pi + arg * (log(arg) - 1.d0 + &
                  (((a4*temp + a3)*temp + a2)*temp + a1)*temp) - log(product)
if (reflect) then
  temp = SIN(pi * x)
  fn_val = log(pi/temp) - fn_val
end if
end function lngamma

function random_binomial2(n, pp, first) result(ival)
!**********************************************************************
!     Translated to Fortran 90 by Alan Miller from:
!                              RANLIB
!
!     Library of Fortran Routines for Random Number Generation
!
!                      Compiled and Written by:
!
!                           Barry W. Brown
!                            James Lovato
!
!               Department of Biomathematics, Box 237
!               The University of Texas, M.D. Anderson Cancer Center
!               1515 Holcombe Boulevard
!               Houston, TX      77030
!
! This work was supported by grant CA-16672 from the National Cancer Institute.

!                    GENerate BINomial random deviate

!                              function

!     Generates a single random deviate from a binomial
!     distribution whose number of trials is N and whose
!     probability of an event in each trial is P.

!                              Arguments

!     N  --> The number of trials in the binomial distribution
!            from which a random deviate is to be generated.
!                              integer N

!     P  --> The probability of an event in each trial of the
!            binomial distribution from which a random deviate
!            is to be generated.
!                              real P

!     FIRST --> Set FIRST = .TRUE. for the first call to perform initialization
!               the set FIRST = .FALSE. for further calls using the same pair
!               of parameter values (N, P).
!                              logical FIRST

!     random_binomial2 <-- A random deviate yielding the number of events
!                from N independent trials, each of which has
!                a probability of event P.
!                              integer random_binomial

!                              Method

!     This is algorithm BTPE from:

!         Kachitvichyanukul, V. and Schmeiser, B. W.
!         Binomial Random Variate Generation.
!         Communications of the ACM, 31, 2 (February, 1988) 216.

!**********************************************************************

!*****DETERMINE APPROPRIATE ALGORITHM AND WHETHER SETUP IS NECESSARY

!     ..
!     .. Scalar Arguments ..
real, intent(IN)    :: pp
integer, intent(IN) :: n
logical, intent(IN) :: first
integer             :: ival
!     ..
!     .. Local Scalars ..
real            :: alv, amaxp, f, f1, f2, u, v, w, w2, x, x1, x2, ynorm, z, z2
real, parameter :: zero = 0.0, half = 0.5, one = 1.0
integer         :: i, ix, ix1, k, mp
integer, save   :: m
real, save      :: p, q, xnp, ffm, fm, xnpq, p1, xm, xl, xr, c, al, xll,  &
                   xlr, p2, p3, p4, qn, r, g

!     ..
!     .. Executable Statements ..

!*****SETUP, PERFORM ONLY WHEN parameterS CHANGE

if (first) then
  p = MIN(pp, one-pp)
  q = one - p
  xnp = n * p
end if

if (xnp > 30.) then
  if (first) then
    ffm = xnp + p
    m = ffm
    fm = m
    xnpq = xnp * q
    p1 = INT(2.195*SQRT(xnpq) - 4.6*q) + half
    xm = fm + half
    xl = xm - p1
    xr = xm + p1
    c = 0.134 + 20.5 / (15.3 + fm)
    al = (ffm-xl) / (ffm - xl*p)
    xll = al * (one + half*al)
    al = (xr - ffm) / (xr*q)
    xlr = al * (one + half*al)
    p2 = p1 * (one + c + c)
    p3 = p2 + c / xll
    p4 = p3 + c / xlr
  end if

!*****GENERATE VARIATE, Binomial mean at least 30.

  20 call random_number(u)
  u = u * p4
  call random_number(v)

!     TRIANGULAR REGION

  if (u <= p1) then
    ix = xm - p1 * v + u
    GO TO 110
  end if

!     PARALLElogRAM REGION

  if (u <= p2) then
    x = xl + (u-p1) / c
    v = v * c + one - abs(xm-x) / p1
    if (v > one .OR. v <= zero) GO TO 20
    ix = x
  ELSE

!     LEFT TAIL

    if (u <= p3) then
      ix = xl + log(v) / xll
      if (ix < 0) GO TO 20
      v = v * (u-p2) * xll
    ELSE

!     RIGHT TAIL

      ix = xr - log(v) / xlr
      if (ix > n) GO TO 20
      v = v * (u-p3) * xlr
    end if
  end if

!*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST

  k = abs(ix-m)
  if (k <= 20 .OR. k >= xnpq/2-1) then

!     expLICIT EVALUATION

    f = one
    r = p / q
    g = (n+1) * r
    if (m < ix) then
      mp = m + 1
      do i = mp, ix
        f = f * (g/i-r)
      end do

    ELSE if (m > ix) then
      ix1 = ix + 1
      do i = ix1, m
        f = f / (g/i-r)
      end do
    end if

    if (v > f) then
      GO TO 20
    ELSE
      GO TO 110
    end if
  end if

!     SQUEEZING USING UPPER AND LOWER BOUNDS ON log(F(X))

  amaxp = (k/xnpq) * ((k*(k/3. + .625) + .1666666666666)/xnpq + half)
  ynorm = -k * k / (2.*xnpq)
  alv = log(v)
  if (alv<ynorm - amaxp) GO TO 110
  if (alv>ynorm + amaxp) GO TO 20

!     STIRLING'S (actually de Moivre's) FORMULA TO MACHINE ACCURACY FOR
!     THE FINAL ACCEPTANCE/REJECTION TEST

  x1 = ix + 1
  f1 = fm + one
  z = n + 1 - fm
  w = n - ix + one
  z2 = z * z
  x2 = x1 * x1
  f2 = f1 * f1
  w2 = w * w
  if (alv - (xm*log(f1/x1) + (n-m+half)*log(z/w) + (ix-m)*log(w*p/(x1*q)) +    &
      (13860.-(462.-(132.-(99.-140./f2)/f2)/f2)/f2)/f1/166320. +               &
      (13860.-(462.-(132.-(99.-140./z2)/z2)/z2)/z2)/z/166320. +                &
      (13860.-(462.-(132.-(99.-140./x2)/x2)/x2)/x2)/x1/166320. +               &
      (13860.-(462.-(132.-(99.-140./w2)/w2)/w2)/w2)/w/166320.) > zero) then
    GO TO 20
  ELSE
    GO TO 110
  end if

ELSE
!     INVERSE CDF logIC FOR MEAN LESS THAN 30
  if (first) then
    qn = q ** n
    r = p / q
    g = r * (n+1)
  end if

  90 ix = 0
  f = qn
  call random_number(u)
  100 if (u >= f) then
    if (ix > 110) GO TO 90
    u = u - f
    ix = ix + 1
    f = f * (g/ix - r)
    GO TO 100
  end if
end if

110 if (pp > half) ix = n - ix
ival = ix
end function random_binomial2

function random_neg_binomial(sk, p) result(ival)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! function GENERATES A RANdoM NEGATIVE BINOMIAL VARIATE USING UNSTORED
! INVERSION AND/OR THE REPRODUCTIVE PROPERTY.

!    SK = NUMBER OF FAILURES REQUIRED (Dagpunar's words!)
!       = the `power' parameter of the negative binomial
!           (0 < real)
!    P = BERNOULLI SUCCESS PROBABILITY
!           (0 < real < 1)

! THE parameter H IS SET SO THAT UNSTORED INVERSION ONLY IS USED WHEN P <= H,
! OTHERWISE A COMBINATION OF UNSTORED INVERSION AND
! THE REPRODUCTIVE PROPERTY IS USED.

real, intent(IN)   :: sk, p
integer            :: ival

!     Local variables
! THE parameter ULN = -log(MACHINE'S SMALLEST real NUMBER).

real, parameter    :: h = 0.7
real               :: q, x, st, uln, v, r, s, y, g
integer            :: k, i, n

if (sk <= zero .OR. p <= zero .OR. p >= one) then
  write(*, *) 'IMPERMISSIBLE DISTRIBUTION parameter VALUES'
  stop
end if

q = one - p
x = zero
st = sk
if (p > h) then
  v = one/log(p)
  k = st
  do i = 1,k
    do
      call random_number(r)
      if (r > zero) exit
    end do
    n = v*log(r)
    x = x + n
  end do
  st = st - k
end if

s = zero
uln = -log(vsmall)
if (st > -uln/log(q)) then
  write(*, *) ' P IS TOO LARGE FOR THIS VALUE OF SK'
  stop
end if

y = q**st
g = st
call random_number(r)
do
  if (y > r) exit
  r = r - y
  s = s + one
  y = y*p*g/s
  g = g + one
end do

ival = x + s + half
end function random_neg_binomial

function random_von_Mises(k, first) result(fn_val)

!     Algorithm VMD from:
!     Dagpunar, J.S. (1990) `Sampling from the von Mises distribution via a
!     comparison of random numbers', J. of Appl. Statist., 17, 165-168.

!     Fortran 90 code by Alan Miller
!     CSIRO Division of Mathematical & Information Sciences

!     Arguments:
!     k (real)        parameter of the von Mises distribution.
!     first (logical) set to .TRUE. the first time that the function
!                     is called, or the first time with a new value
!                     for k.   When first = .TRUE., the function sets
!                     up starting values and may be very much slower.

real, intent(IN)     :: k
logical, intent(IN)  :: first
real                 :: fn_val

!     Local variables

integer          :: j, n
integer, save    :: nk
real, save       :: p(20), theta(0:20)
real             :: sump, r, th, lambda, rlast
real (dp)        :: dk

if (first) then                        ! Initialization, if necessary
  if (k < zero) then
    write(*, *) '** Error: argument k for random_von_Mises = ', k
    return
  end if

  nk = k + k + one
  if (nk > 20) then
    write(*, *) '** Error: argument k for random_von_Mises = ', k
    return
  end if

  dk = k
  theta(0) = zero
  if (k > half) then

!     Set up array p of probabilities.

    sump = zero
    do j = 1, nk
      if (j < nk) then
        theta(j) = acos(one - j/k)
      ELSE
        theta(nk) = pi
      end if

!     Numerical integration of e^[k.cos(x)] from theta(j-1) to theta(j)

      call integral(theta(j-1), theta(j), p(j), dk)
      sump = sump + p(j)
    end do
    p(1:nk) = p(1:nk) / sump
  ELSE
    p(1) = one
    theta(1) = pi
  end if                         ! if k > 0.5
end if                           ! if first

call random_number(r)
do j = 1, nk
  r = r - p(j)
  if (r < zero) exit
end do
r = -r/p(j)

do
  th = theta(j-1) + r*(theta(j) - theta(j-1))
  lambda = k - j + one - k*cos(th)
  n = 1
  rlast = lambda

  do
    call random_number(r)
    if (r > rlast) exit
    n = n + 1
    rlast = r
  end do

  if (n .NE. 2*(n/2)) exit         ! is n even?
  call random_number(r)
end do

fn_val = SIGN(th, (r - rlast)/(one - rlast) - half)
end function random_von_Mises

subroutine integral(a, b, result, dk)
!     Gaussian integration of exp(k.cosx) from a to b.
real (dp), intent(IN) :: dk
real, intent(IN)      :: a, b
real, intent(OUT)     :: result
!     Local variables
real (dp)  :: xmid, range, x1, x2,                                    &
  x(3) = (/0.238619186083197_dp, 0.661209386466265_dp, 0.932469514203152_dp/), &
  w(3) = (/0.467913934572691_dp, 0.360761573048139_dp, 0.171324492379170_dp/)
integer    :: i

xmid = (a + b)/2._dp
range = (b - a)/2._dp

result = 0._dp
do i = 1, 3
  x1 = xmid + x(i)*range
  x2 = xmid - x(i)*range
  result = result + w(i)*(exp(dk*cos(x1)) + exp(dk*cos(x2)))
end do
result = result * range
end subroutine integral

function random_Cauchy() result(fn_val)
!     Generate a random deviate from the standard Cauchy distribution
real     :: fn_val
!     Local variables
real     :: v(2)
do
  call random_number(v)
  v = two*(v - half)
  if (abs(v(2)) < vsmall) cycle               ! Test for zero
  if (v(1)**2 + v(2)**2 < one) exit
end do
fn_val = v(1) / v(2)
end function random_Cauchy

subroutine random_order(order, n)
!     Generate a random ordering of the integers 1 ... n.
integer, intent(IN)  :: n
integer, intent(OUT) :: order(n)
!     Local variables
integer :: i, j, k
real    :: wk

do i = 1, n
  order(i) = i
end do

!     Starting at the end, swap the current last indicator with one
!     randomly chosen from those preceeding it.

do i = n, 2, -1
  call random_number(wk)
  j = 1 + i * wk
  if (j < i) then
    k = order(i)
    order(i) = order(j)
    order(j) = k
  end if
end do
end subroutine random_order

subroutine seed_random_number(iounit)
integer, intent(IN)  :: iounit
! Local variables
integer              :: k
integer, allocatable :: seed(:)
call random_seed(SIZE=k)
allocate( seed(k) )
write(*, '(a, i2, a)')' Enter ', k, ' integers for random no. seeds: '
read(*, *) seed
write(iounit, '(a, (7i10))') ' Random no. seeds: ', seed
call random_seed(put=seed)
deallocate( seed )
end subroutine seed_random_number

  subroutine random_shuffle(x)
!------------------------------------------------------------------------------
! Subroutine: random_shuffle
!
! Purpose:
!   Randomly shuffles the elements of the input array using the Fisher-Yates
!   algorithm. The shuffle is performed in place, ensuring each permutation
!   is equally likely.
!
! Arguments:
!   x - The array to shuffle. Must be a one-dimensional real array with
!       intent(INOUT).
!
! Example:
!     call random_shuffle(data)
!------------------------------------------------------------------------------

    real(kind=dp), intent(inout) :: x(:)
    integer :: i, j
    real(kind=dp) :: temp, rand_val

    do i = size(x), 2, -1
      call random_number(rand_val)
      j = int(rand_val * i) + 1
      temp = x(j)
      x(j) = x(i)
      x(i) = temp
    end do
  end subroutine random_shuffle

subroutine random_seed_init(iseed)
  ! Initializes the random number generator seed by setting each seed element 
  ! to a predefined value plus iseed.
  integer, intent(in) :: iseed
  integer, allocatable :: seed(:)
  integer :: seed_size, i
  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  do i = 1, seed_size
    seed(i) = int(1.0e6_dp * real(i, kind=dp)) + iseed
  end do
  call random_seed(put=seed)
end subroutine random_seed_init

end module random
