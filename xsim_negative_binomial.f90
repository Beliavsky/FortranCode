program xsim_negative_binomial
use binomial_mod, only: nfail, negative_binomial
use kind_mod, only: i64, dp
integer(kind=i64), parameter :: k=3, nsim=10**6, max_fail=5
integer(kind=i64)            :: isim, ifail, isucc, res(k), counts(0:max_fail, k)
real(kind=dp), parameter     :: p=0.6_dp
real(kind=dp)                :: prob_nb(0:max_fail, k), cumul_prob(0:max_fail)
print "('p =',f7.4)", p
print "('#sim = ',i0)", nsim
counts = 0
do isim=1,nsim
   res = nfail(k, p)
   do isucc=1,k
      ifail = res(isucc)
      if (ifail <= max_fail) counts(ifail, isucc) = counts(ifail, isucc) + 1
   end do
end do
print "(/,'simulated probability',/,a10,*(1x,i9))", "k/r", (ifail, ifail=0, max_fail)
do isucc=1,k
   print "(i10,*(1x,f9.6))", isucc,counts(:,isucc)/real(nsim, kind=dp)
end do
print "(/,'theoretical probability',/,a10,*(1x,i9))", "k/r", (ifail, ifail=0, max_fail)
do isucc=1, k
   prob_nb(:, isucc) = [(negative_binomial(ifail,isucc,p), ifail=0, max_fail)]
   print "(i10,*(1x,f9.6))", isucc,prob_nb(:, isucc)
end do
print "(/,'theoretical cumul probability',/,a10,*(1x,i9))", "k/r", (ifail, ifail=0, max_fail)
do isucc=1, k
   cumul_prob(0) = prob_nb(0, isucc)
   do ifail=1,max_fail
      cumul_prob(ifail) = cumul_prob(ifail-1) + prob_nb(ifail, isucc)
   end do
   print "(i10,*(1x,f9.6))", isucc,cumul_prob
end do
end program xsim_negative_binomial
