      double precision function erfc(x)
      implicit none
      double precision x
c
      double precision t
      double precision p,a1,a2,a3,a4,a5
      data p,a1,a2,a3,a4,a5/ 0.327591100d+00, 0.254829592d+00,
     &     -0.284496736d+00, 1.421413741d+00,-1.453152027d+00,
     &      1.061405429d+00/
c
      t=1.d0/(1.d0+p*x)
      erfc=t*(a1+t*(a2+t*(a3+t*(a4+a5*t))))*dexp(-x*x)
c
      return
      end