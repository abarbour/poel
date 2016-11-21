      double precision function erf(x)
      implicit none
      double precision x
c
      double precision erfc
c
      erf=1.d0-erfc(x)
c
      return
      end