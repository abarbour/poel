      subroutine pebsj(nr1,nr2,dk)
      implicit none
      integer nr1,nr2
      double precision dk
c
      include 'peglobal.h'
c
      integer i,ir,ik
      double precision k,x
      double precision bessj0,bessj1
c
      if(r0.gt.sradius)then
        do ik=1,nbsjmax
          k=dble(ik)*dk
          x=k*r0
          disk(ik)=dexp(-0.5d0*x**2)
        enddo
      else
        do ik=1,nbsjmax
          k=dble(ik)*dk
          x=k*r0
          disk(ik)=(2.d0*bessj1(x)-x*bessj0(x))*(2.d0/x)**3
        enddo
      endif
c
      do ik=1,nbsjmax
        k=dble(ik)*dk
        do ir=nr1,nr2
          x=k*r(ir)
          bsj(ik,0,ir)=bessj0(x)
          bsj(ik,1,ir)=bessj1(x)
        enddo
      enddo
c
      return
      end