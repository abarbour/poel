      subroutine pesource(am)
      implicit none
      double precision am
c
      include 'peglobal.h'
c
      integer i
      double precision pi2,a0
c
      pi2=8.d0*datan(1.d0)
      do i=1,5
        sfct(i)=(0.d0,0.d0)
      enddo
      if(pointsource)then
        sfct(6)=dcmplx(-am/pi2,0.d0)
      else
        sfct(6)=dcmplx(-am/pi2/(zsbtm-zstop),0.d0)
      endif
c
      return
      end