      subroutine pewvint(accuracy,nr1,nr2,lf,cs,
     &                   nkc,dk,nklast,converge)
      implicit none
c
      integer nr1,nr2,lf,nkc,nklast
      double precision accuracy,dk
      double complex cs
      logical converge
c
      include 'peglobal.h'
c
c     disp: 1=uz, 2=ur,
c           3=ezz, 4=err, 5=ett, 6=ezr, 7=-dur/dz
c           8=p, 9=vz, 10=vr
c
      integer i,ir,ik,izr,nra,nr0,nk
      double precision k,k0,fac
      double complex af,ck,cdisk
      double complex cfac,ckdk
      double complex cuz(nzrmax),cur(nzrmax),cp(nzrmax),cezz(nzrmax)
      double complex cess(nzrmax),cezr(nzrmax),cdpz(nzrmax)
      double complex preuz(nrmax,nzrmax),preur(nrmax,nzrmax)
      double complex preezz(nrmax,nzrmax),preerr(nrmax,nzrmax)
      double complex preett(nrmax,nzrmax),preezr(nrmax,nzrmax)
      double complex preess(nrmax,nzrmax)
      double complex pretlt(nrmax,nzrmax),prepp(nrmax,nzrmax)
      double complex predvz(nrmax,nzrmax),predvr(nrmax,nzrmax)
      double complex cy(10,nzrmax),cbs(0:1)
      double complex uk(10,nzrmax),preuk(10,nzrmax)
      logical finish,fullspace
c
c     constants
c
      double precision pi,pi2
      double complex c2
      data pi,pi2/3.14159265358979d0,6.28318530717959d0/
      data c2/(2.d0,0.d0)/
c
c     initialization
c
      do izr=1,nzr
        do ir=nr1,nr2
          grnuz(lf,ir,izr)=(0.d0,0.d0)
          grnur(lf,ir,izr)=(0.d0,0.d0)
          grnezz(lf,ir,izr)=(0.d0,0.d0)
          grness(lf,ir,izr)=(0.d0,0.d0)
          grnerr(lf,ir,izr)=(0.d0,0.d0)
          grnett(lf,ir,izr)=(0.d0,0.d0)
          grnezr(lf,ir,izr)=(0.d0,0.d0)
          grntlt(lf,ir,izr)=(0.d0,0.d0)
          grnpp(lf,ir,izr)=(0.d0,0.d0)
          grndvz(lf,ir,izr)=(0.d0,0.d0)
          grndvr(lf,ir,izr)=(0.d0,0.d0)
        enddo
      enddo
      fullspace=n0.eq.1.and.isurfcon.eq.0
      if(fullspace)return
c
      if(r(nr1).le.0.d0)then
        nra=nr1+1
      else
        nra=nr1
      endif
      do izr=1,nzr
        do i=1,10
          uk(i,izr)=(0.d0,0.d0)
          preuk(i,izr)=(0.d0,0.d0)
        enddo
        do ir=nra,nr2
          preuz(ir,izr)=(0.d0,0.d0)
          preur(ir,izr)=(0.d0,0.d0)
          preezz(ir,izr)=(0.d0,0.d0)
          preerr(ir,izr)=(0.d0,0.d0)
          preett(ir,izr)=(0.d0,0.d0)
          preess(ir,izr)=(0.d0,0.d0)
          preezr(ir,izr)=(0.d0,0.d0)
          pretlt(ir,izr)=(0.d0,0.d0)
          prepp(ir,izr)=(0.d0,0.d0)
          predvz(ir,izr)=(0.d0,0.d0)
          predvr(ir,izr)=(0.d0,0.d0)
        enddo
      enddo
c
      nk=1
      do ik=1,nbsjmax
        k=dble(ik)*dk
        ck=dcmplx(k,0.d0)
        cdisk=dcmplx(disk(ik)*k*dk,0.d0)
c
        call pekern(cs,k,cuz(1),cur(1),cp(1),
     &              cezz(1),cess(1),cezr(1),cdpz(1))
        do izr=1,nzr
          cy(1,izr)=cuz(izr)*cdisk
          cy(2,izr)=cur(izr)*cdisk
          cy(3,izr)=cezz(izr)*cdisk
          cy(4,izr)=cess(izr)*cdisk
          cy(5,izr)=(0.d0,0.d0)
          cy(6,izr)=cezr(izr)*cdisk
          cy(7,izr)=ck*cuz(izr)*cdisk
          cy(8,izr)=cp(izr)*cdisk
          cy(9,izr)=cdpz(izr)*cdisk
          cy(10,izr)=ck*cp(izr)*cdisk
c
c         frequency-wavenumber spectra:
c
c         cy(1) = uz
c         cy(2) = ur
c         cy(3) = ezz
c         cy(4) = ess
c         cy(5) = ett (will be derived later)
c         cy(6) = ezr
c         cy(7) = duz/dr
c         cy(8) = p
c         cy(9) = dp/dz
c         cy(10) = dp/dr
c
          if(r(nr1).le.0.d0.and.
     &       (lzr(izr).lt.lstop.or.lzr(izr).gt.lsbtm))then
            uk(1,izr)=uk(1,izr)+cy(1,izr)
            uk(3,izr)=uk(3,izr)+cy(3,izr)
            uk(4,izr)=uk(4,izr)+cy(4,izr)
            uk(8,izr)=uk(8,izr)+cy(8,izr)
            uk(9,izr)=uk(9,izr)+cy(9,izr)
          endif
        enddo
c
        do ir=nra,nr2
          cbs(0)=dcmplx( bsj(ik,0,ir),0.d0)
          cbs(1)=dcmplx(-bsj(ik,1,ir),0.d0)
c
          do izr=1,nzr
            if(r(ir).gt.sradius.or.
     &         lzr(izr).lt.lstop.or.lzr(izr).gt.lsbtm)then
              grnuz(lf,ir,izr)=grnuz(lf,ir,izr)+cy(1,izr)*cbs(0)
              grnur(lf,ir,izr)=grnur(lf,ir,izr)+cy(2,izr)*cbs(1)
              grnezz(lf,ir,izr)=grnezz(lf,ir,izr)+cy(3,izr)*cbs(0)
              grness(lf,ir,izr)=grness(lf,ir,izr)+cy(4,izr)*cbs(0)
              grnezr(lf,ir,izr)=grnezr(lf,ir,izr)+cy(6,izr)*cbs(1)
              grntlt(lf,ir,izr)=grntlt(lf,ir,izr)+cy(7,izr)*cbs(1)
              grnpp(lf,ir,izr)=grnpp(lf,ir,izr)+cy(8,izr)*cbs(0)
              grndvz(lf,ir,izr)=grndvz(lf,ir,izr)+cy(9,izr)*cbs(0)
              grndvr(lf,ir,izr)=grndvr(lf,ir,izr)+cy(10,izr)*cbs(1)
            endif
          enddo
        enddo
c
        if(ik.eq.nk*nkc)then
          nk=nk*2
          if(ik.lt.nklast/2)goto 100
          finish=.true.
c
          do izr=1,nzr
            if(r(nr1).le.0.d0.and.
     &         (lzr(izr).lt.lstop.or.lzr(izr).gt.lsbtm))then
              finish=finish.and.cdabs(uk(1,izr)-preuk(1,izr))
     &              .le.accuracy*cdabs(uk(1,izr))
              preuk(1,izr)=uk(1,izr)
              finish=finish.and.cdabs(uk(3,izr)-preuk(3,izr))
     &              .le.accuracy*cdabs(uk(3,izr))
              preuk(3,izr)=uk(3,izr)
              finish=finish.and.cdabs(uk(4,izr)-preuk(4,izr))
     &              .le.accuracy*cdabs(uk(4,izr))
              preuk(4,izr)=uk(4,izr)
              if(lzr(izr).ne.1.or.isurfcon.ne.1)then
                finish=finish.and.cdabs(uk(8,izr)-preuk(8,izr))
     &                .le.accuracy*cdabs(uk(8,izr))
                preuk(8,izr)=uk(8,izr)
              endif
              if(lzr(izr).ne.1.or.isurfcon.ne.2)then
                finish=finish.and.cdabs(uk(9,izr)-preuk(9,izr))
     &                .le.accuracy*cdabs(uk(9,izr))
                preuk(9,izr)=uk(9,izr)
              endif
            endif
            do ir=nra,nr2
              if(r(ir).gt.sradius.or.
     &           lzr(izr).lt.lstop.or.lzr(izr).gt.lsbtm)then
                finish=finish.and.
     &          cdabs(grnuz(lf,ir,izr)-preuz(ir,izr)).le.
     &          accuracy*cdabs(grnuz(lf,ir,izr))
                preuz(ir,izr)=grnuz(lf,ir,izr)
                finish=finish.and.
     &          cdabs(grnur(lf,ir,izr)-preur(ir,izr)).le.
     &          accuracy*cdabs(grnur(lf,ir,izr))
                preur(ir,izr)=grnur(lf,ir,izr)
                finish=finish.and.
     &          cdabs(grnezz(lf,ir,izr)-preezz(ir,izr)).le.
     &          accuracy*cdabs(grnezz(lf,ir,izr))
                preezz(ir,izr)=grnezz(lf,ir,izr)
                finish=finish.and.
     &          cdabs(grness(lf,ir,izr)-preess(ir,izr)).le.
     &          accuracy*cdabs(grness(lf,ir,izr))
                preess(ir,izr)=grness(lf,ir,izr)
                if(lzr(izr).ne.1.or.isurfcon.ne.1)then
                  finish=finish.and.
     &            cdabs(grnezr(lf,ir,izr)-preezr(ir,izr)).le.
     &            accuracy*cdabs(grnezr(lf,ir,izr))
                  preezr(ir,izr)=grnezr(lf,ir,izr)
                endif
                finish=finish.and.
     &          cdabs(grntlt(lf,ir,izr)-pretlt(ir,izr)).le.
     &          accuracy*cdabs(grntlt(lf,ir,izr))
                pretlt(ir,izr)=grntlt(lf,ir,izr)
                if(lzr(izr).ne.1.or.isurfcon.ne.1)then
                  finish=finish.and.
     &            cdabs(grnpp(lf,ir,izr)-prepp(ir,izr)).le.
     &            accuracy*cdabs(grnpp(lf,ir,izr))
                  prepp(ir,izr)=grnpp(lf,ir,izr)
                endif
                if(lzr(izr).ne.1.or.isurfcon.ne.2)then
                  finish=finish.and.
     &            cdabs(grndvz(lf,ir,izr)-predvz(ir,izr)).le.
     &            accuracy*cdabs(grndvz(lf,ir,izr))
                  predvz(ir,izr)=grndvz(lf,ir,izr)
                endif
                finish=finish.and.
     &          cdabs(grndvr(lf,ir,izr)-predvr(ir,izr)).le.
     &          accuracy*cdabs(grndvr(lf,ir,izr))
                predvr(ir,izr)=grndvr(lf,ir,izr)
              endif
            enddo
          enddo
c
          if(finish)goto 200
        endif
100     continue
      enddo
200   nklast=ik
      if(nklast.ge.nbsjmax)then
        converge=.false.
        nklast=nbsjmax
        print *,' Warning in pewvint: int. accuracy not achieved!'
        print *,' Suggested solution: change the required accuracy'
        print *,'                     or increase the source radius.'
      else
        converge=.true.
      endif
c
      do izr=1,nzr
        if(r(nr1).le.0.d0.and.
     &       (lzr(izr).lt.lstop.or.lzr(izr).gt.lsbtm))then
          grnuz(lf,nr1,izr)=uk(1,izr)
          grnezz(lf,nr1,izr)=uk(3,izr)
          grnerr(lf,nr1,izr)=uk(4,izr)/c2
          grnett(lf,nr1,izr)=uk(4,izr)/c2
          grnpp(lf,nr1,izr)=uk(8,izr)
          grndvz(lf,nr1,izr)=uk(9,izr)
        endif
        do ir=nra,nr2
c
c         ett = ur/r 
c
          grnett(lf,ir,izr)=grnur(lf,ir,izr)/dcmplx(r(ir),0.d0)
c 
c         err = ess - ett
c
          grnerr(lf,ir,izr)=grness(lf,ir,izr)-grnett(lf,ir,izr)
c
c         Tilt = -dur/dz
c
          grntlt(lf,ir,izr)=grntlt(lf,ir,izr)
     &                     -(2.d0,0.d0)*grnezr(lf,ir,izr)
        enddo
      enddo
c
c     end of total wavenumber integral
c
      write(*,'(i6,a,E12.5,a,i7)')lf,'.',dimag(cs)/pi2,
     &     'Hz: wavenumbers used = ',nklast
      return
      end