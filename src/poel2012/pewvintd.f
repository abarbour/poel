      subroutine pewvintd(accuracy,nr1,nr2,lf,cs,
     &                    nkc,dk,nklast,converge)
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
      integer i,id,ir,ik,izr,nnk
      double precision k,k0,fac,kdk,kfdk
      double complex ck,cdisk,cdk,c2dk,cdk2,ckf,ckf2
      double complex cfac
      double complex cuz(nzrmax),cur(nzrmax),cp(nzrmax),cezz(nzrmax)
      double complex cess(nzrmax),cezr(nzrmax),cdpz(nzrmax)
      double complex preuz(nrmax,nzrmax),preur(nrmax,nzrmax)
      double complex preezz(nrmax,nzrmax),preerr(nrmax,nzrmax)
      double complex preett(nrmax,nzrmax),preezr(nrmax,nzrmax)
      double complex preess(nrmax,nzrmax)
      double complex pretlt(nrmax,nzrmax),prepp(nrmax,nzrmax)
      double complex predvz(nrmax,nzrmax),predvr(nrmax,nzrmax)
      double complex uz0(nrmax,nzrmax),ur0(nrmax,nzrmax)
      double complex ezz0(nrmax,nzrmax),err0(nrmax,nzrmax)
      double complex ett0(nrmax,nzrmax),ezr0(nrmax,nzrmax)
      double complex ess0(nrmax,nzrmax)
      double complex tlt0(nrmax,nzrmax),pp0(nrmax,nzrmax)
      double complex dvz0(nrmax,nzrmax),dvr0(nrmax,nzrmax)
      double complex cy(10,nzrmax),cbs(0:1),cm2(10)
      double complex dyp(10,0:ndmax,nzrmax)
      double complex dy0(10,0:ndmax,nzrmax)
      double complex dym(10,0:ndmax,nzrmax)
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
      k0=dmax1(dble(nd)*100.d0*dk,10.d0*pi2/r(nr2))
c
      cm2(1)=(0.d0,0.d0)
      cm2(2)=(1.d0,0.d0)
      cm2(3)=(0.d0,0.d0)
      cm2(4)=(0.d0,0.d0)
      cm2(5)=(0.d0,0.d0)
      cm2(6)=(1.d0,0.d0)
      cm2(7)=(1.d0,0.d0)
      cm2(8)=(0.d0,0.d0)
      cm2(9)=(0.d0,0.d0)
      cm2(10)=(1.d0,0.d0)
      do izr=1,nzr
        do id=0,nd
          do i=1,10
            dym(i,id,izr)=(0.d0,0.d0)
            dy0(i,id,izr)=(0.d0,0.d0)
            dyp(i,id,izr)=(0.d0,0.d0)
          enddo
        enddo
        do ir=nr1,nr2
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
c
          uz0(ir,izr)=(0.d0,0.d0)
          ur0(ir,izr)=(0.d0,0.d0)
          ezz0(ir,izr)=(0.d0,0.d0)
          ess0(ir,izr)=(0.d0,0.d0)
          err0(ir,izr)=(0.d0,0.d0)
          ett0(ir,izr)=(0.d0,0.d0)
          ezr0(ir,izr)=(0.d0,0.d0)
          tlt0(ir,izr)=(0.d0,0.d0)
          pp0(ir,izr)=(0.d0,0.d0)
          dvz0(ir,izr)=(0.d0,0.d0)
          dvr0(ir,izr)=(0.d0,0.d0)
        enddo
      enddo
c
      nnk=1
      cdk=dcmplx(dk,0.d0)
      c2dk=dcmplx(2.d0*dk,0.d0)
      cdk2=dcmplx(dk*dk,0.d0)
      do ik=1,nbsjmax
        k=dble(ik)*dk
        kdk=k*dk
        ck=dcmplx(k,0.d0)
        cdisk=dcmplx(disk(ik),0.d0)
c
        call pekern(cs,k,cuz(1),cur(1),cp(1),
     &              cezz(1),cess(1),cezr(1),cdpz(1))
c
        do izr=1,nzr
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
        enddo
c
        do izr=1,nzr
          do i=1,10
            dym(i,0,izr)=dy0(i,0,izr)
            dy0(i,0,izr)=dyp(i,0,izr)
          enddo
        enddo
c
        if(ik.le.nd)then
          do izr=1,nzr
            do i=1,10
              dyp(i,0,izr)=(0.d0,0.d0)
            enddo
          enddo
        else
          cfac=dcmplx(dexp(-((k-dble(nd)*dk)/k0)**2),0.d0)
          do izr=1,nzr
            do i=1,10
              dyp(i,0,izr)=cy(i,izr)*((1.d0,0.d0)-cfac)
              cy(i,izr)=cy(i,izr)*cfac
            enddo
          enddo
        endif
c
        do ir=nr1,nr2
          cbs(0)=dcmplx( bsj(ik,0,ir)*kdk,0.d0)
          cbs(1)=dcmplx(-bsj(ik,1,ir)*kdk,0.d0)
c
          do izr=1,nzr
            uz0(ir,izr)=uz0(ir,izr)+cy(1,izr)*cbs(0)
            ur0(ir,izr)=ur0(ir,izr)+cy(2,izr)*cbs(1)
            ezz0(ir,izr)=ezz0(ir,izr)+cy(3,izr)*cbs(0)
            ess0(ir,izr)=ess0(ir,izr)+cy(4,izr)*cbs(0)
            ezr0(ir,izr)=ezr0(ir,izr)+cy(6,izr)*cbs(1)
            tlt0(ir,izr)=tlt0(ir,izr)+cy(7,izr)*cbs(1)
            pp0(ir,izr)=pp0(ir,izr)+cy(8,izr)*cbs(0)
            dvz0(ir,izr)=dvz0(ir,izr)+cy(9,izr)*cbs(0)
            dvr0(ir,izr)=dvr0(ir,izr)+cy(10,izr)*cbs(1)
          enddo
        enddo
c
        do id=1,min0(ik-1,nd)
          ckf=dcmplx(k-dble(id)*dk,0.d0)
          ckf2=ckf*ckf
          do izr=1,nzr
            do i=1,10
              dym(i,id,izr)=dy0(i,id,izr)
              dy0(i,id,izr)=dyp(i,id,izr)
              dyp(i,id,izr)=dy0(i,id-1,izr)*(zrs2(izr)+cm2(i)/ckf2)
     &         -(dyp(i,id-1,izr)-dym(i,id-1,izr))/c2dk/ckf
     &         -(dyp(i,id-1,izr)+dym(i,id-1,izr)
     &          -c2*dy0(i,id-1,izr))/cdk2
            enddo
          enddo
        enddo
c
        if(ik.gt.nd.and.ik.lt.nbsjmax-nd)then
          kfdk=(k-dble(nd)*dk)*dk
          do ir=nr1,nr2
            cbs(0)=dcmplx( bsj(ik-nd,0,ir)*kfdk,0.d0)
            cbs(1)=dcmplx(-bsj(ik-nd,1,ir)*kfdk,0.d0)
            do izr=1,nzr
c
              grnuz(lf,ir,izr)=grnuz(lf,ir,izr)
     &                        +dyp(1,nd,izr)*cbs(0)
              grnur(lf,ir,izr)=grnur(lf,ir,izr)
     &                        +dyp(2,nd,izr)*cbs(1)
              grnezz(lf,ir,izr)=grnezz(lf,ir,izr)
     &                         +dyp(3,nd,izr)*cbs(0)
              grness(lf,ir,izr)=grness(lf,ir,izr)
     &                         +dyp(4,nd,izr)*cbs(0)
              grnezr(lf,ir,izr)=grnezr(lf,ir,izr)
     &                         +dyp(6,nd,izr)*cbs(1)
              grntlt(lf,ir,izr)=grntlt(lf,ir,izr)
     &                         +dyp(7,nd,izr)*cbs(1)
              grnpp(lf,ir,izr)=grnpp(lf,ir,izr)
     &                        +dyp(8,nd,izr)*cbs(0)
              grndvz(lf,ir,izr)=grndvz(lf,ir,izr)
     &                         +dyp(9,nd,izr)*cbs(0)
              grndvr(lf,ir,izr)=grndvr(lf,ir,izr)
     &                         +dyp(10,nd,izr)*cbs(1)
            enddo
          enddo
        endif
c
        if(ik.eq.nnk*nkc)then
          nnk=nnk*2
          if(ik.lt.nklast/2)goto 100
          finish=.true.
c
          do izr=1,nzr
            do ir=nr1,nr2
                
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
            enddo
          enddo
c
          if(finish)goto 200
        endif
100     continue
      enddo
200   nklast=ik
      if(nklast.gt.nbsjmax)then
        converge=.false.
        nklast=nbsjmax
        print *,' Warning in pewvint: int. accuracy not achieved!'
        print *,' Suggested solution: change the required accuracy'
        print *,'                     or increase the source radius.'
      else
        converge=.true.
      endif
c
      do ir=nr1,nr2
        do izr=1,nzr
          cfac=dcmplx(1.d0/(zrs2(izr)+r(ir)**2)**nd,0.d0)
          grnuz(lf,ir,izr)=uz0(ir,izr)+grnuz(lf,ir,izr)*cfac
          grnur(lf,ir,izr)=ur0(ir,izr)+grnur(lf,ir,izr)*cfac
          grnezz(lf,ir,izr)=ezz0(ir,izr)+grnezz(lf,ir,izr)*cfac
          grness(lf,ir,izr)=ess0(ir,izr)+grness(lf,ir,izr)*cfac
          grnezr(lf,ir,izr)=ezr0(ir,izr)+grnezr(lf,ir,izr)*cfac
          grntlt(lf,ir,izr)=tlt0(ir,izr)+grntlt(lf,ir,izr)*cfac
          grnpp(lf,ir,izr)=pp0(ir,izr)+grnpp(lf,ir,izr)*cfac
          grndvz(lf,ir,izr)=dvz0(ir,izr)+grndvz(lf,ir,izr)*cfac
          grndvr(lf,ir,izr)=dvr0(ir,izr)+grndvr(lf,ir,izr)*cfac
        enddo
      enddo
c
      do izr=1,nzr
        do ir=nr1,nr2
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