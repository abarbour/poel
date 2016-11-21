      subroutine pespectra(ipr,accuracy,nr1,nr2,wvint,ierr)
      implicit none
c
      integer ipr,nr1,nr2,ierr
      double precision accuracy
      external wvint
c
      include 'peglobal.h'
c
      integer i,ir,izr,l,ls,n,izs,nzs,lf,jf
      integer it,its,itga,itgb,ntg,nfg,nkc,nklast
      double precision f,dfg,t,dtg,fc,fcut,k1,k2,dk,dk0,am
      double precision a,b,zs,dzs,alpha,beta,dismax
      double precision fgrnabs,fgrnmax,delta
      double precision tg,chi,la0,mu0,alf0,chi0,diffus0
      double precision duzfs,durfs,dezzfs,derrfs,dettfs,dezrfs
      double precision dtltfs,dppfs,ddvzfs,ddvrfs
      double precision uzfs(2*nfmax),urfs(2*nfmax)
      double precision ezzfs(2*nfmax),errfs(2*nfmax),ettfs(2*nfmax)
      double precision ezrfs(2*nfmax),tltfs(2*nfmax)
      double precision ppfs(2*nfmax),dvzfs(2*nfmax),dvrfs(2*nfmax)
      double precision tgrn(2*nfmax),obs(4*nfmax)
      double precision maxuz(nrmax,nzrmax),maxur(nrmax,nzrmax)
      double precision maxezz(nrmax,nzrmax),maxerr(nrmax,nzrmax)
      double precision maxett(nrmax,nzrmax),maxezr(nrmax,nzrmax)
      double precision maxtlt(nrmax,nzrmax),maxpp(nrmax,nzrmax)
      double precision maxdvz(nrmax,nzrmax),maxdvr(nrmax,nzrmax)
      double complex preuz(nrmax,nzrmax),preur(nrmax,nzrmax)
      double complex preezz(nrmax,nzrmax),preerr(nrmax,nzrmax)
      double complex preett(nrmax,nzrmax),preezr(nrmax,nzrmax)
      double complex preess(nrmax,nzrmax)
      double complex pretlt(nrmax,nzrmax),prepp(nrmax,nzrmax)
      double complex predvz(nrmax,nzrmax),predvr(nrmax,nzrmax)
      double complex cs,swap,cuz,cur,cezz,cerr,cett,cezr,ctlt
      double complex cpp,cdvz,cdvr
      logical converge,again,dfgfound
c
      double precision pi2
      data pi2/6.28318530717959d0/
c
      ierr=0
      print *,' ==============================================='
      write(*,'(i3,a,F12.4,a,f12.4,a)')ipr,'. sub-profile: ',
     &        r(nr1),' -> ',r(nr2),' m'
c
c     determine wavenumber sampling rate
c
      ls=lstop
      n=nno(ls)
      la0=la(n)
      mu0=mu(n)
      chi0=dm(n)*alf(n)/qa(n)
      alf0=alf(n)
      diffus0=dm(n)*(la(n)+2.d0*mu(n))/(la(n)+2.d0*mu(n)+alf(n)*qa(n))
c
      write(*,'(a)')'  ==============================================='
      write(*,'(a)')'  Finding required sampling rate and range'
c
      delta=0.d0
      do l=1,lp-1
        delta=delta+hp(l)
      enddo
      dismax=dsqrt(delta**2+r(nr2)**2)
c
      dfg=dmin1(0.5d0/twindow,10.d0*diffus0/dismax**2)
      alpha=-dfg*dlog(0.01d0)
c
      do izr=1,nzr
        do ir=nr1,nr2
          preuz(ir,izr)=(0.d0,0.d0)
          preur(ir,izr)=(0.d0,0.d0)
          preezz(ir,izr)=(0.d0,0.d0)
          preerr(ir,izr)=(0.d0,0.d0)
          preett(ir,izr)=(0.d0,0.d0)
          preezr(ir,izr)=(0.d0,0.d0)
          pretlt(ir,izr)=(0.d0,0.d0)
          prepp(ir,izr)=(0.d0,0.d0)
          predvz(ir,izr)=(0.d0,0.d0)
          predvr(ir,izr)=(0.d0,0.d0)
        enddo
      enddo
c
      cs=dcmplx(alpha,pi2*dfg)
      dk0=0.5d0*pi2/dismax
c
      nklast=0
      dk=dk0
      nkc=256
100   call pebsj(nr1,nr2,dk)
      call wvint(accuracy,nr1,nr2,1,cs,
     &           nkc,dk,nklast,converge)
      again=.false.
c
      do izr=1,nzr
        do ir=nr1,nr2
          again=again.or.
     &        cdabs(grnuz(1,ir,izr)-preuz(ir,izr)).gt.
     &        accuracy*cdabs(grnuz(1,ir,izr))
          preuz(ir,izr)=grnuz(1,ir,izr)
          again=again.or.
     &        cdabs(grnur(1,ir,izr)-preur(ir,izr)).gt.
     &        accuracy*cdabs(grnur(1,ir,izr))
          preur(ir,izr)=grnur(1,ir,izr)
          again=again.or.
     &        cdabs(grnezz(1,ir,izr)-preezz(ir,izr)).gt.
     &        accuracy*cdabs(grnezz(1,ir,izr))
          preezz(ir,izr)=grnezz(1,ir,izr)
          again=again.or.
     &        cdabs(grnerr(1,ir,izr)-preerr(ir,izr)).gt.
     &        accuracy*cdabs(grnerr(1,ir,izr))
          preerr(ir,izr)=grnerr(1,ir,izr)
          again=again.or.
     &        cdabs(grnett(1,ir,izr)-preett(ir,izr)).gt.
     &        accuracy*cdabs(grnett(1,ir,izr))
          preett(ir,izr)=grnett(1,ir,izr)
          if(lzr(izr).ne.1.or.isurfcon.ne.1)then
            again=again.or.
     &          cdabs(grnezr(1,ir,izr)-preezr(ir,izr)).gt.
     &          accuracy*cdabs(grnezr(1,ir,izr))
            preezr(ir,izr)=grnezr(1,ir,izr)
          endif
          again=again.or.
     &        cdabs(grntlt(1,ir,izr)-pretlt(ir,izr)).gt.
     &        accuracy*cdabs(grntlt(1,ir,izr))
          pretlt(ir,izr)=grntlt(1,ir,izr)
          if(lzr(izr).ne.1.or.isurfcon.ne.1)then
            again=again.or.
     &          cdabs(grnpp(1,ir,izr)-prepp(ir,izr)).gt.
     &          accuracy*cdabs(grnpp(1,ir,izr))
            prepp(ir,izr)=grnpp(1,ir,izr)
          endif
          if(lzr(izr).ne.1.or.isurfcon.ne.2)then
            again=again.or.
     &          cdabs(grndvz(1,ir,izr)-predvz(ir,izr)).gt.
     &          accuracy*cdabs(grndvz(1,ir,izr))
            predvz(ir,izr)=grndvz(1,ir,izr)
          endif
          again=again.or.
     &        cdabs(grndvr(1,ir,izr)-predvr(ir,izr)).gt.
     &        accuracy*cdabs(grndvr(1,ir,izr))
          predvr(ir,izr)=grndvr(1,ir,izr)
        enddo
      enddo
c
      if(again.and.converge)then
        dk=0.5d0*dk
	  goto 100
      endif
c
c     wavenumber sampling rate dk found
c
      call pebsj(nr1,nr2,dk)
c
c     determine cut-off frequency
c
      write(*,'(a)')'  ==============================================='
      write(*,'(a)')'  Finding necessary cut-off frequency'
c
      do izr=1,nzr
        do ir=nr1,nr2
          maxuz(ir,izr)=cdabs(grnuz(1,ir,izr))
          maxur(ir,izr)=cdabs(grnur(1,ir,izr))
          maxezz(ir,izr)=cdabs(grnezz(1,ir,izr))
          maxerr(ir,izr)=cdabs(grnerr(1,ir,izr))
          maxett(ir,izr)=cdabs(grnett(1,ir,izr))
          maxezr(ir,izr)=cdabs(grnezr(1,ir,izr))
          maxtlt(ir,izr)=cdabs(grntlt(1,ir,izr))
          maxpp(ir,izr)=cdabs(grnpp(1,ir,izr))
          maxdvz(ir,izr)=cdabs(grndvz(1,ir,izr))
          maxdvr(ir,izr)=cdabs(grndvr(1,ir,izr))
        enddo
      enddo
c
      nklast=0
      lf=2
      nfg=4
300   fcut=dble(nfg-1)*dfg
      cs=dcmplx(alpha,pi2*fcut)
      call wvint(accuracy,nr1,nr2,lf,cs,
     &           nkc,dk,nklast,converge)
      again=.false.
c
      do izr=1,nzr
        do ir=nr1,nr2
          am=cdabs(grnuz(lf,ir,izr))
          maxuz(ir,izr)=dmax1(maxuz(ir,izr),am)
          again=again.or.am.gt.accuracy*maxuz(ir,izr)
c
          am=cdabs(grnur(lf,ir,izr))
          maxur(ir,izr)=dmax1(maxur(ir,izr),am)
          again=again.or.am.gt.accuracy*maxur(ir,izr)
c
          am=cdabs(grnezz(lf,ir,izr))
          maxezz(ir,izr)=dmax1(maxezz(ir,izr),am)
          again=again.or.am.gt.accuracy*maxezz(ir,izr)
c
          am=cdabs(grnerr(lf,ir,izr))
          maxerr(ir,izr)=dmax1(maxerr(ir,izr),am)
          again=again.or.am.gt.accuracy*maxerr(ir,izr)
c
          am=cdabs(grnett(lf,ir,izr))
          maxett(ir,izr)=dmax1(maxuz(ir,izr),am)
          again=again.or.am.gt.accuracy*maxett(ir,izr)
c
          if(lzr(izr).ne.1.or.isurfcon.ne.1)then
            am=cdabs(grnezr(lf,ir,izr))
            maxezr(ir,izr)=dmax1(maxezr(ir,izr),am)
            again=again.or.am.gt.accuracy*maxezr(ir,izr)
          endif
c
          am=cdabs(grntlt(lf,ir,izr))
          maxtlt(ir,izr)=dmax1(maxtlt(ir,izr),am)
          again=again.or.am.gt.accuracy*maxtlt(ir,izr)
c
          if(lzr(izr).ne.1.or.isurfcon.ne.1)then
            am=cdabs(grnpp(lf,ir,izr))
            maxpp(ir,izr)=dmax1(maxpp(ir,izr),am)
            again=again.or.am.gt.accuracy*maxpp(ir,izr)
          endif
c
          if(lzr(izr).ne.1.or.isurfcon.ne.2)then
            am=cdabs(grndvz(lf,ir,izr))
            maxdvz(ir,izr)=dmax1(maxdvz(ir,izr),am)
            again=again.or.am.gt.accuracy*maxdvz(ir,izr)
          endif
c
          am=cdabs(grndvr(lf,ir,izr))
          maxdvr(ir,izr)=dmax1(maxdvr(ir,izr),am)
          again=again.or.am.gt.accuracy*maxdvr(ir,izr)
        enddo
      enddo
c
      if(again.and.nfg.lt.nfmax.and.fcut.lt.1.d0/dt)then
        lf=lf+1
        nfg=2*nfg
        goto 300
      endif
c
      if(again)then
        ierr=1
        nfg=min0(nfg,nfmax)
      else if(nfg.lt.nfmin)then
        ierr=0
        nfg=nfmin
      else
        ierr=0
      endif
      fcut=dble(nfg-1)*dfg
      write(*,'(a)')'  ==============================================='
      write(*,'(i6,a,E12.5,a)')nfg,' frequencies used,'
     &                             //' cut-off: ',fcut,'Hz'
      write(*,'(a)')'  ==============================================='
c
      nklast=0
      do lf=1,nfg
        f=dble(lf-1)*dfg
        cs=dcmplx(alpha,pi2*f)
        call wvint(accuracy,nr1,nr2,lf,cs,
     &             nkc,dk,nklast,converge)
      enddo
c
c     convert to time domain by FFT
c
      dtg=1.d0/(dble(2*nfg)*dfg)
      ntg=2*nfg
c
      if(istype.gt.0)then
        call pesinterp(dtg)
      endif
c
c     add Rudnicki's analytical solution
c
      do izr=1,nzr
c
        n=nno(lzr(izr))
        chi=dm(n)*alf(n)/qa(n)
c
        do ir=nr1,nr2
c
          if(pointsource)then
            nzs=1
            dzs=0.d0
            delta=fana(lzr(izr))
          else
            if(rsdis(ir,izr).le.0.d0)then
              nzs=0
            else
              nzs=min0(1000,100*(1+idnint((zsbtm-zstop)/rsdis(ir,izr))))
              dzs=(zsbtm-zstop)/dble(nzs)
              delta=fana(lzr(izr))/dble(nzs)
            endif
          endif
          do it=1,ntg
            tg=dble(it-1)*dtg
            uzfs(it)=0.d0
            urfs(it)=0.d0
            ppfs(it)=0.d0
            ezzfs(it)=0.d0
            errfs(it)=0.d0
            ettfs(it)=0.d0
            ezrfs(it)=0.d0
            tltfs(it)=0.d0
            dvzfs(it)=0.d0
            dvrfs(it)=0.d0
            if(analytic.and.(r(ir).gt.sradius.or.
     &         lzr(izr).lt.lstop.or.lzr(izr).gt.lsbtm))then
              do izs=1,nzs
                zs=zstop+(dble(izs-1)+0.5d0)*dzs
                call rudnicki(zr(izr),r(ir),tg,zs,
     &                    la0,mu0,alf0,chi0,diffus0,
     &                    duzfs,durfs,dppfs,
     &                    dezzfs,derrfs,dettfs,
     &                    dezrfs,dtltfs,ddvzfs,ddvrfs,istype)
                uzfs(it)=uzfs(it)+delta*duzfs
                urfs(it)=urfs(it)+delta*durfs
                ppfs(it)=ppfs(it)+delta*dppfs
                ezzfs(it)=ezzfs(it)+delta*dezzfs
                errfs(it)=errfs(it)+delta*derrfs
                ettfs(it)=ettfs(it)+delta*dettfs
                ezrfs(it)=ezrfs(it)+delta*dezrfs
                tltfs(it)=tltfs(it)+delta*dtltfs
                dvzfs(it)=dvzfs(it)+delta*ddvzfs
                dvrfs(it)=dvrfs(it)+delta*ddvrfs
              enddo
              if(mirror.and.isurfcon.eq.1)then
                do izs=1,nzs
                  zs=-(zstop+(dble(izs-1)+0.5d0)*dzs)
                  call rudnicki(zr(izr),r(ir),tg,zs,
     &                        la0,mu0,alf0,chi0,diffus0,
     &                        duzfs,durfs,dppfs,
     &                        dezzfs,derrfs,dettfs,
     &                        dezrfs,dtltfs,ddvzfs,ddvrfs,istype)
                  uzfs(it)=uzfs(it)-delta*duzfs
                  urfs(it)=urfs(it)-delta*durfs
                  ppfs(it)=ppfs(it)-delta*dppfs
                  ezzfs(it)=ezzfs(it)-delta*dezzfs
                  errfs(it)=errfs(it)-delta*derrfs
                  ettfs(it)=ettfs(it)-delta*dettfs
                  ezrfs(it)=ezrfs(it)-delta*dezrfs
                  tltfs(it)=tltfs(it)-delta*dtltfs
                  dvzfs(it)=dvzfs(it)-delta*ddvzfs
                  dvrfs(it)=dvrfs(it)-delta*ddvrfs
                enddo
              else if(mirror.and.isurfcon.eq.2)then
                do izs=1,nzs
                  zs=-(zstop+(dble(izs-1)+0.5d0)*dzs)
                  call rudnicki(zr(izr),r(ir),tg,zs,
     &                        la0,mu0,alf0,chi0,diffus0,
     &                        duzfs,durfs,dppfs,
     &                        dezzfs,derrfs,dettfs,
     &                        dezrfs,dtltfs,ddvzfs,ddvrfs,istype)
                  uzfs(it)=uzfs(it)+delta*duzfs
                  urfs(it)=urfs(it)+delta*durfs
                  ppfs(it)=ppfs(it)+delta*dppfs
                  ezzfs(it)=ezzfs(it)+delta*dezzfs
                  errfs(it)=errfs(it)+delta*derrfs
                  ettfs(it)=ettfs(it)+delta*dettfs
                  ezrfs(it)=ezrfs(it)+delta*dezrfs
                  tltfs(it)=tltfs(it)+delta*dtltfs
                  dvzfs(it)=dvzfs(it)+delta*ddvzfs
                  dvrfs(it)=dvrfs(it)+delta*ddvrfs
                enddo
              endif
            endif
          enddo
c
          call subfft(grnuz(1,ir,izr),tgrn(1),uzfs(1),
     &          fcut,nfg,nfmax,dfg,dtg,alpha,istype,
     &          ts(1),sinj(1),nts,obs,nt,dt,accuracy)
c
          call subfft(grnur(1,ir,izr),tgrn(1),urfs(1),
     &          fcut,nfg,nfmax,dfg,dtg,alpha,istype,
     &          ts(1),sinj(1),nts,obs,nt,dt,accuracy)
c
          call subfft(grnezz(1,ir,izr),tgrn(1),ezzfs(1),
     &          fcut,nfg,nfmax,dfg,dtg,alpha,istype,
     &          ts(1),sinj(1),nts,obs,nt,dt,accuracy)
c
          call subfft(grnerr(1,ir,izr),tgrn(1),errfs(1),
     &          fcut,nfg,nfmax,dfg,dtg,alpha,istype,
     &          ts(1),sinj(1),nts,obs,nt,dt,accuracy)
c
          call subfft(grnett(1,ir,izr),tgrn(1),ettfs(1),
     &          fcut,nfg,nfmax,dfg,dtg,alpha,istype,
     &          ts(1),sinj(1),nts,obs,nt,dt,accuracy)
c
          call subfft(grnezr(1,ir,izr),tgrn(1),ezrfs(1),
     &          fcut,nfg,nfmax,dfg,dtg,alpha,istype,
     &          ts(1),sinj(1),nts,obs,nt,dt,accuracy)
c
          call subfft(grntlt(1,ir,izr),tgrn(1),tltfs(1),
     &          fcut,nfg,nfmax,dfg,dtg,alpha,istype,
     &          ts(1),sinj(1),nts,obs,nt,dt,accuracy)
c
          call subfft(grnpp(1,ir,izr),tgrn(1),ppfs(1),
     &          fcut,nfg,nfmax,dfg,dtg,alpha,istype,
     &          ts(1),sinj(1),nts,obs,nt,dt,accuracy)
c
          call subfft(grndvz(1,ir,izr),tgrn(1),dvzfs(1),
     &          fcut,nfg,nfmax,dfg,dtg,alpha,istype,
     &          ts(1),sinj(1),nts,obs,nt,dt,accuracy)
c
          call subfft(grndvr(1,ir,izr),tgrn(1),dvrfs(1),
     &          fcut,nfg,nfmax,dfg,dtg,alpha,istype,
     &          ts(1),sinj(1),nts,obs,nt,dt,accuracy)
c
          do it=1,(1+nt)/2
            grndvz(it,ir,izr)=-dcmplx(chi,0.d0)*grndvz(it,ir,izr)
            grndvr(it,ir,izr)=-dcmplx(chi,0.d0)*grndvr(it,ir,izr)
          enddo
        enddo
      enddo
      return
      end