      subroutine peoutput(ierr)
      implicit none
c
      integer ierr
c
      include 'peglobal.h'
c
      integer isn,i,it,ir,izr
      double precision a,b,uz,ur,ezz,err,ett,ezr,tlt,pp,dvz,dvr
c
      if(lzr(1).eq.1.and.isurfcon.eq.1)then
c
c       free surface
c
        do ir=1,nr
          do it=1,(1+nt)/2
            grnpp(it,ir,1)=(0.d0,0.d0)
            grnezr(it,ir,1)=(0.d0,0.d0)
            grnezz(it,ir,1)=(grnerr(it,ir,1)+grnett(it,ir,1))
     &                     *dcmplx(-la(1)/(la(1)+2.d0*mu(1)),0.d0)
          enddo
        enddo
      else if(lzr(1).eq.1.and.isurfcon.eq.2)then
c
c       impermeable surface
c
        do ir=1,nr
          do it=1,(1+nt)/2
            grndvz(it,ir,1)=(0.d0,0.d0)
            grnezr(it,ir,1)=(0.d0,0.d0)
            grnezz(it,ir,1)=(grnerr(it,ir,1)+grnett(it,ir,1))
     &                     *dcmplx(-la(1)/(la(1)+2.d0*mu(1)),0.d0)
     &                     +grnpp(it,ir,1)
     &                     *dcmplx(alf(1)/(la(1)+2.d0*mu(1)),0.d0)
          enddo
        enddo
      endif
c
c     output time series
c
      if(seluz.eq.1)then
        call tsoutput(fileuz,txttime,txtuz,txtdis(1),txtdep(1),
     &       grnuz(1,1,1),nfmax,nrmax,nzrmax,nt,nr,nzr,dt)
      endif
c
      if(selur.eq.1)then
        call tsoutput(fileur,txttime,txtur,txtdis(1),txtdep(1),
     &       grnur(1,1,1),nfmax,nrmax,nzrmax,nt,nr,nzr,dt)
      endif
c
      if(selezz.eq.1)then
        call tsoutput(fileezz,txttime,txtezz,txtdis(1),txtdep(1),
     &       grnezz(1,1,1),nfmax,nrmax,nzrmax,nt,nr,nzr,dt)
      endif
c
      if(selerr.eq.1)then
        call tsoutput(fileerr,txttime,txterr,txtdis(1),txtdep(1),
     &       grnerr(1,1,1),nfmax,nrmax,nzrmax,nt,nr,nzr,dt)
      endif
c
      if(selett.eq.1)then
        call tsoutput(fileett,txttime,txtett,txtdis(1),txtdep(1),
     &       grnett(1,1,1),nfmax,nrmax,nzrmax,nt,nr,nzr,dt)
      endif
c
      if(selezr.eq.1)then
        call tsoutput(fileezr,txttime,txtezr,txtdis(1),txtdep(1),
     &       grnezr(1,1,1),nfmax,nrmax,nzrmax,nt,nr,nzr,dt)
      endif
c
      if(seltlt.eq.1)then
        call tsoutput(filetlt,txttime,txttlt,txtdis(1),txtdep(1),
     &       grntlt(1,1,1),nfmax,nrmax,nzrmax,nt,nr,nzr,dt)
      endif
c
      if(selpp.eq.1)then
        call tsoutput(filepp,txttime,txtpp,txtdis(1),txtdep(1),
     &       grnpp(1,1,1),nfmax,nrmax,nzrmax,nt,nr,nzr,dt)
      endif
c
      if(seldvz.eq.1)then
        call tsoutput(filedvz,txttime,txtdvz,txtdis(1),txtdep(1),
     &       grndvz(1,1,1),nfmax,nrmax,nzrmax,nt,nr,nzr,dt)
      endif
c
      if(seldvr.eq.1)then
        call tsoutput(filedvr,txttime,txtdvr,txtdis(1),txtdep(1),
     &       grndvr(1,1,1),nfmax,nrmax,nzrmax,nt,nr,nzr,dt)
      endif
c
c     output snapshots
c
      do isn=1,nsn
        open(30,file=filesn(isn),status='unknown')
        write(30,'(a24,$)')'    Depth[m] Distance[m]'
        write(30,'(a16,$)')'              Uz'
        write(30,'(a16,$)')'              Ur'
        write(30,'(a16,$)')'             Ezz'
        write(30,'(a16,$)')'             Err'
        write(30,'(a16,$)')'             Ett'
        write(30,'(a16,$)')'             Ezr'
        write(30,'(a16,$)')'             Tlt'
        write(30,'(a16,$)')'              Pp'
        write(30,'(a16,$)')'             Dvz'
        write(30,'(a16)')  '             Dvr'
c
        it=min0(nt,1+idint(timesn(isn)/dt))
        b=dmod(timesn(isn)/dt,1.d0)
        a=1.d0-b
        do izr=1,nzr
          do ir=1,nr
            write(30,'(2f12.3,$)')zr(izr),r(ir)
            if(it.eq.1)then
              uz=dreal(grnuz(1,ir,izr))
              ur=dreal(grnur(1,ir,izr))
              ezz=dreal(grnezz(1,ir,izr))
              err=dreal(grnerr(1,ir,izr))
              ett=dreal(grnett(1,ir,izr))
              ezr=dreal(grnezr(1,ir,izr))
              tlt=dreal(grntlt(1,ir,izr))
              pp=dreal(grnpp(1,ir,izr))
              dvz=dreal(grndvz(1,ir,izr))
              dvr=dreal(grndvr(1,ir,izr))
            else if(it.eq.nt)then
              if(it.eq.2*(it/2))then
                i=it/2
                uz=dimag(grnuz(i,ir,izr))
                ur=dimag(grnur(i,ir,izr))
                ezz=dimag(grnezz(i,ir,izr))
                err=dimag(grnerr(i,ir,izr))
                ett=dimag(grnett(i,ir,izr))
                ezr=dimag(grnezr(i,ir,izr))
                tlt=dimag(grntlt(i,ir,izr))
                pp=dimag(grnpp(i,ir,izr))
                dvz=dimag(grndvz(i,ir,izr))
                dvr=dimag(grndvr(i,ir,izr))
              else
                i=1+it/2
                uz=dreal(grnuz(i,ir,izr))
                ur=dreal(grnur(i,ir,izr))
                ezz=dreal(grnezz(i,ir,izr))
                err=dreal(grnerr(i,ir,izr))
                ett=dreal(grnett(i,ir,izr))
                ezr=dreal(grnezr(i,ir,izr))
                tlt=dreal(grntlt(i,ir,izr))
                pp=dreal(grnpp(i,ir,izr))
                dvz=dreal(grndvz(i,ir,izr))
                dvr=dreal(grndvr(i,ir,izr))
              endif
            else if(it.gt.2*(it/2))then
              i=(1+it)/2
              uz=a*dreal(grnuz(i,ir,izr))+b*dimag(grnuz(i,ir,izr))
              ur=a*dreal(grnur(i,ir,izr))+b*dimag(grnur(i,ir,izr))
              ezz=a*dreal(grnezz(i,ir,izr))+b*dimag(grnezz(i,ir,izr))
              err=a*dreal(grnerr(i,ir,izr))+b*dimag(grnerr(i,ir,izr))
              ett=a*dreal(grnett(i,ir,izr))+b*dimag(grnett(i,ir,izr))
              ezr=a*dreal(grnezr(i,ir,izr))+b*dimag(grnezr(i,ir,izr))
              tlt=a*dreal(grntlt(i,ir,izr))+b*dimag(grntlt(i,ir,izr))
              pp=a*dreal(grnpp(i,ir,izr))+b*dimag(grnpp(i,ir,izr))
              dvz=a*dreal(grndvz(i,ir,izr))+b*dimag(grndvz(i,ir,izr))
              dvr=a*dreal(grndvr(i,ir,izr))+b*dimag(grndvr(i,ir,izr))
            else
              i=it/2
              uz=a*dimag(grnuz(i,ir,izr))+b*dreal(grnuz(i+1,ir,izr))
              ur=a*dimag(grnur(i,ir,izr))+b*dreal(grnur(i+1,ir,izr))
              ezz=a*dimag(grnezz(i,ir,izr))+b*dreal(grnezz(i+1,ir,izr))
              err=a*dimag(grnerr(i,ir,izr))+b*dreal(grnerr(i+1,ir,izr))
              ett=a*dimag(grnett(i,ir,izr))+b*dreal(grnett(i+1,ir,izr))
              ezr=a*dimag(grnezr(i,ir,izr))+b*dreal(grnezr(i+1,ir,izr))
              tlt=a*dimag(grntlt(i,ir,izr))+b*dreal(grntlt(i+1,ir,izr))
              pp=a*dimag(grnpp(i,ir,izr))+b*dreal(grnpp(i+1,ir,izr))
              dvz=a*dimag(grndvz(i,ir,izr))+b*dreal(grndvz(i+1,ir,izr))
              dvr=a*dimag(grndvr(i,ir,izr))+b*dreal(grndvr(i+1,ir,izr))
            endif
            write(30,'(10E16.8)')uz,ur,ezz,err,ett,ezr,tlt,pp,dvz,dvr
          enddo
        enddo
100     close(30)
      enddo
      return
      end