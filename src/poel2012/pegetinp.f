      subroutine pegetinp(accuracy)
      implicit none
      double precision accuracy
c
      include 'peglobal.h'
c
c     work space
c
      integer i,is,ir,izr,j,l,ls,ifs,ierr
      integer ieqdis,ieqdep
      double precision r1,r2,dr,zr1,zr2,dzr
      double precision z,pi,zrs2min,zrs2max
      double precision nu,nuu,skempton,diffus
      character comments*180
      logical sourcelayer,receiverlayer
c
      pi=4.d0*datan(1.d0)
      open(10,file=inputfile,status='old')
c
c     source parameters A. source depth & radius
c     ==========================================
c
      call getdata(10,comments)
      read(comments,*)zstop,zsbtm
      if(zstop.gt.zsbtm)then
        stop ' Error in pegetinp: source top depth > bottom depth!'
      else if(zstop.eq.zsbtm)then
        pointsource=.true.
      else
        pointsource=.false.
      endif
      call getdata(10,comments)
      read(comments,*)sradius
c
c     source parameters B. source type
c     ================================
c
      call getdata(10,comments)
      read(comments,*)istype
c
c     source parameters C. source time function
c     =========================================
c
      if(istype.le.0)then
        call getdata(10,comments)
        read(comments,*)ppini
      else
        call getdata(10,comments)
        read(comments,*)nts0
        if(nts0.gt.ntsmax)then
          stop ' Error in pegetinp: ntsmax defined too small!'
        endif
	  do i=1,nts0
          call getdata(10,comments)
          read(comments,*)j,ts0(i),sinj0(i)
          ts0(i)=ts0(i)
          sinj0(i)=sinj0(i)
        enddo
      endif
c
c     receiver parameters: depth sampling
c     ===================================
c
      call getdata(10,comments)
      read(comments,*)ieqdep
      call getdata(10,comments)
      read(comments,*)nzr
      if(nzr.gt.nzrmax)then
        stop ' Error in pegetinp: nzrmax defined too small!'
      endif
      if(ieqdep.eq.1)then
        read(10,*)zr1,zr2
        if(zr2.lt.zr1.or.zr1.lt.0.d0)then
          stop ' Error in pegetinp: zr2 < zr1 or 0!'
        endif
        if(nzr.gt.1)then
          dzr=(zr2-zr1)/dble(nzr-1)
        else
          dzr=0.d0
        endif
        do i=1,nzr
          zr(i)=zr1+dzr*dble(i-1)
        enddo
      else
        read(10,*)(zr(i),i=1,nzr)
        do i=2,nzr
          if(zr(i).lt.zr(i-1).or.zr(i).lt.0.d0)then
            stop ' Error in pegetinp: zr(i) < zr(i-1) or 0!'
          endif
        enddo
      endif
c
c     receiver parameters: distance sampling
c     ======================================
c
      call getdata(10,comments)
      read(comments,*)ieqdis
      call getdata(10,comments)
      read(comments,*)nr
      if(nr.gt.nrmax)then
        stop ' Error in pegetinp: nrmax defined too small!'
      endif
      if(ieqdis.eq.1)then
        read(10,*)r1,r2
        if(r2.lt.r1.or.r1.lt.0.d0)then
          stop ' Error in pegetinp: r2 < r1 or 0!'
        endif
        if(nr.gt.1)then
          dr=(r2-r1)/dble(nr-1)
        else
          dr=0.d0
        endif
        do i=1,nr
          r(i)=r1+dr*dble(i-1)
        enddo
      else
        read(10,*)(r(i),i=1,nr)
        do i=2,nr
          if(r(i).lt.r(i-1).or.r(i).lt.0.d0)then
            stop ' Error in pegetinp: r(i) < r(i-1) or 0!'
          endif
        enddo
      endif
c
c     receiver parameters: time sampling
c     ==================================
c
      call getdata(10,comments)
      read(comments,*)twindow
      call getdata(10,comments)
      read(comments,*)nt
      if(twindow.le.0.d0.or.nt.le.1)then
        stop ' Error in pegetinp: time window or sampling no <= 0!'
      else if(nt.gt.2*nfmax)then
        stop ' Error in pegetinp: no of time samples to large!'
      endif
      dt=twindow/dble(nt-1)
c
c     wavenumber integration parameters
c     =================================
c
      call getdata(10,comments)
      read(comments,*)accuracy
      if(accuracy.le.0.d0.or.accuracy.ge.1.d0)then
        print *,' Warning in pegetinp: unsuitable accuracy ',accuracy
        print *,' has been changed to 0.05 (5% integration accuracy)'
        accuracy=0.05
      endif
c
c     output A: displacement
c      =====================
c
      call getdata(10,comments)
      read(comments,*)seluz,selur
      call getdata(10,comments)
      read(comments,*)fileuz,fileur
c
c     output B: strain tensor & tilt
c     ==============================
c
      call getdata(10,comments)
      read(comments,*)selezz,selerr,selett,selezr,seltlt
      call getdata(10,comments)
      read(comments,*)fileezz,fileerr,fileett,fileezr,filetlt
c
c     output C: pore pressure & Darcy velocity
c     ========================================
c
      call getdata(10,comments)
      read(comments,*)selpp,seldvz,seldvr
      call getdata(10,comments)
      read(comments,*)filepp,filedvz,filedvr
c
c     output D: snapshots
c     ===================
c
      call getdata(10,comments)
      read(comments,*)nsn
      if(nsn.gt.nsnmax)then
        stop ' Error in pegetinp: too large no of snapshots!'
      endif
      do i=1,nsn
        call getdata(10,comments)
        read(comments,*)timesn(i),filesn(i)
        if(timesn(i).lt.0.d0.or.timesn(i).gt.twindow)then
          stop ' Error in pegetinp: bad snapshot time!'
        endif
      enddo
c
c     global model parameters
c     =======================
c
      call getdata(10,comments)
      read(comments,*)isurfcon
      if(isurfcon.lt.0.or.isurfcon.gt.2)then
        stop ' Error in pegetinp: isurfcon is wrong!'
      endif
      call getdata(10,comments)
      read(comments,*)l
      if(l.gt.lmax)then
        stop ' Error in pegetinp: lmax defined too small'
      endif
c
c     multilayered model parameters
c     =============================
c
      do i=1,l
        call getdata(10,comments)
        read(comments,*)j,h(i),mu(i),nu,nuu,skempton,diffus
        if(h(i).lt.0.d0.or.mu(i).le.0.d0.or.nu.le.0.d0.or.nuu.le.0.d0
     &    .or.skempton.le.0.d0.or.diffus.le.0.d0)then
          stop ' Error in pegetinp: inconsistent model parameter(s)!'
        endif
        la(i)=mu(i)*nu/(0.5d0-nu)
c
c       alf = alpha (effective stress coefficient)
c       qa = Q(biot) * alpha
c       dm = Q(biot) * kappa / eta
c
        alf(i)=1.5d0*(nuu-nu)/((0.5d0-nu)*(1+nuu)*skempton)
        qa(i)=(1.d0+nuu)*mu(i)*skempton/(0.5d0-nuu)/3.d0
        dm(i)=diffus*(1.d0-nuu)*(0.5d0-nu)/((1.d0-nu)*(0.5d0-nuu))
      enddo
      do i=2,l
        if(h(i).lt.h(i-1))then
          stop ' Error in pegetinp: inconsistent layer depth(s)!'
        endif
      enddo
c
c     end of inputs
c     =============
c
      close(10)
c
c     determine upper und lower parameter values of each layer
c
      l0=1
      z1(l0)=0.d0
      do i=2,l
        if(h(i).gt.h(i-1))then
          z1(l0)=h(i-1)
          la1(l0)=la(i-1)
          mu1(l0)=mu(i-1)
          alf1(l0)=alf(i-1)
          qa1(l0)=qa(i-1)
          dm1(l0)=dm(i-1)
c
          z2(l0)=h(i)
          la2(l0)=la(i)
          mu2(l0)=mu(i)
          alf2(l0)=alf(i)
          qa2(l0)=qa(i)
          dm2(l0)=dm(i)
          l0=l0+1
        else
          z1(l0)=h(i)
          la1(l0)=la(i)
          mu1(l0)=mu(i)
          alf1(l0)=alf(i)
          qa1(l0)=qa(i)
          dm1(l0)=dm(i)
        endif
      enddo
      z1(l0)=h(l)
      la1(l0)=la(l)
      mu1(l0)=mu(l)
      alf1(l0)=alf(l)
      qa1(l0)=qa(l)
      dm1(l0)=dm(l)
c
c     construction of sublayers
c
      ierr=0
      call pesublay(ierr)
c
      call pelayer(ierr)
c
      if(lstop.lt.lsbtm)then
        if(nno(lstop).ne.nno(lsbtm-1))then
          stop ' Error in pegetinp: line source not in same layer!'
        endif
      endif
      if(lstop.eq.1)then
        stop ' Error in pegetinp: source coincides the surface!'
      else if(nno(lstop).ne.nno(lstop-1).or.
     &        nno(lsbtm).ne.nno(lsbtm-1))then
        stop ' Error in pegetinp: source coincides material interface!'
      endif
c
c     check necessity of analytical model
c
      mirror=nno(lstop).eq.1.and.isurfcon.gt.0
c
      if(istype.le.0)then
c
c       convert initial pore pressure to equivalent impulsive fluid injection
c
        nts=1
        sinj(1)=ppini*alf(nno(lstop))/qa(nno(lstop))
      endif
c
      do izr=1,nzr
        zrs2(izr)=dmin1((zr(izr)-zstop)**2,(zr(izr)-zsbtm)**2)
        do ir=1,nr
          rsdis(ir,izr)=dsqrt(zrs2(izr)+r(ir)**2)
        enddo
      enddo
c
      rsdismax=0.d0
      do izr=1,nzr
        do ir=1,nr
          rsdismax=dmax1(rsdismax,rsdis(ir,izr))
        enddo
      enddo
c
      rsdismin=rsdismax
      do izr=1,nzr
        do ir=1,nr
          if(rsdis(ir,izr).gt.0.d0)then
            rsdismin=dmin1(rsdismin,rsdis(ir,izr))
          endif
        enddo
      enddo
      if(rsdismin.le.0.d0)then
        stop ' Error in pegetinp: zero source-receiver distance!'
      endif
c
      zrsmin=rsdismax
      do izr=1,nzr
        if(zrs2(izr).gt.0.d0)then
          zrsmin=dmin1(zrsmin,dsqrt(zrs2(izr)))
        endif
      enddo
c
      call pesource(1.d0)
c
      write(*,*)'  no  layer-no       z[m]'
      z=0.d0
      do j=1,lp
        sourcelayer=.false.
        do ls=lstop,lsbtm
          sourcelayer=sourcelayer.or.j.eq.ls
        enddo
        receiverlayer=.false.
        do izr=1,nzr
          receiverlayer=receiverlayer.or.j.eq.lzr(izr)
        enddo
c
        if(sourcelayer.and.receiverlayer)then
          write(*,1001)j,nno(j),z,
     &         ' Here is source and one of receivers level'
        else if(receiverlayer)then
          write(*,1001)j,nno(j),z,' Here is one of receivers level'
        else if(sourcelayer)then
          write(*,1001)j,nno(j),z,' Here is source level'
        else
          write(*,1001)j,nno(j),z
        endif
        z=z+hp(j)
      enddo
c
      write(*,'(a)')' Receiver depth sampling [m]:'
      l=nzr/8
      do i=1,l
        write(*,'(8f10.3)')(zr(j),j=8*(i-1)+1,8*i)
      enddo
      if(8*l.lt.nzr)then
        write(*,'(8f10.3)')(zr(j),j=8*l+1,nzr)
      endif
c
      write(*,'(a)')' Receiver distance sampling [m]:'
      l=nr/8
      do i=1,l
        write(*,'(8f10.3)')(r(j),j=8*(i-1)+1,8*i)
      enddo
      if(8*l.lt.nr)then
        write(*,'(8f10.3)')(r(j),j=8*l+1,nr)
      endif
c
      txttime='     Time[s]'
c
	txtuz='  Uz'
	txtur='  Ur'
	txtezz=' Ezz'
	txterr=' Err'
	txtett=' Ett'
	txtezr=' Ezr'
	txttlt=' Tlt'
	txtpp='  Pp'
	txtdvr=' Dvz'
	txtdvr=' Dvr'
      do j=1,nzr
        txtdep(j)(1:2)='_Z'
        i=j/1000
        txtdep(j)(3:3)=char(ichar('0')+i)
        i=mod(j,1000)/100
        txtdep(j)(4:4)=char(ichar('0')+i)
        i=mod(j,100)/10
        txtdep(j)(5:5)=char(ichar('0')+i)
        i=mod(j,10)
        txtdep(j)(6:6)=char(ichar('0')+i)
      enddo
      do j=1,nr
        txtdis(j)(1:2)='_R'
        i=j/1000
        txtdis(j)(3:3)=char(ichar('0')+i)
        i=mod(j,1000)/100
        txtdis(j)(4:4)=char(ichar('0')+i)
        i=mod(j,100)/10
        txtdis(j)(5:5)=char(ichar('0')+i)
        i=mod(j,10)
        txtdis(j)(6:6)=char(ichar('0')+i)
      enddo
c
1001  format(2i7,f12.3,a)
      return
      end