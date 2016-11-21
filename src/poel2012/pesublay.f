      subroutine pesublay(ierr)
      implicit none
c
      integer ierr
c
      include 'peglobal.h'
c
c     work space
c
      integer i,i0,l
      double precision dh,dla,dmu,dalf,dqa,ddm,z,dz
c
      n0=0
c
      do l=1,l0-1
        dz=z2(l)-z1(l)
        dla=2.d0*dabs(la2(l)-la1(l))/(la2(l)+la1(l))
        dmu=2.d0*dabs(mu2(l)-mu1(l))/(mu2(l)+mu1(l))
        dalf=2.d0*dabs(alf2(l)-alf1(l))/(alf2(l)+alf1(l))
        dqa=2.d0*dabs(qa2(l)-qa1(l))/(qa2(l)+qa1(l))
        ddm=2.d0*dabs(dm2(l)-dm1(l))/(dm2(l)+dm1(l))
        i0=idnint(dmax1(1.d0,dla/reslm,dmu/reslm,dqa/reslm,
     &                         dalf/resla,ddm/resld))
        dla=(la2(l)-la1(l))/dz
        dmu=(mu2(l)-mu1(l))/dz
        dalf=(alf2(l)-alf1(l))/dz
        dqa=(qa2(l)-qa1(l))/dz
        ddm=(dm2(l)-dm1(l))/dz
        dh=dz/dble(i0)
        do i=1,i0
          n0=n0+1
          if(n0.gt.lmax)then
            print *,' Error in pesublay: lmax defined too small!'
            stop
          endif
          h(n0)=dh
          z=(dble(i)-0.5d0)*dh
          la(n0)=la1(l)+dla*z
          mu(n0)=mu1(l)+dmu*z
          alf(n0)=alf1(l)+dalf*z
          qa(n0)=qa1(l)+dqa*z
          dm(n0)=dm1(l)+ddm*z
        enddo
      enddo
c
c     last layer is halfspace
c
      n0=n0+1
      if(n0.gt.lmax)then
        print *,' Error in pesublay: lmax defined too small!'
        stop
      endif
      h(n0)=0.d0
      la(n0)=la1(l0)
      mu(n0)=mu1(l0)
      alf(n0)=alf1(l0)
      qa(n0)=qa1(l0)
      dm(n0)=dm1(l0)
c
      write(*,'(7a)')'  no',' thick(m)    ','  la(Pa)    ',
     &    '  mu(Pa)    ','  alpha    ','   Qa(Pa)  ',
     &    '   dm(m^2/s)  '
      do i=1,n0
        write(*,1001)i,h(i),la(i),mu(i),
     &               alf(i),qa(i),dm(i)
      enddo
1001  format(i4,f11.4,5E12.4)
      ierr=0
      return
      end