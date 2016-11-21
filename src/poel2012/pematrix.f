      subroutine pematrix(a,cs,k,z,n)
      implicit none
c
      include 'peglobal.h'
c
      integer n
      double precision k,z
      double complex cs
      double complex a(6,6)
c
      integer i,j
      double complex ck,cx,cz,cla,cmu,c2bk2
      double complex ca,cb,cqa,cdm,cp,cm,cqaa,et,etu,xi,xiu
      double complex beti,delt,cphs,cpha,cfac
c
      double precision eps
      double complex c1,c2,c3,c4
      data eps/1.0d-03/
      data c1,c2,c3,c4/(1.d0,0.d0),(2.d0,0.d0),(3.d0,0.d0),
     &                 (4.d0,0.d0)/
c
      ck=dcmplx(k,0.d0)
      cz=dcmplx(z,0.d0)
      cx=dcmplx(k*z,0.d0)
      cla=dcmplx(la(n),0.d0)
      cmu=dcmplx(mu(n),0.d0)
      ca=dcmplx(alf(n),0.d0)
      cqa=dcmplx(qa(n),0.d0)
      cdm=dcmplx(dm(n),0.d0)
      cp=c1+cx
      cm=c1-cx
      cqaa=cqa*ca
      et=cla+cmu
      etu=et+cqaa
      xi=cla+c2*cmu
      xiu=xi+cqaa
      cb=cdm*xi/xiu
      c2bk2=c2*cb*ck*ck
c
      a(1,1)=c1
      a(2,1)=c2*cmu*ck
      a(3,1)=c1
      a(4,1)=c2*cmu*ck
      a(5,1)=(0.d0,0.d0)
      a(6,1)=(0.d0,0.d0)
c
      a(1,2)=c1
      a(2,2)=-c2*cmu*ck
      a(3,2)=-c1
      a(4,2)=c2*cmu*ck
      a(5,2)=(0.d0,0.d0)
      a(6,2)=(0.d0,0.d0)
c
      a(1,3)=c1+etu*cm/cmu
      a(2,3)=c2*etu*cm*ck
      a(3,3)=-c1-etu*cx/cmu
      a(4,3)=-c2*etu*cx*ck
      a(5,3)=-c2*cqa*ck
      a(6,3)=c2*ca*cdm*ck*ck
c
      a(1,4)=c1+etu*cp/cmu
      a(2,4)=-c2*etu*cp*ck
      a(3,4)=c1-etu*cx/cmu
      a(4,4)=c2*etu*cx*ck
      a(5,4)=c2*cqa*ck
      a(6,4)=c2*ca*cdm*ck*ck
c
      a(1,5)=ca
      a(2,5)=c2*ca*cmu*ck*ck/ka(n)
      a(3,5)=ca*ck/ka(n)
      a(4,5)=c2*ca*cmu*ck
      a(5,5)=xi*cs/(cb*ka(n))
      a(6,5)=-cs*ca*xiu/cqa
c
      a(1,6)=ca
      a(2,6)=-c2*ca*cmu*ck*ck/ka(n)
      a(3,6)=-ca*ck/ka(n)
      a(4,6)=c2*ca*cmu*ck
      a(5,6)=-xi*cs/(cb*ka(n))
      a(6,6)=-cs*ca*xiu/cqa
c
      if(smalls(n))then
c
c       a(i,5)=(a(i,5)*cpha-ca*a(i,1))/cs
c       a(i,6)=(a(i,6)/cpha-ca*a(i,2))/cs
c
c       delt = (ka - k)/s
c       cphs = (exp((ka-k)*z)-1)/s
c
        beti=dcmplx((1.d0+alf(n)*qa(n)/(la(n)+2.d0*mu(n)))/dm(n),0.d0)
        delt=(0.5d0,0.d0)*beti/ck
        cfac=(0.5d0,0.d0)*beti/ck
        do i=2,10
          cfac=cfac*(beti/ck**2)
     &        *cs*dcmplx((0.5d0-dble(i-1))/dble(i),0.d0)
          delt=delt+cfac
        enddo
c
        cpha=cdexp((ka(n)-ck)*cz)
        if(cdabs((ka(n)-ck)*cz).gt.eps)then
          cphs=(cpha-(1.d0,0.d0))/cs
        else
          cphs=delt*cz
          cfac=delt*cz
          do i=2,10
            cfac=cfac*delt*cz*cs/dcmplx(dble(i),0.d0)
            cphs=cphs+cfac
          enddo
        endif
c
        a(1,5)=a(1,5)*cphs
        a(2,5)=-c2*ca*cmu*ck*delt/ka(n)+a(2,5)*cphs
        a(3,5)=-ca*delt/ka(n)+a(3,5)*cphs
        a(4,5)=a(4,5)*cphs
        a(5,5)=xi/(cb*ka(n))+a(5,5)*cphs
        a(6,5)=-ca*xiu/cqa+a(6,5)*cphs
c
        a(1,6)=-a(1,6)*cphs/cpha
        a(2,6)=c2*ca*cmu*ck*delt/ka(n)-a(2,6)*cphs/cpha
        a(3,6)=ca*delt/ka(n)-a(3,6)*cphs/cpha
        a(4,6)=-a(4,6)*cphs/cpha
        a(5,6)=-xi/(cb*ka(n))-a(5,6)*cphs/cpha
        a(6,6)=-ca*xiu/cqa-a(6,6)*cphs/cpha
      endif
c
c     convert y6 =-xi*dp/dz to y6 = -xi*dp/dz + alpha*s*y1
c
      do j=1,6
        a(6,j)=a(6,j)+ca*cs*a(1,j)
      enddo
c
      return
      end