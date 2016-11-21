      subroutine pematinv(a,cs,k,n)
      implicit none
c
      include 'peglobal.h'
c
      integer n
      double precision k
      double complex cs
      double complex a(6,6)
c
      integer i,j,l,m
      double complex ck,ck2,cla,cmu,c2bk2
      double complex cv1,cv2,cv3,cv4,cv5,b(6,6),c(6,6)
      double complex ca,cb,cqa,cdm,cp,cm,cqaa,et,et1,xi,xi1
      double complex beti,delt,cfac
c
      double complex c1,c2,c3,c4
      data c1,c2,c3,c4/(1.d0,0.d0),(2.d0,0.d0),(3.d0,0.d0),
     &                 (4.d0,0.d0)/
c
      ck=dcmplx(k,0.d0)
      ck2=dcmplx(k*k,0.d0)
      cla=dcmplx(la(n),0.d0)
      cmu=dcmplx(mu(n),0.d0)
      ca=dcmplx(alf(n),0.d0)
      cqa=dcmplx(qa(n),0.d0)
      cdm=dcmplx(dm(n),0.d0)
      cqaa=cqa*ca
      et=cla+cmu
      et1=et+cqaa
      xi=cla+c2*cmu
      xi1=xi+cqaa
      cb=cdm*xi/xi1
      c2bk2=c2*cb*ck2
      cv1=cqaa*cmu*c2bk2
      cv2=c2*xi*xi1
      cv3=c2*cmu*ck*cv2
      cv4=cv2*cs
      cv5=cv3*cs
c
      a(3,1)=cmu/(c2*xi1)
      a(3,2)=c1/(c4*xi1*ck)
      a(3,3)=-a(3,1)
      a(3,4)=-a(3,2)
      a(3,5)=(0.d0,0.d0)
      a(3,6)=(0.d0,0.d0)
c
      a(4,1)=a(3,1)
      a(4,2)=-a(3,2)
      a(4,3)=a(3,1)
      a(4,4)=-a(3,2)
      a(4,5)=(0.d0,0.d0)
      a(4,6)=(0.d0,0.d0)
c
      if(smalls(n))then
c
c       a(1,j)=a(1,j)+ca*a(5,j)
c       a(2,j)=a(2,j)+ca*a(6,j)
c       a(5,j)=a(5,j)*cs
c       a(6,j)=a(6,j)*cs
c
c       delt = (ka - k)/s
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
        a(1,1)=(0.d0,0.d0)
        a(1,2)=(cv1*delt/ck+xi*cmu)/cv3
        a(1,3)=(-cv1*delt/ck+xi*et1)/cv2
        a(1,4)=xi*xi1/cv3
        a(1,5)=ca*cb*delt/(c2*xi)
        a(1,6)=(0.d0,0.d0)
c
        a(5,1)=cqa*cmu*c2bk2/cv2
        a(5,2)=cb*cqa*ka(n)/cv2
        a(5,3)=-c2*cmu*ck*a(5,2)
        a(5,4)=-cb*cqa*ck/cv2
        a(5,5)=cb*ka(n)/(c2*xi)
        a(5,6)=-cqa/(c2*ca*xi1)
      else
        a(1,1)=-cv1/cv4
        a(1,2)=(-cv1+cs*xi*cmu)/cv5
        a(1,3)=(cv1+cs*xi*et1)/cv4
        a(1,4)=(cv1+cs*xi*xi1)/cv5
        a(1,5)=-ca*cb*ck/(c2*xi*cs)
        a(1,6)=cqa/(c2*xi1*cs)
c
        a(5,1)=cqa*cmu*c2bk2/cv4
        a(5,2)=cb*cqa*ka(n)/cv4
        a(5,3)=-c2*cmu*ck*a(5,2)
        a(5,4)=-cb*cqa*ck/cv4
        a(5,5)=cb*ka(n)/(c2*xi*cs)
        a(5,6)=-cqa/(c2*ca*xi1*cs)
      endif
c
      a(2,1)=a(1,1)
      a(2,2)=-a(1,2)
      a(2,3)=-a(1,3)
      a(2,4)=a(1,4)
      a(2,5)=-a(1,5)
      a(2,6)=a(1,6)
c
      a(6,1)=a(5,1)
      a(6,2)=-a(5,2)
      a(6,3)=-a(5,3)
      a(6,4)=a(5,4)
      a(6,5)=-a(5,5)
      a(6,6)=a(5,6)
c
c     corresponding transform y6 =-xi*dp/dz to y6 = -xi*dp/dz + alpha*s*y1
c
      do i=1,6
        a(i,1)=a(i,1)-ca*cs*a(i,6)
      enddo
c
      return
      end