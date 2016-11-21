      subroutine pedifmai(difmai,cs,k,n)
      implicit none
c
      integer n
      double precision k
      double complex cs
      double complex difmai(6,6)
c
      include 'peglobal.h'
c
      integer i,j
      double complex ck,ck2,cla,clau,cmu,cv
      double complex ca,cqa,cdm,et,etu,xi,xiu
c
      double complex c0,c1,c2,c3,c4
      data c0,c1,c2,c3,c4/(0.d0,0.d0),(1.d0,0.d0),(2.d0,0.d0),
     &                    (3.d0,0.d0),(4.d0,0.d0)/
c
      ck=dcmplx(k,0.d0)
      ck2=dcmplx(k*k,0.d0)
      cla=dcmplx(la(n),0.d0)
      clau=dcmplx(la(n)+alf(n)*qa(n),0.d0)
      cmu=dcmplx(mu(n),0.d0)
      ca=dcmplx(alf(n),0.d0)
      cqa=dcmplx(qa(n),0.d0)
      cdm=dcmplx(dm(n),0.d0)
      et=cla+cmu
      etu=clau+cmu
      xi=cla+c2*cmu
      xiu=clau+c2*cmu
      cv=xiu*cs+xi*cdm*ck2
c
      difmai(1,1)=c0
      difmai(1,2)=c1/(cmu*ck2)
      difmai(1,3)=-c1/ck
      difmai(1,4)=c0
      difmai(1,5)=c0
      difmai(1,6)=c0
c
      difmai(2,1)=c4*cmu*(etu*cs+et*cdm*ck2)/cv
      difmai(2,2)=c0
      difmai(2,3)=c0
      difmai(2,4)=-(clau*cs+cla*cdm*ck2)/(ck*cv)
      difmai(2,5)=c0
      difmai(2,6)=c2*cmu*cqa/cv
c
      difmai(3,1)=(clau*cs+cla*cdm*ck2)/(ck*cv)
      difmai(3,2)=c0
      difmai(3,3)=c0
      difmai(3,4)=(cs+cdm*ck2)/(ck2*cv)
      difmai(3,5)=c0
      difmai(3,6)=cqa/(ck*cv)
c
      difmai(4,1)=c0
      difmai(4,2)=c1/ck
      difmai(4,3)=c0
      difmai(4,4)=c0
      difmai(4,5)=c0
      difmai(4,6)=c0
c
      difmai(5,1)=-c2*cmu*cqa*cs/cv
      difmai(5,2)=c0
      difmai(5,3)=c0
      difmai(5,4)=cqa*cs/(ck*cv)
      difmai(5,5)=c0
      difmai(5,6)=-xi*cqa/(ca*cv)
c
      difmai(6,1)=c0
      difmai(6,2)=c0
      difmai(6,3)=c0
      difmai(6,4)=c0
      difmai(6,5)=-cdm*ca/cqa
      difmai(6,6)=c0
c
c     transform y6 -> y6+a*s*y1
c
      do j=1,6
        difmai(6,j)=difmai(6,j)+ca*cs*difmai(1,j)
      enddo
      do i=1,6
        difmai(i,1)=difmai(i,1)-ca*cs*difmai(i,6)
      enddo
c
      return
      end