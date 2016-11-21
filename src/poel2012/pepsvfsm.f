      subroutine pepsvfsm(y,cs,k,n)
      implicit none
c
c     calculation of response to p-sv source
c     y(6): solution vector (complex)
c     k: wave number
c     cs: complex Laplace variable
c     n: model layer whose parameter is adopted for the full space
c
      include 'peglobal.h'
c
      integer n
      double precision k
      double complex cs
      double complex y(6,nzmax)
c
c     work space
c
      integer i,j,l,key
      double complex ck
      double complex yup(6,3,-1:nzmax),ylw(6,3,-1:nzmax)
      double complex ysup(6,3,-1:nzmax),yslw(6,3,-1:nzmax)
      double complex ma0(6,6),dmai(6,6),c0(6,3)
      double complex coef(6,6),b(6),coefl(12,12),bl(12)
c
      ck=dcmplx(k,0.d0)
c
c===============================================================================
c
      do l=-1,0
        call pematrix(maup(1,1,l),cs,k,hp(l),n)
        call pematinv(maiup(1,1,l),cs,k,n)
      enddo
      do l=-1,lp-1
        call pematrix(malw(1,1,l),cs,k,-hp(l),n)
        call pematinv(mailw(1,1,l),cs,k,n)
      enddo
c
c     matrix propagation from surface to source
c
      do j=1,3
        do i=1,6
          c0(i,j)=(0.d0,0.d0)
        enddo
      enddo
c
c     without free surface: full-space theory
c
      c0(1,1)=(1.d0,0.d0)
      c0(3,2)=(1.d0,0.d0)
      c0(5,3)=(1.d0,0.d0)
      call pematrix(ma0,cs,k,0.d0,n)
      call caxcb(ma0,c0,6,6,3,yup(1,1,-1))
c
c===============================================================================
c
c     matrix propagation from half-space to source
c
      do j=1,3
        do i=1,6
          c0(i,j)=(0.d0,0.d0)
        enddo
      enddo
c
c     coefficient vectors in the half-space
c
      c0(2,1)=(1.d0,0.d0)
      c0(4,2)=(1.d0,0.d0)
      c0(6,3)=(1.d0,0.d0)
      call pematrix(ma0,cs,k,0.d0,n)
      call caxcb(ma0,c0,6,6,3,ylw(1,1,lp))
c
      call pepropfs(lp,0,ylw,cs,k,n)
c
c===============================================================================
c
c     conditions on the source surface
c
c
c     source function
c
      if(pointsource)then
        do i=1,6
          b(i)=sfct(i)
        enddo
c
        do i=1,6
          do j=1,3
            coef(i,j)=yup(i,j,-1)
            coef(i,j+3)=-ylw(i,j,0)
          enddo
        enddo
        key=0
        call cdsvd500(coef,b,6,1,1.d-30,key)
        if(key.eq.0)then
          print *,' Warning in pepsvfsm: anormal exit from cdsvd500!'
          return
        endif
        do l=1,lp
          do i=1,6
            y(i,l)=(0.d0,0.d0)
            do j=1,3
              y(i,l)=y(i,l)+b(j+3)*ylw(i,j,l)
            enddo
          enddo
        enddo
      else
c
c       vertical line source
c
        do j=1,3
          do i=1,6
            ysup(i,j,-1)=yup(i,j,-1)
          enddo
        enddo
c
        call pepropfs(-1,0,ysup,cs,k,n)
c
        do j=1,3
          do i=1,6
            yslw(i,j,0)=ylw(i,j,0)
          enddo
        enddo
c
        call pepropfs(0,-1,yslw,cs,k,n)
c
        call pedifmai(dmai,cs,k,n)
        call caxcb(dmai,sfct,6,6,1,b)
c
        do i=1,6
c
          bl(i  )=b(i)
          bl(i+6)=b(i)
c
          do j=1,3
            coefl(i,j  )= yup(i,j,-1)
            coefl(i,j+3)=-ysup(i,j,-1)
            coefl(i,j+6)=-yslw(i,j,-1)
            coefl(i,j+9)= (0.d0,0.d0)
          enddo
c
          do j=1,3
            coefl(i+6,j  )= (0.d0,0.d0)
            coefl(i+6,j+3)=-ysup(i,j,0)
            coefl(i+6,j+6)=-yslw(i,j,0)
            coefl(i+6,j+9)= ylw(i,j,0)
          enddo
        enddo
        key=0
        call cdsvd500(coefl,bl,12,1,1.d-30,key)
        if(key.eq.0)then
          print *,' Warning in pepsvsm: anormal exit from cdsvd500!'
          return
        endif
        do l=1,lp
          do i=1,6
            y(i,l)=(0.d0,0.d0)
            do j=1,3
              y(i,l)=y(i,l)+bl(j+9)*ylw(i,j,l)
            enddo
          enddo
        enddo
      endif
c
c     convert poroelastic normal tress to matrix stress
c
      do l=1,lp
        y(2,l)=y(2,l)+dcmplx(alf(n),0.d0)*y(5,l)
        y(6,l)=y(6,l)-dcmplx(alf(n),0.d0)*cs*y(1,l)
      enddo
c
      return
      end