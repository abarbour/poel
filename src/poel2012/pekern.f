      subroutine pekern(cs,k,cuz,cur,cp,cezz,cess,cezr,cdpz)
      implicit none
c
c     calculation of response function in Laplace domain
c     k: wave number
c     cs: complex Laplace variable
c     uz,ur = displacement
c     p = pore pressure
c     ezz, ess, ezr = normal, surface, shear strain
c     dpz = pore pressure vertical gradient
c
      include 'peglobal.h'
c
      double precision k
      double complex cs
      double complex cuz(nzrmax),cur(nzrmax),cp(nzrmax)
      double complex cezz(nzrmax),cess(nzrmax),cezr(nzrmax)
      double complex cdpz(nzrmax)
c
      integer i,n,ly,izr
      double precision beti
      double complex ck,cla,cmu,chi,cla0,cmu0,chi0
      double complex y(6,nzmax),y0(6,nzmax),ym(6,nzmax)
c
      double precision eps
      data eps/1.0d-03/
c
      ck=dcmplx(k,0.d0)
c
      do n=1,n0
        beti=(1.d0+alf(n)*qa(n)/(la(n)+2.d0*mu(n)))/dm(n)
        ka(n)=cdsqrt(ck*ck+cs*dcmplx(beti,0.d0))
        smalls(n)=cdabs(ka(n)-ck)/k.lt.eps
        if(smalls(n))then
          ck3(n)=ck
        else
          ck3(n)=ka(n)
        endif
      enddo
c
      call pepsv(y,cs,k)
      do izr=1,nzr
        ly=lzr(izr)
        n=nno(ly)
        cla=dcmplx(la(n),0.d0)
        cmu=dcmplx(mu(n),0.d0)
        chi=dcmplx(dm(n)*alf(n)/qa(n),0.d0)
c
        cuz(izr)=y(1,ly)
        cur(izr)=y(3,ly)
        cp(izr)=y(5,ly)
c
        cezz(izr)=(y(2,ly)+cla*ck*y(3,ly))/(cla+(2.d0,0.d0)*cmu)
        cess(izr)=-ck*y(3,ly)
        cezr(izr)=(0.5d0,0.d0)*y(4,ly)/cmu
        cdpz(izr)=-y(6,ly)/chi
      enddo
      if(.not.analytic)return
c
      n=nno(lstop)
      cla0=dcmplx(la(n),0.d0)
      cmu0=dcmplx(mu(n),0.d0)
      chi0=dcmplx(dm(n)*alf(n)/qa(n),0.d0)
c
      call pepsvfs(y0,cs,k,n)
c
      if(mirror.and.isurfcon.eq.1)then
        call pepsvfsm(ym,cs,k,n)
        do ly=1,lp
          do i=1,6
            y0(i,ly)=y0(i,ly)-ym(i,ly)
          enddo
        enddo
      else if(mirror.and.isurfcon.eq.2)then
        call pepsvfsm(ym,cs,k,n)
        do ly=1,lp
          do i=1,6
            y0(i,ly)=y0(i,ly)+ym(i,ly)
          enddo
        enddo
      endif
c
      do ly=1,lp
        do i=1,6
          y0(i,ly)=dcmplx(fana(ly),0.d0)*y0(i,ly)
        enddo
      enddo
c
      do izr=1,nzr
        ly=lzr(izr)
        cuz(izr)=cuz(izr)-y0(1,ly)
        cur(izr)=cur(izr)-y0(3,ly)
        cp(izr)=cp(izr)-y0(5,ly)
c
        cezz(izr)=cezz(izr)
     &     -(y0(2,ly)+cla0*ck*y0(3,ly))/(cla0+(2.d0,0.d0)*cmu0)
        cess(izr)=cess(izr)+ck*y0(3,ly)
        cezr(izr)=cezr(izr)-(0.5d0,0.d0)*y0(4,ly)/cmu0
        cdpz(izr)=cdpz(izr)+y0(6,ly)/chi0
      enddo
c
      return
      end