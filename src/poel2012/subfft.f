      subroutine subfft(fgrn,tgrn,tgcorr,fcut,nfg,nfmax,dfg,dtg,alpha,
     &                  istype,ts,sinj,nts,obs,nt,dt,eps)
      implicit none
c
      integer nfg,nfmax,nts,nt,istype
      double precision fcut,dfg,dtg,dt,alpha,eps
      double precision ts(nts),sinj(nts)
      double precision tgrn(2*nfmax),tgcorr(2*nfmax)
      double precision obs(4*nfmax)
      double complex fgrn(2*nfmax)
c
c     local memory
c
      integer i,it,itga,itgb,lf,jf,ntg,its
      double precision a,b,t,fc,fgrnabs,fgrnmax,beta,pi2
c
      pi2=8.d0*datan(1.d0)
      ntg=2*nfg
c
c     low-pass filter has to be used
c
c     find a proper corner frequency as high as possible
c
      fc=fcut/eps
      beta=1.d0+alpha/(pi2*fc)
300   fgrnmax=0.d0
      do lf=1,nfg
        fgrnabs=cdabs(fgrn(lf))/(beta**2
     &         +(dble(lf-1)*dfg/fc)**2)**1.5d0
        fgrnmax=dmax1(fgrnmax,fgrnabs)
      enddo
      if(fgrnabs.gt.eps*fgrnmax)then
        fc=0.75d0*fc
        beta=1.d0+alpha/(pi2*fc)
        goto 300
      endif
c
c     use low-pass filter with the corner frequency found
c
      do lf=1,nfg
        fgrn(lf)=fgrn(lf)/(dcmplx(beta,dble(lf-1)*dfg/fc))**3
      enddo
c
	jf=1
	do lf=2*nfg,nfg+2,-1
	  jf=jf+1
	  fgrn(lf)=dconjg(fgrn(jf))
	enddo
	fgrn(nfg+1)=(0.d0,0.d0)
c
c	convention for Fourier transform:
c	f(t)=\int F(f) exp(i2\pi f t) df
c
      do lf=1,2*nfg
        obs(2*lf-1)=dreal(fgrn(lf))
        obs(2*lf  )=dimag(fgrn(lf))
      enddo
c
      call four1(obs(1),2*nfg,+1)
c
      if(istype.eq.0)then
c
c       calculate response to the impulsive fluid injection
c
        tgrn(1)=0.d0
        do it=2,ntg
          tgrn(it)=tgcorr(it)
     &            +obs(2*it-1)*dfg*dexp(alpha*dble(it-1)*dtg)
        enddo
c
        do it=1,2*((1+nt)/2)
          t=dble(it-1)*dt
          obs(it)=0.d0
          b=dmod(t/dtg,1.d0)
          a=1.d0-b
          itga=min0(1+idint(t/dtg),ntg)
          itgb=min0(2+idint(t/dtg),ntg)
          obs(it)=obs(it)+sinj(1)*(tgrn(itga)*a+tgrn(itgb)*b)
        enddo
      else
c
c       calculate response to the Heaviside source history
c
        tgrn(1)=0.d0
        do it=2,ntg
          tgrn(it)=tgrn(it-1)+obs(2*it-1)
     &            *dtg*dfg*dexp(alpha*dble(it-1)*dtg)
        enddo
c
        do it=1,ntg
          tgrn(it)=tgrn(it)+tgcorr(it)
        enddo
c
c       convolution with source time function
c
        do it=1,2*((1+nt)/2)
          t=dble(it-1)*dt
          obs(it)=0.d0
          do its=1,nts
            if(ts(its).le.t)then
              b=dmod((t-ts(its))/dtg,1.d0)
              a=1.d0-b
              itga=min0(1+idint((t-ts(its))/dtg),ntg)
              itgb=min0(2+idint((t-ts(its))/dtg),ntg)
              obs(it)=obs(it)+sinj(its)*(tgrn(itga)*a+tgrn(itgb)*b)
            endif
          enddo
        enddo
      endif
c
      do it=1,(1+nt)/2
        fgrn(it)=dcmplx(obs(2*it-1),obs(2*it))
      enddo
c
      return
      end