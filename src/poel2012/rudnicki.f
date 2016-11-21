      subroutine rudnicki(z,r,t,zs,la,mu,alpha,chi,d,
     &                    uz,ur,p,
     &                    ezz,err,ett,ezr,tlt,dpz,dpr,istype)
      implicit none
c
c     response of a whole-space poroelastic medium
c     to a point Heaviside injection source
c
c     input:
c     z,r = receiver position
c     t = time
c     zs = injection depth
c     la,mu = elasticity parameters
c     alpha = effective stress
c     chi = Darcy conductivity
c     d = hydraulic diffusivity
c
      integer istype
      double precision z,r,t,zs,la,mu,alpha,chi,d
c
c     output:
c     uz,ur,p = displacement vector and pore pressure increase
c     ezz,err,ett,ezr,tlt = strain tensor and tilt of borehole
c     daz,dar = Darcy flux
c     
      double precision uz,ur,p,ezz,err,ett,ezr,tlt,dpz,dpr
c
c     local work memory
c
      double precision u0,p0,dis,eta,doteta,fct,pfct,ppfct,dp,perfc
      double precision erf,erfc
c
      double precision pi
      data pi/3.14159265358979d0/
c
      dis=dsqrt((z-zs)**2+r*r)
      if(t.le.0.d0.or.dis.le.0.d0)then
        uz=0.d0
        ur=0.d0
        p=0.d0
        ezz=0.d0
        err=0.d0
        ett=0.d0
        ezr=0.d0
        tlt=0.d0
        dpz=0.d0
        dpr=0.d0
      else
        pi=4.d0*datan(1.d0)
        p0=1.d0/(4.d0*pi*chi)
        u0=p0*alpha/(2.d0*(la+2.d0*mu))
        eta=0.5d0*dis/dsqrt(d*t)
        fct=erfc(eta)
     &     +0.5d0*erf(eta)/eta**2
     &     -dexp(-eta**2)/(dsqrt(pi)*eta)
        pfct=-erf(eta)/eta**3+2.d0*dexp(-eta**2)/(dsqrt(pi)*eta**2)
c
        if(istype.eq.0)then
c
c         response to impulsive source
c
          perfc=-2.d0*dexp(-eta**2)/dsqrt(pi)
          ppfct=perfc/eta**3+3.d0*erf(eta)/eta**4
     &         -4.d0*dexp(-eta**2)*(1.d0+eta**2)/(dsqrt(pi)*eta**3)
          doteta=-0.5d0*eta/t
          
          uz=u0*pfct*doteta*(z-zs)/dis
          ur=u0*pfct*doteta*r/dis
          p=p0*perfc*doteta/dis
c
          ezz=u0*(pfct*r**2+(ppfct*eta+pfct)*(z-zs)**2)*doteta/dis**3
          err=u0*(pfct*(z-zs)**2+(ppfct*eta+pfct)*r**2)*doteta/dis**3
          ett=u0*pfct*doteta/dis
          ezr=u0*eta*ppfct*doteta*(z-zs)*r/dis**3
          tlt=-ezr
c
          dp=p0*(perfc+2.d0*(1.d0-2.d0*eta**2)*dexp(-eta**2)
     &                /dsqrt(pi))*doteta
          dpz=-dp*(z-zs)/dis**3
          dpr=-dp*r/dis**3
        else if(istype.eq.1)then
c
c         response to heaviside source
c
          uz=u0*fct*(z-zs)/dis
          ur=u0*fct*r/dis
          p=p0*erfc(eta)/dis
c
          ezz=u0*(fct*r**2+pfct*eta*(z-zs)**2)/dis**3
          err=u0*(fct*(z-zs)**2+pfct*eta*r**2)/dis**3
          ett=u0*fct/dis
          ezr=u0*(eta*pfct-fct)*(z-zs)*r/dis**3
          tlt=-ezr
c
          dp=p0*(erfc(eta)+2.d0*eta*dexp(-eta**2)/dsqrt(pi))
          dpz=-dp*(z-zs)/dis**3
          dpr=-dp*r/dis**3
        endif
      endif 
c
      return
      end