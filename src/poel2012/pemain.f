      program poel2012
      implicit none
c
      include 'peglobal.h'
c
c     work space
c
      integer i,ir,ipr,izr,nr1,nr2,ierr
      double precision accuracy,d1,d2,dbeta
      external pewvint,pewvintd
c
c     read input file file
c
      print *,'######################################################'
      print *,'#                                                    #'
      print *,'#                Welcome to Program                  #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#         PPPP     OOO     EEEEE    L                #'
      print *,'#         P   P   O   O    E        L                #'
      print *,'#         PPPP    O   O    EEEE     L                #'
      print *,'#         P       O   O    E        L                #'
      print *,'#         P        OOO     EEEEE    LLLLLL           #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#                   Version-2012                     #'
      print *,'#                   ============                     #'
      print *,'#             Last modified: July 2012               #'
      print *,'#                                                    #'
      print *,'#                   * * * * *                        #'
      print *,'#         A semi-analytical code for simulating      #'
      print *,'#    fully coupled deformation-diffusion processes   #'
      print *,'#          in a layered poroelastic half-space       #'
      print *,'#          induced by an injection (pump) test       #'
      print *,'#          or an initial pore pressure anomaly       #'
      print *,'#                   * * * * *                        #'
      print *,'#                                                    #'
      print *,'#                   Copyright                        #'
      print *,'#                      by                            #'
      print *,'#                 Rongjiang Wang                     #'
      print *,'#              (wang@gfz-potsdam.de)                 #'
      print *,'#           GeoForschungsZentrum Potsdam             #'
      print *,'#                                                    #'
      print *,'######################################################'
      print *,'                                                      '
      write(*,'(a,$)')' Please type the file name of input data: '
      read(*,'(a)')inputfile
c
c     get input data and construct layered model
c
      call pegetinp(accuracy)
c
      write(*,'(a)')'  Wavenumber integration ...'
c
      r0=dmax1(sradius,0.2d0*rsdismin,0.2d0*r(1))
c
      d1=dmax1(rsdismin,r(1))
      d2=r(nr)
      i=1+idint(dlog(d2/d1)/dlog(5.d0))
      dbeta=dexp(dlog(d2/d1)/dble(i))
c
      ipr=1
      nr1=1
100   nr2=nr1
      do ir=nr1+1,nr
        if(r(ir).le.dmax1(sradius,dbeta*r(nr1),dbeta*rsdismin))then
          nr2=ir
        endif
      enddo
c
      if(r(nr1).le.r0)then
        nd=0
      else
        nd=min0(max0(0,
     &   idnint(dlog10(r(nr2)/dsqrt(r0**2+zrsmin**2)))),ndmax)
      endif
c
      analytic=r(nr1).ge.5.d0*sradius
c
      if(nd.le.0)then
        call pespectra(ipr,accuracy,nr1,nr2,pewvint,ierr)
      else
        call pespectra(ipr,accuracy,nr1,nr2,pewvintd,ierr)
      endif
c
      if(nr2.lt.nr)then
        ipr=ipr+1
        nr1=nr2+1
        goto 100
      endif
c
c     output time series and snapshots
c
      write(*,'(a)')'  Output ...'
      call peoutput(ierr)
c
      print *,'######################################################'
      print *,'#                                                    #'
      print *,'#          End of computations with POEL12           #'
      print *,'#                                                    #'
      print *,'######################################################'
      stop
      end
