c     GLOBAL INDEX PARAMETERS FOR DEFINING ARRAYS
c     ===========================================
c     lmax: max. no of total homogeneous layers
c     nzrmax: max. no of receiver depths
c     nzmax: max. interface index (= lmax+nzrmax+1)
c     nrmax: max. no of distances
c     nfmin: min. no of frequency samples for Green's functions (power of 2)
c     nfmax: max. no of frequency samples  for Green's functions (power of 2)
c     nbsjmax: max. no of wavenumber samples
c     ntsmax: max. no of pump rate samples
c     nsnmax = max. no of snapshot outputs
c
      integer lmax,nzrmax,nzmax,nrmax
      integer nfmin,nfmax,nbsjmax,ntsmax
      integer nsnmax,ndmax
      parameter(lmax=50)
      parameter(nzrmax=51)
      parameter(nzmax=lmax+nzrmax+2)
      parameter(nrmax=51)
      parameter(nfmin=128)
      parameter(nfmax=4096)
      parameter(nbsjmax=65536)
      parameter(ntsmax=1000)
      parameter(nsnmax=100)
      parameter(ndmax=4)
c
c     DISCRETISATION ACCURACY FOR LAYERS WITH CONSTANT GRADIENT
c     =========================================================
c     reslm: for moduls
c     resla: for alpha
c     resld: for diffusivity
c
      double precision reslm,resla,resld
      parameter(reslm=0.05d0,resla=0.05d0,resld=0.250d0)
c
c     original model parameters
c
      integer l0
      double precision z1(lmax),z2(lmax)
      double precision la1(lmax),la2(lmax),mu1(lmax),mu2(lmax)
      double precision alf1(lmax),alf2(lmax),qa1(lmax),qa2(lmax)
      double precision dm1(lmax),dm2(lmax)
      common /imodel0/ l0
      common /dmodel0/ z1,z2,la1,la2,mu1,mu2,alf1,alf2,qa1,qa2,
     &             dm1,dm2
c       
c     model parameter:
c     n0: number of homogeneous layers
c
      integer n0
      double precision h(lmax),la(lmax),mu(lmax)
      double precision alf(lmax),qa(lmax),dm(lmax)
      common /imodel/ n0
      common /dmodel/ h,la,mu,alf,qa,dm
c
      double complex ka(lmax),ck3(lmax)
      logical smalls(lmax)
      common /dphases/ ka,ck3
      common /lphases/ smalls
c
      integer lp,nno(nzmax)
      double precision hp(-1:nzmax)
      common /isublay/ lp,nno
      common /dsublay/ hp
c
c     zrec: receiver depth
c     lzrec: sublayer no of receiver
c
      integer nzr
      integer lzr(nzrmax)
      double precision zr(nzrmax)
      common /ireceiver/ nzr,lzr
      common /dreceiver/ zr
c
      integer nr,nd
      double precision r(nrmax)
      common /irprofile/ nr,nd
      common /drprofile/ r
c
      integer nt
      double precision twindow,dt
      common /itimes/ nt
      common /dtimes/ twindow,dt
c
c     source parameters
c
      integer istype,lstop,lsbtm
      double precision zstop,zsbtm,sradius,ppini
      double complex sfct(6)
      logical pointsource,analytic,mirror
      common /isource/ istype,lstop,lsbtm
      common /dsource/ zstop,zsbtm,sradius,ppini,sfct
      common /lsource/ pointsource,analytic,mirror
c
c     table of J_n(x), n = 0, 1
c
      double precision rsdismax,rsdismin,zrsmin,r0
      double precision disk(nbsjmax),bsj(nbsjmax,0:1,nrmax)
      double precision zrs2(nzrmax),fana(nzrmax)
      double precision rsdis(nrmax,nzrmax)
      common /bessels/ rsdismax,rsdismin,zrsmin,
     &                 r0,disk,bsj,zrs2,fana,rsdis
c
      integer isurfcon
      common /surface/ isurfcon
c
c     source functions
c
      integer nts0,nts
      double precision ts0(ntsmax),sinj0(ntsmax)
      double precision ts(ntsmax),sinj(ntsmax)
      common /iinjection/ nts0,nts
      common /dinjection/ ts0,sinj0,ts,sinj
c
c     psv layer matrics
c
      double complex maup(6,6,-1:nzmax),maiup(6,6,-1:nzmax)
      double complex malw(6,6,-1:nzmax),mailw(6,6,-1:nzmax)
      common /psvlayma/ maup,maiup,malw,mailw
c
c     title text
c
      character*12 txttime
      character*4 txtuz,txtur,txtpp,txtezz,txterr,
     &        txtett,txtezr,txttlt,txtdvz,txtdvr
      character*6 txtdep(nzrmax),txtdis(nrmax)
      common /title/txttime,txtuz,txtur,txtpp,txtezz,txterr,
     &        txtett,txtezr,txttlt,txtdvz,txtdvr,
     &        txtdep,txtdis
c
c     memories for time series
c
      integer seluz,selur,selpp,selezz,selerr,
     &        selett,selezr,seltlt,seldvz,seldvr
      character*80 inputfile,fileuz,fileur,filepp,
     &        fileezz,fileerr,fileett,fileezr,filetlt,
     &        filedvz,filedvr
      double complex grnuz(nfmax,nrmax,nzrmax)
      double complex grnur(nfmax,nrmax,nzrmax)
      double complex grnpp(nfmax,nrmax,nzrmax)
      double complex grnezz(nfmax,nrmax,nzrmax)
      double complex grness(nfmax,nrmax,nzrmax)
      double complex grnerr(nfmax,nrmax,nzrmax)
      double complex grnett(nfmax,nrmax,nzrmax)
      double complex grnezr(nfmax,nrmax,nzrmax)
      double complex grntlt(nfmax,nrmax,nzrmax)
      double complex grndvz(nfmax,nrmax,nzrmax)
      double complex grndvr(nfmax,nrmax,nzrmax)
      common /iselection/ seluz,selur,selpp,selezz,selerr,
     &        selett,selezr,seltlt,seldvz,seldvr
      common /filenames/ inputfile,fileuz,fileur,filepp,
     &        fileezz,fileerr,fileett,fileezr,filetlt,
     &        filedvz,filedvr
      common /dispfield/ grnuz,grnur,grnpp,grnezz,grness,
     &                   grnerr,grnett,grnezr,grntlt,
     &                   grndvz,grndvr
c
c     memories for snapshots
c
      integer nsn
      double precision timesn(nsnmax)
      character*80 filesn(nsnmax)
      common /isnapshots/ nsn
      common /dsnapshots/ timesn
      common /fsnapshots/ filesn