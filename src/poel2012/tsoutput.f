      subroutine tsoutput(outfile,txttime,txtcmp,txtdis,txtdep,
     &  cdata,nfmax,nrmax,nzrmax,nt,nr,nzr,dt)
      implicit none
      integer nfmax,nrmax,nzrmax,nt,nr,nzr
      double precision dt
      double complex cdata(nfmax,nrmax,nzrmax)
      character*80 outfile
      character*12 txttime
      character*4 txtcmp
      character*6 txtdis(nrmax),txtdep(nzrmax)
c
      integer it,ir,izr
c
      open(30,file=outfile,status='unknown')
      write(30,'(a12,$)')txttime
      do izr=1,nzr-1
        do ir=1,nr
	    write(30,'(a16,$)')txtcmp//txtdep(izr)//txtdis(ir)
	  enddo
	enddo
      do ir=1,nr-1
        write(30,'(a16,$)')txtcmp//txtdep(nzr)//txtdis(ir)
      enddo
      write(30,'(a16)')txtcmp//txtdep(nzr)//txtdis(nr)
      do it=1,(1+nt)/2
        if(2*it-1.gt.nt)goto 100
        write(30,'(E12.5,$)')dble(2*it-2)*dt
        do izr=1,nzr-1
          do ir=1,nr
            write(30,'(E16.8,$)')dreal(cdata(it,ir,izr))
          enddo
        enddo
        do ir=1,nr-1
          write(30,'(E16.8,$)')dreal(cdata(it,ir,nzr))
        enddo
        write(30,'(E16.8)')dreal(cdata(it,nr,nzr))
        if(2*it.gt.nt)goto 100
        write(30,'(E12.5,$)')dble(2*it-1)*dt
        do izr=1,nzr-1
          do ir=1,nr
            write(30,'(E16.8,$)')dimag(cdata(it,ir,izr))
          enddo
        enddo
        do ir=1,nr-1
          write(30,'(E16.8,$)')dimag(cdata(it,ir,nzr))
        enddo
        write(30,'(E16.8)')dimag(cdata(it,nr,nzr))
      enddo
100   close(30)
      return
      end