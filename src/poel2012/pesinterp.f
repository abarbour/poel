      subroutine pesinterp(dtg)
      implicit none
      double precision dtg
c
      include 'peglobal.h'
c
      integer its,its0,dnts
      double precision ddtg,dsinj
c
      nts=0
      do its0=2,nts0
        if(sinj0(its0).eq.sinj0(its0-1))then
c
c         ignored if no change in injection
c
        else if(ts0(its0).eq.ts0(its0-1))then		
          nts=nts+1
          if(nts.gt.ntsmax)then
            stop ' Error in pesinterp: ntsmax defined too small!'
          endif
          ts(nts)=ts0(its0)
          sinj(nts)=sinj0(its0)-sinj0(its0-1)
        else if(ts0(its0).gt.ts0(its0-1))then
          dnts=1+idint((ts0(its0)-ts0(its0-1))/dmin1(dt,dtg))
          ddtg=(ts0(its0)-ts0(its0-1))/dble(dnts)
          dsinj=(sinj0(its0)-sinj0(its0-1))/dble(dnts)
          do its=1,dnts
            nts=nts+1
            if(nts.gt.ntsmax)then
              stop ' Error in pesinterp: ntsmax defined too small!'
            endif
            ts(nts)=ts0(its0-1)+dble(its)*ddtg
            sinj(nts)=dsinj
            if(ts(nts).ge.twindow)return
          enddo
        else
          stop ' Error in sinterp: wrong injection series!'
        endif
      enddo
	return
      end