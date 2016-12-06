MAIN=pemain.f
H=peglobal.h
SUBS=bessj0.f \
 bessj1.f \
 caxcb.f \
 cdsvd500.f \
 cmemcpy.f \
 erf.f \
 erfc.f \
 four1.f \
 getdata.f \
 pebsj.f \
 pedifmai.f \
 pegetinp.f \
 pekern.f \
 pelayer.f \
 pematinv.f \
 pematrix.f \
 peoutput.f \
 pepropfs.f \
 peproppsv.f \
 pepsv.f \
 pepsvfs.f \
 pepsvfsm.f \
 pesinterp.f \
 pesnaps.f \
 pesource.f \
 pespectra.f \
 pesublay.f \
 pewvint.f \
 pewvintd.f \
 rudnicki.f \
 subfft.f \
 tsoutput.f

PO=poel2012

all: clean $(PO)

$(PO):
	#gfortran -o $@ $(MAIN) $(H) $(SUBS)
	gfortran -o $@ $(MAIN) $(SUBS)

clean:
	rm -f *.o
	rm -f poel

