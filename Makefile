COMPILER = ifort  # gfortran -w       #  -ffixed-line-length-none

sncf.x : sncf.o nm.o bsec.o cbcs_new.o  bcsrad2.o test4e.o  libra_new.o from_nag.o timelog.o
	$(COMPILER) -o sncf.x sncf.o nm.o bsec.o cbcs_new.o   bcsrad2.o test4e.o  libra_new.o from_nag.o timelog.o
sncf.o : sncf.f 
	$(COMPILER) -c sncf.f
nm.o : nm.f
	$(COMPILER) -c nm.f
bsec.o : bsec.f
	$(COMPILER) -c bsec.f
cbcs_new.o : cbcs_new.f
	$(COMPILER) -c cbcs_new.f 	
bcsrad2.o : bcsrad2.f
	$(COMPILER) -c bcsrad2.f
test4e.o : test4e.f 
	$(COMPILER) -c test4e.f 
# gr0_j.o : gr0_j.f 
#	$(COMPILER) -c gr0_j.f 
libra_new.o : libra_new.f
	$(COMPILER) -c libra_new.f
from_nag.o : from_nag.f
	$(COMPILER) -c from_nag.f
timelog.o : timelog.f
	$(COMPILER) -c timelog.f

clean:
	rm -f *\.mod *genmod* *\.o *~
