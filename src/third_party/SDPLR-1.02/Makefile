sdplr:
	cd source; make

mexsdplr:
	cd source; make -f Makefile.mex

gsl:
	cd gsl-1.5; make

all: clean gsl sdplr

mex: clean gsl mexsdplr

clean:
	cd gsl-1.5; make clean
	cd source; make clean

cleanall: clean
	rm -rf sdplr
	rm -rf mexsdplr.*
	rm -rf lib/libgsl.a

