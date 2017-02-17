S = ../src_3.38b/
M = $(S)math/
COMP = ifort -o $@ -w -O2 -I$(M) 
OPTI = opti2ltr.f90 $(M)math_lib.o $(M)matr_lib.o $(M)eigen_r8.o $(S)matr_aux.o
DIST = distrwgh.f90 $(M)math_lib.o
EMIT = emitgfit.f90 $(M)math_lib.o $(M)zip.o
SLIC = slice.f90    $(M)math_lib.o

all        : opti2ltr distrwgh emitgfit slice
opti2ltr   : $(OPTI)
	$(COMP) $(OPTI)
distrwgh   : $(DIST)
	$(COMP) $(DIST)
emitgfit   : $(EMIT)
	$(COMP) $(EMIT)
slice      : $(SLIC)
	$(COMP) $(SLIC)
