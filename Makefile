#-----------setup server env ------------------
Platform=Medivh
#Platform=COSMA
#Platform=Bright
ifeq ($(shell uname -n), Medivh)
	Platform=Medivh
endif
ifeq ($(USER), jvbq85)
	Platform=COSMA
endif


ifeq ($(Platform), COSMA)
	HDFLIB= -lhdf5_hl -lhdf5
	GSLLIB=-lgsl -lgslcblas
	CC = icc
	#CC=h5cc
endif

ifeq ($(Platform), Bright)
	HDFLIB= -lhdf5_hl -lhdf5 -lz  #shao cluster
	GSLINC=-I/data/raid3/jxhan/opt/include #shao
	GSLLIB=-L/data/raid3/jxhan/opt/lib -lgsl -lgslcblas #shao
	CC = icc
endif

ifeq ($(Platform), Medivh)
	HDFLIB= -lhdf5_hl -lhdf5
	HDFINC= -I/usr/include/mpi
	GSLLIB=-lgsl -lgslcblas
	CC = icc
	#CC=h5cc
endif
#--------other flags----------------------------
DEBUG=on

ifeq ($(CC),icc)
LDFLAGS+= -limf
OMPLIB= -openmp
endif
ifeq ($(CC),gcc)
LDFLAGS+= -lm 
OMPLIB= -fopenmp
endif

#SRCFLAG= -DSRCFILE='"shearmap_MODLE2Z0X0M0.seed1024"'

CFLAGS += $(HDFINC) $(GSLINC) -DROOTDIR=$(ROOTDIR) $(OMPLIB)
LDFLAGS += $(HDFLIB) $(GSLLIB) $(OMPLIB)

#~ ifeq ($(USER),kam)      #my laptop has hdf v1.6
ifeq ($(shell uname -n),kam-laptop)
CFLAGS+= -DHDF_V16
endif

ifeq ($(DEBUG),on)
CFLAGS+= -g -Wall
LDFLAGS+= -g
endif

#-----File Dependencies----------------------
SRC_COMM = hdf_util.c globals.c cosmology.c mymath.c models.c tracer.c halo.c template.c nfw.c
OBJS_COMM= $(patsubst %.f90,%.f.o,$(SRC_COMM:%.c=%.o))

#-----targets and common rules--------------------------------
default: lib

#the default rule will handle the rest
# %.o : %.c
# 	$(CC) $< $(CFLAGS) -c -o $@
# 
%.f.o : %.f90
	$(FC) $< $(FFLAGS) -c -o $@
	
lib: CFLAGS+=-fPIC
lib: FFLAGS+=-fPIC
lib: liboPDF.so

liboPDF.so:$(OBJS_COMM)
	$(CC) -shared -Wl,-soname,liboPDF.so -o liboPDF.so $^ $(LDFLAGS)
#----------------------
submit: tmpdir=exe/mockfitFmin_mc$(ESTIMATOR)
submit: FORCE
	mkdir -p $(tmpdir)
	@$(MAKE) lib -B
	cp dynio.py libdyn.so mockFmin.py job.bsub $(tmpdir)
	cd $(tmpdir);\
	bsub <job.bsub #do something immediately after the previous command; keep in same line so that the same subshell is used.
# 	@$(MAKE) -C $(tmpdir) job #to use this, a makefile has to be in the subdir first.
	
#-----Other stuff----------------------------
.PHONY : clean depend distclean

synccosma:
	rsync -avz $(shell pwd)/../ jvbq85@login.cosma.dur.ac.uk:data/DynDistr/code

synccosmalocal:
	rsync -e "ssh -p 4800" -avz $(shell pwd)/../ jvbq85@localhost:data/DynDistr/code


depend:
	makedepend --$(CFLAGS)-- -Y $(SRC_COMM) $(SRC_MAIN) $(SRC)
	
clean:
	rm -f $(OBJS) $(OBJS_COMM) $(OBJS_MAIN)
	
distclean: clean
	rm -f $(EXEC) liboPDF.so *~ *.pyc py/*.pyc
#-----end--- auto dependencies below---------	
# DO NOT DELETE

hdf_util.o: hdf_util.h mymath.h
globals.o: globals.h
cosmology.o: cosmology.h globals.h
mymath.o: mymath.h
models.o: mymath.h globals.h cosmology.h tracer.h halo.h models.h
tracer.o: mymath.h globals.h hdf_util.h tracer.h halo.h models.h
halo.o: mymath.h globals.h cosmology.h halo.h template.h nfw.h tracer.h
template.o: mymath.h globals.h cosmology.h halo.h template.h TemplateData.h
nfw.o: globals.h cosmology.h halo.h tracer.h
