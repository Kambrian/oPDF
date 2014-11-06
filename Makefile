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
	FC = ifort
	#CC=h5cc
	MPICC=mpicc 	
	ROOTDIR='"/cosma/home/jvbq85/data/DynDistr"'
# 	LAPACKLIB=-mkl #no longer needed.
# 	LAPACKINC=-I$(MKLINCLUDE)
endif

ifeq ($(Platform), Bright)
	HDFLIB= -lhdf5_hl -lhdf5 -lz  #shao cluster
	GSLINC=-I/data/raid3/jxhan/opt/include #shao
	GSLLIB=-L/data/raid3/jxhan/opt/lib -lgsl -lgslcblas #shao
	CC = icc
	FC = ifort
	MPICC:=mpicc -cc=$(CC)
	ROOTDIR='"/data/raid3/jxhan/Lensing"'
endif

ifeq ($(Platform), Medivh)
	HDFLIB= -lhdf5_hl -lhdf5
	HDFINC= -I/usr/include/mpi
	GSLLIB=-lgsl -lgslcblas
	CC = icc
	FC = ifort
	#CC=h5cc
	MPICC=icc -lmpi
# 	LAPACKLIB=-mkl
# 	LAPACKINC=-I$(MKLINCLUDE)
	ROOTDIR='"/work/Projects/DynDistr"' 	
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
SRC_MAIN=scan.c NumericalConvergence.c saveHDFstars.c mock_stars.c
EXEC=$(basename $(SRC_MAIN))
OBJS_MAIN=$(addsuffix .o, $(EXEC))
SRC_COMM = hdf_util.c globals.c cosmology.c mymath.c models.c tracer.c halo.c template.c nfw.c
OBJS_COMM= $(patsubst %.f90,%.f.o,$(SRC_COMM:%.c=%.o))
# OBJS_COMM = $(SRC_COMM:%.c=%.o)
#any additional OBJS below:
# SRC = wenting.c wenting.f90
# OBJS= $(patsubst %.f90,%.f.o,$(SRC:%.c=%.o))

#-----targets and common rules--------------------------------
default: scan
all: $(EXEC)

#the default rule will handle the rest
$(EXEC): $(OBJS_COMM) 

# %.o : %.c
# 	$(CC) $< $(CFLAGS) -c -o $@
# 
%.f.o : %.f90
	$(FC) $< $(FFLAGS) -c -o $@
	
#mpi flags
#gama_WL_rand: CC:=$(MPICC)

#Force recompile
#shear.o: FORCE

lib: CFLAGS+=-fPIC
lib: FFLAGS+=-fPIC
lib: liboPDF.so

liboPDF.so:models.o $(OBJS_COMM)
	$(CC) -shared -Wl,-soname,liboPDF.so -o libdyn.so $^ $(LDFLAGS)
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
.PHONY : clean depend distclean FORCE

FORCE:

synccosma:
	rsync -avz $(shell pwd)/../ jvbq85@login.cosma.dur.ac.uk:data/DynDistr/code

synccosmalocal:
	rsync -e "ssh -p 4800" -avz $(shell pwd)/../ jvbq85@localhost:data/DynDistr/code

syncuv2:
	rsync -e "ssh -p 4702" -avz $(shell pwd)/../ jxhan@localhost:data/DynDistr/code
	
syncbright:
	rsync -e "ssh -p 4700" -avz $(shell pwd)/../ jxhan@localhost:DynDistr/code

module:
	module add intel_comp/2012.0.032
	module add hdf5/intel_2012.0.032/1.8.9

depend:
	makedepend --$(CFLAGS)-- -Y $(SRC_COMM) $(SRC_MAIN) $(SRC)
	
clean:
	rm -f $(OBJS) $(OBJS_COMM) $(OBJS_MAIN)
	
distclean: clean
	rm -f $(EXEC) libdyn.so *~ *.pyc py/*.pyc
#-----end--- auto dependencies below---------	
# DO NOT DELETE

hdf_util.o: hdf_util.h mymath.h
cosmology.o: cosmology.h globals.h
mymath.o: mymath.h
models.o: mymath.h globals.h cosmology.h tracer.h halo.h models.h
tracer.o: hdf_util.h tracer.h halo.h models.h
halo.o: globals.h cosmology.h halo.h template.h nfw.h
template.o: globals.h cosmology.h halo.h template.h TemplateData.h
nfw.o: globals.h cosmology.h halo.h tracer.h
scan.o: mymath.h cosmology.h globals.h models.h
NumericalConvergence.o: mymath.h cosmology.h globals.h models.h
saveHDFstars.o: mymath.h hdf_util.h
mock_stars.o: mymath.h cosmology.h globals.h tracer.h halo.h models.h
