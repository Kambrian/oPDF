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
	MPICC:=mpicc -cc=$(CC)
	ROOTDIR='"/data/raid3/jxhan/Lensing"'
endif

ifeq ($(Platform), Medivh)
	HDFLIB= -lhdf5_hl -lhdf5
	HDFINC= -I/usr/include/mpi
	GSLLIB=-lgsl -lgslcblas
	CC = icc
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

CFLAGS += $(HDFINC) $(GSLINC) -DROOTDIR=$(ROOTDIR)
LDFLAGS += $(HDFLIB) $(GSLLIB)

#~ ifeq ($(USER),kam)      #my laptop has hdf v1.6
ifeq ($(shell uname -n),kam-laptop)
CFLAGS+= -DHDF_V16
endif

ifeq ($(DEBUG),on)
CFLAGS+= -g -Wall
LDFLAGS+= -g
endif

#-----File Dependencies----------------------
SRC_MAIN=models.c
EXEC=$(basename $(SRC_MAIN))
OBJS_MAIN=$(addsuffix .o, $(EXEC))
SRC_COMM = io.c hdf_util.c cosmology.c mymath.c 
OBJS_COMM = $(SRC_COMM:%.c=%.o)
#any additional OBJS below:
#SRC = 
#OBJS= $(SRC:%.c=%.o)

#-----targets and common rules--------------------------------
default: models
all: $(EXEC)

% : %.o $(OBJS_COMM)
	$(CC) $^ $(LDFLAGS) -o $@

%.o : %.c
	$(CC) $< $(CFLAGS) -c -o $@

#Additional Dependencies

#Additional Flags
models lib: CFLAGS+=$(OMPLIB)	
	LDFLAGS+=$(OMPLIB)
#mpi flags
#gama_WL_rand: CC:=$(MPICC)

#Force recompile
#shear.o: FORCE

lib: CFLAGS+=-fPIC
lib: libdyn.so

libdyn.so:models.o $(OBJS_COMM)
	$(CC) -shared -Wl,-soname,libdyn.so -o libdyn.so $^ $(LDFLAGS)
	
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
	rm -f $(EXEC) libdyn.so *~ *.pyc
#-----end--- auto dependencies below---------	
# DO NOT DELETE

io.o: mymath.h hdf_util.h io.h
hdf_util.o: hdf_util.h mymath.h
cosmology.o: mymath.h cosmology.h
mymath.o: mymath.h
models.o: mymath.h cosmology.h io.h models.h
