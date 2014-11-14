#-----------setup server env ------------------
# Platform=Medivh
Platform=COSMA
#Platform=Bright
ifeq ($(shell uname -n), Medivh)
	Platform=Medivh
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
	#HDFINC+=-DHDF_V16 #uncomment this if you have HDF version 1.6 or less
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

ifeq ($(DEBUG),on)
CFLAGS+= -g -Wall
LDFLAGS+= -g
endif

VPATH=C $(DOC_SRC_DIR)
#-----File Dependencies----------------------
SRC_COMM = hdf_util.c globals.c cosmology.c mymath.c models.c tracer.c halo.c template.c nfw.c
OBJS_COMM= $(patsubst %.f90,%.f.o,$(SRC_COMM:%.c=%.o))
DOC= html latexpdf

OPDF_LIB=liboPDF.so
#-----targets and common rules--------------------------------
default: lib

lib: $(OPDF_LIB)


doc: $(DOC)

#the default rule will handle the rest
# %.o : %.c
# 	$(CC) $< $(CFLAGS) -c -o $@
# 
%.f.o : %.f90
	$(FC) $< $(FFLAGS) -c -o $@
	
$(OPDF_LIB): CFLAGS+=-fPIC
$(OPDF_LIB): FFLAGS+=-fPIC
$(OPDF_LIB):$(OBJS_COMM)
	$(CC) -shared -Wl,-soname,$(OPDF_LIB) -o $(OPDF_LIB) $^ $(LDFLAGS)
#------documentation----------------------

DOC_DIR=Doc
DOC_SRC_DIR=$(DOC_DIR)/source
$(DOC): doc.rst tutorial.rst api.rst $(OPDF_LIB)
	$(MAKE) $@ -C $(DOC_DIR)

$(DOC_SRC_DIR)/tutorial.rst: tutorial.ipynb docclean
	ipython nbconvert --to rst $^
	mv tutorial.rst tutorial_files $(DOC_SRC_DIR)/

html_tutorial: $(DOC_DIR)/tutorial.html
$(DOC_DIR)/tutorial.html: tutorial.ipynb
	ipython nbconvert --to html $^
	mv tutorial.html $@

docclean:
	cd $(DOC_SRC_DIR) && rm -rf tutorial.rst tutorial_files
#-----Other stuff----------------------------
.PHONY : clean depend distclean docclean

depend:
	makedepend --$(CFLAGS)-- -Y $(SRC_COMM) $(SRC_MAIN) $(SRC)
	
clean:
	rm -f $(OBJS) $(OBJS_COMM) $(OBJS_MAIN)
	
distclean: clean docclean
	rm -f $(EXEC) $(OPDF_LIB) *~ *.pyc py/*.pyc
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
