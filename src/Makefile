PROG = gfn0
MYLIB = libgfn0.a
OBJDIR = build
#--------------------------------------------------------------------------
 OSTYPE=LINUXI
 #MKLROOT=/opt/intel/mkl
#--------------------------------------------------------------------------

SRC1 += gfn0_math_wrapper.o
SRC1 += slater_module.o
SRC1 += wfn_module.o
SRC1 += cn_module.o
SRC1 += eeq_module.o
SRC1 += dftd4param.o
SRC1 += overlap_module.o

SRC2 += gfn0_types.o
SRC2 += gfn0_dftd4.o
SRC2 += gfn0_srb.o
SRC2 += gfn0_param.o
SRC2 += gfn0_basis.o
SRC2 += gfn0_qm.o
SRC2 += gfn0_prints.o
SRC2 += gfn0_module.o

OBJS := $(SRC1) $(SRC2)
#--------------------------------------------------------------------------

SUBDIRS += .

vpath % $(SUBDIRS)
OBJS := $(addprefix $(OBJDIR)/, $(OBJS))
MYLIB := $(addprefix $(OBJDIR)/, $(MYLIB))
#--------------------------------------------------------------------------

ifeq ($(OSTYPE),LINUXI)
  PREOPTS :=
  FC := ifort
  CC := icc
  LINKER = ifort -g -O -static -fopenmp  -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
  LIBS := -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a -Wl,--end-group 
  PREFLAG := -E -P
  FFLAGS := -g -O -qopenmp -align array64byte -module $(OBJDIR)
#  FFLAGS += -check all -g -debug all -traceback -fpe0
  FFLAGS += -traceback
  CCFLAGS += -D_Float128=__float128
  CCFLAGS := -O -DLINUX
endif

ifeq ($(OSTYPE),LINUXG)
  PREOPTS :=
  FC := gfortran
  CC := gcc
  LINKER := gfortran -static-libgfortran -fopenmp 
  LIBS :=
  PREFLAG :=
  FFLAGS :=  -J$(MODDIR)
  CCFLAGS := -O -DLINUX
endif


ifeq ($(OSTYPE),SOMETHINGELSE)
  PREOPTS :=
  FC := 
  CC := 
  LINKER := 
  LIBS :=
  PREFLAG := 
  FFLAGS := 
  CCFLAGS := 
endif


.PHONY: all
.PHONY: setup
.PHONY: clean
.PHONY: mklib

all: setup
all: $(PROG)
mklib: $(MYLIB)
setup: $(OBJDIR)

#--------------------------------------------------------------------------
# Setting the rules to make the object files:
#--------------------------------------------------------------------------
# implizite Regel zur Erzeugung von *.o aus *.F ausschalten
%.o: %.F
fortran.o : $(CUDA)/src/fortran.c
	$(CC) $(CCFLAGS) -I$(CUDA)/src/ $(INC) -c  $(CUDA)/src/fortran.c

# aus *.F mache ein *.f
%.f: %.F
	@echo "making $@ from $<"
	$(FC) $(PREFLAG) $(PREOPTS) $< -o $@

# aus *.f mache ein *.o
$(OBJDIR)/%.o: %.f
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

$(OBJDIR)/%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/%.o: %.F90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

# aus *.c mache ein *.o
$(OBJDIR)/%.o: %.c
	@echo "making $@ from $<"
	$(CC) $(CCFLAGS) -c $< -o $@

# linken
$(PROG): $(OBJS)
#	$(LINKER) $(OBJS) $(LIBS) -o $(PROG)

$(MYLIB): $(OBJS)
	ar rc $(MYLIB) $(OBJS)

$(OBJDIR):
	mkdir -p $@

#aufraeumen
clean:
	rm -f $(filter %.o,$(OBJS)) $(PROG) 
	rm -f $(OBJDIR)/*.mod *.mod $(PROG) 
	rm -f $(patsubst %.F, %.f, $(wildcard *.F))
