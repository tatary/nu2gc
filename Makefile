#F77 = ifort
#F77 = mpiifort
F77 = mpif90
#FFLAGS = -O3 -xP -axP -C -CA -CB -CS -CU -CV -extend_source
#FFLAGS = -fast -fPIE
#FFLAGS = -xW -no-ipo -extend_source -warn # astro6
#FFLAGS = -O1 -CB -traceback -g -ftrapuv -extend_source -warn all -check all #-i-dynamic -mcmodel=large
#FFLAGS = -extend_source -g
#FFLAGS = -O1 -xHOST -mcmodel=large -shared-intel -g -traceback -CB -ftrapuv -extend_source -warn all
#FFLAGS = -O1 -xHOST -mcmodel=large -shared-intel
#FFLAGS = -fast -extend_source -warn all
#FFLAGS = -O3 -lm -fbacktrace -Wall -m64 -lmpi -fimplicit-none -fbounds-check 
#FFLAGS = -O3 -lm -lmpi -fallow-argument-mismatch
FFLAGS = -O3 -lm -lmpich
COSMO = lcdm
#PPS = 071219
#PPS = 080507
PPS = nugc
OPT = cst2
TARGET = ${COSMO}_${OPT}_141a.out

MODULE_LAE = moduleLAE.o
LAE = LAE.o

MODULE_SFH = moduleSFH.o
SFH = SFH.o

MODULE_SN = moduleSN.o
SN = SN.o

MODULE_CLCOOL = moduleCLCOOL.o

SRCS = com.f90 \
	cos_${COSMO}_${OPT}.f90 star_formation.f90 \
	optdepth.f90 mcmc.f90 read_clcool.f90
OBJS = $(SRCS:.f90=.o)
MAIN = main_${PPS}.o
MODULE = global_var.o


# Primary Target
.PHONY: all
all: $(MODULE) $(MODULE_LAE) $(MODULE_SFH) $(MODULE_SN) $(MODULE_CLCOOL) $(OBJS) $(MAIN) $(LAE) $(SFH) $(SN) $(TARGET) clean_unnec


# Rules of Making Targets
$(MODULE): $(MODULE:.o=.f90) $(MODULE_LAE:.o=.f90) $(MODULE_SFH:.o=.f90) $(MODULE_SN:.o=.f90) $(MODULE_CLCOOL:.o=.f90) Makefile
	$(F77) $(FFLAGS) -c $<

$(MODULE_LAE): $(MODULE_LAE:.o=.f90) $(MODULE)

$(MODULE_SFH): $(MODULE_SFH:.o=.f90) $(MODULE)

$(MODULE_SN): $(MODULE_SN:.o=.f90) $(MODULE)

$(MODULE_CLCOOL): $(MODULE_CLCOOL:.o=.f90) $(MODULE)

$(MAIN): $(MAIN:.o=.f90) out_norm_l9a.f90 out_norm_l9b.f90 $(MODULE)
	$(F77) $(FFLAGS) -c $<

$(LAE): $(MODULE_LAE)

$(SFH): $(MODULE_SFH)

$(SN): $(MODULE_SN)

$(TARGET): $(OBJS) $(MAIN) $(MODULE) $(LAE)  $(MODULE_LAE) $(SFH) $(MODULE_SFH) $(SN) $(MODULE_SN)
	$(F77) $(FFLAGS) $(OBJS) $(MAIN) $(MODULE) $(LAE) $(MODULE_LAE) $(SFH) $(MODULE_SFH) $(SN) $(MODULE_SN) $(MODULE_CLCOOL) -o $@


# File Dependence
%.o: %.f90 $(MODULE)
	$(F77) $(FFLAGS) -c $<


# Target for Removing Files
.PHONY: clean
clean:
	rm -f $(TARGET) $(OBJS) $(MODULE) $(MODULE_LAE) $(LAE) $(SFH) $(MODULE_SFH) $(SN) $(MODULE_SN) $(MODULE_CLCOOL) core* *~ *.o *.mod

.PHONY: clean_unnec
clean_unnec:
	rm -f *__genmod*
