#-----------------------------------------------------------------------------
# Replace this with the location of LHAPDF on your system (if desired)
#LHAPDFLIB   = /home/ferrera/lhapdf-5.8.5/lib/
LHAPDFLIB   = /mnt/shared/lhapdf_install/lib
#/Users/grazzini/lhapdf/lib

DYNNLOHOME        = $(PWD)
SOURCEDIR       = $(PWD)/src
VPATH		= $(DIRS)
BIN		= $(DYNNLOHOME)/bin
INCPATH  	= $(SOURCEDIR)/Inc
OUTPUT_OPTION	= -o $(DYNNLOHOME)/obj/$@

#   NATIVE -- DYNNLO internal routines
#   LHAPDF -- Les Houches library
PDFROUTINES = LHAPDF
#PDFROUTINES = NATIVE

#FC = g77 -fno-automatic -fno-f2c -O1 -I$(INCPATH)
#FC = ifort -O3 -I$(INCPATH)
#FC = ifort -g
FC = gfortran -O2 -I$(INCPATH)
#FC=gfortran -O2 -I$(INCPATH) `lhapdf-config --ldflags` -lstdc++

DIRS	=	$(DYNNLOHOME):\
		$(DYNNLOHOME)/obj:\
		$(SOURCEDIR)/User:$(SOURCEDIR)/Matrix:\
		$(SOURCEDIR)/Need:$(SOURCEDIR)/Phase:\
		$(SOURCEDIR)/Integrate:\
               
# -----------------------------------------------------------------------------
# Specify the object files. 


MATRIXFILES = \
hjetfill.o \
qqb_w.o \
qqb_w_g.o \
qqb_z.o \
qqb_z_g.o \
qqb_w2jet.o \
qqb_w1jet_gs.o \
qqb_w_gvec.o \
qqb_z2jet.o \
qqb_z1jet_gs.o \
qqb_z1jet.o \
qqb_z_gvec.o \
aqqb_zbb.o \
ampqqb_qqb.o \
w2jetsq.o \
z2jetsq.o \
subqcd.o \
storecsz.o \
qqb_w1jet_v.o \
qqb_z1jet_z.o \
qqb_z1jet_v.o \
qqb_w1jet_z.o \
virt5.o \
A5NLO.o \
A51.o \
A52.o \
i3m.o \
qqb_QQb.o \
qqb_QQb_g.o \



INTEGRATEFILES = \
vegas.o \
mbook.o \
ran0.o \
ran1.o \




NEEDFILES = \
boostx.o \
boost.o \
branch.o \
ckmfill.o \
coupling.o \
couplz.o \
dipoles.o \
dipoles_fac.o \
dipolesub.o \
dot.o \
dotem.o \
getptildejet.o \
gtperp.o \
higgsp.o \
higgsw.o \
histofinLH.o \
includedipole.o \
lowintHstnew.o \
lowint_incldip.o \
masscuts.o \
dynnlo.o \
dyinit.o \
hexit.o \
integrate.o \
ptyrap.o \
r.o \
setup.o \
realint.o \
realvirt2.o \
count_ttb.o \
myli2.o \
myli3new.o \
sethparams.o \
smalls.o \
spinork.o \
spinoru.o \
storedip.o \
storeptilde.o \
strcat.o \
swapjet.o \
transform.o \
virtint.o \
writeinfo.o \
zeromsq.o \
besselkfast.o \
integration.o \
isolation.o \
ddilog.o \
newton1.o \
lnrat.o \
lfunctions.o \
dclaus.o \
a06.o \
a09.o \
scaleset.o \
FinVirt.o \
born_amp_project.o \

PHASEFILES = \
gen2.o \
gen3.o \
gen4MIO.o \
gen4h.o \
genBORN2.o \
genBORN4.o \
phase3.o \
phase4.o \
phi1_2new.o \
phi1_2m.o \
phi3m0.o \
breitw.o \
gen6hp.o \
phase6h.o \
phi1_2h.o \
genCT2.o \
genCT2m.o \
genBORN2m.o \
phase3m.o \
phi3m.o \
gen3m.o \
gen2m.o \



USERFILES = \
deltarj.o \
genclust2.o \
genclust_kt.o \
genclust_cone.o \
cuts.o \
getet.o \
mdata.o \
miscclust.o \
plotter.o \

LIBDIR=.
LIBFLAGS=

ifeq ($(PDFROUTINES),LHAPDF)
   PDFFILES += \
   pdf_lhapdf.o \
   pdfset_lhapdf.o
   NEEDFILES += \
   alfamz_lhapdf.o
   LIBDIR += -L$(LHAPDFLIB)
   LIBFLAGS += -lLHAPDF
   PDFMSG='   ----> DYNNLO compiled with LHAPDF routines <----'
else
ifeq ($(PDFROUTINES),NATIVE)
   PDFFILES += \
   pdf.o \
   pdfset.o
   NEEDFILES += \
   alfamz.o
   PDFMSG='   ----> DYNNLO compiled with its own PDFs <----'
endif
endif

CODE = $(NEEDFILES) $(PHASEFILES) \
          $(USERFILES) $(MATRIXFILES) $(INTEGRATEFILES) $(PDFFILES) \
          

dynnlo: $(CODE)
	$(FC) $(FFLAGS) -L$(LIBDIR) -o $@ \
	$(patsubst %,obj/%,$(CODE)) $(LIBFLAGS)
	mv dynnlo bin/
	@echo $(PDFMSG)

# -----------------------------------------------------------------------------


