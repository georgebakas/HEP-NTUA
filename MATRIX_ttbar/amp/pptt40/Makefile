#AMPTT_CXXFLAGS = -fPIC -O2 -I$(GINAC_DIR)/include -I$(CLN_DIR)/include
#AMPTT_CXXFLAGS += -I$(shell $(LHAPDF_CONFIG) --incdir)
#AMPTT_CXXFLAGS +=  $(shell $(GSL_CONFIG) --cflags)

#AMPTT_SOFLAGS = -shared -L$(CHAPLIN_DIR)/lib -lchaplin
AMPTT_SOFLAGS = -shared #-Wl,-rpath,$(CHAPLIN_DIR)/lib -L$(CHAPLIN_DIR)/lib -lchaplin

$(OBJDIR)/pptt40/pptt40.%.o: $(AMPDIR)/pptt40/pptt40.%.cpp
	$(CXX) $(CXXFLAGS) $(AMPTT_CXXFLAGS) -I$(SRCDIR)/include -I$(SRCDIR)/header -I$(AMPDIR)/pptt40 -c -o $@ $<

Opptt40amplitudedoublevirtual = \
	$(OBJDIR)/pptt40/pptt40.amplitude.doublevirtual.o

LIBpptt40amplitudedoublevirtual = $(LIBDIR)/pptt40/libpptt40amplitudedoublevirtual.a

$(LIBpptt40amplitudedoublevirtual): $(Opptt40amplitudedoublevirtual)
	ar rs $@ $(Opptt40amplitudedoublevirtual)


$(OBJDIR)/pptt40/softfunction/%.o: $(AMPDIR)/pptt40/softfunction/%.f
	$(FC) -fPIC -O2 $(AMPTT_SOFLAGS) -c -o $@ $<

Ottampgrid = \
	$(OBJDIR)/pptt40/softfunction/virtgg.o \
	$(OBJDIR)/pptt40/softfunction/h_function.o \
	$(OBJDIR)/pptt40/softfunction/gghe.o \
	$(OBJDIR)/pptt40/softfunction/ggth.o \
	$(OBJDIR)/pptt40/softfunction/virtqq.o \
	$(OBJDIR)/pptt40/softfunction/qqhe.o \
	$(OBJDIR)/pptt40/softfunction/qqth.o \
	$(OBJDIR)/pptt40/softfunction/splines.o \
	$(OBJDIR)/pptt40/softfunction/polylogs.o \
	$(OBJDIR)/pptt40/softfunction/bispev.o \
	$(OBJDIR)/pptt40/softfunction/fpback.o \
	$(OBJDIR)/pptt40/softfunction/fpdisc.o \
	$(OBJDIR)/pptt40/softfunction/fprank.o \
	$(OBJDIR)/pptt40/softfunction/fpsurf.o \
	$(OBJDIR)/pptt40/softfunction/polylogs_cplx.o \
	$(OBJDIR)/pptt40/softfunction/T3T4.o \
	$(OBJDIR)/pptt40/softfunction/Fex1.o \
	$(OBJDIR)/pptt40/softfunction/fpbisp.o \
	$(OBJDIR)/pptt40/softfunction/fpgivs.o \
	$(OBJDIR)/pptt40/softfunction/fprati.o \
	$(OBJDIR)/pptt40/softfunction/sq_Fex1.o \
	$(OBJDIR)/pptt40/softfunction/TITj.o \
	$(OBJDIR)/pptt40/softfunction/Fex1_sq.o \
	$(OBJDIR)/pptt40/softfunction/fpbspl.o \
	$(OBJDIR)/pptt40/softfunction/fporde.o \
	$(OBJDIR)/pptt40/softfunction/fprota.o \
	$(OBJDIR)/pptt40/softfunction/surfit.o

LIBttampgrid = $(LIBDIR)/pptt40/libttampgrid.so
#LIBttampgrid = $(LIBDIR)/pptt40/libttampgrid.a

$(LIBttampgrid): $(Ottampgrid)
	$(FC) -o $@ $+ $(AMPTT_SOFLAGS)
#	ar rs $@ $(Ottampgrid) 

pptt40-clean:
	$(RM) \
	$(Ottampgrid) \
	$(LIBttampgrid)

.PHONY: clean
