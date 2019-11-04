# determines wether LHAPDF 5 or 6 are used. This is needed to know in the C++ code (via the set preprocessor) whether to set the correct N_decoupled values in LHAPDF 5.
LHAPDF_VERSION = $(shell $(LHAPDF_CONFIG) --version | head -c 1)
ifeq ($(LHAPDF_VERSION), 5)
  $(info IMPORTANT NOTE: Using $$LHAPDF_VERSION = $(LHAPDF_VERSION) and adding preprocessor variable -DLHAPDF5 to compiler Flags to correctly choose N_decoupled in certain PDF sets.)
  CXXFLAGS += -DLHAPDF5
endif

LIBLHAPDF = -Wl,-rpath,$(shell $(LHAPDF_CONFIG) --libdir) -L$(shell $(LHAPDF_CONFIG) --libdir) -lLHAPDF $(LIBMORE)
CXXFLAGS += -I$(shell $(LHAPDF_CONFIG) --incdir)

LIBGSL = -Wl,-rpath,$(shell $(GSL_CONFIG) --prefix)/lib $(shell $(GSL_CONFIG) --libs) -lgsl
CXXFLAGS +=  $(shell $(GSL_CONFIG) --cflags)

LIBOpenLoops = -Wl,-rpath,$(shell $(OpenLoops_CONFIG) --libdir) -L$(shell $(OpenLoops_CONFIG) --libdir) $(shell $(OpenLoops_CONFIG) --ldflags)
CXXFLAGS += -I$(shell $(OpenLoops_CONFIG) --incdir)

#LIBGiNaC = -L$(GINAC_DIR)/lib -L$(CLN_DIR)/lib -lginac -lcln
LIBGiNaC = -L$(GINAC_DIR)/lib -L$(CLN_DIR)/lib -lginac -lcln -lquadmath
#LIBGiNaC = -L$(GINAC_CLNDIR)/lib -lginac -lcln -lquadmath
LDFLAGS += -L$(GINAC_DIR)/lib -Wl,-rpath,$(GINAC_DIR)/lib -lginac 
LDFLAGS += -L$(CLN_DIR)/lib -Wl,-rpath,$(CLN_DIR)/lib -lcln 

LIBchaplin = -Wl,-rpath,$(CHAPLIN_DIR)/lib64 -L$(CHAPLIN_DIR)/lib64 -lchaplin

CXXFLAGS += -I$(SRCDIR)/include -I$(SRCDIR)/header -I$(SRCDIR)/classes/header -I$(SRCDIR)/phasespace 
CXXFLAGS += -I$(SRCDIR)/NJsubtraction
CXXFLAGS += -I$(SRCDIR)/qTsubtraction/header


GENERICFILES = \
	$(SRCDIR)/include/const.cxx \
	$(SRCDIR)/include/definitions.cxx \
	$(SRCDIR)/include/definitions.observable.set.cxx \
	$(SRCDIR)/include/definitions.phasespace.set.cxx
#$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(GENERICFILES)
#	$(CXX) $(CXXFLAGS) -c -o $@ $< 


$(OBJDIR)/routines/%.o: $(SRCDIR)/routines/%.cpp
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

OJroutines = \
	$(OBJDIR)/routines/routines.generic.o \
	$(OBJDIR)/routines/routines.output.o \
	$(OBJDIR)/routines/routines.particle.code.o \
	$(OBJDIR)/routines/routines.distributions.o \
	$(OBJDIR)/routines/oldversion.routines.pdf.o 

LIBroutines = $(LIBDIR)/libroutines.a
$(LIBroutines): $(OJroutines)
	ar rs $@ $(OJroutines)
#$(LIBroutines): $(LIBroutines)($(OJroutines))
#	ar ruv $(LIBroutines) $(OJroutines)
#LIBroutines = $(LIBDIR)/libroutines.so
#$(LIBroutines): $(OJroutines)
#	$(CXX) -o $(LIBroutines) $(CXXFLAGS) $(LINKFLAGS) $(OJroutines)


$(OBJDIR)/classes/%.o: $(SRCDIR)/classes/%.cpp
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

OJclasses = \
	$(OBJDIR)/classes/contribution.set.o \
	$(OBJDIR)/classes/define.particle.set.o \
	$(OBJDIR)/classes/event.set.o \
	$(OBJDIR)/classes/fiducialcut.o \
	$(OBJDIR)/classes/inputparameter.set.o \
	$(OBJDIR)/classes/model.set.o \
	$(OBJDIR)/classes/user.defined.o \
	$(OBJDIR)/classes/observable.set.o \
	$(OBJDIR)/NJsubtraction/observable.NJsubtraction.definition.o \
	$(OBJDIR)/classes/observable.pspweight.o \
	$(OBJDIR)/classes/observable.pdf.o \
	$(OBJDIR)/classes/observable.event.selection.o \
	$(OBJDIR)/classes/observable.ncollinear.o \
	$(OBJDIR)/classes/observable.initialization.o \
	$(OBJDIR)/classes/observable.initialization.scale.o \
	$(OBJDIR)/classes/observable.initialization.contribution.o \
	$(OBJDIR)/classes/observable.initialization.OpenLoops.o \
	$(OBJDIR)/classes/observable.distribution.o \
	$(OBJDIR)/classes/observable.object.o \
	$(OBJDIR)/classes/observable.output.o \
	$(OBJDIR)/classes/observable.proceeding.o \
	$(OBJDIR)/classes/observable.ME2.o \
	$(OBJDIR)/classes/runresumption.set.o \
	$(OBJDIR)/classes/fourvector.o \
	$(OBJDIR)/classes/particle.o \
	$(OBJDIR)/classes/xdistribution.o \
	$(OBJDIR)/classes/dddistribution.o \
	$(OBJDIR)/classes/logger.o \
	$(OBJDIR)/classes/dipole.set.o \
	$(OBJDIR)/classes/ioperator.set.o \
	$(OBJDIR)/classes/collinear.set.o \
	$(OBJDIR)/classes/multicollinear.set.o \
	$(OBJDIR)/classes/correlationoperator.set.o \
	$(OBJDIR)/classes/call.generic.o \
	$(OBJDIR)/classes/topwidth.o \
	$(OBJDIR)/classes/oldversion.observable.initialization.contribution.o

LIBclasses = $(LIBDIR)/libclasses.a
$(LIBclasses): $(OJclasses)
	ar rs $@ $(OJclasses)
#$(LIBclasses): $(LIBclasses)($(OJclasses))
#	ar ruv $(LIBclasses) $(OJclasses)
#LIBclasses = $(LIBDIR)/libclasses.so
#$(LIBclasses): $(OJclasses)
#	$(CXX) -o $(LIBclasses) $(CXXFLAGS) $(LINKFLAGS) $(OJclasses)


$(OBJDIR)/summary/%.o: $(SRCDIR)/summary/%.cpp
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

OJsummary = \
	$(OBJDIR)/summary/summary.generic.o \
	$(OBJDIR)/summary/summary.generic.output.distribution.o \
	$(OBJDIR)/summary/summary.order.o \
	$(OBJDIR)/summary/summary.order.result.o \
	$(OBJDIR)/summary/summary.order.distribution.o \
	$(OBJDIR)/summary/summary.order.output.result.o \
	$(OBJDIR)/summary/summary.order.output.distribution.o \
	$(OBJDIR)/summary/summary.list.o \
	$(OBJDIR)/summary/summary.list.result.o \
	$(OBJDIR)/summary/summary.list.distribution.o \
	$(OBJDIR)/summary/summary.list.output.result.o \
	$(OBJDIR)/summary/summary.list.output.distribution.o \
	$(OBJDIR)/summary/summary.contribution.o \
	$(OBJDIR)/summary/summary.contribution.result.o \
	$(OBJDIR)/summary/summary.contribution.distribution.o \
	$(OBJDIR)/summary/summary.contribution.output.result.o \
	$(OBJDIR)/summary/summary.contribution.output.distribution.o \
	$(OBJDIR)/summary/summary.subprocess.o \
	$(OBJDIR)/summary/summary.subprocess.result.o \
	$(OBJDIR)/summary/summary.subprocess.distribution.o

LIBsummary = $(LIBDIR)/libsummary.a
$(LIBsummary): $(OJsummary)
	ar rs $@ $(OJsummary)
#$(LIBsummary): $(LIBsummary)($(OJsummary))
#	ar ruv $(LIBsummary) $(OJsummary)
#LIBsummary = $(LIBDIR)/libsummary.so
#$(LIBsummary): $(OJsummary)
#	$(CXX) -o $(LIBsummary) $(CXXFLAGS) $(LINKFLAGS) $(OJsummary)


$(OBJDIR)/munich/%.o: $(SRCDIR)/munich/%.cpp
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

OJmunich = \
	$(OBJDIR)/munich/munich.o \
	$(OBJDIR)/munich/munich.integration.born.o \
	$(OBJDIR)/munich/munich.integration.VA.QCD.o \
	$(OBJDIR)/munich/munich.integration.CA.QCD.o \
	$(OBJDIR)/munich/munich.integration.RA.QCD.o \
	$(OBJDIR)/munich/munich.integration.VA.QEW.o \
	$(OBJDIR)/munich/munich.integration.CA.QEW.o \
	$(OBJDIR)/munich/munich.integration.RA.QEW.o \
	$(OBJDIR)/munich/munich.integration.VA.MIX.o \
	$(OBJDIR)/munich/munich.integration.RA.MIX.o \
	$(OBJDIR)/munich/munich.integration.VT.QCD.o \
	$(OBJDIR)/munich/munich.integration.CT.QCD.o \
	$(OBJDIR)/munich/munich.integration.VT2.QCD.o \
	$(OBJDIR)/munich/munich.integration.CT2.QCD.o \
	$(OBJDIR)/munich/munich.integration.VJ.QCD.o \
	$(OBJDIR)/munich/munich.integration.CJ.QCD.o \
	$(OBJDIR)/munich/munich.errorhandling.o

LIBmunich = $(LIBDIR)/libmunich.a
$(LIBmunich): $(OJmunich)
	ar rs $@ $(OJmunich)
#$(LIBmunich): $(LIBmunich)($(OJmunich))
#	ar ruv $(LIBmunich) $(OJmunich)
#LIBmunich = $(LIBDIR)/libmunich.so
#$(LIBmunich): $(OJmunich)
#	$(CXX) -o $(LIBmunich) $(CXXFLAGS) $(LINKFLAGS) $(OJmunich)


$(OBJDIR)/phasespace/%.o: $(SRCDIR)/phasespace/%.cpp
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

OJphasespace = \
	$(OBJDIR)/phasespace/phasespace.set.o \
	$(OBJDIR)/phasespace/phasespace.initialization.o \
	$(OBJDIR)/phasespace/phasespace.mapping.o \
	$(OBJDIR)/phasespace/phasespace.mapping.initial.o \
	$(OBJDIR)/phasespace/phasespace.optimization.o \
	$(OBJDIR)/phasespace/phasespace.correction.o \
	$(OBJDIR)/phasespace/phasespace.initial.o \
	$(OBJDIR)/phasespace/phasespace.generator.o \
	$(OBJDIR)/phasespace/multichannel.set.o \
	$(OBJDIR)/phasespace/importancesampling.set.o \
	$(OBJDIR)/phasespace/randomvariable.o \
	$(OBJDIR)/phasespace/randommanager.o \
	$(OBJDIR)/phasespace/oldversion.phasespace.set.o \
	$(OBJDIR)/phasespace/oldversion.phasespace.initialization.o

LIBphasespace = $(LIBDIR)/libphasespace.a
$(LIBphasespace): $(OJphasespace)
	ar rs $@ $(OJphasespace)
#$(LIBphasespace): $(LIBphasespace)($(OJphasespace))
#	ar ruv $(LIBphasespace) $(OJphasespace)
#LIBphasespace = $(LIBDIR)/libphasespace.so
#$(LIBphasespace): $(OJphasespace)
#	$(CXX) -o $(LIBphasespace) $(CXXFLAGS) $(LINKFLAGS) $(OJphasespace)


$(OBJDIR)/dipolesubtraction/%.o: $(SRCDIR)/dipolesubtraction/%.cpp
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

OJdipolesubtraction = \
	$(OBJDIR)/dipolesubtraction/CS.phasespace.dipoles.o \
	$(OBJDIR)/dipolesubtraction/CDST.phasespace.dipoles.o \
	$(OBJDIR)/dipolesubtraction/CS.phasespace.group.dipoles.o \
	$(OBJDIR)/dipolesubtraction/CDST.phasespace.group.dipoles.o \
	$(OBJDIR)/dipolesubtraction/QEW.CS.dipoles.o \
	$(OBJDIR)/dipolesubtraction/QCD.CS.dipoles.o \
	$(OBJDIR)/dipolesubtraction/dipolesubtraction.ME2.VA.o \
	$(OBJDIR)/dipolesubtraction/dipolesubtraction.ME2.CA.o \
	$(OBJDIR)/dipolesubtraction/dipolesubtraction.ME2.RA.o \
	$(OBJDIR)/dipolesubtraction/splitting.o \
	$(OBJDIR)/dipolesubtraction/observable.dipolesubtraction.VA.o \
	$(OBJDIR)/dipolesubtraction/observable.dipolesubtraction.CA.o \
	$(OBJDIR)/dipolesubtraction/observable.dipolesubtraction.RA.o 

LIBdipolesubtraction = $(LIBDIR)/libdipolesubtraction.a
$(LIBdipolesubtraction): $(OJdipolesubtraction)
	ar rs $@ $(OJdipolesubtraction)
#$(LIBdipolesubtraction): $(LIBdipolesubtraction)($(OJdipolesubtraction))
#	ar ruv $(LIBdipolesubtraction) $(OJdipolesubtraction)
#LIBdipolesubtraction = $(LIBDIR)/libdipolesubtraction.so
#$(LIBdipolesubtraction): $(OJdipolesubtraction)
#	$(CXX) -o $(LIBdipolesubtraction) $(CXXFLAGS) $(LINKFLAGS) $(OJdipolesubtraction)


$(OBJDIR)/qTsubtraction/%.o: $(SRCDIR)/qTsubtraction/%.cpp
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

$(OBJDIR)/qTsubtraction/%.o : $(SRCDIR)/qTsubtraction/%.f
	$(F90) -c -o $@ $(FFLAGS) -J $(SRCDIR)/qTsubtraction $<

$(OBJDIR)/qTsubtraction/external/%.o : $(SRCDIR)/qTsubtraction/external/%.f
	$(F90) -c -o $@ $(FFLAGS) -J $(SRCDIR)/qTsubtraction/external $<

OJqTsubtraction = \
	$(OBJDIR)/qTsubtraction/qTsubtraction.ME2.VT.o \
	$(OBJDIR)/qTsubtraction/external/testcall.o \
	$(OBJDIR)/qTsubtraction/external/testcall_gamma1.o \
	$(OBJDIR)/qTsubtraction/external/setup.o \
	$(OBJDIR)/qTsubtraction/external/qqb_QQb.o \
	$(OBJDIR)/qTsubtraction/external/FinVirt.o \
	$(OBJDIR)/qTsubtraction/external/born_amp_project.o \
	$(OBJDIR)/qTsubtraction/external/dot.o \
	$(OBJDIR)/qTsubtraction/external/dotem.o \
	$(OBJDIR)/qTsubtraction/external/myli2.o \
	$(OBJDIR)/qTsubtraction/external/myli3new.o \
	$(OBJDIR)/qTsubtraction/external/ptyrap.o \
	$(OBJDIR)/qTsubtraction/specialfunctions.o \
	$(OBJDIR)/qTsubtraction/auxiliaryfunctions.o \
	$(OBJDIR)/qTsubtraction/spinorproducts.o \
	$(OBJDIR)/qTsubtraction/hpl.o \
	$(OBJDIR)/qTsubtraction/hplog.o \
	$(OBJDIR)/qTsubtraction/tdhpl.o \
	$(OBJDIR)/qTsubtraction/observable.qTsubtraction.o \
	$(OBJDIR)/qTsubtraction/observable.qTsubtraction.coefficients.o \
	$(OBJDIR)/qTsubtraction/observable.qTsubtraction.massive.o \
	$(OBJDIR)/qTsubtraction/observable.qTsubtraction.CT.o \
	$(OBJDIR)/qTsubtraction/observable.qTsubtraction.CT2.o \
	$(OBJDIR)/qTsubtraction/observable.qTsubtraction.VT.o \
	$(OBJDIR)/qTsubtraction/observable.qTsubtraction.VT2.o \
	$(OBJDIR)/qTsubtraction/more.o \
	$(OBJDIR)/qTsubtraction/oldversion.counterterm_qqbar.o \
	$(OBJDIR)/qTsubtraction/oldversion.counterterm_gg.o \
	$(OBJDIR)/qTsubtraction/oldversion.observable.qTsubtraction.o \
	$(OBJDIR)/qTsubtraction/oldversion.observable.qTsubtraction.CT.o \
	$(OBJDIR)/qTsubtraction/oldversion.observable.qTsubtraction.CT2.o \
	$(OBJDIR)/qTsubtraction/oldversion.observable.qTsubtraction.VT.o \
	$(OBJDIR)/qTsubtraction/oldversion.observable.qTsubtraction.VT2.o

LIBqTsubtraction = $(LIBDIR)/libqTsubtraction.a
$(LIBqTsubtraction): $(OJqTsubtraction)
	ar rs $@ $(OJqTsubtraction)
#$(LIBqTsubtraction): $(LIBqTsubtraction)($(OJqTsubtraction))
#	ar ruv $(LIBqTsubtraction) $(OJqTsubtraction)
#LIBqTsubtractionsubtraction = $(LIBDIR)/libqTsubtraction.so
#$(LIBqTsubtractionsubtraction): $(OJqTsubtractionsubtraction)
#	$(CXX) -o $(LIBqTsubtractionsubtraction) $(CXXFLAGS) $(LINKFLAGS) $(OJqTsubtractionsubtraction)


$(OBJDIR)/NJsubtraction/%.o: $(SRCDIR)/NJsubtraction/%.cpp
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

OJNJsubtraction = \
	$(OBJDIR)/NJsubtraction/NJsubtraction.ME2.VJ.o \
	$(OBJDIR)/NJsubtraction/observable.NJsubtraction.CJ.o \
	$(OBJDIR)/NJsubtraction/observable.NJsubtraction.VJ.o \
	$(OBJDIR)/NJsubtraction/NJsubtraction.o

#	$(OBJDIR)/NJsubtraction/observable.NJsubtraction.o \
	$(OBJDIR)/NJsubtraction/observable.NJsubtraction.VJ2.o \
	$(OBJDIR)/NJsubtraction/NJsubtraction.o

LIBNJsubtraction = $(LIBDIR)/libNJsubtraction.a
$(LIBNJsubtraction): $(OJNJsubtraction)
	ar rs $@ $(OJNJsubtraction)


$(OBJDIR)/varia/%.o: $(SRCDIR)/varia/%.cpp
	$(CXX) -c -o $@ $^ $(CXXFLAGS)
#$(OBJDIR)/varia/%.o : $(SRCDIR)/varia/%.f90
#	$(F90) -c -o $@ $(FFLAGS) -J$(SRCDIR)/varia -I$(shell $(OpenLoops_CONFIG) --moddir) $<

OJvaria = \
	$(OBJDIR)/varia/determination.ckm.factor.o \
	$(OBJDIR)/varia/event.selection.o \
	$(OBJDIR)/varia/photon.algorithm.o \
	$(OBJDIR)/varia/jet.algorithm.o

LIBvaria = $(LIBDIR)/libvaria.a
$(LIBvaria): $(OJvaria)
	ar rs $@ $(OJvaria)
#$(LIBvaria): $(LIBvaria)($(OJvaria))
#	ar ruv $(LIBvaria) $(OJvaria)
#LIBvaria = $(LIBDIR)/libvaria.so
#$(LIBvaria): $(OJvaria)
#	$(CXX) -o $(LIBvaria) $(CXXFLAGS) $(LINKFLAGS) $(OJvaria)

LIBgeneralC = $(LIBroutines) $(LIBmunich) $(LIBtools) $(LIBphasespace) $(LIBdipolesubtraction) $(LIBqTsubtraction) $(LIBNJsubtraction) $(LIBclasses) $(LIBsummary) $(LIBvaria)

LIBgeneralL = $(LIBmunich) $(LIBdipolesubtraction) $(LIBqTsubtraction) $(LIBNJsubtraction) $(LIBphasespace) $(LIBsummary) $(LIBclasses) $(LIBroutines) $(LIBtools) $(LIBvaria) $(LIBOpenLoops) $(LIBLHAPDF) $(LIBGSL) $(STDLIBS) $(LDFLAGS)


lib: $(LIBgeneralC)
	$(CXX) -v -c -I$(SRCDIR)/include -I$(SRCDIR)/header $(LIBgeneralL) $(LIBgeneralL)

# amplitude code
-include $(AMPDIR)/Makefile

# process code
-include $(PRCDIR)/Makefile



$(OBJDIR)/tools/integrate.topwidth.o: $(SRCDIR)/tools/integrate.topwidth.cpp
	$(CXX) $(CXXFLAGS) -I$(SRCDIR)/include -I$(SRCDIR)/header -c -o $@ $<

integrate.topwidth: $(OBJDIR)/tools/integrate.topwidth.o $(LIBroutines) $(LIBtools) $(LIBphasespace) $(LIBdipolesubtraction) $(LIBclasses) $(LIBvaria) $(LIBppeexjjphasespaceborn) $(LIBppeexjjphasespacereal) $(LIBppeexjj) 
	$(CXX) -v -o bin/tools/topwidth -I$(SRCDIR)/include -I$(SRCDIR)/header $(OBJDIR)/tools/integrate.topwidth.o $(LIBppeexjj) $(LIBppeexjjphasespaceborn) $(LIBppeexjjphasespacereal) $(LIBroutines) $(LIBtools) $(LIBphasespace) $(LIBdipolesubtraction) $(LIBclasses) $(LIBvaria) $(LIBOpenLoops) $(LIBLHAPDF) $(LIBGSL) $(STDLIBS) $(LDFLAGS)

$(OBJDIR)/OLtest/OL%.o: $(SRCDIR)/OLtest/OL%.cpp
	$(CXX) $(CXXFLAGS) -I$(AMPDIR)/pptt40 -c -o $@ $<

OLtest: $(OBJDIR)/OLtest/OL_cpp.o $(LIBgeneralC) $(LIBppttx20phasespaceborn) $(LIBppttx20phasespacereal) $(LIBppttx20phasespacedoublereal) $(LIBppttx20) $(LIBpptt40amplitudedoublevirtual) 
	$(CXX) -v -o OLtest -I$(SRCDIR)/include -I$(SRCDIR)/header -I$(SRCDIR)/classes/header -I$(SRCDIR)/phasespace -I$(PRCDIR)/ppttx20 -I$(AMPDIR)/pptt40 $(OBJDIR)/OLtest/OL_cpp.o $(LIBppttx20) $(LIBpptt40amplitudedoublevirtual) $(LIBppttx20phasespaceborn) $(LIBppttx20phasespacereal) $(LIBppttx20phasespacedoublereal) $(LIBgeneralL)

clean-obj:
	$(RM) \
	$(OJclasses) \
	$(OJmunich) \
	$(OJsummary) \
	$(OJroutines) \
	$(OJphasespace) \
	$(OJdipolesubtraction) \
	$(OJqTsubtraction) \
	$(OJNJsubtraction) \
	$(OJvaria)
obj-clean: clean-obj

clean-lib:
	$(RM) \
	$(LIBclasses) \
	$(LIBmunich) \
	$(LIBsummary) \
	$(LIBroutines) \
	$(LIBphasespace) \
	$(LIBdipolesubtraction) \
	$(LIBqTsubtraction) \
	$(LIBNJsubtraction) \
	$(LIBvaria)
lib-clean: clean-lib


clean:
	make obj-clean
	make lib-clean
