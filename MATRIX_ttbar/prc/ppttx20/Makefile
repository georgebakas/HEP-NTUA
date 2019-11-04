$(OBJDIR)/ppttx20/ppttx20.%.o: $(PRCDIR)/ppttx20/ppttx20.%.cpp
	$(CXX) $(CXXFLAGS) -I$(SRCDIR)/include -I$(SRCDIR)/header -I$(PRCDIR)/ppttx20 -I$(AMPDIR)/pptt40 -c -o $@ $<

ppttx20SPECIFYUSERMODIFYFILES = \
	$(PRCDIR)/ppttx20/user/specify.particles.cxx \
	$(PRCDIR)/ppttx20/user/specify.cuts.cxx \
	$(PRCDIR)/ppttx20/user/specify.cuts.test.cxx \
	$(PRCDIR)/ppttx20/user/specify.prepare.scales.cxx \
	$(PRCDIR)/ppttx20/user/specify.scales.cxx \
	$(PRCDIR)/ppttx20/user/specify.moments.cxx \
	$(PRCDIR)/ppttx20/user/testpoint.born.cxx \
	$(PRCDIR)/ppttx20/user/testpoint.collinear.cxx \
	$(PRCDIR)/ppttx20/user/testpoint.real.cxx \
	$(PRCDIR)/ppttx20/user/testpoint.realcollinear.cxx \
	$(PRCDIR)/ppttx20/user/testpoint.doublereal.cxx
$(OBJDIR)/ppttx20/ppttx20.user.modify.o: $(PRCDIR)/ppttx20/ppttx20.user.modify.cpp $(ppttx20SPECIFYUSERMODIFYFILES)
	$(CXX) $(CXXFLAGS) -I$(SRCDIR)/include -I$(SRCDIR)/header -I$(PRCDIR)/ppttx20 -c -o $@ $<

$(OBJDIR)/ppttx20/ppttx20.phasespace.%.o: $(PSGDIR)/ppttx20/ppttx20.phasespace.%.cpp
	$(CXX) $(CXXFLAGS) -I$(SRCDIR)/include -I$(SRCDIR)/header -I$(PRCDIR)/ppttx20 -c -o $@ $<

Oppttx20phasespaceborn = \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.born.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.born.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.born.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.born.gg_ttx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.born.gg_ttx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.born.gg_ttx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.born.ddx_ttx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.born.ddx_ttx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.born.ddx_ttx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.born.uux_ttx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.born.uux_ttx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.born.uux_ttx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.born.bbx_ttx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.born.bbx_ttx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.born.bbx_ttx.o 

Oppttx20phasespacereal = \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.real.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.real.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.real.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.real.gg_ttxg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.real.gg_ttxg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.real.gg_ttxg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.real.gd_ttxd.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.real.gd_ttxd.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.real.gd_ttxd.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.real.gu_ttxu.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.real.gu_ttxu.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.real.gu_ttxu.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.real.gb_ttxb.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.real.gb_ttxb.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.real.gb_ttxb.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.real.gdx_ttxdx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.real.gdx_ttxdx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.real.gdx_ttxdx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.real.gux_ttxux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.real.gux_ttxux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.real.gux_ttxux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.real.gbx_ttxbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.real.gbx_ttxbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.real.gbx_ttxbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.real.ddx_ttxg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.real.ddx_ttxg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.real.ddx_ttxg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.real.uux_ttxg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.real.uux_ttxg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.real.uux_ttxg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.real.bbx_ttxg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.real.bbx_ttxg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.real.bbx_ttxg.o 

Oppttx20phasespacedoublereal = \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.gg_ttxgg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.gg_ttxgg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.gg_ttxgg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.gg_ttxddx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.gg_ttxddx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.gg_ttxddx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.gg_ttxuux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.gg_ttxuux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.gg_ttxuux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.gg_ttxbbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.gg_ttxbbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.gg_ttxbbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.gd_ttxgd.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.gd_ttxgd.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.gd_ttxgd.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.gu_ttxgu.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.gu_ttxgu.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.gu_ttxgu.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.gb_ttxgb.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.gb_ttxgb.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.gb_ttxgb.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.gdx_ttxgdx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.gdx_ttxgdx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.gdx_ttxgdx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.gux_ttxgux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.gux_ttxgux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.gux_ttxgux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.gbx_ttxgbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.gbx_ttxgbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.gbx_ttxgbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.dd_ttxdd.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.dd_ttxdd.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.dd_ttxdd.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.du_ttxdu.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.du_ttxdu.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.du_ttxdu.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.ds_ttxds.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.ds_ttxds.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.ds_ttxds.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.dc_ttxdc.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.dc_ttxdc.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.dc_ttxdc.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.db_ttxdb.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.db_ttxdb.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.db_ttxdb.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.ddx_ttxgg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.ddx_ttxgg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.ddx_ttxgg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.ddx_ttxddx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.ddx_ttxddx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.ddx_ttxddx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.ddx_ttxuux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.ddx_ttxuux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.ddx_ttxuux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.ddx_ttxssx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.ddx_ttxssx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.ddx_ttxssx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.ddx_ttxccx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.ddx_ttxccx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.ddx_ttxccx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.ddx_ttxbbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.ddx_ttxbbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.ddx_ttxbbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.dux_ttxdux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.dux_ttxdux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.dux_ttxdux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.dsx_ttxdsx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.dsx_ttxdsx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.dsx_ttxdsx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.dcx_ttxdcx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.dcx_ttxdcx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.dcx_ttxdcx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.dbx_ttxdbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.dbx_ttxdbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.dbx_ttxdbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.uu_ttxuu.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.uu_ttxuu.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.uu_ttxuu.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.uc_ttxuc.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.uc_ttxuc.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.uc_ttxuc.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.ub_ttxub.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.ub_ttxub.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.ub_ttxub.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.udx_ttxudx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.udx_ttxudx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.udx_ttxudx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.uux_ttxgg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.uux_ttxgg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.uux_ttxgg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.uux_ttxddx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.uux_ttxddx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.uux_ttxddx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.uux_ttxuux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.uux_ttxuux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.uux_ttxuux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.uux_ttxssx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.uux_ttxssx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.uux_ttxssx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.uux_ttxccx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.uux_ttxccx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.uux_ttxccx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.uux_ttxbbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.uux_ttxbbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.uux_ttxbbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.usx_ttxusx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.usx_ttxusx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.usx_ttxusx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.ucx_ttxucx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.ucx_ttxucx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.ucx_ttxucx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.ubx_ttxubx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.ubx_ttxubx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.ubx_ttxubx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.bb_ttxbb.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.bb_ttxbb.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.bb_ttxbb.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.bdx_ttxbdx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.bdx_ttxbdx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.bdx_ttxbdx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.bux_ttxbux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.bux_ttxbux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.bux_ttxbux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.bbx_ttxgg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.bbx_ttxgg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.bbx_ttxgg.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.bbx_ttxddx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.bbx_ttxddx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.bbx_ttxddx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.bbx_ttxuux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.bbx_ttxuux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.bbx_ttxuux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.bbx_ttxbbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.bbx_ttxbbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.bbx_ttxbbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.dxdx_ttxdxdx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.dxdx_ttxdxdx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.dxdx_ttxdxdx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.dxux_ttxdxux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.dxux_ttxdxux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.dxux_ttxdxux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.dxsx_ttxdxsx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.dxsx_ttxdxsx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.dxsx_ttxdxsx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.dxcx_ttxdxcx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.dxcx_ttxdxcx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.dxcx_ttxdxcx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.dxbx_ttxdxbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.dxbx_ttxdxbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.dxbx_ttxdxbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.uxux_ttxuxux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.uxux_ttxuxux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.uxux_ttxuxux.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.uxcx_ttxuxcx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.uxcx_ttxuxcx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.uxcx_ttxuxcx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.uxbx_ttxuxbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.uxbx_ttxuxbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.uxbx_ttxuxbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.x.doublereal.bxbx_ttxbxbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.c.doublereal.bxbx_ttxbxbx.o \
	$(OBJDIR)/ppttx20/ppttx20.phasespace.g.doublereal.bxbx_ttxbxbx.o 

Oppttx20 = \
	$(OBJDIR)/ppttx20/ppttx20.user.modify.o \
	$(OBJDIR)/ppttx20/ppttx20.subprocesses.o \
	$(OBJDIR)/ppttx20/ppttx20.process.determination.born.o \
	$(OBJDIR)/ppttx20/ppttx20.process.determination.real.o \
	$(OBJDIR)/ppttx20/ppttx20.process.determination.doublereal.o 

LIBppttx20phasespaceborn = $(LIBDIR)/ppttx20/libppttx20phasespaceborn.a
LIBppttx20phasespacereal = $(LIBDIR)/ppttx20/libppttx20phasespacereal.a
LIBppttx20phasespacedoublereal = $(LIBDIR)/ppttx20/libppttx20phasespacedoublereal.a
LIBppttx20 = $(LIBDIR)/ppttx20/libppttx20.a

$(LIBppttx20phasespaceborn): $(Oppttx20phasespaceborn)
	ar rs $@ $(Oppttx20phasespaceborn)
$(LIBppttx20phasespacereal): $(Oppttx20phasespacereal)
	ar rs $@ $(Oppttx20phasespacereal)
$(LIBppttx20phasespacedoublereal): $(Oppttx20phasespacedoublereal)
	ar rs $@ $(Oppttx20phasespacedoublereal)
$(LIBppttx20): $(Oppttx20)
	ar rs $@ $(Oppttx20)


ppttx20: $(OBJDIR)/ppttx20/ppttx20.main.o $(LIBgeneralC) $(LIBppttx20phasespaceborn) $(LIBppttx20phasespacereal) $(LIBppttx20phasespacedoublereal) $(LIBppttx20) $(LIBpptt40amplitudedoublevirtual) $(LIBttampgrid) 
	$(CXX) -v -o $(BINDIR)/ppttx20 -I$(SRCDIR)/include -I$(SRCDIR)/header -I$(PRCDIR)/ppttx20 -I$(AMPDIR)/pptt40 $(OBJDIR)/ppttx20/ppttx20.main.o $(LIBppttx20) $(LIBpptt40amplitudedoublevirtual) $(LIBttampgrid) $(LIBchaplin) $(LIBppttx20phasespaceborn) $(LIBppttx20phasespacereal) $(LIBppttx20phasespacedoublereal) $(LIBgeneralL) 
ppttx20-clean-obj-prc:
	$(RM) \
	$(OBJDIR)/ppttx20/ppttx20.main.o \
	$(Oppttx20)

ppttx20-clean-obj-psg:
	$(RM) \
	$(Oppttx20phasespaceborn) \
	$(Oppttx20phasespacereal) \
	$(Oppttx20phasespacedoublereal)

ppttx20-clean-obj-amp:
	$(RM) \
	$(Opptt40amplitudedoublevirtual)

ppttx20-clean-lib-prc:
	$(RM) \
	$(LIBppttx20)

ppttx20-clean-lib-psg:
	$(RM) \
	$(LIBppttx20phasespaceborn) \
	$(LIBppttx20phasespacereal) \
	$(LIBppttx20phasespacedoublereal)

ppttx20-clean-lib-amp:
	$(RM) \
	$(LIBpptt40amplitudedoublevirtual)

ppttx20-clean-obj-complete:
	make ppttx20-clean-obj-prc
	make ppttx20-clean-obj-psg
	make ppttx20-clean-obj-amp

ppttx20-clean-lib-complete:
	make ppttx20-clean-lib-prc
	make ppttx20-clean-lib-psg
	make ppttx20-clean-lib-amp

ppttx20-clean-prc:
	make ppttx20-clean-obj-prc
	make ppttx20-clean-lib-prc

ppttx20-clean-psg:
	make ppttx20-clean-obj-psg
	make ppttx20-clean-lib-psg

ppttx20-clean-amp:
	make ppttx20-clean-obj-amp
	make ppttx20-clean-lib-amp

ppttx20-clean-complete:
	make ppttx20-clean-obj-complete
	make ppttx20-clean-lib-complete

ppttx20-distclean:
	make clean
	make ppttx20-clean-complete
