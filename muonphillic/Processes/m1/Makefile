CALCHEP:=/home/duncan/UChicago/MuonBeam/calchep_3.8.10
WORK:=/home/duncan/UChicago/MuonBeam/muonphillic
include $(CALCHEP)/FlagsForMake
include EXTLIBmake
ifneq ($(SONAME),)
   soname_VandP := $(SONAME)../../../VandP.so
else
   soname_VandP :=
endif

ifneq ($(LHAPDFPATH),)
  EXTLIB := -L$(LHAPDFPATH) -lLHAPDF $(EXTLIB)
endif

VandP.so:VandP.c $(DEPEND) models/extlib1.mdl
	$(CC) $(CFLAGS)  -shared -o VandP.so $(soname_VandP) VandP.c $(CALCHEP)/include/VandPgate.c  $(EXTLIB) $(CALCHEP)/lib/dynamic_me.a  $(CALCHEP)/lib/libSLHAplus.a $(CALCHEP)/lib/ntools.a $(CALCHEP)/lib/serv.a $(CALCHEP)/lib/dummy.a  -lm $(lQuad) -lpthread
