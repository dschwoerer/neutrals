
BOUT_TOP	= ../..

SOURCEH		= neutrals.hxx neutrals_diffusion.hxx
SOURCEC		= test.cxx neutrals.cxx neutrals_diffusion.cxx cross_section.cxx radiation.cxx helper.cxx
TARGET		= test
CXXFLAGS = -g
include $(BOUT_TOP)/make.config

FORCE:
git_version.hxx: FORCE
	@sh gen_version_header.sh

neutrals.o: git_version.hxx neutrals.cxx neutrals.hxx

lib:
	@$(AR) $(ARFLAGS) $(SOUREC%.cxx=.o) libneutrals.o
