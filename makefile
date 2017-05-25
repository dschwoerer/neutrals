
BOUT_TOP	= ../..

SOURCEH		= neutrals.hxx neutrals_diffusion.hxx
SOURCEC		= test.cxx neutrals.cxx neutrals_diffusion.cxx cross_section.cxx radiation.cxx helper.cxx
TARGET		= test
CXXFLAGS = -g
include $(BOUT_TOP)/make.config
