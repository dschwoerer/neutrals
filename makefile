
BOUT_TOP ?= ../..

SOURCEH	= neutrals.hxx neutrals_diffusion.hxx
SOURCEC	= neutrals.cxx neutrals_diffusion.cxx cross_section.cxx \
	  radiation.cxx neutrals_parallel.cxx helper.cxx \
	  radiation_factory.cxx

ifdef MODULE_DIR
SUB_NAME        = neutrals
TARGET		= sub
else
SOURCEC	       += calc.cxx
TARGET		= calc

all: rates.pdf

rates.pdf: rates.gpl rates.dat
	gpl rates.gpl
	cp rates.pdf /home/dave/Documents/phd/slides/bremen17/graphic/neut_rates.pdf

rates.dat: calc
	./calc > rates.dat
endif
CXXFLAGS = -g
include $(BOUT_TOP)/make.config

FORCE:
git_version.hxx: FORCE
	@sh gen_version_header.sh

neutrals.o: git_version.hxx neutrals.cxx neutrals.hxx

lib:
	@$(AR) $(ARFLAGS) $(SOUREC%.cxx=.o) libneutrals.o

generated=cross_section_factory.cxx

%xx: %xx.in.py
	@echo "  Generating $@"
	python3 $< > $@

doc:
	./doc.py mem > members.rst
	./doc.py tree > tree.rst
	./doc.py def > defaults.rst
