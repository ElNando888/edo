AM_CXXFLAGS = -I$(top_srcdir)/src
AM_CXXFLAGS += $(OPENMP_CFLAGS)
AM_CXXFLAGS += $(VRNA_CFLAGS)

GENGETOPT_CMDL = edo_cmdl.c edo_cmdl.h
GENGETOPT_FILES = edo.ggo

EXTRA_DIST = $(GENGETOPT_FILES) $(GENGETOPT_CMDL)

bin_PROGRAMS = edo

edo_LDADD = $(AM_LDFLAGS)
edo_LDADD += $(VRNA_LIBS)


edo_cmdl.h edo_cmdl.c: edo.ggo
	gengetopt -i edo.ggo

edo_SOURCES = edo.cpp edo_fold.h edo_fold.cpp edo_cmdl.h edo_cmdl.c
