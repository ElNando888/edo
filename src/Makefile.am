AM_CXXFLAGS = -I$(top_srcdir)/src
AM_CXXFLAGS += $(OPENMP_CFLAGS)
AM_CXXFLAGS += $(VRNA_CFLAGS)

bin_PROGRAMS = edo
edo_SOURCES = edo.cpp edo_fold.h edo_fold.cpp

edo_LDADD = $(AM_LDFLAGS)
edo_LDADD += $(VRNA_LIBS)

