# compiler flags
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/matlab $(matlab_CPPFLAGS)

# matlab wrapper directory
rnfftmatlabdir = $(datadir)/nfft/matlab/rnfft

# library
lib_LTLIBRARIES = librnfft.la
librnfft_la_SOURCES = rnfftmex.c

librnfft_la_LIBADD = $(top_builddir)/libnfft3@PREC_SUFFIX@_matlab.la @matlab_fftw3_LIBS@ $(top_builddir)/matlab/libmatlab.la $(matlab_LIBS)

librnfft_la_LDFLAGS = -no-undefined -module -shared -shrext $(matlab_mexext) -avoid-version @matlab_fftw3_LDFLAGS@ $(matlab_LDFLAGS)

if HAVE_MATLAB_THREADS
librnfft_la_CFLAGS = $(OPENMP_CFLAGS)
endif


dist_rnfftmatlab_DATA = rnfftmex.m

# target all-am builds .libs/librnfft@matlab_mexext@
rnfftmex@matlab_mexext@: all-am
	rm -f ./rnfftmex@matlab_mexext@
	$(LN_S) .libs/librnfft@matlab_mexext@ ./rnfftmex@matlab_mexext@

all: all-am rnfftmex@matlab_mexext@

clean: clean-am
	rm -f ./rnfftmex@matlab_mexext@

install-exec-hook:
	rm -f $(DESTDIR)$(libdir)/rnfftmex@matlab_mexext@
	$(LN_S) $(DESTDIR)$(libdir)/librnfft@matlab_mexext@ $(DESTDIR)$(libdir)/rnfftmex@matlab_mexext@

uninstall-hook:
	rm -f $(DESTDIR)$(libdir)/rnfftmex@matlab_mexext@

