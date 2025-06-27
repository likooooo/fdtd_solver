/* compile-time configuration of meep that affects user linkage */

// Enable single precision builds.
#ifndef MEEP_SINGLE
#define MEEP_SINGLE 0
#endif

// Enable SWIG threads for Python. When enabled, the GIL is released for
// all SWIG wrapped calls by defaults.
#ifndef MEEP_SWIG_PYTHON_THREADS
#define MEEP_SWIG_PYTHON_THREADS 0
#endif

// Enable debugging for SWIG wrapped calls for Python. When enabled, all
// SWIG wrapped python calls are logged.
#ifndef MEEP_SWIG_PYTHON_DEBUG
#define MEEP_SWIG_PYTHON_DEBUG 0
#endif


/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* define to enable debugging code */
/* #undef DEBUG */

/* Define to dummy `main' function (if any) required to link to the Fortran
libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to a macro mangling the given C identifier (in lower and upper
case), which must not contain underscores, for linking with Fortran. */
#define F77_FUNC(name,NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define F77_FUNC_(name,NAME) name ## _

/* Define if fenv.h declares this. */
#define HAVE_DECL_FEENABLEEXCEPT 1
/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the `feenableexcept' function. */
#define HAVE_FEENABLEEXCEPT 1

/* Define if you have LAPACK library. */
#define HAVE_LAPACK 1


/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define if you have a BLAS library. */
#define HAVE_BLAS 1

/* Define to 1 if you have the `BSDgettimeofday' function. */
/* #undef HAVE_BSDGETTIMEOFDAY */

/* Define to 1 if you have the `cblas_daxpy' function. */
#define HAVE_CBLAS_DAXPY 1

/* Define to 1 if you have the `cblas_ddot' function. */
#define HAVE_CBLAS_DDOT 1

/* If we have the ctl_printf_callback variable */
#define HAVE_CTL_PRINTF_CALLBACK 1

/* define if the compiler supports basic C++11 syntax */
#define HAVE_CXX11 1



/* If we have libGDSII::GetLayers */
#define HAVE_GDSII_GETLAYERS 1

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the <guile/gh.h> header file. */
/* #undef HAVE_GUILE_GH_H */

/* Define to 1 if you have the `H5Pset_fapl_mpio' function. */
#define HAVE_H5PSET_FAPL_MPIO 1

/* Define to 1 if you have the `H5Pset_mpi' function. */
/* #undef HAVE_H5PSET_MPI */

// /* Define if you have libharminv */
// #define HAVE_HARMINV 0

// /* Define if we have & link HDF5 */
#define HAVE_HDF5 1

/* Define to 1 if you have the <immintrin.h> header file. */
#define HAVE_IMMINTRIN_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `jn' function. */
#define HAVE_JN 1

/* Define to 1 if you have the `ctl' library (-lctl). */
/* #undef HAVE_LIBCTL */

/* If we have the libctl_quiet variable */
/* #undef HAVE_LIBCTL_QUIET */

/* Define to 1 if you have the `dfftw' library (-ldfftw). */
/* #undef HAVE_LIBDFFTW */

/* Define to 1 if you have the `dl' library (-ldl). */
/* #undef HAVE_LIBDL */

/* Define to 1 if you have the `fftw' library (-lfftw). */
/* #undef HAVE_LIBFFTW */

/* Define to 1 if you have the `fftw3' library (-lfftw3). */
#define HAVE_LIBFFTW3 1

// /* Define to 1 if you have the `GDSII' library (-lGDSII). */
// #define HAVE_LIBGDSII 0

/* Define to 1 if you have the `gsl' library (-lgsl). */
/* #undef HAVE_LIBGSL */

/* Define to 1 if you have the `gslcblas' library (-lgslcblas). */
/* #undef HAVE_LIBGSLCBLAS */

/* Define to 1 if you have the `guile' library (-lguile). */
/* #undef HAVE_LIBGUILE */

/* Define to 1 if you have the <libguile.h> header file. */
/* #undef HAVE_LIBGUILE_H */

/* Define to 1 if you have the `guile-ltdl' library (-lguile-ltdl). */
/* #undef HAVE_LIBGUILE_LTDL */

/* Define to 1 if you have the `hdf5' library (-lhdf5). */
#define HAVE_LIBHDF5 1

/* Define to 1 if you have the `ltdl' library (-lltdl). */
/* #undef HAVE_LIBLTDL */

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the `readline' library (-lreadline). */
/* #undef HAVE_LIBREADLINE */

/* Define to 1 if you have the `z' library (-lz). */
#define HAVE_LIBZ 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define if you have libmpb */
// #define HAVE_MPB 1

/* If we have the mpb_printf_callback variable */
#define HAVE_MPB_PRINTF_CALLBACK 1

/* Define if you have the MPI library. */
/* #undef HAVE_MPI */

/* Define to enable OpenMP */
/* #undef HAVE_OPENMP */

/* Define to 1 if you have the `scm_make_smob_type' function. */
/* #undef HAVE_SCM_MAKE_SMOB_TYPE */

/* define if we have SCM_NEWSMOB */
/* #undef HAVE_SCM_NEWSMOB */

/* define if we have SCM_SMOB_DATA */
/* #undef HAVE_SCM_SMOB_DATA */

/* define if we have SCM_SMOB_PREDICATE */
/* #undef HAVE_SCM_SMOB_PREDICATE */

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to catch and ignore SIGFPE signals */
/* #undef IGNORE_SIGFPE */

/* Define to the sub-directory in which libtool stores uninstalled libraries.
    */
#define LT_OBJDIR ".libs/"

/* Define if mpi.h needs SEEK macros to be undefined */
/* #undef NEED_UNDEF_SEEK_FOR_MPI */

/* Name of package */
#define PACKAGE "meep"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT ""

/* Define to the full name of this package. */
#define PACKAGE_NAME "meep"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "meep 1.0.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "meep"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.0.0"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "1.0.0"

/* define to nothing if C99 _Pragma is not supported */
/* #undef _Pragma */

/* Define to the equivalent of the C99 'restrict' keyword, or to
    nothing if this is not supported.  Do not define if restrict is
    supported directly.  */
#define restrict 
/* Work around a bug in Sun C++: it does not support _Restrict or
    __restrict__, even though the corresponding Sun C compiler ends up with
    "#define restrict _Restrict" or "#define restrict __restrict__" in the
    previous line.  Perhaps some future version of Sun C++ will work with
    restrict; if so, hopefully it defines __RESTRICT like Sun C does.  */
#if defined __SUNPRO_CC && !defined __RESTRICT
# define _Restrict
# define __restrict__
#endif

