#ifndef ADJUST_VERBOSITY_H
#define ADJUST_VERBOSITY_H

// NOTE: This header assumes it has been #included *after* headers that declare
// verbosity, HAVE_MPB, and mpb_verbosity.

namespace meep {

// This is a RAII (Resource Allocation Is Initialization) class which adjusts
// the mpb_verbosity level when created, and restores it when deleted.
class adjust_mpb_verbosity {
public:
  adjust_mpb_verbosity() {
#if defined(HAVE_MPB) &&                                                                           \
    (MPB_VERSION_MAJOR > 1 || (MPB_VERSION_MAJOR == 1 && MPB_VERSION_MINOR >= 11))
    old_level = mpb_verbosity;
    mpb_verbosity = verbosity - 1;
    if (mpb_verbosity < 0) mpb_verbosity = 0;
    if (mpb_verbosity > 3) mpb_verbosity = 3;
#else
    // avoid warnings
    (void)old_level;
#endif
  }

  ~adjust_mpb_verbosity() {
#if defined(HAVE_MPB) &&                                                                           \
    (MPB_VERSION_MAJOR > 1 || (MPB_VERSION_MAJOR == 1 && MPB_VERSION_MINOR >= 11))
    mpb_verbosity = old_level;
#endif
  }

private:
  int old_level;
};

} // namespace meep

#endif // ADJUST_VERBOSITY_H
