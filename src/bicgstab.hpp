
#ifndef BICGSTAB_H
#define BICGSTAB_H

#include "meep.hpp"

namespace meep {

/* signature of a user callback function, passed to bicgstabL, that computes y = Ax.  data is
   an opaque pointer to any user data that you need, passed through from the Adata pointer to
   bicgstabL (e.g. this could be a data structure that records information about A, such as
   its size n). */
typedef void (*bicgstab_op)(const realnum *x, realnum *y, void *data);

/* solve Ax=b by biCGSTAB(L).   L is the dimension of the subspace (typically < 10) — larger L
   gives more expensive iterations but faster convergence; n is the dimension of the vectors and
   matrices (x and b are length n, A is n by n); x is the initial solution guess on input (often
   zero) and the approximate solution on output; b is the right hand side; A is a function that
   implements A-times-vector (see above) and Adata is an opaque user data pointer passed through to
   A; tol is the stopping tolerance |Ax-b| < tol*|b|; on *iters is the maximum allowed number of
   iterations and on output it is the actual number; quiet=1 to suppress status messages during the
   solver.

   You must also allocate a workspace, usually by calling the routine twice: on the first call, pass
   NULL for the work pointer, and the return value nwork is the required size of the workspace. Then
   allocate work = new realnum[nwork], and call it again.

   For non-NULL nwork, returns 0 on success, 1 if the maximum number of iterations was reached, and
   -1 if a breakdown in convergence was detected. */
ptrdiff_t bicgstabL(const int L, const size_t n, realnum *x, bicgstab_op A, void *Adata,
                    const realnum *b, const double tol,
                    int *iters,    // input *iters = max iters, output = actual iters
                    realnum *work, // if you pass work=NULL, bicgstab returns nwork
                    const bool quiet);

} // namespace meep

#endif /* BICGSTAB_H */
