
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "meep.hpp"
#include "meep_internals.hpp"
// #include "config.h"

using namespace std;

namespace meep {

void fields::update_pols(field_type ft) {
  for (int i = 0; i < num_chunks; i++)
    if (chunks[i]->is_mine())
      if (chunks[i]->update_pols(ft)) {
        chunk_connections_valid = false;
        assert(changed_materials);
      }
}

bool fields_chunk::update_pols(field_type ft) {
  bool allocated_fields = false;

  realnum *w[NUM_FIELD_COMPONENTS][2];
  FOR_COMPONENTS(c) DOCMP2 { w[c][cmp] = f_w[c][cmp] ? f_w[c][cmp] : f[c][cmp]; }

  for (polarization_state *p = pol[ft]; p; p = p->next) {

    // Lazily allocate internal polarization data:
    if (!p->data) {
      p->data = p->s->new_internal_data(f, gv);
      if (p->data) {
        p->s->init_internal_data(f, dt, gv, p->data);
        allocated_fields = true;
      }
    }

    // Finally, timestep the polarizations:
    p->s->update_P(w, f_w_prev, dt, gv, p->data);
  }

  return allocated_fields;
}

} // namespace meep

