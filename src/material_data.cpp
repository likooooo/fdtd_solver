
#include "material_data.hpp"

#include <algorithm>

#include "meep/mympi.hpp"

namespace meep_geom {

bool transition::operator==(const transition &other) const {
  return (from_level == other.from_level && to_level == other.to_level &&
          transition_rate == other.transition_rate && frequency == other.frequency &&
          vector3_equal(sigma_diag, other.sigma_diag) && gamma == other.gamma &&
          pumping_rate == other.pumping_rate);
}

bool transition::operator!=(const transition &other) const { return !(*this == other); }

medium_struct::medium_struct(double epsilon)
    : epsilon_diag{epsilon, epsilon, epsilon}, epsilon_offdiag{}, mu_diag{1, 1, 1}, mu_offdiag{},
      E_susceptibilities(), H_susceptibilities(), E_chi2_diag{}, E_chi3_diag{}, H_chi2_diag{},
      H_chi3_diag{}, D_conductivity_diag{}, B_conductivity_diag{} {}

void medium_struct::check_offdiag_im_zero_or_abort() const {
  if (epsilon_offdiag.x.im != 0 || epsilon_offdiag.y.im != 0 || epsilon_offdiag.z.im != 0 ||
      mu_offdiag.x.im != 0 || mu_offdiag.y.im != 0 || mu_offdiag.z.im != 0) {
    meep::abort("Found non-zero imaginary part of epsilon or mu offdiag.\n");
  }
}

material_data::material_data()
    : which_subclass(MEDIUM), medium(), user_func(NULL), user_data(NULL), do_averaging(false),
      epsilon_data(NULL), epsilon_dims{}, grid_size{}, weights(NULL), medium_1(),
      medium_2(), material_grid_kinds{U_DEFAULT} {}

void material_data::copy_from(const material_data &from) {
  which_subclass = from.which_subclass;
  medium = from.medium;

  user_func = from.user_func;
  // NOTE: the user_data field here opaque/void - so this is the best we can do.
  user_data = from.user_data;
  do_averaging = from.do_averaging;

  std::copy(std::begin(from.epsilon_dims), std::end(from.epsilon_dims), std::begin(epsilon_dims));
  if (from.epsilon_data) {
    size_t N = from.epsilon_dims[0] * from.epsilon_dims[1] * from.epsilon_dims[2];
    epsilon_data = new double[N];
    std::copy_n(from.epsilon_data, N, epsilon_data);
  }

  grid_size = from.grid_size;
  if (from.weights) {
    size_t N = from.grid_size.x * from.grid_size.y * from.grid_size.z;
    weights = new double[N];
    std::copy_n(from.weights, N, weights);
  }

  medium_1 = from.medium_1;
  medium_2 = from.medium_2;
  beta = from.beta;
  eta = from.eta;
  damping = from.damping;
  material_grid_kinds = from.material_grid_kinds;
}

material_type_list::material_type_list() : items(NULL), num_items(0) {}

} // namespace meep_geom
