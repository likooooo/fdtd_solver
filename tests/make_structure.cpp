#include <meep.hpp>
#include <material_data.hpp>
#include <meepgeom.hpp>
#include <memory>
#include <assert.h>
#include "material_lib.hpp"
#include <cstdarg>
#include <filesystem>

using namespace meep;
using complexnum = std::complex<realnum>;

struct material_struct
{
    meep_geom::material_data d;
    realnum z{0};
    realnum depth{0};
};
namespace meep_geom
{
    extern void GDSIIError(const char *Routine);
    struct geometric_object_list_raii final : public geometric_object_list 
    {
        geometric_object_list_raii(){
            geometric_object_list::num_items = 0;
            geometric_object_list::items = nullptr;
        };
        geometric_object_list_raii(const geometric_object_list_raii&) = delete;
        geometric_object_list_raii(geometric_object_list_raii&& other){
            this->num_items = other.num_items;
            this->items = other.items;
            other.num_items = 0;
            other.items = nullptr;
        }
        ~geometric_object_list_raii(){
            if(geometric_object_list::num_items){
                assert(nullptr != this->items);
                delete [] this->items;
                return;
            }
            assert(nullptr == this->items);
    
            // int length = geometry.num_items;
            // for (int i = 0; i < length; i++) {
            //   geometric_object_destroy(geometry.items[i]);
            // }
            // delete[] geometry.items;
        }
    };
}


struct user_config
{
    realnum freq;
    realnum PML_thickness;
    vec low, high;
    //== 传播方向
    vec k0{vec(0, 0, 1)};
    //== 震动方向
    std::array<std::complex<realnum>, 3> Exyz{0, 1, 0}; // total Y-polarization

    std::vector<material_struct> materials;
    std::vector<meep_geom::geometric_object_list_raii> geo_objects;
    void set_lambda(realnum lambda){
        freq = realnum(1)/lambda;
    }
    void set_roi(realnum x0, realnum y0, realnum x1, realnum y1)
    {
        low = vec(x0, y0, 0);
        high = vec(x1, y1, 0);
    }
    material_struct& append_material(complexnum nk, realnum depth)
    {
        materials.push_back(material_struct());
        material_struct& md = materials.back();
        if(materials.size() > 1){
            material_struct& prev = materials.rbegin()[1];
            md.z = prev.z + prev.depth;
        }
        else{
            md.z = PML_thickness;
        }
        md.depth = depth;
        md.d.which_subclass = meep_geom::material_data::MEDIUM;
        auto &m = md.d.medium;
        complexnum eps = nk * nk;
        m.epsilon_diag.x = eps.real();
        m.epsilon_diag.y = m.epsilon_diag.x;
        m.epsilon_diag.z = m.epsilon_diag.x;  
        m.D_conductivity_diag.x = 2 * M_PI * freq * eps.imag() / eps.real();
        m.D_conductivity_diag.y = m.D_conductivity_diag.x;
        m.D_conductivity_diag.z = m.D_conductivity_diag.x;
        return md;
    }
    void append_materials(complexnum* eps, realnum* depth, size_t n)
    {
        for(size_t i = 0; i < n; i++) append_material(eps[i], depth[i]);
    }
    void append_repeat_materials(complexnum* eps, realnum* depth, size_t n, size_t repeat_count){
        for(size_t i =0; i < repeat_count; i++) append_materials(eps, depth, n);
    }
    meep_geom::geometric_object_list_raii& append_geo_objects(const std::string& path, int layer_id, const material_struct& m)
    {
#ifdef HAVE_LIBGDSII
        vec center = (low + high) * realnum(0.5);
        assert(0 == center.z());
        geo_objects.push_back(meep_geom::geometric_object_list_raii());
        meep_geom::geometric_object_list_raii& geo_list = geo_objects.back();
        geo_list = meep_geom::get_GDSII_prisms(
            m.d, path.c_str(), layer_id, 
            center.x(), center.y(),
            m.z, m.z + m.depth
        );
        return geo_list;
#else
    meep_geom::GDSIIError(__PRETTY_FUNCTION__);
    return geo_objects.front();
#endif
    }
    grid_volume create_grid_volume(realnum physical_step = 0.5) const
    {
        vec physical_span = high - low; 
        physical_span.set_direction(direction::Z, materials.back().z + materials.back().depth + PML_thickness);
        return vol3d(physical_span.x(), physical_span.y(), physical_span.z(), 1/ physical_step); 
    }
    boundary_region create_boundary_conditions(){
        return meep::pml(PML_thickness, meep::direction::Z,  meep::boundary_side::Low) + 
            meep::pml(PML_thickness, meep::direction::Z,  meep::boundary_side::High);
    }
}config;

void verb_printf(const char *fmt, ...) 
{
    if(0 == meep::verbosity) return;
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stdout, fmt, ap);
    fflush(stdout);
    va_end(ap);
}

class material_struct_function : public material_function
{
    const std::vector<material_struct>* p;
    std::vector<double> eps_cutline;
    std::vector<double> conductivity_cutline;
    double step;
    double from;
    double to;
public:
    material_struct_function(const std::vector<material_struct>* material_structs, double inva) : p(material_structs), step(inva){
        make_row_data(step);
    }
    void make_row_data(double step)
    {
        if(p->empty()){
            eps_cutline = std::vector<double>();
            conductivity_cutline = std::vector<double>(); 
            return;
        }
        from = p->front().z;
        to = p->back().z + p->back().depth;
        double resolution = 1.0 / step;
        size_t n = std::ceil((to - from) * resolution);
        eps_cutline = std::vector<double>(n);
        conductivity_cutline = std::vector<double>(n);
        auto update = [&](size_t idx, const material_struct& m){
            eps_cutline.at(idx) = m.d.medium.epsilon_diag.z;
            conductivity_cutline.at(idx) = m.d.medium.D_conductivity_diag.z;
        };
        size_t i = 0;
        for(auto it = p->begin();i < (n - 1) && it != p->end(); i++){
            if(from + (i + 1) * step >= (it->z + it->depth)){
                double coef_next =  ((from + (i + 1) * step - (it->z + it->depth)))/step;
                double coef_current = 1 - coef_next;
                eps_cutline.at(i) = it->d.medium.epsilon_diag.z * coef_current;
                conductivity_cutline.at(i) = it->d.medium.D_conductivity_diag.z * coef_current;
                it++;
                eps_cutline.at(i) += it->d.medium.epsilon_diag.z * coef_next;
                conductivity_cutline.at(i) += it->d.medium.D_conductivity_diag.z * coef_next;
            }
            else{
                update(i, *it);
            }
        }
        for(; i<n; i++) update(i, p->back());
        verb_printf("make_row_data i=%ld, n=%ld\n", i, n);
    }
    bool has_conductivity(component c) override {return true;}
    double chi1p1(field_type ft, const vec &r) override
    { 
        return __set_grid_material_impl(r, eps_cutline);
    }
    double conductivity(component c, const vec &r) override
    {
        return __set_grid_material_impl(r, conductivity_cutline);
    }
private:
    double __set_grid_material_impl(const vec &r, const std::vector<double>& cutline){
        // printf("%f\t%f\t%f\n", r.x(), r.y(), r.z());
        //== PML
        if(r.z() < from) return cutline.front();
        if(r.z() >= to)  return cutline.back();
        
        double idx = (r.z() - from) / step;
        int i = std::floor(idx);
        if(i + 1 < cutline.size()){
            double coef_next = idx - i;
            double coef_current = 1 - coef_next;
            return cutline.at(i) * coef_current + cutline.at(i + 1) * coef_next;
        }
        return cutline.back();
    }
};
std::complex<double> component_fun_abs_arg(const std::complex<realnum> *fields, const vec &loc, void *data_)
{
    return std::complex(std::abs(fields[0]), std::arg(fields[0]));
}
std::complex<double> src_phase_function(const vec &p, direction dir)
{
    // printf("%f\t%f\t%f\n", p.x(), p.y(), p.z());
    double theta = 2.0 * M_PI * config.freq * 1 * config.k0.in_direction(dir);
    return std::exp(std::complex<double>(0.0, theta));
}
int meep::verbosity = 1;
int main(int argc, char** argv)
{
    constexpr int static_args_size = 4;
    assert((argc - static_args_size)%2 == 0 && argc >= (static_args_size + 2));
    double lambda = atof(argv[1]);
    config.set_lambda(lambda);
    config.PML_thickness = atof(argv[2]);
    double physical_step = atof(argv[3]);
    config.set_roi(0, 0, lambda * 2, lambda * 2);
    verb_printf("Struct\t Epsilon\t Z(from, to)\n");
    for(int i = static_args_size; i < argc; i+= 2){
        const auto& ms =config.append_material(material_lib<realnum>::get(argv[i]), atof(argv[i + 1]));
        verb_printf("%s\t(%.3f,%.3f)\t(%.2f, %.2f)\n", argv[i], ms.d.medium.epsilon_diag.x, ms.d.medium.D_conductivity_diag.x, ms.z, ms.depth);
    }
    grid_volume v = config.create_grid_volume(physical_step);
    material_struct_function func(&config.materials, physical_step);
    structure s(v, func, config.create_boundary_conditions());
    s.outdir = "/tmp";

    fields f(&s);
    bool enable_bloch = true;
    if(enable_bloch){
        double alpha = 0, beta = 0;
        double kx = config.freq * alpha;
        double ky = config.freq * beta;
        vec k(kx, ky, std::sqrt(config.freq * config.freq - kx * kx - ky * ky));
        f.use_bloch(direction::X, kx);
        f.use_bloch(direction::Y, ky);
    }

    gaussian_src_time src(config.freq, config.freq);
    realnum intergral = std::sqrt(2 * M_PI) * lambda;
    vec src_pos = (config.low + config.high)/2;
    double src_pos_z = config.materials.front().z + config.materials.front().depth / 2;
    src_pos.set_direction(direction::Z, src_pos_z);
    component src_component = Ey;
    vec from = config.low;
    from.set_direction(direction::Z, src_pos_z);
    vec to = config.high;
    to.set_direction(direction::Z, src_pos_z);
    {
        f.add_volume_source(Ey, src, volume(from, to), 
            [](const vec & p){return src_phase_function(p, direction::Z);}, config.Exyz.at(Y) / intergral
        );
    }
    std::array<component,3> k_components {meep::component::Ex, meep::component::Ey, meep::component::Ez};
    from.set_direction(direction::Z, src_pos_z + 1 * lambda);
    to.set_direction(direction::Z, src_pos_z +   1 * lambda);
    dft_fields monitor = f.add_dft_fields(
        k_components.data(), k_components.size(), 
        volume(from, to), std::vector<double>{config.freq}, true
    );
    while (f.time() < (f.last_source_time())) f.step();
    size_t dims[3];
    meep::direction directions[3];
    int rank = f.get_array_slice_dimensions(monitor.where, dims, directions, true);
    std::array<complexnum*, 3> final_exyz;
    for(size_t i = 0; i < k_components.size(); i++){
        // f.output_hdf5(k_components.at(i), v.surroundings());
        // f.output_hdf5((std::string(component_name(k_components.at(i))) + std::string("_abs_ang")).c_str(), 1, &k_components.at(i), 
        //     component_fun_abs_arg, nullptr, v.surroundings()
        // );
        final_exyz.at(i) = f.get_dft_array(monitor, k_components.at(i), 0, &rank, dims);
    }
    if(meep::verbosity){
        f.output_hdf5(Dielectric, v.surroundings());
    }
    for(size_t i = 0; i < 100; i +=10)
        printf("%ld %f, %f(degree)\n", i, std::abs(final_exyz.at(1)[i]), 180.0 / M_PI * std::arg(final_exyz.at(1)[i]));
}