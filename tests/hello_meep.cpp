
#include <meep.hpp>
#include <material_data.hpp>
#include <meepgeom.hpp>
#include <iostream>

int meep::verbosity = 1;


void config_fdtd();

int main(int argc, char** argv)
{
    meep::initialize instance(argc, argv);
    meep::volume illum_source_pos_ = meep::volume(meep::vec(0, 0, 0));            // um
    meep::volume field_monitor_pos_ = meep::volume(meep::vec(0, 0, 0));         // um
    meep::volume energy_monitor_volume_ = meep::volume(meep::vec(0, 0, 0));       // um
    std::unique_ptr<meep::boundary_region> boundary_regions_;
    geometric_object_list geo_object_list_;
    std::unique_ptr<meep::dft_fields> dft_field_; 
    std::unique_ptr<meep::fields> field_;               // without mask
    std::unique_ptr<meep::fields> filled_field_;        // with mask
    
    using OLreal = float;
    OLreal PI = 3.1415926;
    OLreal theta = 20 * PI / 180;
    OLreal phi = 30 * PI / 180;
    OLreal pol = 40 * PI / 180;

    OLreal kx = -std::sin(theta) * std::cos(phi);
    OLreal ky = -std::sin(theta) * std::sin(phi);

    OLreal Ex0 = std::cos(pol) * std::cos(theta) * std::cos(phi) - std::sin(pol) * std::sin(phi);
    OLreal Ey0 = std::cos(pol) * std::cos(theta) * std::sin(phi) + std::sin(pol) * std::cos(phi);
    OLreal Ez0 = -std::cos(pol) * std::sin(theta);
    config_fdtd();
    // using OLMedium = std::complex<double>;
    // OLMedium material1(2.25, 0);
    // WaferSlice slice1;
    // slice1.gds_layer_number = 1;
    // slice1.medium = material1;
    // slice1.slice_bottom_position = -400.0;
    // slice1.thickness = 500.0;
}
double dummy_eps(const meep::vec &v) { std::cout << v.str() << std::endl;return 1.0; }
void config_fdtd()
{
    double PML_thickness = 0.5;
    double mesh_step = 0.01;
    //-----------------Parameter definition-------------------
    double wvl = 0.365;    // um
    double dx = mesh_step;  // grid size, um
    double eps_Substrate = 3 * 3;
    double eps_SiO2 = 1.5 * 1.5;
    double eps_SiON = (std::complex<double>(2.3, 0.3) * std::complex<double>(2.3, 0.3)).real(); 

    // simulation domain size (without PML)
    double domain_size_x = 0.5;
    double domain_size_y = 1.4;
    double domain_size_z = dx;
    // simulation domain size (with PML)
    double cell_size_x = domain_size_x;
    double cell_size_y = domain_size_y + 2*PML_thickness;
    double cell_size_z = domain_size_z;

    // source size and center point position
    double source_size_x = cell_size_x;
    double source_size_y = 0;
    double source_size_z = 0;
    double source_center_x = 0;
    double source_center_y = 0;
    double source_center_z = 0;

    // monitor size and center point position [used if monitor_type = "plane"]
    double monitor_size_x = domain_size_x;
    double monitor_size_y = 0;
    double monitor_size_z = 0;
    double ref_monitor_center_x = 0;
    double ref_monitor_center_y = 0.1;
    double ref_monitor_center_z = 0;
    double trans_monitor_center_x = 0;
    double trans_monitor_center_y = 0.6;
    double trans_monitor_center_z = 0;
    double courant_factor = 0.5;

    // incident angle
    double theta_inc = 30; // degree [corresponding to -theta in Lumerical!]
    auto get_simualtion_cell_volume = [&](){
        meep::grid_volume gv = meep::vol3d(cell_size_x/dx, cell_size_y/dx, cell_size_z/dx, 1.0); 
        // gv.center_origin();
        // gv.shift_origin(meep::vec( 0,0, cell_size_z/dx/2) );
        return gv;
    };
    auto get_boundary_conditions = [&](){
        meep::boundary_region bound = 
            meep::pml(PML_thickness, meep::direction::Z,  meep::boundary_side::Low) + 
            meep::pml(PML_thickness, meep::direction::Z,  meep::boundary_side::High);
        return bound;
    };
    auto gv = get_simualtion_cell_volume();
    auto boundary_regions = get_boundary_conditions();
    meep::structure filled_structure(gv, dummy_eps, boundary_regions,  meep::identity(), 0, courant_factor);

    meep::fields(&filled_structure).output_hdf5(meep::Dielectric, gv.surroundings());

    // #-----------------Define simulation domain and resolution--------------
    // pml_layers = [mp.PML(thickness=PML_thickness, direction=mp.Y)]
    // cell_size = mp.Vector3(cell_size_x, cell_size_y, cell_size_z) 
    // resolution = 1/dx  # pixels/um

    // #-----------------Define soruce-------------
    // fcen = 1/wvl
    // source_center = mp.Vector3(source_center_x, source_center_y, source_center_z)
    // source_size = mp.Vector3(source_size_x, source_size_y, source_size_z)

    // rot_angle = np.radians(theta_inc)  # If error was induced in run phase, use math.radians()
        
    // n_source = 1
    // k = mp.Vector3(0,fcen*n_source,0).rotate(mp.Vector3(z=1), rot_angle) 
    // print(f"k={k}")
    // print(f"rot_angle={rot_angle}")
    // def pw_amp(k, x0):
    //     def _pw_amp(x):
    //         return cmath.exp(1j * 2 * math.pi * k.dot(x + x0))  # x+x0 is coorinates of mesh points on source
    //     return _pw_amp
    // sources = [mp.Source(src = mp.GaussianSource(fcen, fwidth=0.2*fcen, is_integrated=True),
    //                         center=source_center,
    //                         size=source_size,
    //                         component=mp.Ez,
    //                         amp_func=pw_amp(k, source_center))
    // ] 

    // cell geometry [with PML]
    // cell_geometry_.xmin = domain_xmin - boundary_padding_[0];
    // cell_geometry_.xmax = domain_xmax + boundary_padding_[1];
    // cell_geometry_.ymin = domain_ymin - boundary_padding_[2];
    // cell_geometry_.ymax = domain_ymax + boundary_padding_[3];
    // cell_geometry_.zmin = domain_zmin - boundary_padding_[4];
    // cell_geometry_.zmax = domain_zmax + boundary_padding_[5];
    // LOG(MID_ENGINE, LTYPE_INFO, "FDTD: domain volume WITH pml region, min: %lfum, %lfum, %lfum, max: %lfum, %lfum, %lfum.", 
    //     cell_geometry_.xmin, cell_geometry_.ymin, cell_geometry_.zmin,
    //     cell_geometry_.xmax, cell_geometry_.ymax, cell_geometry_.zmax); 

    // // used for auto shutoff control
    // energy_monitor_volume_ = meep::volume(
    //     meep::vec( cell_geometry_.xmin, 
    //             cell_geometry_.ymin, 
    //             cell_geometry_.zmin),
    //     meep::vec( cell_geometry_.xmax, 
    //             cell_geometry_.ymax, 
    //             cell_geometry_.zmax)
    // );

    // illum_source_pos_ = meep::volume(
    //     meep::vec( cell_geometry_.xmin, 
    //                cell_geometry_.ymin, 
    //                nm_to_um(nm_illum_z_)),
    //     meep::vec( cell_geometry_.xmax, 
    //                cell_geometry_.ymax, 
    //                nm_to_um(nm_illum_z_))
    //   );

    // field_monitor_pos_ = meep::volume(
    //     meep::vec(nm_field_monitor_from_.X(), 
    //               nm_field_monitor_from_.Y(),
    //               nm_field_monitor_from_.Z()),
    //     meep::vec(nm_field_monitor_to_.X(), 
    //               nm_field_monitor_to_.Y(),
    //               nm_field_monitor_to_.Z())
    //   );
}