#include <vector>
#include <valarray>
#include <algorithm>
#include <stdexcept>
#include <netcdf.h>
#include "qsc.hpp"

using namespace qsc;

#ifdef SINGLE
#define nc_get_var_qscfloat nc_get_var_float
#else
#define nc_get_var_qscfloat nc_get_var_double
#endif

namespace qsc {
  
  /** A class to streamline the process of reading a NetCDF file.
   */
  class NetCDFReader {
  private:
    int ncid, ndims, nvars, ngatts, unlimdimid;
    static void ERR(int);
    
  public:
    NetCDFReader(std::string);
    
    // Scalars:
    void get(std::string, int&);
    void get(std::string, qscfloat&);
    // Vectors
    void get(std::string, Vector&);
	     
    void close();
  };
}

qsc::NetCDFReader::NetCDFReader(std::string filename) {
  int retval;
  if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid)))
    ERR(retval);
  // Get the number of dimensions, variables, etc:
  if ((retval = nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid)))
      ERR(retval);
}

void qsc::NetCDFReader::ERR(int e) {
  throw std::runtime_error(nc_strerror(e));
}

void qsc::NetCDFReader::get(std::string varname, int& var) {
  int var_id, retval;
  if ((retval = nc_inq_varid(ncid, varname.c_str(), &var_id)))
      ERR(retval);
  if ((retval = nc_get_var_int(ncid, var_id, &var)))
      ERR(retval);
}

void qsc::NetCDFReader::get(std::string varname, qscfloat& var) {
  int var_id, retval;
  if ((retval = nc_inq_varid(ncid, varname.c_str(), &var_id)))
      ERR(retval);
  if ((retval = nc_get_var_qscfloat(ncid, var_id, &var)))
      ERR(retval);
}

void qsc::NetCDFReader::get(std::string varname, Vector& var) {
  int var_id, retval;
  if ((retval = nc_inq_varid(ncid, varname.c_str(), &var_id)))
      ERR(retval);
  if ((retval = nc_get_var_qscfloat(ncid, var_id, &var[0])))
      ERR(retval);
}

void qsc::NetCDFReader::close() {
  int retval;
  if ((retval = nc_close(ncid))) ERR(retval);
}

//////////////////////////////////////////////////////////////////////
// End of definition of the qsc::NetCDFReader class.
//////////////////////////////////////////////////////////////////////

/** Read in quantities for a Qsc object from a NetCDF file. This
    function is used for testing.
 */
void Qsc::read_netcdf(std::string filename, char C_or_F) {
  if (C_or_F != 'C' && C_or_F != 'F') throw std::runtime_error("C_or_F must be C or F");
  int fortran = (C_or_F == 'F');
  
  if (verbose > 0) std::cout << "About to try reading netcdf file " << filename << std::endl;
  qsc::NetCDFReader nc(filename);

  if (fortran) {
    nc.get("N_phi", nphi);
  } else {
    nc.get("nphi", nphi);
  }
  allocate();

  // Scalars
  nc.get("nfp", nfp);
  nc.get("eta_bar", eta_bar);
  nc.get("B2c", B2c);
  nc.get("B2s", B2s);
  nc.get("p2", p2);
  nc.get("d_phi", d_phi);
  nc.get("B0", B0);
  nc.get("G0", G0);
  nc.get("sG", sG);
  nc.get("spsi", spsi);
  nc.get("axis_length", axis_length);
  nc.get("d_l_d_varphi", d_l_d_varphi);
  nc.get("B0_over_abs_G0", B0_over_abs_G0);
  nc.get("abs_G0_over_B0", abs_G0_over_B0);
  nc.get("rms_curvature", rms_curvature);
  nc.get("mean_elongation", mean_elongation);
  nc.get("mean_R", mean_R);
  nc.get("mean_Z", mean_Z);
  nc.get("standard_deviation_of_R", standard_deviation_of_R);
  nc.get("standard_deviation_of_Z", standard_deviation_of_Z);
  nc.get("iota", iota);
  if (fortran) {
    nc.get("axis_helicity", helicity);
    nc.get("sigma_initial", sigma0);
  } else {
    // Quantities that are not in the fortran code
    nc.get("helicity", helicity);
    nc.get("I2", I2);
    nc.get("sigma0", sigma0);
    nc.get("grid_max_curvature", grid_max_curvature);
    nc.get("grid_max_elongation", grid_max_elongation);
    nc.get("grid_min_R0", grid_min_R0);
    nc.get("max_newton_iterations", max_newton_iterations);
    nc.get("max_linesearch_iterations", max_linesearch_iterations);
    nc.get("newton_tolerance", newton_tolerance);
    nc.get("iota_N", iota_N);
  }
  // nc.get("", );

  // Vectors
  nc.get("phi", phi);
  nc.get("curvature", curvature);
  nc.get("torsion", torsion);
  nc.get("sigma", sigma);
  nc.get("X1c", X1c);
  nc.get("Y1c", Y1c);
  nc.get("Y1s", Y1s);
  nc.get("R0", R0);
  nc.get("Z0", Z0);
  nc.get("d_l_d_phi", d_l_d_phi);
  nc.get("elongation", elongation);
  nc.get("Boozer_toroidal_angle", Boozer_toroidal_angle);
  if (fortran) {
    nc.get("modBinv_sqrt_half_grad_B_colon_grad_B", L_grad_B_inverse);
    L_grad_B = ((qscfloat)1.0) / L_grad_B_inverse;
  } else {
    nc.get("d2_l_d_phi2", d2_l_d_phi2);
    nc.get("L_grad_B", L_grad_B);
    nc.get("L_grad_B_inverse", L_grad_B_inverse);
    nc.get("d_X1c_d_varphi", d_X1c_d_varphi);
    nc.get("d_Y1c_d_varphi", d_Y1c_d_varphi);
    nc.get("d_Y1s_d_varphi", d_Y1s_d_varphi);
  }

  nc.close();
  
}