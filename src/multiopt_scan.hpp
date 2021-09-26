#ifndef QSC_MULTIOPT_SCAN_H
#define QSC_MULTIOPT_SCAN_H

#include <vector>
#include <valarray>
#include <chrono>
#include <mpi.h>
#include "multiopt.hpp"

#ifdef SINGLE
#define MPI_QSCFLOAT MPI_FLOAT
#else
#define MPI_QSCFLOAT MPI_DOUBLE
#endif

namespace qsc {

  class MultiOptScan {
  private:
    enum {ATTEMPTS,
      KEPT,
      REJECTED_DUE_TO_R0,
      REJECTED_DUE_TO_IOTA,
      REJECTED_DUE_TO_ELONGATION,
      REJECTED_DUE_TO_L_GRAD_B,
      REJECTED_DUE_TO_B20_VARIATION,
      REJECTED_DUE_TO_L_GRAD_GRAD_B,
      REJECTED_DUE_TO_D2_VOLUME_D_PSI2,
      REJECTED_DUE_TO_DMERC,
      REJECTED_DUE_TO_R_SINGULARITY,
      N_FILTERS};

    std::chrono::time_point<std::chrono::steady_clock> start_time;
    int filters_local[N_FILTERS];
    void defaults();
    
  public:
    MultiOpt mo_ref, mo;
    MPI_Comm mpi_comm;
    int mpi_rank, n_procs;
    bool proc0;
    qscfloat max_seconds, print_status_period, save_period;
    int n_scan, n_scan_all, filters[N_FILTERS];
    qscfloat filter_fractions[N_FILTERS];
    qscfloat min_R0_to_keep, min_iota_to_keep, max_elongation_to_keep;
    qscfloat min_L_grad_B_to_keep, min_L_grad_grad_B_to_keep;
    qscfloat max_B20_variation_to_keep, min_r_singularity_to_keep;
    qscfloat max_d2_volume_d_psi2_to_keep, min_DMerc_to_keep;
    bool keep_all;
    int verbose;
    std::string outfilename;
    int ndim;
    std::vector<std::string> params;
    Vector params_max, params_min;
    std::valarray<bool> params_log;
    std::valarray<int> params_n, params_stage;
    std::vector<Vector> params_vals;
    int axis_nmax_plus_1;
    bool quit_after_init;
    
    const int n_parameters_base = 38;
    const int n_int_parameters_base = 3;
    const int n_int_parameters = n_int_parameters_base + N_FILTERS;
    int n_parameters;
    Matrix parameters;
    std::valarray<int> int_parameters;
    Vector parameters_single;
    std::valarray<int> int_parameters_single;
    std::valarray<int> n_solves_kept, attempts_per_proc;
    qscfloat total_cpu_seconds;
    
    Vector scan_eta_bar, scan_sigma0, scan_B2s, scan_B2c;
    Matrix scan_R0c, scan_R0s, scan_Z0c, scan_Z0s;
    Vector scan_min_R0, scan_max_curvature;
    Vector scan_iota, scan_max_elongation;
    Vector scan_min_L_grad_B, scan_min_L_grad_grad_B;
    Vector scan_r_singularity, scan_B20_variation, scan_B20_residual;
    Vector scan_d2_volume_d_psi2, scan_DMerc_times_r2;
    Vector scan_standard_deviation_of_R, scan_standard_deviation_of_Z;
    std::valarray<int> scan_helicity;
    
    Vector scan_weight_B20, scan_weight_iota, scan_target_iota;
    Vector scan_weight_elongation, scan_weight_curvature;
    Vector scan_weight_R0, scan_target_min_R0;
    Vector scan_weight_d2_volume_d_psi2, scan_max_d2_volume_d_psi2;
    Vector scan_weight_XY2, scan_weight_XY2Prime;
    Vector scan_weight_Z2, scan_weight_Z2Prime;
    Vector scan_weight_XY3, scan_weight_XY3Prime;
    Vector scan_weight_grad_B, scan_weight_grad_grad_B, scan_weight_r_singularity;
    Vector scan_weight_axis_length, scan_weight_standard_deviation_of_R;
    
    MultiOptScan();
    void run(std::string);
    void input(std::string);
    void init();
    void scan();
    void eval_scan_index(int);
    int proc0_recv();
    void print_status();
    void filter_global_arrays();
    void write_netcdf();
  };
}

#endif
