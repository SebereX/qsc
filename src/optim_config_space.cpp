#include <stdio.h>
#include <cmath>
#include <chrono>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h> 
#include <gsl/gsl_min.h> 
#include <gsl/gsl_multimin.h> 
#include <gsl/gsl_vector.h>
#include "qsc.hpp"

using namespace qsc;
using namespace std::chrono;


// Function definition
double eta_res(double x, void *p);
double B20_res(double x, void *p);
int choose_eta(Qsc * q, qscfloat * etabar_out, int verbose=0);
int choose_B22c(Qsc * q, qscfloat * B22c_out, int verbose=0);
double findZOpt(const gsl_vector *z, void *p);
double choose_Z_axis(Qsc* q, double* z_vals, int verbose=0);
int find_Z_axis(Qsc* q, int nfp, int num_harm, double r0c[], int verbose=0);

int main(int argc,char **argv) {
    /* Main function that reads the specified inputs to construct 2D spaces of NAE-QS configurations
   following the prescription in arXiv:2204.10234, and outputs them onto a .txt file. 
   The code takes as an input:
    %d - number of field periods, N
    %d - number of harmonics describing the magnetic axis (>=2)
    %d %d - number of NAE configurations in each direction of the 2D space (corresponding to the
            first two harmonics)
    %f %f %f %f - R_N in (%f,%f) and R_2N in (%f,%f), limits of the 2D space (values of the coefficients)
    %f .... - R_3N, R_4N... vales of the coefficients for the higher order harmonics (whose number follows
                is given by # of harm. (2nd input) - 2). 
    %s %s - output .txt files. The first will include properties of the configurations in the space from NAE other
            than the complete list of Z/R, which is saved to the second file.              */
            
    Qsc q;
    // Read inputs to the code
    int nfp = std::stoi(argv[1]);
    int num_harm = std::stoi(argv[2]);
    int num_axes[2];
    for (int i=0;i<2;i++)
        num_axes[i] = std::stoi(argv[i+3]);
    double lim_axes[2][2];
    for (int i=0;i<2;i++){
        lim_axes[i][0] = std::stod(argv[2*i+5]);
        lim_axes[i][1] = std::stod(argv[2*i+6]);
    }
    double high_axes[num_harm-2];
    for(int i=0;i<num_harm-2;i++)
        high_axes[i] = std::stod(argv[9+i]);

    // Create 2D grid for the first two harmonics
    double a_vec[num_axes[0]];
    for (int i=0;i<num_axes[0];i++)
        a_vec[i] = lim_axes[0][0]+i*(lim_axes[0][1]-lim_axes[0][0])/(num_axes[0]+1.0);
    double b_vec[num_axes[1]];
    for (int i=0;i<num_axes[1];i++)
        b_vec[i] = lim_axes[1][0]+i*(lim_axes[1][1]-lim_axes[1][0])/(num_axes[1]+1.0);

    // Get starting timepoint
    int iter = 1;
    auto start = high_resolution_clock::now();

    // Open the files to write on
    FILE * pFile, zFile;
    pFile = fopen (argv[num_harm+7],"w");
    zFile = fopen (argv[num_harm+8],"w");
    fprintf (pFile, "a_val b_val z_a/a z_b/b hel iota eta_bar elong B2c delB20 avB20 rSing shear \n");
    for (int i=0;i<num_harm;i++)
        fprintf(zFile,"z_%d ",nfp*(i+1));
    fprintf(zFile,"\n");
    
    // Construct the optimised QS-NAE configurations over the whole space
    for(int i=0;i<num_axes[0];i++){
        for(int j=0;j<num_axes[1];j++){
            // Set the axis harmonics
            double r0c[num_harm];
            r0c[0] = a_vec[i];
            r0c[1] = b_vec[j];
            for (int k=2;k<num_harm;k++)
                r0c[k] = high_axes[k-2];

            // Optimise Z harmonics to minimise the QS deviation
            int status = find_Z_axis(&q, nfp, num_harm, r0c);
            q.calculate_shear(0);   // Calculate magnetic shear
            std::cout << iter << " " << q.B20_grid_variation << std::endl; // Output to the console iteration # and delB20

            // Output space properties to the first file
            fprintf (pFile, "%f %f %f %f %d %f %f %f %f %f %f %f %f \n",a_vec[i], b_vec[j], q.Z0s[1]/q.R0c[1], q.Z0s[2]/q.R0c[2],
                q.helicity, q.iota, q.eta_bar, q.mean_elongation, q.B2c, q.B20_grid_variation, q.B20_mean, q.r_singularity_robust, q.iota2);
            // Output Z harmonic (normalised to R harmonics) found for each point in space
            for(k=0;k<num_harm;k++)
                fprintf(zFile,"%f ",q.Z0s[k+1]/q.R0c[k+1]);
            fprintf(zFile,"\n");
            iter++;
        } 
    }
    // Close output files
    fclose (pFile);
    fclose (zFile);

    // Get ending timepoint
    auto stop = high_resolution_clock::now();
 
    // Get duration. Substart timepoints to get duration. To cast it to proper unit
    // use duration cast method
    auto duration = duration_cast<microseconds>(stop - start);
 
    std::cout << "Time taken by function: "
         << duration.count()/1000000.0 << " seconds" << std::endl;

    return 0;
}

double eta_res(double x, void *p) {
    // Function that evaluates -|iotabar_0| as a function of eta, to find its extremum
    Qsc * q = (Qsc *) p;
    (q->eta_bar) = x;
    q->init_axis();
    q->solve_sigma_equation();
    double val = -std::abs(q->iota + q->helicity * q->nfp);
    return val;
    } 

double B20_res(double x, void *p) {
    // Function that evaluates delB20 as a function of B2c, to find its minimum
    Qsc * q = (Qsc *) p;
    (q->at_least_order_r2) = true;
    (q->B2c) = x;
    q->calculate_r2();
    double val = (q->B20_grid_variation);
    return val;
}

double findZOpt(const gsl_vector *z, void *p){
    // GSL class function that returns optimised delB20 for given Z harmonic/R harmonic choice
    Qsc * q = (Qsc *) p;
    int status;
    int num_harm = (q->R0c.size())-1;
    // Put together axis harmonics
    double zs[num_harm+1];
    double x;
    for(int i=0; i<num_harm; i++){
        x = gsl_vector_get(z, i);
        (q->Z0s)[i+1] = x*(q->R0c)[i+1];
    }
    // First find etabar by extremising iotabar_0
    q->init_axis();
    qscfloat eta_bar_opt;
    status = choose_eta(q, &eta_bar_opt);
    if(status){
        return GSL_NAN; // Standard failed return to GSL optimiser
    }
    // Evaluate the first order opt. NAE
    q->eta_bar = eta_bar_opt;
    q->at_least_order_r2 = false;
    q->calculate();

    // Second order search for B2c to minimise delB20
    qscfloat B22c_opt;
    status = choose_B22c(q, &B22c_opt);
    if(status)
        return GSL_NAN; // Standard failed return to GSL optimiser
    // Evaluate the second order opt. NAE
    q->B2c = B22c_opt;
    q->at_least_order_r2 = true;
    q->calculate_r2();
    
    return q->B20_grid_variation;
}

int find_Z_axis(Qsc* q, int nfp, int num_harm, double r0_in[], int verbose) {
    // Function that places in q the Z-optimised QS NAE configuration for input R harmonics
    //  q - handle to QSC object where optimised NAE configuration is left
    //  nfp - number of field periods
    //  num_harm - number of harmonics describing the magnetic axis
    //  r0_in - R harmonic array
    //  verbose - default is set to 0, but a non-zero value prints the evolution of the optimisation and times elapsed

    // Set basic properties of QSC object
    q->verbose = verbose;
    q->nfp = nfp;
    q->nphi = 61;

    q->B0 = 1;

    // Set axis harmonics if these are provided in r0_in
    if (num_harm == 0){
        q->R0c.resize(3, 0);
        q->R0s.resize(3, 0);
        q->Z0c.resize(3, 0);
        q->Z0s.resize(3, 0);

        q->R0c[0] = 1;
        q->R0c[1] = 0.23;
        q->R0c[2] = 0.0012;

        q->Z0s[1] = 0.23;
        q->Z0s[2] = 0.0012;

        num_harm = 2;
    } else{
        q->R0c.resize(num_harm+1,1);
        q->R0s.resize(num_harm+1, 0);
        q->Z0c.resize(num_harm+1, 0);
        q->Z0s.resize(num_harm+1, 0);
        q->R0c[0] = 1;
        for (int i=0;i<num_harm;i++){
            q->R0c[i+1] = r0_in[i];
            q->Z0s[i+1] = r0_in[i];
        }
    }

    q->at_least_order_r2 = true;

    // Get starting timepoint
    auto start = high_resolution_clock::now();
    
    // Initialise QSC object
    q->allocate();
    q->init_axis();

    // Run Z harmonic optimisation
    double Z_vals[num_harm];
    choose_Z_axis(q, Z_vals,verbose);

    if (verbose!=0){
        // Get ending timepoint
        auto stop = high_resolution_clock::now();
    
        std::cout << "Z values:" ;
        for(int i=0;i<num_harm;i++)
            std::cout << Z_vals[i] << " ";
        std::cout << std::endl;
        std::printf("etabar: %f, B2c: %f, delB20: %f\n",q->eta_bar,q->B2c,q->B20_grid_variation);
        // Get duration. Substart timepoints to get duration. To cast it to proper unit
        // use duration cast method
        auto duration = duration_cast<microseconds>(stop - start);
    
        std::cout << "Time taken by function: "
            << duration.count()/1000000.0 << " seconds" << std::endl;
    }

    return 0;
}

int choose_eta(Qsc * q, qscfloat * etabar_out, int verbose) {
    // Etabar optimisation using eta_res as cost function

    int status; 
    int iter=0, max_iter=100; 
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;

    /* Set the bracket for Brent search. Use the origin for a, the curvature maximum for b (upper bound of etabar),
     and a nearby point from b for m. This will fail if the extrema is too close to b. Generally works well */
    double m = 0.99*(q->curvature.max()), m_expected=M_PI; 
    double a=0.0, b = (q->curvature.max()); 

    // Construct GSL function to perform optimisation
    gsl_function F; 
    F.function=&eta_res; 
    F.params=q; 

    // Use Brent scalar optimisation 
    T = gsl_min_fminimizer_brent; 
    s = gsl_min_fminimizer_alloc(T); 

    // Avoid optimisation stopping if the search somehow fails
    gsl_set_error_handler_off();
    status = gsl_min_fminimizer_set(s,&F,m,a,b); 
    if (status){
        // The initial bracket does not work
        std::cout << "Error with initialising eta_bar optimisation!!" << std::endl;
        *etabar_out = GSL_NAN; // Standard error return to the upper level optimisation
        return status;
    }
    if (verbose!=0){
        // Print the optimisation evolution
        std::printf("using %s method\n", gsl_min_fminimizer_name(s)); 
        std::printf("%5s [%9s,%9s] %9s %9s\n",
            "iter","lower","upper","min",
                "err(est)");
        std::printf("%5d [%.7f,%.7f] %.7f %.7f\n",
            iter,a,b,
            m,b-a); 
        }

    // Optimisation iteration
    do {
        iter++; 
        status=gsl_min_fminimizer_iterate(s); 
        m=gsl_min_fminimizer_x_minimum(s); 
        a=gsl_min_fminimizer_x_lower(s); 
        b=gsl_min_fminimizer_x_upper(s); 
        status =gsl_min_test_interval(a,b,0.001,0.0); // Stop criterion when smaller than 0.001 change
        if(verbose!=0){
            if(status==GSL_SUCCESS) 
                printf("Converged:\n"); 
            printf("%5d [%.7f,%.7f] "
                "%.7f %.7f\n", 
                iter,a,b, 
                m,b-a);
            }
        }
    while(status==GSL_CONTINUE&&iter<max_iter); // Stop if exceeded max number of iter
    gsl_min_fminimizer_free(s); 
    
    // Output optimal etabar
    *etabar_out = m;
    return status; 
    }

int choose_B22c(Qsc * q, qscfloat * B22c_out, int verbose) {
    // B2c optimisation using B20_res as cost function

    int status; 
    int iter=0, max_iter=100; 
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;

    // Set bracket for Brent scalar search
    double m;
    // Bounds of the bracket are taken to be B22c=+-20, to avoid extreme shaping
    double a = 3.0*(q->B0)*(q->eta_bar)*(q->eta_bar)/4.0-20.0/4.0*(q->B0)*(q->B0)*(q->B0)*(q->B0);
    double b = 3.0*(q->B0)*(q->eta_bar)*(q->eta_bar)/4.0+20.0/4.0*(q->B0)*(q->B0)*(q->B0)*(q->B0);
    // Try to find an appropriate m
    double lb = B20_res(a, q);
    double rb = B20_res(b, q);
    double lrb = B20_res(a+0.001, q); // Evaluate delB20 close to the boundary
    double rlb = B20_res(b-0.001, q);
    if(verbose!=0){
        std::printf("%f %f \n",a,b);
        std::printf("%f %f %f %f \n",lb,lrb,rlb,rb);
    }

    /* The delB20(B2c) cost function can be regarded, roughly as a V shaped function, each side
     of the valley being monotonic (true in the larger scale). So, in the considerations of the bracket
     if the left bound is lower than a nearby point to its right, then we may assume that a is 
     the minimmum (similarly if the right bound is lower than a point immediatly to its left, then b
     will be the minimum). Otherwise, we may take as a middle point of the bracket one of the nearby values
     to the boundary as shown */
    if (lb<lrb) {
        *B22c_out = a;
        return 0;
    } else if (lb>lrb) {
        if (rb<rlb){
            *B22c_out = b;
            return 0;
        } else {
            if (rlb>lrb)
                m = a+0.001;
            else
                m = b-0.001;
        } 
    }
    // Construct optimisation function
    gsl_function F; 
    F.function = &B20_res; 
    F.params = q; 
    T = gsl_min_fminimizer_brent; 
    s = gsl_min_fminimizer_alloc(T); 
    // Get starting timepoint
    // auto start = high_resolution_clock::now();
    status = gsl_min_fminimizer_set(s,&F,m,a,b); 
    if (status){
        // If bracket fails (usually due to the assumption of monotonicity failing)
        std::cout << "Error with initialising B22c optimisation!!" << std::endl;
        *B22c_out = GSL_NAN;
        return status;
    }
    if(verbose!=0){
        std::printf("using %s method\n", gsl_min_fminimizer_name(s)); 
        std::printf("%5s [%9s,%9s] %9s %9s\n",
            "iter","lower","upper","min",
                "err(est)");
        std::printf("%5d [%.7f,%.7f] %.7f %.7f\n",
            iter,a,b,
            m,b-a); 
        }

    // Optimisation iteration
    do {
        iter++; 
        status=gsl_min_fminimizer_iterate(s); 
        m=gsl_min_fminimizer_x_minimum(s); 
        a=gsl_min_fminimizer_x_lower(s); 
        b=gsl_min_fminimizer_x_upper(s); 
        status =gsl_min_test_interval(a,b,0.001,0.0); // Stop criterion when smaller than 0.001 change
        if(verbose!=0){
            if(status==GSL_SUCCESS) 
                printf("Converged:\n"); 
            printf("%5d [%.7f,%.7f] "
                "%.7f %.7f\n", 
                iter,a,b, 
                m,b-a);
            }
        }
    while(status==GSL_CONTINUE&&iter<max_iter); // Stop if exceeded the number of max iterations
    gsl_min_fminimizer_free(s); 

    // Output B2c value
    *B22c_out = m;
    return status; 
    }

double choose_Z_axis(Qsc* q, double* z_vals, int verbose){
    // Z harmonics optimisation 
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;

    size_t iter = 0;
    int max_iter = 200;
        
    int status;
    double size;
    // Starting point default at Z harm/R harm=0.8 
    int num_harm = (q->R0c.size())-1;
    x = gsl_vector_alloc (num_harm);
    for (int i=0;i<num_harm;i++)
        gsl_vector_set (x, i, 0.8);

    // Set initial step sizes to 0.1 
    ss = gsl_vector_alloc (num_harm);
    gsl_vector_set_all (ss, 0.1);

    // Initialize method and iterate
    minex_func.n = num_harm;
    minex_func.f = findZOpt;
    minex_func.params = q;

    s = gsl_multimin_fminimizer_alloc (T, num_harm);
    status = gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    if (status){
        // Error initialising (probably as a result of a GSL_NAN from one of the underlying optimisations)
        std::cout << "Error with initialising Z axis optimisation!!" << std::endl;
        return status;
    }

    // Optimisation iterations
    do
        {
            iter++;
            status = gsl_multimin_fminimizer_iterate(s);
            if (status)
                break;
            size = gsl_multimin_fminimizer_size (s);
            status = gsl_multimin_test_size (size, 3e-3); // Stop criterion size change
            if (verbose!=0){
                if (status == GSL_SUCCESS)
                    printf ("converged to minimum at\n");
                printf ("%5ld %10.3e %10.3e f() = %7.3f size = %.3f\n",
                iter,
                gsl_vector_get (s->x, 0),
                gsl_vector_get (s->x, 1),
                s->fval, size);
            }
        }
    while (status == GSL_CONTINUE && iter < max_iter);
    // Indicate if max. number of iterations has been reached
    if (iter==max_iter)
        std::cout << "Max. number of iterations reached!" << std::endl;
    gsl_vector_free(x);
    gsl_vector_free(ss);

    // Output Z values to z_vals
    for (int i=0;i<num_harm;i++)
        z_vals[i] = gsl_vector_get (s->x, i);

    // Evaluate optimised NAE configuration
    findZOpt(s->x,q);
    q->r2_diagnostics();

    gsl_multimin_fminimizer_free (s);
    return status;
}