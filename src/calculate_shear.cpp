#include <cmath>
#include "qsc.hpp"

using namespace qsc;

void qsc::Qsc::calculate_shear(qscfloat B31c=0) {
    // Shorthand introduced: we have to ransform to 1/B**2 expansion parameters, taking into account the 
    // difference in the definition of the radial coordinate. In the work of Rodriguez et al.,
    // Phys. Plasmas, (2021), epsilon=sqrt(psi) while in the work of Landreman et al.,
    // J. Plasma Physics (2019) it is defined r=\sqrt(2*psi/B0). Need to transform between the
    // two.

    qscfloat eps_scale = std::sqrt(2.0/B0); 

    // sign_psi = self.spsi
    // sign_G   = self.sG  // Sign is taken to be positive for simplicity. To include this, need to track expressions
    qscfloat dldp = abs_G0_over_B0;

    // Transformation to 1/B**2 parameters 
    qscfloat B0e, B1ce, eta, B31s, I4;
    Vector B20e;
    B0e = 1/B0/B0;
    eta = eta_bar*std::sqrt(2.0/B0);
    B1ce = -2*B0e*eta_bar*std::sqrt(2.0)*std::pow(B0e,0.25);
    B20e = (0.75*eta_bar*eta_bar*B0 - B20)*4*B0e*B0e;
    B31s = 0; // To preserve stellarator symmetry
    I4 = 0; // Take current variations at this order to be 0
            
    // Compute Z31c and Z31s from Cp2: we assume standard equilibria, meaning that we may
    // pick Bpsi0=0 and Bpsi1=0
    
    Vector Z31s(nphi); 
    Vector Z31c = -1.0/3.0/G0/X1c/Y1s*eps_scale*(2*iota_N*(X1c*X2s - Y2c*Y1s + Y1c*Y2s) + eps_scale*eps_scale*(-2*G0*X2s*Y1c*Z20 +
        2*G0* X2c*Y1s*Z20 + 2*G0*X1c*Y2s*Z20 - 4*G0*X2s*Y1c*Z2c - 2*G0* X20*Y1s*Z2c +
        4*G0*X1c*Y2s*Z2c) - abs_G0_over_B0*(-torsion*(2*X20*Y1c + X2c*Y1c - 2*X1c*Y20 - X1c*Y2c +
        X2s*Y1s) + I2*(2*X20*Y1c + X2c*Y1c - 2*X1c*Y20 - X1c*Y2c + X2s*Y1s)*eps_scale*eps_scale - 
        2*curvature*X1c*Z20 - curvature*X1c*Z2c) + eps_scale*eps_scale*(2*G0*X20*Y1c*Z2s + 4*G0*X2c*Y1c*Z2s - 
        2*G0*X1c*Y20*Z2s - 4*G0*X1c*Y2c*Z2s) + 2*X1c*d_X20_d_varphi + X1c*d_X2c_d_varphi+2*Y1c*d_Y20_d_varphi +
        Y1c*d_Y2c_d_varphi + Y1s*d_Y2s_d_varphi);

    // printf("%e\n", Z31c[30]);
         
    Vector dZ31cdp(nphi);

    matrix_vector_product(d_d_varphi, Z31c, dZ31cdp);
            
    Z31s = 1.0/3.0/G0/X1c/Y1s*eps_scale*(2*iota_N*(X1c*X2c + Y1c*Y2c + Y1s*Y2s) + eps_scale*eps_scale*(-2*G0*X2c*Y1c*Z20 + 
        2*G0*X1c*Y2c*Z20 - 2*G0*X2s*Y1s*Z20 + 2*G0*X20*Y1c*Z2c - 2*G0*X1c*Y20*Z2c +
        4*G0*X2s*Y1s*Z2c + 2*G0*X20*Y1s*Z2s - 4*G0*X2c*Y1s*Z2s) + abs_G0_over_B0*(eps_scale*eps_scale*(I2*X2s*Y1c + 
        2*I2*X20*Y1s - I2*X2c*Y1s - I2*X1c*Y2s) - torsion*(X2s*Y1c + 2*X20*Y1s - X2c*Y1s -
        X1c*Y2s) - curvature*X1c*Z2s) - X1c*d_X2s_d_varphi - 2*Y1s*d_Y20_d_varphi + Y1s*d_Y2c_d_varphi - Y1c*d_Y2s_d_varphi);
            
    Vector dZ31sdp(nphi);
    matrix_vector_product(d_d_varphi, Z31s, dZ31sdp);
            
    // Equation J3: expression for X31c/s
    Vector X31c, X31s;
    X31c = 0.5/abs_G0_over_B0/abs_G0_over_B0/curvature*(-2*G0*eps_scale*eps_scale*(G2 + iota_N*I2)*B1ce - G0*G0*B31c +
        eps_scale*eps_scale*eps_scale*(2*dldp*dldp*torsion*torsion*X1c*X20 +
        2*iota_N*iota_N*X1c*X2c + dldp*dldp*torsion*torsion*X1c*X2c + dldp*dldp*curvature*curvature*X1c*(2*X20 + X2c) - 
        3*dldp*iota_N*torsion*X2s*Y1c + 2*dldp*dldp*torsion*torsion*Y1c*Y20 + 2*iota_N*iota_N*Y1c*Y2c +
        dldp*dldp*torsion*torsion*Y1c*Y2c + 2*dldp*iota_N*torsion*X20*Y1s + 3*dldp*iota_N*torsion*X2c*Y1s +
        3*dldp*iota_N*torsion*X1c*Y2s + 2*iota_N*iota_N*Y1s*Y2s + dldp*dldp*torsion*torsion*Y1s*Y2s) + 
        2*dldp*iota_N*Z31s + eps_scale*eps_scale*eps_scale*(2*iota_N*X2s*d_X1c_d_varphi - 2*dldp*torsion*Y20*d_X1c_d_varphi -
        dldp*torsion*Y2c*d_X1c_d_varphi - 2*dldp*torsion*Y1c*d_X20_d_varphi + 2*d_X1c_d_varphi*d_X20_d_varphi - 
        dldp*torsion*Y1c*d_X2c_d_varphi + d_X1c_d_varphi*d_X2c_d_varphi - iota_N*X1c*d_X2s_d_varphi - dldp*torsion*Y1s*d_X2s_d_varphi + 
        2*dldp*torsion*X20*d_Y1c_d_varphi + dldp*torsion*X2c*d_Y1c_d_varphi +
        2*iota_N*Y2s*d_Y1c_d_varphi + 2*dldp*torsion*X1c*d_Y20_d_varphi + 2*iota_N*Y1s*d_Y20_d_varphi + 
        2*d_Y1c_d_varphi*d_Y20_d_varphi + dldp*torsion*X1c*d_Y2c_d_varphi + iota_N*Y1s*d_Y2c_d_varphi +
        d_Y1c_d_varphi*d_Y2c_d_varphi + dldp*torsion*X2s*d_Y1s_d_varphi - 2*iota_N*Y2c*d_Y1s_d_varphi - 
        iota_N*Y1c*d_Y2s_d_varphi + d_Y1s_d_varphi*d_Y2s_d_varphi + dldp*curvature*(-3*iota_N*X1c*Z2s - 
        dldp*torsion*(Y1c*(2*Z20 + Z2c) + Y1s*Z2s) + 2*Z20*d_X1c_d_varphi + Z2c*d_X1c_d_varphi - 2*X1c*d_Z20_d_varphi - 
        X1c*d_Z2c_d_varphi)) + 2*dldp*dZ31cdp);
                
    X31s = 0.5/dldp/dldp/curvature*(-G0*G0*B31s + eps_scale*eps_scale*eps_scale*(dldp*dldp*curvature*curvature*X1c*X2s +
        dldp*dldp*torsion*torsion*X1c*X2s + 2*dldp*dldp*torsion*torsion*Y20*Y1s - dldp*dldp*torsion*torsion*Y2c*Y1s + 
        dldp*dldp*torsion*torsion*Y1c*Y2s + 2*iota_N*iota_N*(X1c*X2s - Y2c*Y1s + Y1c*Y2s) - 2*dldp*dldp*curvature*torsion*Y1s*Z20 + 
        dldp*dldp*curvature*torsion*Y1s*Z2c - dldp*dldp*curvature*torsion*Y1c*Z2s - dldp*torsion*Y2s*d_X1c_d_varphi +
        dldp*curvature*Z2s*d_X1c_d_varphi - 2*dldp*torsion*Y1s*d_X20_d_varphi + dldp*torsion*Y1s*d_X2c_d_varphi - 
        dldp*torsion*Y1c*d_X2s_d_varphi + d_X1c_d_varphi*d_X2s_d_varphi + dldp*torsion*X2s*d_Y1c_d_varphi +
        2*dldp*torsion*X20*d_Y1s_d_varphi - dldp*torsion*X2c*d_Y1s_d_varphi + 2*d_Y20_d_varphi*d_Y1s_d_varphi - 
        d_Y2c_d_varphi*d_Y1s_d_varphi + dldp*torsion*X1c*d_Y2s_d_varphi + d_Y1c_d_varphi*d_Y2s_d_varphi +
        iota_N*(-dldp*torsion*(2*X20*Y1c - 3*X2c*Y1c - 2*X1c*Y20 + 3*X1c*Y2c - 3*X2s*Y1s) + dldp*curvature*X1c*
        (-2*Z20 + 3*Z2c)- 2*X2c*d_X1c_d_varphi - 2*X1c*d_X20_d_varphi + X1c*d_X2c_d_varphi - 2*Y2c*d_Y1c_d_varphi -
        2*Y1c*d_Y20_d_varphi + Y1c*d_Y2c_d_varphi - 2*Y2s*d_Y1s_d_varphi + Y1s*d_Y2s_d_varphi) -
        dldp*curvature*X1c*d_Z2s_d_varphi) - 2*dldp*Z31c*iota_N + 2*dldp*dZ31sdp);
    
    Vector dX31sdp(nphi);
    matrix_vector_product(d_d_varphi, X31s, dX31sdp);

                        
    // Equation Cb2
    Vector Y31s;
    Y31s = 0.25/G0/X1c*eps_scale*(-2*(G2 + iota_N*I2)*X1c*Y1s*eps_scale*eps_scale + 2*iota_N*I2*X1c*Y1s*eps_scale*eps_scale -
        dldp*(4*curvature*X20 - torsion*I2*(X1c*X1c + Y1c*Y1c + Y1s*Y1s)*eps_scale*eps_scale) + 4*G0*(X31s*Y1c/eps_scale + 
        2*X2s*Y2c*eps_scale*eps_scale - X31c*Y1s/eps_scale - 2*X2c*Y2s*eps_scale*eps_scale) - eps_scale*eps_scale*(
        I2*Y1c*d_X1c_d_varphi + I2*X1c*d_Y1c_d_varphi) + 4*d_Z20_d_varphi); 

    Vector dY31sdp(nphi);
    matrix_vector_product(d_d_varphi, Y31s, dY31sdp);
    
    // From the equation for Bt to order n=4, and looking at m=0
    Vector LamTilde;
    LamTilde = 2.0/Y1s/Y1s*(G0*B0e*I4 + ((G2 + iota_N*I2)*B0e*eps_scale*eps_scale + G0*B20e)*I2) + 
        1/Y1s/Y1s*eps_scale*eps_scale*(-2*iota_N*(2*X2c*X2c + X1c*X31c/eps_scale/eps_scale/eps_scale +
        2*X2s*X2s + 2*Y2c*Y2c + 2*Y2s*Y2s + Y1s*Y31s/eps_scale/eps_scale/eps_scale +
        2*Z2c*Z2c + 2*Z2s*Z2s) + 2*dldp*(-torsion*(-X31s*Y1c/eps_scale/eps_scale/eps_scale -
        2*X2s*Y2c + X31c*Y1s/eps_scale/eps_scale/eps_scale + 2*X2c*Y2s + X1c*Y31s/eps_scale/eps_scale/eps_scale) + 
        curvature*(-2*X2s*Z2c + 2*X2c*Z2s + X1c*Z31s/eps_scale/eps_scale/eps_scale)) -
        X31s*d_X1c_d_varphi/eps_scale/eps_scale/eps_scale - 2*X2s*d_X2c_d_varphi + 2*X2c*d_X2s_d_varphi + 
        X1c*dX31sdp/eps_scale/eps_scale/eps_scale - Y31s*d_Y1c_d_varphi/eps_scale/eps_scale/eps_scale -
        2*Y2s*d_Y2c_d_varphi + 2*Y2c*d_Y2s_d_varphi + Y1c*dY31sdp/eps_scale/eps_scale/eps_scale - 2*Z2s*d_Z2c_d_varphi + 2*Z2c*d_Z2s_d_varphi);

    // Need to compute the integration factor necessary for computing the shear
    Matrix DMred(nphi-1,nphi-1);
    for (int i=0;i<nphi-1;i++){
        for(int j=0;j<nphi-1;j++)
            DMred(i,j) = d_d_varphi(1+i,1+j);   // The differentiation matrix has a linearly dependent row, focus on submatrix
    }

    Vector temp(nphi-1);
    for(int i=0;i<nphi-1;i++)
        temp[i] = sigma[i+1];

    std::valarray<int> shear_ipiv(nphi-1);
    // Distinguish between the stellarator symmetric case and the non-symmetric one at order r^1.
    // Distinction leads to the expSig function being periodic (stell. sym.) or not.
    if (sigma0 == 0 && R0s.max() == 0 && R0s.min() == 0 && Z0c.max() == 0 && Z0c.max() == 0){
        // Case in which sigma is stellarator-symmetric:
        linear_solve(DMred,temp,shear_ipiv);   // Invert differentiation matrix: as if first entry a zero, need to add it later
        Vector integSig(nphi);
        integSig[0] = 1;
        for (int i=0;i<nphi-1;i++)
            integSig[i+1] = std::exp(2*iota_N*temp[i]);  // Add the first entry 0
        LamTilde = integSig*LamTilde*d_varphi_d_phi;
        integSig = integSig*(X1c*X1c + Y1c*Y1c + Y1s*Y1s)/Y1s/Y1s*d_varphi_d_phi;
        iota2 = B0/2.0*LamTilde.sum()/integSig.sum();
    } else{
        // Case in which sigma is not stellarator-symmetric:
        // d_phi_d_varphi = 1 + np.matmul(d_d_varphi,self.phi-self.varphi)
        qscfloat avSig;
        X31c = sigma*d_varphi_d_phi;
        avSig = X31c.sum()/nphi;     // Separate the piece that gives secular part, so all things periodic
        temp = temp-avSig;
        linear_solve(DMred,temp,shear_ipiv);   // Invert differentiation matrix: as if first entry a zero, need to add it later
        Vector integSig(nphi+1);
        integSig[0] = 1;
        for (int i=0;i<nphi-1;i++)
            integSig[i+1] = std::exp(2*iota_N*(temp[i] + avSig*Boozer_toroidal_angle[i+1]));  // Include the secular piece
        integSig[nphi] = std::exp(2*iota_N*(avSig*2*pi/nfp));   // Add endpoint at 2*pi for better integration
        qscfloat numer = 0, denom = 0;
        for (int i=0;i<nphi-1;i++){
            numer = numer + (integSig[i]*LamTilde[i]+integSig[i+1]*LamTilde[i+1])*(Boozer_toroidal_angle[i+1]-
                Boozer_toroidal_angle[i])/2.0;
        }
        numer = numer + (integSig[nphi-1]*LamTilde[nphi-1]+integSig[nphi]*LamTilde[0])*(2*pi/nfp -
            Boozer_toroidal_angle[nphi-1])/2.0;
        for (int i=0;i<nphi-1;i++){
            denom = denom + (integSig[i]*(X1c[i]*X1c[i] + Y1c[i]*Y1c[i] + Y1s[i]*Y1s[i])/ Y1s[i]/Y1s[i]+
            integSig[i+1]*(X1c[i+1]*X1c[i+1] + Y1c[i+1]*Y1c[i+1] + Y1s[i+1]*Y1s[i+1])/ Y1s[i+1]/Y1s[i+1])*(Boozer_toroidal_angle[i+1]-
                Boozer_toroidal_angle[i])/2.0;
        }
        denom = denom + (integSig[nphi-1]*(X1c[nphi-1]*X1c[nphi-1] + Y1c[nphi-1]*Y1c[nphi-1] + Y1s[nphi-1]*Y1s[nphi-1])/ 
            Y1s[nphi-1]/Y1s[nphi-1]+integSig[nphi]*(X1c[0]*X1c[0] + Y1c[0]*Y1c[0] + Y1s[0]*Y1s[0])/ Y1s[0]/Y1s[0])*(2*pi/nfp -
                Boozer_toroidal_angle[nphi-1])/2.0;
        iota2 = B0 / 2.0 * numer / denom;
    }


}