#ifndef LIGHTNUCLEICLASS_H
#define LIGHTNUCLEICLASS_H

#define PI 3.14159265 //[cite: 10]
#define hbarc 0.197   // Conversion constant[cite: 10]

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <boost/math/special_functions/bessel.hpp> // Required for modified Bessel functions in full_yield[cite: 10]

class LightNuclei {
    // Thermodynamic tables for Equation of State (EoS) mapping[cite: 10]
    std::vector<double> Ts;
    std::vector<double> es;
    std::vector<double> mStar;
    std::vector<std::vector<double>> particles;
    std::vector<std::vector<double>> mus;
    int Nextd = 40; // Number of predefined particle species[cite: 10]
    
public:
    // Loads equation of state and particle data from .dat files[cite: 10]
    void init();

    // Fetches parameters once per cell to avoid redundant interpolations across pT/phi bins[cite: 10]
    bool get_cell_properties(double e, double& T, std::vector<double>& mu, double& m_star);

    // Computes the Cooper-Frye emission yield for all particle species from a single fluid cell[cite: 10]
    // aniso_inv (optional, 6 numbers xx yy zz xy xz yz of Sigma^-1, det Sigma = 1, see aniso_tensor): local momentum anisotropy of the cell.
    // The distribution in the cell rest frame is f(sqrt(m^2 + (p*_x/s_x)^2 + (p*_y/s_y)^2 + (p*_z/s_z)^2)); since the
    // Jacobian of the rescaling is s_x*s_y*s_z = 1 the particle number of the cell (and the total yield) is unchanged.
    std::vector<double> emission(double pt, double phi, double y, const std::vector<double>& v, const std::vector<double>& dsigma, double theta, double T, const std::vector<double>& mu_for_e, double mStar_val, const double* aniso_inv = nullptr);

    // Converts the local momentum anisotropy e_xx, e_yy, e_zz (= 3<p_i p_i>/sum, e_xx+e_yy+e_zz = 3) into the scales
    // s_i with s_i^2 proportional to e_ii and s_x*s_y*s_z = 1. strength = 1 uses the input as is, 0 switches the anisotropy off.
    static void aniso_scales(double exx, double eyy, double ezz, double strength, double* s);

    // Full local momentum tensor: e = (exx, eyy, ezz, exy, exz, eyz) with e_ij = 3<p_i p_j>/(<p_x p_x>+<p_y p_y>+<p_z p_z>)
    // (isotropic: e_ii = 1, e_ij = 0). Builds C = 1 + strength (e - 1), normalises it to det = 1 (Sigma = C / det^(1/3)) and returns
    // inv = Sigma^-1 as (xx, yy, zz, xy, xz, yz). If C is not positive definite the off-diagonal terms are dropped.
    // The rest-frame distribution is f(sqrt(m^2 + p*^T Sigma^-1 p*)); since det Sigma = 1 the yield per cell is unchanged.
    static void aniso_tensor(const double* e6, double strength, double* inv6);
    
    // Calculates the center-of-mass momentum relative to the fluid cell velocity[cite: 10]
    double qCM(double p0, double pz, double pt_A, double phi, const std::vector<double>& v, double theta);
    
    // Boltzmann factor calculation[cite: 10]
    double fB(double E, double T);
    
    // Interpolation routines mapping energy density (e) to chemical potential, effective mass, and temperature[cite: 10]
    std::vector<double> interpolate_mu(double e, const std::vector<double>& es_vec);
    double interpolate_mStar(double e, const std::vector<double>& es_vec);
    double interpolate_T(double e, const std::vector<double>& es_vec);
    std::vector<double> interpolate_mu_for_T(double T, const std::vector<double>& Ts_vec);

    // Computes the thermal yield for a specific momentum state[cite: 10]
    double yield(double T, double E, const std::vector<double>& tmu, const std::vector<double>& particle);
    
    // Computes the integrated thermal yield[cite: 10]
    double full_yield(double T, const std::vector<double>& tmu, const std::vector<double>& particle);
};

#endif
