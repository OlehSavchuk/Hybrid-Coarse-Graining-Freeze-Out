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
    std::vector<double> emission(double pt, double phi, double y, const std::vector<double>& v, const std::vector<double>& dsigma, double theta, double T, const std::vector<double>& mu_for_e, double mStar_val);
    
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
