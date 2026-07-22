#include "lightnuclei.h"

// Initialize the class by reading EoS and particle tables from external data files[cite: 9]
void LightNuclei::init() {
    std::ifstream file("Ts.dat");
    double tvalue;
    for(int i = 0; i < 100; i++) { file >> tvalue; this->Ts.push_back(tvalue); }
    file.close();

    file.open("es.dat");
    for(int i = 0; i < 100; i++) { file >> tvalue; this->es.push_back(tvalue); }
    file.close();

    file.open("particles.dat");
    for(int i = 0; i < Nextd; i++) {
        std::vector<double> particle(6);
        for(int j = 0; j < 6; j++) file >> particle[j];
        this->particles.push_back(particle);
    }
    file.close();

    file.open("mus.dat");
    for(int i = 0; i < 100; i++) {
        std::vector<double> particle(3);
        for(int j = 0; j < 3; j++) file >> particle[j];
        this->mus.push_back(particle);
    }
    file.close();

    file.open("mStar.dat");
    for(int i = 0; i < 100; i++) { file >> tvalue; this->mStar.push_back(tvalue); }
    file.close();
}

// Retrieves all thermodynamic properties for a given cell's rest-frame energy density[cite: 9]
bool LightNuclei::get_cell_properties(double e, double& T, std::vector<double>& mu, double& m_star) {
    T = interpolate_T(e, es);
    mu = interpolate_mu(e, es);
    m_star = interpolate_mStar(e, es);
    return true;
}

// Calculates the effective momentum Q in the center of mass of the fluid cell[cite: 9]
double LightNuclei::qCM(double p0, double pz, double pt_A, double phi, const std::vector<double>& v, double theta) {
    return v[0] * (p0 - v[1] * pt_A * cos(phi - theta) - v[2] * pz);
}

// --- 1D Linear Interpolation Functions for EoS Variables --- //

std::vector<double> LightNuclei::interpolate_mu(double e, const std::vector<double>& es_vec) {
    int i = 0;
    while(e < es_vec[i] && i < 99) i++;
    std::vector<double> mu_for_e(3);
    for(int j = 0; j < 3; j++) mu_for_e[j] = this->mus[i][j];
    return mu_for_e;
}

std::vector<double> LightNuclei::interpolate_mu_for_T(double T, const std::vector<double>& Ts_vec) {
    int i = 0;
    while(T < Ts_vec[i] && i < 99) i++;
    std::vector<double> mu_for_T(3);
    for(int j = 0; j < 3; j++) mu_for_T[j] = this->mus[i][j];
    return mu_for_T;
}

double LightNuclei::interpolate_T(double e, const std::vector<double>& es_vec) {
    int i = 0;
    while(e < es_vec[i] && i < 99) i++;
    return this->Ts[i];
}

double LightNuclei::interpolate_mStar(double e, const std::vector<double>& es_vec) {
    int i = 0;
    while(e < es_vec[i] && i < 99) i++;
    return this->mStar[i];
}

// ----------------------------------------------------------- //

// Calculates the thermal emission array for all particles from a single hypersurface element[cite: 9]
std::vector<double> LightNuclei::emission(double pt, double phi, double y, const std::vector<double>& v, const std::vector<double>& dsigma, double theta, double T, const std::vector<double>& mu_for_e, double mStar_val) {
    std::vector<double> produced(Nextd, 0.0);
    
    for(int i = 0; i < Nextd; i++) {
        // Adjust bare mass with the effective mass shift[cite: 9]
        double m = particles[i][1] - particles[i][2] * mStar_val;
        double A = particles[i][2]; // Atomic mass number[cite: 9]
        double pt_A = A * (pt / 1000.0); // Transverse momentum scaled by mass number[cite: 9]
        
        // Kinematic variables[cite: 9]
        double mT = sqrt(m * m + pt_A * pt_A);
        double p0 = mT * cosh(y);
        double pz = mT * sinh(y);
        
        double QCM = qCM(p0, pz, pt_A, phi, v, theta);
        
        // Contraction of the normal surface vector with the four-momentum: p^\mu d\sigma_\mu[cite: 9]
        double dsigmadotp = (dsigma[0] * p0 + pt_A * (cos(phi) * dsigma[1] + sin(phi) * dsigma[2]) + dsigma[3] * pz);
        
        if(dsigmadotp != 0) {
            // Compute final yield based on the statistical distribution[cite: 9]
            produced[i] = dsigmadotp * yield(T, QCM, mu_for_e, particles[i]);
        }
    }
    return produced;
}

double LightNuclei::fB(double E, double T) {
    return exp(-E / T); // Standard Boltzmann factor[cite: 9]
}

// Statistical thermal yield for a specific species at a specific energy[cite: 9]
double LightNuclei::yield(double T, double E, const std::vector<double>& tmu, const std::vector<double>& particle) {
    // Phase space prefactor 1/(2\pi)^3[cite: 9]
    static const double prefactor = 1.0 / pow(2 * PI, 3);
    double tyield = particle[0] * prefactor;
    
    // Sum the chemical potentials[cite: 9]
    double tsum = 0;
    for(int i = 2; i < 5; i++) {
        tsum += particle[i] * tmu[i-2];
    }
    
    // Exponential term in the Fermi-Dirac / Bose-Einstein denominator[cite: 9]
    double exp_term = exp(-(tsum - E) / T);
    double denom = exp_term + particle[5];

    if(denom > 0) {
        return tyield / denom;
    }
    return 0;
}

// Integrated full yield utilizing modified Bessel functions of the second kind[cite: 9]
double LightNuclei::full_yield(double T, const std::vector<double>& tmu, const std::vector<double>& particle) {
    double tyield = (particle[0] * particle[1] * particle[1] * T) * boost::math::cyl_bessel_k(2, particle[1] / T) / (2 * PI * PI * pow(0.197, 3));
    double tsum = 0;
    for(int i = 2; i < 5; i++) {
        tsum += particle[i] * tmu[i-2];
    }
    tyield *= exp(tsum / T);
    return tyield;
}
