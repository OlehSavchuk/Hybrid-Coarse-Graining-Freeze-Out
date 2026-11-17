// Checks the local momentum anisotropy of LightNuclei::emission():
//  (1) the total yield of a surface element (all 40 species) is unchanged by the anisotropy (moving cell, space-like dsigma),
//  (2) for a static cell the second moments <p_i^2> reproduce the requested e_xx, e_yy, e_zz,
//  (3) s = (1,1,1) gives the isotropic emission.
// Build (repository root):  g++ -std=c++17 -O2 -fopenmp -Iinclude tests/test_conservation.cpp src/lightnuclei.cpp -o test_conservation
// Run in the repository root (needs Ts.dat, es.dat, mus.dat, mStar.dat, particles.dat).
#include "lightnuclei.h"
#include <algorithm>
#include <cstdio>

struct Result { std::vector<double> N, pxx, pyy, pzz; };

// Integrates (E dN/d^3p) with the measure pt_A dpt_A dphi dy over phase space; also the weighted p_x^2, p_y^2, p_z^2.
static Result integrate(LightNuclei& LN, const std::vector<double>& A, const std::vector<double>& mbare, double mStar,
                        const std::vector<double>& v, double theta, const std::vector<double>& dsigma,
                        double T, const std::vector<double>& mu, const double* s)
{
    const int npt = 240, nphi = 48, ny = 120;
    const double pmax = 6000.0;      // MeV per nucleon (emission() scales it by A)
    const double ymax = 3.5;
    Result r{std::vector<double>(40, 0.0), std::vector<double>(40, 0.0), std::vector<double>(40, 0.0), std::vector<double>(40, 0.0)};
    for (int ip = 0; ip < npt; ip++) {
        const double pn = (ip + 0.5) * pmax / npt;
        for (int jp = 0; jp < nphi; jp++) {
            const double phi = -PI + 2 * PI / nphi * (jp + 0.5);
            for (int iy = 0; iy < ny; iy++) {
                const double y = -ymax + 2 * ymax / ny * (iy + 0.5);
                const std::vector<double> e = LN.emission(pn, phi, y, v, dsigma, theta, T, mu, mStar, s);
                for (int i = 0; i < 40; i++) {
                    const double ptA = A[i] * pn / 1000.0;
                    const double m = mbare[i] - A[i] * mStar;
                    const double pz = sqrt(m * m + ptA * ptA) * sinh(y);
                    const double w = A[i] * A[i] * pn * (pmax / npt) / 1e6 * (2 * PI / nphi) * (2 * ymax / ny) * e[i];
                    r.N[i] += w;
                    r.pxx[i] += w * ptA * ptA * cos(phi) * cos(phi);
                    r.pyy[i] += w * ptA * ptA * sin(phi) * sin(phi);
                    r.pzz[i] += w * pz * pz;
                }
            }
        }
    }
    return r;
}

int main()
{
    LightNuclei LN;
    LN.init();
    std::vector<double> A(40), mbare(40);
    { std::ifstream f("particles.dat"); for (int i = 0; i < 40; i++) { double x[6]; for (int j = 0; j < 6; j++) f >> x[j]; mbare[i] = x[1]; A[i] = x[2]; } }

    double T, mStar; std::vector<double> mu;
    LN.get_cell_properties(1.0, T, mu, mStar);
    std::printf("cell: T = %.4f GeV, effective-mass shift %.4f\n", T, mStar);

    const double ex[3] = {1.6, 0.9, 0.5};            // e_xx, e_yy, e_zz (sum 3)
    double s[3], s1[3] = {1, 1, 1};
    LightNuclei::aniso_scales(ex[0], ex[1], ex[2], 1.0, s);
    std::printf("scales s = %.4f %.4f %.4f, product = %.12f\n", s[0], s[1], s[2], s[0] * s[1] * s[2]);

    // (1) moving cell with time-like and space-like surface vector components
    {
        const double vx = 0.30, vy = 0.20, vz = 0.10;
        const double g = 1.0 / sqrt(1 - vx * vx - vy * vy - vz * vz);
        std::vector<double> v = {g, sqrt(vx * vx + vy * vy), vz};
        const double theta = atan2(vy, vx);
        std::vector<double> dsigma = {1.0, 0.15, -0.10, 0.20};
        Result iso = integrate(LN, A, mbare, mStar, v, theta, dsigma, T, mu, nullptr);
        Result one = integrate(LN, A, mbare, mStar, v, theta, dsigma, T, mu, s1);
        Result an  = integrate(LN, A, mbare, mStar, v, theta, dsigma, T, mu, s);
        double maxdiff = 0, maxone = 0;
        std::printf("moving cell: total yield, isotropic / anisotropic / ratio (selected species)\n");
        for (int i = 0; i < 40; i++) {
            if (iso.N[i] <= 0) continue;
            maxdiff = std::max(maxdiff, std::abs(an.N[i] / iso.N[i] - 1.0));
            maxone = std::max(maxone, std::abs(one.N[i] / iso.N[i] - 1.0));
            if (i < 6 || i == 20) std::printf("  species %2d (A=%g): %.6e  %.6e  %.6f\n", i, A[i], iso.N[i], an.N[i], an.N[i] / iso.N[i]);
        }
        std::printf("max |N_aniso/N_iso - 1| over all species: %.2e;  s = 1 versus isotropic: %.2e\n", maxdiff, maxone);
    }

    // (2) static cell, dsigma = (1,0,0,0): second moments in the rest frame
    {
        std::vector<double> v = {1.0, 0.0, 0.0};
        std::vector<double> dsigma = {1.0, 0.0, 0.0, 0.0};
        Result iso = integrate(LN, A, mbare, mStar, v, 0.0, dsigma, T, mu, nullptr);
        Result an  = integrate(LN, A, mbare, mStar, v, 0.0, dsigma, T, mu, s);
        std::printf("static cell: 3<p_i^2>/sum  (requested %.3f %.3f %.3f; isotropic gives 1 1 1)\n", ex[0], ex[1], ex[2]);
        for (int i : {0, 1, 3, 6}) {
            const double tr = iso.pxx[i] + iso.pyy[i] + iso.pzz[i], ta = an.pxx[i] + an.pyy[i] + an.pzz[i];
            std::printf("  species %d (A=%g)  isotropic: %.3f %.3f %.3f   anisotropic: %.3f %.3f %.3f   N ratio %.5f\n", i, A[i],
                        3 * iso.pxx[i] / tr, 3 * iso.pyy[i] / tr, 3 * iso.pzz[i] / tr,
                        3 * an.pxx[i] / ta, 3 * an.pyy[i] / ta, 3 * an.pzz[i] / ta, an.N[i] / iso.N[i]);
        }
    }
    return 0;
}
