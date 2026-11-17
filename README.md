# Light Nuclei Thermal Emission

A high-performance C++ implementation for computing the thermal emission of light nuclei from relativistic heavy-ion collision freeze-out surfaces. 

This code processes 4D fluid hypersurfaces (e.g., extracted via Cornelius algorithms), interpolates the Equation of State (EoS) thermodynamic fields (Temperature, Chemical Potential, Effective Mass), and computes the $p_T$, $\phi$, and $y$ differential spectra utilizing the Cooper-Frye formalism. The computationally expensive integrations are parallelized across fluid cells using OpenMP.

For more information please read https://arxiv.org/abs/2605.11275
or presentation in docs/HCGF.pdf

For citations:
```bibtex
@article{Savchuk:2026dij,
    author = "Savchuk, Oleh and Danielewicz, Pawel and Lynch, William and Margueron, J{\'e}r{\^o}me",
    title = "{Freeze-out model of light nuclei formation in heavy-ion collision transport}",
    eprint = "2605.11275",
    archivePrefix = "arXiv",
    primaryClass = "nucl-th",
    month = "5",
    year = "2026"
}
```
## ⚙️ Prerequisites

To compile and run this code, you will need:
* A compiler supporting **C++17** 
* **CMake** 3.15 or higher
* **OpenMP** (For multi-threaded cell integration)
* **Boost C++ Libraries** (Specifically `<boost/math/special_functions/bessel.hpp>` for modified Bessel functions of the second kind)

## Local momentum anisotropy (branches `local-momentum-anisotropy`, `full-momentum-anisotropy`)
HCGF normally emits particles isotropically in the rest frame of every freeze-out cell. With `anisotropy = 1` in `config.txt`
the cell can carry the local momentum anisotropy `exx, eyy, ezz` (`e_ii = 3<p_i p_i>/(<p_x p_x> + <p_y p_y> + <p_z p_z>)`,
computed in the cell rest frame; `e_xx + e_yy + e_zz = 3`), which the freeze-out surface generator interpolates to the surface
(`surface.dat` then has 18 columns: the three extra columns 16-18 are `exx eyy ezz`).

`anisotropy = 2` additionally uses the off-diagonal elements `exy, exz, eyz` (`e_ij = 3<p_i p_j>/sum_k <p_k^2>`, columns 19-21 of `surface.dat`).

In the cell rest frame the distribution is `f(sqrt(m^2 + p*^T Sigma^-1 p*))` instead of `f(E*)`, where `p*` is the particle momentum boosted
into the cell rest frame (pure boost, axes parallel to the lab axes) and `Sigma = C / det(C)^(1/3)` with `C = 1 + s (e - 1)`
(`s` = `aniso_strength`). For `anisotropy = 1` the off-diagonal elements of `C` are set to zero. If `C` is not positive definite the
off-diagonal terms of that cell are dropped; cells with non-finite `e` are emitted isotropically.
* **The yield is conserved:** `det Sigma = 1`, so the Jacobian of the momentum-space transformation is 1 and the number of particles
  of each species emitted by a cell, and therefore the total yield, is unchanged for any statistics (Boltzmann, Fermi, Bose).
  Space-like parts of `d sigma` do not matter because the distribution stays symmetric under `p -> -p`.
* The second moments of the emitted momentum are `<p_i p_j> ~ Sigma_ij`, i.e. `e_ij` is reproduced for every species (the same `Sigma` for all species).
* The off-diagonal elements (in particular `e_xz`) correlate `p_x` with `p_z` and therefore change the rapidity dependence of the directed flow,
  which a diagonal tensor cannot do.
* `aniso_strength` scales `(e - 1)`: `0` reproduces the isotropic code, `1` uses the extracted values.
* If the surface file has fewer columns than expected the program stops with an error.
* `phi_half = 1` computes only `phi` in `[0, pi]` and mirrors the output to `-phi`, which halves the run time.

`tests/test_conservation.cpp` checks, on a moving cell with a space-like `d sigma`, that the total yield of all 40 species is unchanged
(to a few 1e-5, the accuracy of the phase-space integral), that the rest-frame second moments of a static cell reproduce `e_xx, e_yy, e_zz` and, for the full tensor, `e_xy, e_xz, e_yz`,
and that `Sigma = 1` is identical to the isotropic emission. Build with `cmake -DBUILD_TESTS=ON ..` and run `./test_conservation` in the repository root.

## 📂 Expected Data Inputs

The solver expects the following EoS and parameter data files in the working directory to execute `LightNuclei::init()`:
* `Ts.dat`: Temperature mapping
* `es.dat`: Energy density mapping
* `mStar.dat`: Effective mass shift mappings
* `mus.dat`: Chemical potentials
* `particles.dat`: Defined properties of the 40 light nuclei species
* `surface.dat`: The output of your freeze-out surface generator

## 🚀 Build Instructions

We recommend an out-of-source build using CMake to cleanly link Boost and OpenMP. From the root repository directory:

```bash
mkdir build && cd build
cmake ..
make -j4

./lightnuclei ../config.txt
```

