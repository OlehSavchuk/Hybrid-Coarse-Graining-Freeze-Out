# Light Nuclei Thermal Emission

A high-performance C++ implementation for computing the thermal emission of light nuclei from relativistic heavy-ion collision freeze-out surfaces. 

This code processes 4D fluid hypersurfaces (e.g., extracted via Cornelius algorithms), interpolates the Equation of State (EoS) thermodynamic fields (Temperature, Chemical Potential, Effective Mass), and computes the $p_T$, $\phi$, and $y$ differential spectra utilizing the Cooper-Frye formalism. The computationally expensive integrations are parallelized across fluid cells using OpenMP.

For more information please read https://arxiv.org/abs/2605.11275
or presentation in docs/HCGF.pdf
## ⚙️ Prerequisites

To compile and run this code, you will need:
* A compiler supporting **C++17** 
* **CMake** 3.15 or higher
* **OpenMP** (For multi-threaded cell integration)
* **Boost C++ Libraries** (Specifically `<boost/math/special_functions/bessel.hpp>` for modified Bessel functions of the second kind)

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

./light_nuclei_emission ../config.txt
```

