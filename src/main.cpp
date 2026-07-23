#include "lightnuclei.h"
#include <algorithm>
#include <string>
#include <sstream>
#include <map>
#include <filesystem>
#include <omp.h> // Included for multi-threading calculations[cite: 11]

// Struct to hold precomputed fluid cell data mapping to the freeze-out surface[cite: 11]
struct CellData {
    std::vector<double> dsigma; // Surface normal vector[cite: 11]
    std::vector<double> v;      // Fluid velocity[cite: 11]
    double theta;               // Azimuthal angle[cite: 11]
    double T;                   // Temperature[cite: 11]
    std::vector<double> mu;     // Chemical potentials[cite: 11]
    double mStar;               // Effective mass shift[cite: 11]
    bool valid;                 // Flag for physical validity[cite: 11]
};

// Helper function to read run parameters from config.txt[cite: 11]
std::map<std::string, std::string> readConfig(const std::string& filename) {
    std::map<std::string, std::string> config;
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue; // Ignore comments[cite: 11]
        std::istringstream is_line(line);
        std::string key;
        if (std::getline(is_line, key, '=')) {
            std::string value;
            if (std::getline(is_line, value)) {
                // Trim whitespaces[cite: 11]
                key.erase(0, key.find_first_not_of(" \t"));
                key.erase(key.find_last_not_of(" \t") + 1);
                value.erase(0, value.find_first_not_of(" \t"));
                value.erase(value.find_last_not_of(" \t") + 1);
                config[key] = value;
            }
        }
    }
    return config;
}

int main(int argc, char** argv) {
    std::string config_path = "config.txt";
    if (argc > 1) config_path = argv[1];
    
    // Parse simulation parameters[cite: 11]
    auto config = readConfig(config_path);
    
    std::string path = config.count("path_to_surface") ? config["path_to_surface"] : ".";
    double n0 = config.count("n0") ? std::stod(config["n0"]) : 0.01;
    int bins_phi = config.count("bins_phi") ? std::stoi(config["bins_phi"]) : 20;
    int bins_pt = config.count("bins_pt") ? std::stoi(config["bins_pt"]) : 40;
    double t_y_start = config.count("t_y_start") ? std::stod(config["t_y_start"]) : 0.05;
    double t_y_step = config.count("t_y_step") ? std::stod(config["t_y_step"]) : 0.10;
    int t_y_steps_count = config.count("t_y_steps_count") ? std::stoi(config["t_y_steps_count"]) : 10;

    // Define transverse momentum bin centers[cite: 11]
    double ptbins[40];
    for(int i = 0; i < 40; i++) {
        ptbins[i] = 25 + 50 * i;
    }

    LightNuclei LN;
    LN.init(); // Load EoS data[cite: 11]
    
    double volume_fr = 0;
    std::vector<CellData> valid_cells; // Stores verified surface elements[cite: 11]
    
    std::filesystem::path base_dir = path;
    std::filesystem::path surface_file = base_dir / "surface.dat";
    std::ifstream surface("../surface.dat"); // Note: hardcoded fallback, ensure path consistency[cite: 11]

    std::cout << "Started reading fo hypersurface..." << std::endl;
    std::string line;
    
    // Process the hypersurface geometry and thermodynamic fields line by line[cite: 11]
    while (std::getline(surface, line)) {
        if (line.find("nan") != std::string::npos || line.find("NaN") != std::string::npos) continue;

        std::istringstream iss(line);
        std::vector<double> dsigma(4), coordinate_fr(4), velocity_fr(4);
        double enden_fr, temperature_fr, density_fr;

        // Parse d^3\sigma_\mu, coordinates, energy density, temperature, net-baryon density, and 4-velocity[cite: 11]
        for (int i = 0; i < 4; i++) iss >> dsigma[i];
        for (int i = 0; i < 4; i++) iss >> coordinate_fr[i];
        iss >> enden_fr >> temperature_fr >> density_fr;
        for (int i = 0; i < 4; i++) iss >> velocity_fr[i];
        if (!iss) continue;

        // Normalize surface elements by the threshold density[cite: 11]
        dsigma[0] = dsigma[0] * density_fr / n0;
        dsigma[1] = dsigma[1] * density_fr / n0;
        dsigma[2] = dsigma[2] * density_fr / n0;
        dsigma[3] = dsigma[3] * density_fr / n0;
        
        // Scalar product of fluid velocity and normal vector determines contributing volume[cite: 11]
        double t_volume_fr = velocity_fr[0] * (dsigma[0] + dsigma[1] * velocity_fr[1] + dsigma[2] * velocity_fr[2] + dsigma[3] * velocity_fr[3]);
        
        // Accept cells close to the particlization criteria with positive outward flow[cite: 11]
        if(t_volume_fr > 0.0 && (std::abs(density_fr - n0) < 0.1 * n0)) {
            volume_fr += t_volume_fr;
            
            CellData cell;
            cell.dsigma = dsigma;
            cell.v = {velocity_fr[0], sqrt(velocity_fr[1]*velocity_fr[1] + velocity_fr[2]*velocity_fr[2]), velocity_fr[3]};
            cell.theta = atan2(velocity_fr[2], velocity_fr[1]);
            
            double e_fr = enden_fr / density_fr;
            cell.valid = LN.get_cell_properties(e_fr, cell.T, cell.mu, cell.mStar);
            
            if (cell.valid) {
                valid_cells.push_back(cell);
            }
        }
    }
    
    std::cout << "Finished reading. Total valid cells: " << valid_cells.size() << std::endl;
    std::cout << "Running calculations on " << omp_get_max_threads() << " CPU threads..." << std::endl;

    // Loop over the rapidity dimension[cite: 11]
    for(int ii = 0; ii < t_y_steps_count; ii++) {
        double t_y = t_y_start + t_y_step * ii;
        
        std::string out_filename = "spectra-pt-" + std::to_string(int(t_y * 100)) + ".dat";
        std::filesystem::path out_file = base_dir / out_filename;
        std::ofstream file_outs(out_filename);

        // Nested loops traversing the transverse momentum (pT) and azimuthal angle (phi) grid[cite: 11]
        for(int i = 0; i < bins_pt; i++) {
            double pt = ptbins[i];
            for(int j = 0; j < bins_phi; j++) {
                double phi = -PI + 2 * PI / bins_phi * (j + 0.5);
                double y = t_y;
                    
                    std::cout << pt << " " << phi << " " << y << "\r" << std::flush;
                    
                    // Master array for this specific momentum/angle state[cite: 11]
                    std::vector<double> t_yield(40, 0.0);
                    
                    // --- OPENMP PARALLEL REGION ---
                    // Drastically accelerates the integration by dividing fluid cells among available CPU cores[cite: 11]
                    #pragma omp parallel
                    {
                        // Thread-local accumulator to avoid race conditions[cite: 11]
                        std::vector<double> local_yield(40, 0.0);
                        
                        // Divide iteration workload across the valid surface elements[cite: 11]
                        #pragma omp for nowait
                        for (size_t c = 0; c < valid_cells.size(); ++c) {
                            const auto& cell = valid_cells[c];
                            std::vector<double> t_emission = LN.emission(pt, phi, y, cell.v, cell.dsigma, cell.theta, cell.T, cell.mu, cell.mStar);
                            for(size_t m = 0; m < local_yield.size(); m++) {
                                local_yield[m] += t_emission[m];
                            }
                        }
                        
                        // Merge thread-local sub-totals into the master yield array safely[cite: 11]
                        #pragma omp critical
                        {
                            for(size_t m = 0; m < t_yield.size(); m++) {
                                t_yield[m] += local_yield[m];
                            }
                        }
                    }
                    // --- END OPENMP REGION ---
                    
                    // Export the calculated phase-space spectra[cite: 11]
                    file_outs << pt << " " << phi << " " << y;
                    for(size_t l = 0; l < t_yield.size(); l++) {
                        file_outs << " " << t_yield[l];
                    }
                    file_outs << std::endl;
                
            }
        }
        file_outs.close();
    }
    
    std::cout << "\nVOLUME: " << volume_fr << std::endl;
    std::cout << "IT`S OVER" << std::endl;
    return 0;
}
