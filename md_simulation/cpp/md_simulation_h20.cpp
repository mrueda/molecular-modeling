#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <random>
#include <iomanip>
#include <cassert>

// Constants
const double PI = 3.14159265358979323846;

// Boltzmann constant in kcal/mol/K
const double KB = 0.001987;

// Simulation Parameters (matching Perl)
std::string integration_method = "verlet";    // Choose "verlet" or "leapfrog"
double dt = 0.002;                            // Time step (fs)
int steps = 1000;                             // Number of steps (match Perl: 1,000 steps)
double epsilon_LJ = 0.1;                      // Lennard-Jones depth (kcal/mol)
double sigma_LJ = 3.0;                        // Lennard-Jones sigma (Å)
double k_angle = 55.0;                        // Angle stiffness (kcal/mol/rad^2)
double theta_eq = 104.5 * PI / 180.0;         // Equilibrium angle in radians
double r_eq_bond = 0.96;                      // Equilibrium bond length (Å)
double k_bond = 450.0;                        // Bond stiffness (kcal/mol/Å^2)
double temperature = 3000.0;                  // Target temperature (K) (Perl code uses 3000 K)
double tolerance = 1e-6;                      // Tolerance for SHAKE
int max_iterations = 100;                     // Max iterations for SHAKE
double offset_distance = 10.0;                // Distance between molecules

// Random number generator for Gaussian distribution
std::mt19937 rng(12345); // Initialize Mersenne Twister with a fixed seed
std::normal_distribution<double> gaussian_dist(0.0, 1.0);

// Helper Functions
double rand_gaussian(double std_dev) {
    return gaussian_dist(rng) * std_dev;
}

double dot_product(const std::vector<double>& v1, const std::vector<double>& v2) {
    assert(v1.size() == v2.size());
    double dp = 0.0;
    for (size_t i = 0; i < v1.size(); ++i) {
        dp += v1[i] * v2[i];
    }
    return dp;
}

// Structures
struct Atom {
    double mass;
    double charge;
    std::vector<double> pos; // [x, y, z]
    std::vector<double> vel; // [vx, vy, vz]

    Atom(double m, double q, double x, double y, double z)
        : mass(m), charge(q), pos{ x, y, z }, vel{ 0.0, 0.0, 0.0 } {}
};

struct Bond {
    int i;
    int j;

    Bond(int atom1, int atom2) : i(atom1), j(atom2) {}
};

struct Angle {
    int i;
    int j;
    int k;

    Angle(int atom1, int atom2, int atom3) : i(atom1), j(atom2), k(atom3) {}
};

// Function Prototypes
void add_molecules(std::vector<Atom>& atoms, std::vector<Bond>& bonds, std::vector<Angle>& angles, int num_molecules);
void initialize_velocities(std::vector<Atom>& atoms, double temperature);
void apply_thermostat(std::vector<Atom>& atoms, double target_temperature);
void compute_bond_forces(const std::vector<Atom>& atoms, const std::vector<Bond>& bonds, std::vector<std::vector<double>>& forces);
void compute_angle_forces(const std::vector<Atom>& atoms, const std::vector<Angle>& angles, std::vector<std::vector<double>>& forces);
void compute_LJ_forces(const std::vector<Atom>& atoms, std::vector<std::vector<double>>& forces);
void apply_shake_constraints(std::vector<Atom>& atoms, const std::vector<Bond>& bonds, double r_eq, double tolerance, int max_iterations);
void verlet_update(std::vector<Atom>& atoms, const std::vector<std::vector<double>>& forces, double dt);
void leapfrog_update(std::vector<Atom>& atoms, const std::vector<std::vector<double>>& forces, double dt);
void write_simulation_data(int step, const std::vector<Atom>& atoms);

int main() {
    // Initialize atoms, bonds, and angles
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<Angle> angles;

    // Add molecules (e.g., water molecules)
    add_molecules(atoms, bonds, angles, 5); // Adding 5 water molecules

    // Initialize velocities once before main loop
    initialize_velocities(atoms, temperature);

    // Initialize forces array
    std::vector<std::vector<double>> forces(atoms.size(), std::vector<double>(3, 0.0));

    // Main MD Loop
    for (int step = 0; step < steps; ++step) {
        // Reset forces
        for (auto& force : forces) {
            std::fill(force.begin(), force.end(), 0.0);
        }

        // Compute forces
        compute_bond_forces(atoms, bonds, forces);
        compute_angle_forces(atoms, angles, forces);
        compute_LJ_forces(atoms, forces);

        // Update positions and velocities based on the selected method
        if (integration_method == "verlet") {
            verlet_update(atoms, forces, dt);
        }
        else if (integration_method == "leapfrog") {
            leapfrog_update(atoms, forces, dt);
        }
        else {
            std::cerr << "Unknown integration method: " << integration_method << std::endl;
            return EXIT_FAILURE;
        }

        // Apply the thermostat every 10 steps (like Perl code)
        if (step % 10 == 0) {
            apply_thermostat(atoms, temperature);
        }

        // Apply SHAKE constraints to fix bond lengths
        apply_shake_constraints(atoms, bonds, r_eq_bond, tolerance, max_iterations);

        // Output positions every 100 steps
        if (step % 100 == 0) {
            write_simulation_data(step, atoms);
        }
    }

    return 0;
}

// Function Definitions

void add_molecules(std::vector<Atom>& atoms, std::vector<Bond>& bonds, std::vector<Angle>& angles, int num_molecules) {
    for (int molecule_id = 0; molecule_id < num_molecules; ++molecule_id) {
        int base_index = (int)atoms.size(); // Starting index for this molecule
        double x_offset = molecule_id * offset_distance;

        // Add atoms (Oxygen, Hydrogen 1, Hydrogen 2)
        atoms.emplace_back(16.0, -0.8, x_offset, 0.0, 0.0);                       // Oxygen
        atoms.emplace_back(1.0, 0.4, x_offset + r_eq_bond, 0.0, 0.0);             // H1
        atoms.emplace_back(1.0, 0.4, x_offset - 0.48, 0.83, 0.0);                 // H2

        // Add bonds
        bonds.emplace_back(base_index, base_index + 1); // O-H1
        bonds.emplace_back(base_index, base_index + 2); // O-H2

        // Add angle
        angles.emplace_back(base_index + 1, base_index, base_index + 2); // H1-O-H2
    }
}

void initialize_velocities(std::vector<Atom>& atoms, double temperature) {
    for (auto& atom : atoms) {
        double std_dev = std::sqrt(KB * temperature / atom.mass); // Standard deviation for velocity
        for (int dim = 0; dim < 3; ++dim) {
            atom.vel[dim] = rand_gaussian(std_dev);
        }
    }
}

void apply_thermostat(std::vector<Atom>& atoms, double target_temperature) {
    // Calculate current kinetic energy
    double total_kinetic_energy = 0.0;
    for (const auto& atom : atoms) {
        double velocity_squared = 0.0;
        for (auto v : atom.vel) {
            velocity_squared += v * v;
        }
        total_kinetic_energy += 0.5 * atom.mass * velocity_squared;
    }

    // Target kinetic energy
    double target_kinetic_energy = 1.5 * KB * target_temperature * atoms.size();

    // Scaling factor for velocities
    double scaling_factor = std::sqrt(target_kinetic_energy / total_kinetic_energy);

    // Scale velocities
    for (auto& atom : atoms) {
        for (int dim = 0; dim < 3; ++dim) {
            atom.vel[dim] *= scaling_factor;
        }
    }
}

void compute_bond_forces(const std::vector<Atom>& atoms, const std::vector<Bond>& bonds, std::vector<std::vector<double>>& forces) {
    for (const auto& bond : bonds) {
        int i = bond.i;
        int j = bond.j;
        std::vector<double> r(3);
        for (int dim = 0; dim < 3; ++dim) {
            r[dim] = atoms[j].pos[dim] - atoms[i].pos[dim];
        }
        double r_mag = std::sqrt(dot_product(r, r));
        double f_mag = -k_bond * (r_mag - r_eq_bond); // Harmonic bond force
        for (int dim = 0; dim < 3; ++dim) {
            double f = f_mag * r[dim] / r_mag;
            forces[i][dim] -= f; // Apply to atom i
            forces[j][dim] += f; // Apply to atom j
        }
    }
}

void compute_angle_forces(const std::vector<Atom>& atoms, const std::vector<Angle>& angles, std::vector<std::vector<double>>& forces) {
    for (const auto& angle : angles) {
        int i = angle.i;
        int j = angle.j;
        int k = angle.k;

        std::vector<double> r1(3), r2(3);
        for (int dim = 0; dim < 3; ++dim) {
            r1[dim] = atoms[i].pos[dim] - atoms[j].pos[dim];
            r2[dim] = atoms[k].pos[dim] - atoms[j].pos[dim];
        }

        double r1_mag = std::sqrt(dot_product(r1, r1));
        double r2_mag = std::sqrt(dot_product(r2, r2));

        if (r1_mag == 0 || r2_mag == 0) continue;

        double cos_theta = dot_product(r1, r2) / (r1_mag * r2_mag);
        cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
        double theta = std::acos(cos_theta);
        double torque = -k_angle * (theta - theta_eq);

        // Placeholder: actual angle forces not fully implemented
        // Implementing angle forces would require distributing the torque
        // onto the three atoms in a physically correct manner.
    }
}

void compute_LJ_forces(const std::vector<Atom>& atoms, std::vector<std::vector<double>>& forces) {
    int num_atoms = (int)atoms.size();
    for (int i = 0; i < num_atoms - 1; ++i) {
        for (int j = i + 1; j < num_atoms; ++j) {
            std::vector<double> r(3);
            for (int dim = 0; dim < 3; ++dim) {
                r[dim] = atoms[j].pos[dim] - atoms[i].pos[dim];
            }
            double r_mag = std::sqrt(dot_product(r, r));
            if (r_mag == 0) continue;
            double r6 = std::pow(sigma_LJ / r_mag, 6);
            double f_mag = 24 * epsilon_LJ * r6 * (2 * r6 - 1) / (r_mag * r_mag);
            for (int dim = 0; dim < 3; ++dim) {
                double f = f_mag * r[dim] / r_mag;
                forces[i][dim] -= f; // Apply to atom i
                forces[j][dim] += f; // Apply to atom j
            }
        }
    }
}

void apply_shake_constraints(std::vector<Atom>& atoms, const std::vector<Bond>& bonds, double r_eq, double tolerance, int max_iterations) {
    for (int iter = 0; iter < max_iterations; ++iter) {
        double max_error = 0.0;
        for (const auto& bond : bonds) {
            int i = bond.i;
            int j = bond.j;
            std::vector<double> r_ij(3);
            for (int dim = 0; dim < 3; ++dim) {
                r_ij[dim] = atoms[j].pos[dim] - atoms[i].pos[dim];
            }
            double r_mag = std::sqrt(dot_product(r_ij, r_ij));
            double error = r_mag - r_eq;
            if (std::abs(error) > std::abs(max_error)) {
                max_error = error;
            }
            if (std::abs(error) > tolerance) {
                double correction = 0.5 * error / r_mag;
                for (int dim = 0; dim < 3; ++dim) {
                    double delta = correction * r_ij[dim];
                    atoms[i].pos[dim] += delta / atoms[i].mass;
                    atoms[j].pos[dim] -= delta / atoms[j].mass;
                }
            }
        }
        if (std::abs(max_error) < tolerance) {
            break;
        }
    }
}

void verlet_update(std::vector<Atom>& atoms, const std::vector<std::vector<double>>& forces, double dt) {
    for (size_t i = 0; i < atoms.size(); ++i) {
        double mass = atoms[i].mass;
        for (int dim = 0; dim < 3; ++dim) {
            atoms[i].pos[dim] += atoms[i].vel[dim] * dt + 0.5 * forces[i][dim] * dt * dt / mass;
            atoms[i].vel[dim] += 0.5 * forces[i][dim] * dt / mass;
        }
    }
}

void leapfrog_update(std::vector<Atom>& atoms, const std::vector<std::vector<double>>& forces, double dt) {
    // Note: A full leapfrog step typically requires forces at half steps.
    // This is a simplified version. The Perl code is similar (also simplified).
    for (size_t i = 0; i < atoms.size(); ++i) {
        double mass = atoms[i].mass;
        for (int dim = 0; dim < 3; ++dim) {
            // Half step velocity update
            atoms[i].vel[dim] += 0.5 * forces[i][dim] * dt / mass;
            // Position update
            atoms[i].pos[dim] += atoms[i].vel[dim] * dt;
            // Half step velocity update again would be done after recomputing forces
        }
    }
}

void write_simulation_data(int step, const std::vector<Atom>& atoms) {
    for (size_t i = 0; i < atoms.size(); ++i) {
        double x = atoms[i].pos[0];
        double y = atoms[i].pos[1];
        double z = atoms[i].pos[2];
        std::cout << step << " " << i << " "
                  << std::fixed << std::setprecision(3)
                  << x << " " << y << " " << z << "\n";
    }
}

