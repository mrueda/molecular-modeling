#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <limits>
#include <sstream> // Added to use std::ostringstream

// Constants
const double PI = 3.14159265358979323846;

// Structure to represent an Atom (Amino Acid)
struct Atom {
    std::string name;
    double x;
    double y;
    double z;

    Atom(std::string n, double xpos, double ypos, double zpos)
        : name(n), x(xpos), y(ypos), z(zpos) {}
};

// Structure to represent a Docking Pose
struct Pose {
    std::vector<Atom> docked_peptide2;
    double rotation_deg;
    double translation_x;
    double translation_y;
    double score;

    Pose() : rotation_deg(0.0), translation_x(0.0), translation_y(0.0), score(std::numeric_limits<double>::max()) {}
};

// Function to calculate Euclidean distance between two atoms
double calculate_distance(const Atom& a1, const Atom& a2) {
    double dx = a1.x - a2.x;
    double dy = a1.y - a2.y;
    double dz = a1.z - a2.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Function to calculate docking score (minimum distance between any pair of amino acids)
double calculate_score(const std::vector<Atom>& peptide1, const std::vector<Atom>& docked_peptide2) {
    double min_dist = std::numeric_limits<double>::max();
    for (const auto& aa1 : peptide1) {
        for (const auto& aa2 : docked_peptide2) {
            double dist = calculate_distance(aa1, aa2);
            if (dist < min_dist) {
                min_dist = dist;
            }
        }
    }
    return min_dist;
}

// Function to apply rotation around Z-axis and translate
Atom rotate_atom(const Atom& aa, double cos_a, double sin_a, double tx, double ty) {
    double x_rot = aa.x * cos_a - aa.y * sin_a;
    double y_rot = aa.x * sin_a + aa.y * cos_a;
    double z_rot = aa.z; // No change in Z for rotation around Z-axis

    // Apply translation
    double x_new = x_rot + tx;
    double y_new = y_rot + ty;
    double z_new = z_rot;

    return Atom(aa.name, x_new, y_new, z_new);
}

// Function to print peptide coordinates to console
void print_peptide(const std::vector<Atom>& peptide) {
    for (const auto& aa : peptide) {
        std::cout << std::left << std::setw(5) << aa.name << ": ("
                  << std::fixed << std::setprecision(2)
                  << aa.x << ", " << aa.y << ", " << aa.z << ")\n";
    }
}

int main() {
    // ============================
    // Define Peptides
    // ============================

    // Peptide 1: Simplified coordinates for 5 amino acids (e.g., Alanine)
    std::vector<Atom> peptide1 = {
        Atom("Ala1", 0.0, 0.0, 0.0),
        Atom("Ala2", 1.5, 0.0, 0.0),
        Atom("Ala3", 3.0, 0.0, 0.0),
        Atom("Ala4", 4.5, 0.0, 0.0),
        Atom("Ala5", 6.0, 0.0, 0.0)
    };

    // Peptide 2: Simplified coordinates for 5 amino acids (e.g., Glycine)
    std::vector<Atom> peptide2 = {
        Atom("Gly1", 0.0, 0.0, 5.0),
        Atom("Gly2", 1.5, 0.0, 5.0),
        Atom("Gly3", 3.0, 0.0, 5.0),
        Atom("Gly4", 4.5, 0.0, 5.0),
        Atom("Gly5", 6.0, 0.0, 5.0)
    };

    // ============================
    // Docking Parameters
    // ============================

    std::vector<int> rot_angles = {0, 45, 90, 135, 180, 225, 270, 315}; // Degrees around Z-axis
    double translation_step = 1.0;  // Ångströms
    double max_translation = 5.0;    // Maximum translation in each direction (X and Y)

    // Output file for trajectory
    std::string trajectory_file = "docking_trajectory.xyz";

    // Open the trajectory file for writing
    std::ofstream traj_fh(trajectory_file);
    if (!traj_fh.is_open()) {
        std::cerr << "Cannot open " << trajectory_file << " for writing.\n";
        return 1;
    }

    // ============================
    // Docking Process
    // ============================

    Pose best_pose;
    int frame_counter = 0;

    for (const auto& angle_deg : rot_angles) {
        double angle_rad = angle_deg * PI / 180.0;
        double cos_a = std::cos(angle_rad);
        double sin_a = std::sin(angle_rad);

        for (double tx = -max_translation; tx <= max_translation; tx += translation_step) {
            for (double ty = -max_translation; ty <= max_translation; ty += translation_step) {

                // Apply rotation and translation to Peptide 2
                std::vector<Atom> docked_peptide2;
                for (const auto& aa : peptide2) {
                    Atom rotated_translated_aa = rotate_atom(aa, cos_a, sin_a, tx, ty);
                    docked_peptide2.push_back(rotated_translated_aa);
                }

                // Calculate score for this pose
                double score = calculate_score(peptide1, docked_peptide2);

                // Update best score and pose if necessary
                if (score < best_pose.score) {
                    best_pose.score = score;
                    best_pose.docked_peptide2 = docked_peptide2;
                    best_pose.rotation_deg = angle_deg;
                    best_pose.translation_x = tx;
                    best_pose.translation_y = ty;
                }

                // Log this pose to the trajectory file
                // Each frame in XYZ requires:
                // 1. Number of atoms
                // 2. Comment line
                // 3. Atom coordinates

                int num_atoms = peptide1.size() + docked_peptide2.size();
                std::ostringstream comment;
                comment << "Frame " << frame_counter << ": Rotation=" << angle_deg << "°, "
                        << "Translation=(" << tx << ", " << ty << ") Å, Score=" << score;

                traj_fh << num_atoms << "\n" << comment.str() << "\n";

                // Print Peptide1 atoms
                for (const auto& aa1 : peptide1) {
                    traj_fh << std::left << std::setw(5) << aa1.name << " "
                            << std::fixed << std::setprecision(3)
                            << std::setw(8) << aa1.x << " "
                            << std::setw(8) << aa1.y << " "
                            << std::setw(8) << aa1.z << "\n";
                }

                // Print Docked Peptide2 atoms
                for (const auto& aa2 : docked_peptide2) {
                    traj_fh << std::left << std::setw(5) << aa2.name << " "
                            << std::fixed << std::setprecision(3)
                            << std::setw(8) << aa2.x << " "
                            << std::setw(8) << aa2.y << " "
                            << std::setw(8) << aa2.z << "\n";
                }

                frame_counter++;
            }
        }
    }

    // Close the trajectory file
    traj_fh.close();
    std::cout << "Trajectory saved to " << trajectory_file << "\n";

    // ============================
    // Output Best Docking Pose
    // ============================

    std::cout << "Best Docking Score: " << best_pose.score << " Å\n";
    std::cout << "Best Docking Rotation: " << best_pose.rotation_deg << " degrees\n";
    std::cout << "Best Docking Translation: X=" << best_pose.translation_x << " Å, Y=" << best_pose.translation_y << " Å\n\n";

    std::cout << "Peptide 1 Coordinates:\n";
    print_peptide(peptide1);
    std::cout << "\nBest Docked Peptide 2 Coordinates:\n";
    print_peptide(best_pose.docked_peptide2);

    return 0;
}
