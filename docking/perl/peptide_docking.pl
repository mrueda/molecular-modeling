#!/usr/bin/env perl
use strict;
use warnings;
use Math::Trig;      # For trigonometric functions
use Math::Complex;  # For complex numbers (optional)

# ============================
# Define Peptides
# ============================

# Peptide 1: Simplified coordinates for 5 amino acids (e.g., Alanine)
my @peptide1 = (
    { name => 'Ala1', x => 0.0, y => 0.0, z => 0.0 },
    { name => 'Ala2', x => 1.5, y => 0.0, z => 0.0 },
    { name => 'Ala3', x => 3.0, y => 0.0, z => 0.0 },
    { name => 'Ala4', x => 4.5, y => 0.0, z => 0.0 },
    { name => 'Ala5', x => 6.0, y => 0.0, z => 0.0 },
);

# Peptide 2: Simplified coordinates for 5 amino acids (e.g., Glycine)
my @peptide2 = (
    { name => 'Gly1', x => 0.0, y => 0.0, z => 5.0 },
    { name => 'Gly2', x => 1.5, y => 0.0, z => 5.0 },
    { name => 'Gly3', x => 3.0, y => 0.0, z => 5.0 },
    { name => 'Gly4', x => 4.5, y => 0.0, z => 5.0 },
    { name => 'Gly5', x => 6.0, y => 0.0, z => 5.0 },
);

# ============================
# Docking Parameters
# ============================

my @rot_angles = (0, 45, 90, 135, 180, 225, 270, 315); # Degrees around Z-axis
my $translation_step = 1.0;  # Ångströms
my $max_translation = 5.0;   # Maximum translation in each direction (X and Y)
my $translation_range = (-$max_translation .. $max_translation); # Not directly usable in Perl, will handle in loops

# Output file for trajectory
my $trajectory_file = 'docking_trajectory.xyz';

# Open the trajectory file for writing
open(my $traj_fh, '>', $trajectory_file) or die "Cannot open $trajectory_file for writing: $!";

# ============================
# Docking Process
# ============================

my $best_score = 1e10;  # Initialize with a large number
my @best_pose = ();
my $best_rotation = 0;
my ($best_tx, $best_ty) = (0, 0);

# Initialize a frame counter for XYZ file
my $frame_counter = 0;

# Iterate over rotation angles
foreach my $angle_deg (@rot_angles) {
    my $angle_rad = deg2rad($angle_deg);
    my ($cos_a, $sin_a) = (cos($angle_rad), sin($angle_rad));

    # Iterate over translations in X direction
    for (my $tx = -$max_translation; $tx <= $max_translation; $tx += $translation_step) {
        # Iterate over translations in Y direction
        for (my $ty = -$max_translation; $ty <= $max_translation; $ty += $translation_step) {

            # Apply rotation and translation to Peptide 2
            my @docked_peptide2 = ();
            foreach my $aa (@peptide2) {
                # Rotate around Z-axis
                my $x_rot = $aa->{x} * $cos_a - $aa->{y} * $sin_a;
                my $y_rot = $aa->{x} * $sin_a + $aa->{y} * $cos_a;
                my $z_rot = $aa->{z};  # No change in Z for rotation around Z-axis

                # Translate
                my $x_new = $x_rot + $tx;
                my $y_new = $y_rot + $ty;
                my $z_new = $z_rot;

                push @docked_peptide2, { name => $aa->{name}, x => $x_new, y => $y_new, z => $z_new };
            }

            # Calculate score for this pose
            my $score = calculate_score(\@peptide1, \@docked_peptide2);

            # Update best score and pose if necessary
            if ($score < $best_score) {
                $best_score = $score;
                @best_pose = @docked_peptide2;
                $best_rotation = $angle_deg;
                ($best_tx, $best_ty) = ($tx, $ty);
            }

            # Log this pose to the trajectory file
            # Each frame in XYZ requires:
            # 1. Number of atoms
            # 2. Comment line
            # 3. Atom coordinates

            # Total number of atoms: 5 (Peptide1) + 5 (Docked Peptide2) = 10
            my $num_atoms = scalar(@peptide1) + scalar(@docked_peptide2);
            my $comment = "Frame $frame_counter: Rotation=$angle_deg°, Translation=($tx, $ty) Å, Score=$score";
            print $traj_fh "$num_atoms\n$comment\n";

            # Print Peptide1 atoms
            foreach my $aa1 (@peptide1) {
                printf $traj_fh "%-4s %8.3f %8.3f %8.3f\n", $aa1->{name}, $aa1->{x}, $aa1->{y}, $aa1->{z};
            }

            # Print Docked Peptide2 atoms
            foreach my $aa2 (@docked_peptide2) {
                printf $traj_fh "%-4s %8.3f %8.3f %8.3f\n", $aa2->{name}, $aa2->{x}, $aa2->{y}, $aa2->{z};
            }

            $frame_counter++;
        }
    }
}

# Close the trajectory file
close($traj_fh);
print "Trajectory saved to $trajectory_file\n";

# ============================
# Output Best Docking Pose
# ============================

print "Best Docking Score: $best_score Å\n";
print "Best Docking Rotation: $best_rotation degrees\n";
print "Best Docking Translation: X=$best_tx Å, Y=$best_ty Å\n";
print "\nPeptide 1 Coordinates:\n";
print_peptide(\@peptide1);
print "\nBest Docked Peptide 2 Coordinates:\n";
print_peptide(\@best_pose);

# ============================
# Subroutines
# ============================

# Calculate the docking score based on minimum distance between any pair of amino acids
sub calculate_score {
    my ($pep1_ref, $pep2_ref) = @_;
    my $min_dist = 1e10;

    foreach my $aa1 (@$pep1_ref) {
        foreach my $aa2 (@$pep2_ref) {
            my $dist = euclidean_distance($aa1->{x}, $aa1->{y}, $aa1->{z},
                                         $aa2->{x}, $aa2->{y}, $aa2->{z});
            if ($dist < $min_dist) {
                $min_dist = $dist;
            }
        }
    }
    return $min_dist;
}

# Calculate Euclidean distance between two points
sub euclidean_distance {
    my ($x1, $y1, $z1, $x2, $y2, $z2) = @_;
    return sqrt( ($x1 - $x2)**2 + ($y1 - $y2)**2 + ($z1 - $z2)**2 );
}

# Print peptide coordinates
sub print_peptide {
    my ($pep_ref) = @_;
    foreach my $aa (@$pep_ref) {
        printf "%-5s: (%.2f, %.2f, %.2f)\n", $aa->{name}, $aa->{x}, $aa->{y}, $aa->{z};
    }
}

