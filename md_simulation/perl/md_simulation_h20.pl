#!/usr/bin/env perl
use strict;
use warnings;
use Math::Complex;

# Integration method: Choose 'verlet' or 'leapfrog'
my $integration_method = 'verlet';    # Change to 'verlet' for Verlet integration

# Simulation parameters
my $dt          = 0.002;               # Time step (fs)
my $steps       = 1_000;               # Number of steps
my $epsilon_LJ  = 0.1;                 # Lennard-Jones depth (kcal/mol)
my $sigma_LJ    = 3.0;                 # Lennard-Jones sigma (Å)
my $k_angle     = 55.0;                # Angle stiffness (kcal/mol/rad^2)
my $theta_eq    = 104.5 * pi / 180;    # Equilibrium angle in radians
my $r_eq_bond   = 0.96;                # Equilibrium bond length (Å)
my $k_bond      = 450.0;               # Bond stiffness (kcal/mol/Å^2)
my $temperature = 3000;                 # Target temperature (K)
my $kb          = 0.001987;            # Boltzmann constant (kcal/mol/K)

# Set a fixed seed
srand(12345);

# Initialize atom positions, velocities, and charges for N molecules
# Declare global variables
my @atoms;
my @bonds;
my @angles;

# Add 5 water molecules
add_molecules( \@atoms, \@bonds, \@angles, 5 );

# Forces array
my @forces = ();

# Assign random initial velocities based on temperature
initialize_velocities( \@atoms, $temperature );

# Main MD loop
for my $step ( 0 .. $steps - 1 ) {

    # Reset forces
    @forces = map { [ 0, 0, 0 ] } @atoms;

    # Compute forces
    compute_bond_forces( \@atoms, \@bonds, \@forces );
    compute_angle_forces( \@atoms, \@angles, \@forces );
    compute_LJ_forces( \@atoms, \@forces );    # Non-bonded interactions across all atoms

    # Update positions and velocities based on the selected method
    if ( $integration_method eq 'verlet' ) {
        verlet_update( \@atoms, \@forces, $dt );
    }
    elsif ( $integration_method eq 'leapfrog' ) {
        leapfrog_update( \@atoms, \@forces, $dt );
    }
    else {
        die "Unknown integration method: $integration_method\n";
    }

    # Apply the thermostat every 10 steps
    if ($step % 10 == 0) {
        apply_thermostat(\@atoms, $temperature);
    }

    # Apply SHAKE constraints to fix bond lengths
    apply_shake_constraints( \@atoms, \@bonds, $r_eq_bond, 1e-6, 100 );

    # Output positions every 100 steps
    if ( $step % 100 == 0 ) {
        write_simulation_data( $step, \@atoms );
    }
}

# Subroutine: Initialize random velocities based on temperature
sub initialize_velocities {
    my ( $atoms, $temperature ) = @_;
    foreach my $atom (@$atoms) {
        my $mass    = $atom->{mass};
        my $std_dev = sqrt( $kb * $temperature / $mass );    # Standard deviation for velocity
        foreach my $dim ( 0 .. 2 ) {
            $atom->{vel}[$dim] = rand_gaussian($std_dev);
        }
    }
}

# Subroutine: Gaussian random number generator
sub rand_gaussian {
    my ($std_dev) = @_;
    my $u1        = rand();
    my $u2        = rand();
    my $z         = sqrt( -2 * log($u1) ) * cos( 2 * pi * $u2 );    # Box-Muller transform
    return $z * $std_dev;
}

# Subroutine: Apply a velocity rescaling thermostat
sub apply_thermostat {
    my ( $atoms, $target_temperature ) = @_;

    # Calculate current kinetic energy
    my $total_kinetic_energy = 0;
    foreach my $atom (@$atoms) {
        my $mass             = $atom->{mass};
        my $velocity_squared = sum( map { $_**2 } @{ $atom->{vel} } );
        $total_kinetic_energy += 0.5 * $mass * $velocity_squared;
    }

    # Target kinetic energy
    my $target_kinetic_energy =
      1.5 * $kb * $target_temperature * scalar(@$atoms);

    # Scaling factor for velocities
    my $scaling_factor = sqrt( $target_kinetic_energy / $total_kinetic_energy );

    # Scale velocities
    foreach my $atom (@$atoms) {
        foreach my $dim ( 0 .. 2 ) {
            $atom->{vel}[$dim] *= $scaling_factor;
        }
    }
}

# Subroutine: Compute bond forces (harmonic potential)
sub compute_bond_forces {
    my ( $atoms, $bonds, $forces ) = @_;
    for my $bond (@$bonds) {
        my ( $i, $j ) = @$bond;
        my @r =
          map { $atoms->[$j]->{pos}[$_] - $atoms->[$i]->{pos}[$_] } ( 0 .. 2 );
        my $r_mag = sqrt( sum( map { $_**2 } @r ) );
        my $f_mag = -$k_bond * ( $r_mag - $r_eq_bond );    # Harmonic bond force
        for my $dim ( 0 .. 2 ) {
            my $f = $f_mag * $r[$dim] / $r_mag;
            $forces->[$i][$dim] -= $f;                     # Apply to atom i
            $forces->[$j][$dim] += $f;                     # Apply to atom j
        }
    }
}

# Subroutine: Compute angle forces (harmonic potential)
sub compute_angle_forces {
    my ( $atoms, $angles, $forces ) = @_;
    for my $angle (@$angles) {
        my ( $i, $j, $k ) = @$angle;
        my @r1 =
          map { $atoms->[$i]->{pos}[$_] - $atoms->[$j]->{pos}[$_] } ( 0 .. 2 );
        my @r2 =
          map { $atoms->[$k]->{pos}[$_] - $atoms->[$j]->{pos}[$_] } ( 0 .. 2 );
        my $r1_mag    = sqrt( sum( map { $_**2 } @r1 ) );
        my $r2_mag    = sqrt( sum( map { $_**2 } @r2 ) );
        my $cos_theta = dot_product( \@r1, \@r2 ) / ( $r1_mag * $r2_mag );
        my $theta     = Math::Complex::acos($cos_theta);
        my $torque    = -$k_angle * ( $theta - $theta_eq );

        # Force calculations here (not fully implemented in this snippet)
    }
}

# Subroutine: Compute Lennard-Jones forces
sub compute_LJ_forces {
    my ( $atoms, $forces ) = @_;
    for my $i ( 0 .. $#$atoms - 1 ) {
        for my $j ( $i + 1 .. $#$atoms ) {
            my @r = map { $atoms->[$j]->{pos}[$_] - $atoms->[$i]->{pos}[$_] }
              ( 0 .. 2 );
            my $r_mag = sqrt( sum( map { $_**2 } @r ) );
            my $r6    = ( $sigma_LJ / $r_mag )**6;
            my $f_mag = 24 * $epsilon_LJ * $r6 * ( 2 * $r6 - 1 ) / $r_mag**2;
            for my $dim ( 0 .. 2 ) {
                my $f = $f_mag * $r[$dim] / $r_mag;
                $forces->[$i][$dim] -= $f;    # Apply to atom i
                $forces->[$j][$dim] += $f;    # Apply to atom j
            }
        }
    }
}

# Subroutine: Apply SHAKE constraints
sub apply_shake_constraints {
    my ( $atoms, $bonds, $r_eq, $tolerance, $max_iterations ) = @_;
    for my $iteration ( 1 .. $max_iterations ) {
        my $max_error = 0;
        for my $bond (@$bonds) {
            my ( $i, $j ) = @$bond;
            my @r_ij = map { $atoms->[$j]->{pos}[$_] - $atoms->[$i]->{pos}[$_] }
              ( 0 .. 2 );
            my $r_mag = sqrt( sum( map { $_**2 } @r_ij ) );
            my $error = $r_mag - $r_eq;
            $max_error = $error if abs($error) > abs($max_error);
            if ( abs($error) > $tolerance ) {
                my $correction = 0.5 * $error / $r_mag;
                for my $dim ( 0 .. 2 ) {
                    my $delta = $correction * $r_ij[$dim];
                    $atoms->[$i]->{pos}[$dim] += $delta / $atoms->[$i]->{mass};
                    $atoms->[$j]->{pos}[$dim] -= $delta / $atoms->[$j]->{mass};
                }
            }
        }
        last if abs($max_error) < $tolerance;
    }
}

# Verlet integration
sub verlet_update {
    my ( $atoms, $forces, $dt ) = @_;
    for my $i ( 0 .. $#$atoms ) {
        my $mass = $atoms->[$i]->{mass};
        for my $dim ( 0 .. 2 ) {
            $atoms->[$i]->{pos}[$dim] += $atoms->[$i]->{vel}[$dim] * $dt +
              0.5 * $forces->[$i][$dim] * $dt**2 / $mass;
            $atoms->[$i]->{vel}[$dim] +=
              0.5 * $forces->[$i][$dim] * $dt / $mass;
        }
    }
}

# Subroutine: Leapfrog update
sub leapfrog_update {
    my ($atoms, $forces, $dt) = @_;

    for my $i (0 .. $#$atoms) {
        my $atom = $atoms->[$i];
        my $mass = $atom->{mass};

        # Update velocities (half step)
        for my $dim (0..2) {
            $atom->{vel}[$dim] += 0.5 * $forces->[$i][$dim] * $dt / $mass;
        }

        # Update positions (full step)
        for my $dim (0..2) {
            $atom->{pos}[$dim] += $atom->{vel}[$dim] * $dt;
        }

        # Update velocities (another half step)
        for my $dim (0..2) {
            $atom->{vel}[$dim] += 0.5 * $forces->[$i][$dim] * $dt / $mass;
        }
    }
}

# Subroutine: Write simulation data
sub write_simulation_data {
    my ( $step, $atoms ) = @_;

    # Loop through all atoms in the @atoms array
    for my $i ( 0 .. $#$atoms ) {
        my $x = $atoms->[$i]->{pos}[0];
        my $y = $atoms->[$i]->{pos}[1];
        my $z = $atoms->[$i]->{pos}[2];

        # Print step, atom ID, and coordinates
        printf "%d %d %.3f %.3f %.3f\n", $step, $i, $x, $y, $z;
    }
}

sub add_molecules {
    my ( $atoms_ref, $bonds_ref, $angles_ref, $num_molecules ) = @_;
    my $offset = 10;    # Distance between molecules

    for my $molecule_id ( 0 .. $num_molecules - 1 ) {
        my $base_index = @$atoms_ref;              # Starting index for this molecule
        my $x_offset   = $molecule_id * $offset;

        # Add atoms (Oxygen, Hydrogen 1, Hydrogen 2)
        push @$atoms_ref,
          {
            mass   => 16.0,
            charge => -0.8,
            pos    => [ $x_offset, 0, 0 ],
            vel    => [ 0,         0, 0 ]
          };                                       # Oxygen
        push @$atoms_ref,
          {
            mass   => 1.0,
            charge => 0.4,
            pos    => [ $x_offset + 0.96, 0, 0 ],
            vel    => [ 0,                0, 0 ]
          };                                       # H1
        push @$atoms_ref,
          {
            mass   => 1.0,
            charge => 0.4,
            pos    => [ $x_offset - 0.48, 0.83, 0 ],
            vel    => [ 0,                0,    0 ]
          };                                       # H2

        # Add bonds
        push @$bonds_ref, [ $base_index, $base_index + 1 ];    # O-H1
        push @$bonds_ref, [ $base_index, $base_index + 2 ];    # O-H2

        # Add angle
        push @$angles_ref, [ $base_index + 1, $base_index, $base_index + 2 ];  # H1-O-H2
    }
}

# Helper functions
sub sum { my $sum = 0; $sum += $_ for @_; return $sum; }

sub dot_product {
    my ( $v1, $v2 ) = @_;
    return sum( map { $v1->[$_] * $v2->[$_] } 0 .. $#$v1 );
}
