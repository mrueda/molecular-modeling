# Molecular Modeling Examples

Small educational examples for molecular simulation, peptide docking, plotting, and Fourier transforms in Perl, C++, and Python.

These scripts are intended for learning and experimentation. They are not a production molecular dynamics or docking engine, and the simplified force fields and scoring functions should not be used for scientific conclusions without validation.

## Contents

- `md_simulation/`: simplified water molecular dynamics in Perl and C++.
- `docking/`: simplified peptide docking search in Perl and C++ plus a Python trajectory viewer.
- `misc/discrete_fourier_transform/`: Perl DFT and recursive FFT examples.

## Requirements

- `g++` with C++11 support.
- Perl 5 with core math modules.
- Python 3 for syntax checks and plotting.
- Optional Python packages for plotting: `matplotlib`, `pandas`, and `numpy`.

## Build and Check

```sh
make build
make check
```

`make check` compiles the C++ examples with warnings enabled, checks Perl syntax, and byte-compiles the Python plotting scripts.

## Run Examples

```sh
make run-md-cpp
make run-md-perl
make run-docking-cpp
make run-docking-perl
```

Generated binaries, trajectories, and simulation outputs are ignored by git. Remove them with:

```sh
make clean
```

## Plot Outputs

Water trajectory:

```sh
python3 md_simulation/plots/plot_3d.py md_simulation/cpp/simulation_output.txt
```

Docking trajectory:

```sh
python3 docking/plots/plot_docking.py docking/cpp/docking_trajectory.xyz
```

## Notes on the MD Example

The MD example includes harmonic bond forces, harmonic angle forces, Lennard-Jones interactions, a simple velocity-rescaling thermostat, and SHAKE-like bond constraints. The integration and constraints are intentionally compact for readability.
