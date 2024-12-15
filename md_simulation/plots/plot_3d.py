import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import itertools
import sys

def plot_water_molecules(input_file):
    # Load the data with manually assigned column names
    data = pd.read_csv(
        input_file,
        sep=r'\s+',  # Use regex to split whitespace
        names=['Step', 'Atom', 'X', 'Y', 'Z']  # Assign column names
    )

    # Check the structure of the data
    print("Data preview:")
    print(data.head())

    # Assume water molecules with 3 atoms each
    atoms_per_molecule = 3
    num_atoms = data['Atom'].nunique()
    num_molecules = num_atoms // atoms_per_molecule

    # Create molecule mapping: Assign each atom to a molecule
    molecule_mapping = {
        atom: f"Molecule {atom // atoms_per_molecule + 1}"
        for atom in range(num_atoms)
    }

    # Generate distinct colors for each molecule
    colors = itertools.cycle(plt.cm.tab20.colors)  # Use colormap with 20 unique colors
    molecule_colors = {f"Molecule {i + 1}": next(colors) for i in range(num_molecules)}

    # Initialize 3D plot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Plot each molecule with its unique color
    for molecule, color in molecule_colors.items():
        molecule_atoms = [
            atom for atom, mol in molecule_mapping.items() if mol == molecule
        ]
        for atom in molecule_atoms:
            atom_data = data[data['Atom'] == atom]
            ax.plot(
                atom_data['X'], atom_data['Y'], atom_data['Z'],
                label=f"{molecule} Atom {atom}",
                color=color,
                alpha=0.7
            )

    # Add labels, title, legend, and grid
    ax.set_xlabel('X Coordinate (Å)')
    ax.set_ylabel('Y Coordinate (Å)')
    ax.set_zlabel('Z Coordinate (Å)')
    ax.set_title(f"3D Trajectories of {num_molecules} Water Molecules")
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1), fontsize='small')
    plt.tight_layout()
    plt.show()


# Example usage: python script.py simulation_output.txt
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_water_molecules.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    plot_water_molecules(input_file)
