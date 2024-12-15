import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.animation as animation

# ============================
# Configuration Parameters
# ============================

# Path to the XYZ trajectory file
TRAJECTORY_FILE = 'docking_trajectory.xyz'

# Number of atoms per peptide
ATOMS_PEPTIDE1 = 5  # e.g., 5 Alanine residues
ATOMS_PEPTIDE2 = 5  # e.g., 5 Glycine residues

# Total number of atoms per frame
TOTAL_ATOMS = ATOMS_PEPTIDE1 + ATOMS_PEPTIDE2

# Colors for peptides
COLOR_PEPTIDE1 = 'blue'
COLOR_PEPTIDE2 = 'red'

# Markers for amino acids
MARKER_PEPTIDE1 = 'o'  # Circle
MARKER_PEPTIDE2 = '^'  # Triangle

# ============================
# Function to Read XYZ File
# ============================

def read_xyz(filename):
    """
    Reads an XYZ file and returns a list of frames.
    Each frame is a list of atoms, where each atom is a dictionary with 'name', 'x', 'y', 'z'.
    """
    frames = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        i = 0
        while i < len(lines):
            try:
                num_atoms = int(lines[i].strip())
            except ValueError:
                print(f"Error parsing number of atoms at line {i+1}")
                break
            comment = lines[i+1].strip()
            frame = []
            for j in range(i+2, i+2+num_atoms):
                parts = lines[j].strip().split()
                if len(parts) < 4:
                    print(f"Error parsing atom coordinates at line {j+1}")
                    continue
                atom = {
                    'name': parts[0],
                    'x': float(parts[1]),
                    'y': float(parts[2]),
                    'z': float(parts[3])
                }
                frame.append(atom)
            frames.append(frame)
            i += 2 + num_atoms
    return frames

# ============================
# Function to Separate Peptides
# ============================

def separate_peptides(frame):
    """
    Separates atoms into two peptides based on their indices.
    Returns two lists: peptide1_atoms and peptide2_atoms.
    """
    peptide1_atoms = frame[:ATOMS_PEPTIDE1]
    peptide2_atoms = frame[ATOMS_PEPTIDE1:TOTAL_ATOMS]
    return peptide1_atoms, peptide2_atoms

# ============================
# Function to Initialize Plot
# ============================

def init_plot(ax):
    """
    Initializes the 3D plot.
    """
    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    ax.set_title('Peptide Docking Simulation')
    ax.set_xlim(-10, 20)
    ax.set_ylim(-10, 20)
    ax.set_zlim(-10, 20)
    return

# ============================
# Function to Update Plot
# ============================

def update_plot(num, frames, scat_pep1, scat_pep2, ax):
    """
    Updates the scatter plots for each frame.
    """
    frame = frames[num]
    peptide1, peptide2 = separate_peptides(frame)
    
    # Extract coordinates for Peptide 1
    x1 = [atom['x'] for atom in peptide1]
    y1 = [atom['y'] for atom in peptide1]
    z1 = [atom['z'] for atom in peptide1]
    
    # Extract coordinates for Peptide 2
    x2 = [atom['x'] for atom in peptide2]
    y2 = [atom['y'] for atom in peptide2]
    z2 = [atom['z'] for atom in peptide2]
    
    # Update scatter plots
    scat_pep1._offsets3d = (x1, y1, z1)
    scat_pep2._offsets3d = (x2, y2, z2)
    
    # Optional: Update the title with frame information
    ax.set_title(f'Peptide Docking Simulation - Frame {num+1}/{len(frames)}')
    
    return scat_pep1, scat_pep2

# ============================
# Main Function
# ============================

def main():
    # Read the trajectory frames
    frames = read_xyz(TRAJECTORY_FILE)
    if not frames:
        print("No frames found in the trajectory file.")
        return
    
    print(f"Total frames loaded: {len(frames)}")
    
    # Set up the figure and 3D axis
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    init_plot(ax)
    
    # Initialize scatter plots for both peptides
    initial_pep1, initial_pep2 = separate_peptides(frames[0])
    x1 = [atom['x'] for atom in initial_pep1]
    y1 = [atom['y'] for atom in initial_pep1]
    z1 = [atom['z'] for atom in initial_pep1]
    
    x2 = [atom['x'] for atom in initial_pep2]
    y2 = [atom['y'] for atom in initial_pep2]
    z2 = [atom['z'] for atom in initial_pep2]
    
    scat_pep1 = ax.scatter(x1, y1, z1, c=COLOR_PEPTIDE1, marker=MARKER_PEPTIDE1, label='Peptide 1')
    scat_pep2 = ax.scatter(x2, y2, z2, c=COLOR_PEPTIDE2, marker=MARKER_PEPTIDE2, label='Peptide 2')
    
    ax.legend()
    
    # Create the animation
    ani = animation.FuncAnimation(fig, update_plot, frames=len(frames),
                                  fargs=(frames, scat_pep1, scat_pep2, ax),
                                  interval=200, blit=False)
    
    # Save the animation as a GIF or MP4 (optional)
    # Uncomment the following lines if you want to save the animation
    # ani.save('docking_animation.gif', writer='imagemagick', fps=5)
    # ani.save('docking_animation.mp4', writer='ffmpeg', fps=5)
    
    # Display the plot
    plt.show()

if __name__ == "__main__":
    main()

