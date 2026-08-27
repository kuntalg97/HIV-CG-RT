import MDAnalysis as mda
import numpy as np

# --- Input ---
traj_file = "traj.lammpstrj"   # your LAMMPS dump file
trimer_type = 40                     # atom type to analyze
dimer_type = 37
cutoff = 20.0                        # in Angstrom
output_1 = "cn_trimer.dat"
output_2 = "cn_dimer.dat"

# --- Load trajectory ---
u = mda.Universe(traj_file, format="LAMMPSDUMP")

# --- Select target atoms ---
trimer_atoms = u.select_atoms(f"type {trimer_type}")
dimer_atoms = u.select_atoms(f"type {dimer_type}")

cn_tri_per_frame = []
cn_di_per_frame = []

# --- Loop over all frames ---
for ts in u.trajectory:
    # Compute pairwise distances among target atoms
    dists_tri = np.array(trimer_atoms.positions)
    dists_di  = np.array(dimer_atoms.positions)
    diff_tri = dists_tri[:, np.newaxis, :] - dists_tri[np.newaxis, :, :]
    diff_di = dists_di[:, np.newaxis, :] - dists_di[np.newaxis, :, :]
    dist_matrix_tri = np.linalg.norm(diff_tri, axis=-1)
    dist_matrix_di  = np.linalg.norm(diff_di, axis=-1)
    
    # Remove self-distances
    np.fill_diagonal(dist_matrix_tri, np.inf)
    np.fill_diagonal(dist_matrix_di, np.inf)
    
    # Count neighbors within cutoff
    cn_each_atom_tri = np.sum(dist_matrix_tri < cutoff, axis=1)
    cn_each_atom_di  = np.sum(dist_matrix_di < cutoff, axis=1)
    avg_cn_tri = np.mean(cn_each_atom_tri)
    avg_cn_di  = np.mean(cn_each_atom_di)
    
    cn_tri_per_frame.append(avg_cn_tri)
    cn_di_per_frame.append(avg_cn_di)

# --- Save to file ---
frames = np.arange(len(cn_tri_per_frame))

np.savetxt(output_1, np.column_stack((frames, cn_tri_per_frame)),
           fmt="%d %.6f", header="Frame  CN")

np.savetxt(output_2, np.column_stack((frames, cn_di_per_frame)),
           fmt="%d %.6f", header="Frame  CN")
