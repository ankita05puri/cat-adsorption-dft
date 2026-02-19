from ase.io import read
import numpy as np

def main():
    atoms = read("outputs/co_on_top.xyz")

    # Indices
    c_indices = [i for i, a in enumerate(atoms) if a.symbol == "C"]
    pt_indices = [i for i, a in enumerate(atoms) if a.symbol == "Pt"]

    if not c_indices:
        raise RuntimeError("No carbon atoms found in co_on_top.xyz")
    if not pt_indices:
        raise RuntimeError("No Pt atoms found in co_on_top.xyz")

    # Use the last carbon (CO you added)
    c_idx = c_indices[-1]

    # --- Min C–Pt distance (IMPORTANT: mic=False to avoid vacuum wrap artifacts) ---
    dmin = min(atoms.get_distance(c_idx, j, mic=False) for j in pt_indices)
    print(f"Min C–Pt distance (Å) [mic=False]: {dmin:.3f}")

    # --- Height of C above the top Pt layer (z-direction sanity check) ---
    z_pt_max = atoms.positions[pt_indices, 2].max()
    dz = atoms.positions[c_idx, 2] - z_pt_max
    print(f"C height above top Pt layer (Å): {dz:.3f}")

if __name__ == "__main__":
    main()