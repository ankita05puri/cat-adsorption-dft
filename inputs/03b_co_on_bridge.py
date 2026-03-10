# 03b_co_on_bridge.py

from pathlib import Path

import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from gpaw import GPAW, PW, FermiDirac


BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "outputs"
OUT.mkdir(exist_ok=True)

VACUUM_Z = 15.0
FMAX = 0.05
STEPS = 200


def fix_bottom_n_layers(atoms, n_layers=2, tol=1e-3):
    z = atoms.positions[:, 2]
    z_sorted = np.sort(z)

    layers = []
    current = [z_sorted[0]]
    for val in z_sorted[1:]:
        if abs(val - current[-1]) < tol:
            current.append(val)
        else:
            layers.append(float(np.mean(current)))
            current = [val]
    layers.append(float(np.mean(current)))

    z_fix_max = layers[n_layers - 1] + tol
    mask = z <= z_fix_max
    atoms.set_constraint(FixAtoms(mask=mask))
    return mask


def build_calc(txt_path: Path):
    return GPAW(
        mode=PW(400),
        xc="PBE",
        kpts={"size": (4, 4, 1), "gamma": True},
        occupations=FermiDirac(0.1),
        symmetry={"point_group": False, "time_reversal": False},
        txt=str(txt_path),
    )


def find_top_layer_indices(atoms, tol=1e-3):
    z = atoms.positions[:, 2]
    z_top = z.max()
    top_idx = [i for i, zi in enumerate(z) if abs(zi - z_top) < tol]
    return top_idx


def choose_bridge_pair(atoms):
    """
    Pick a nearest-neighbor pair of Pt atoms in the top layer.
    """
    top_idx = find_top_layer_indices(atoms)
    top_positions = atoms.positions[top_idx]

    best_pair = None
    best_dist = 1e9

    for a_i, i in enumerate(top_idx):
        for a_j, j in enumerate(top_idx):
            if j <= i:
                continue
            d = atoms.get_distance(i, j, mic=True)
            # Pt-Pt nearest neighbor on surface should be the shortest top-layer pair
            if d < best_dist:
                best_dist = d
                best_pair = (i, j)

    return best_pair, best_dist


def build_co_on_bridge():
    slab_path = OUT / "pt111_relaxed.xyz"
    if not slab_path.exists():
        raise FileNotFoundError(
            "Need outputs/pt111_relaxed.xyz from 01_slab_relax.py"
        )

    slab = read(slab_path)
    slab.pbc = (True, True, False)
    slab.wrap()
    slab.center(axis=2, vacuum=VACUUM_Z)

    mask = fix_bottom_n_layers(slab, n_layers=2)
    print(f"[constraint] Fixed atoms: {int(mask.sum())} / {len(slab)}")

    (i, j), d_ptpt = choose_bridge_pair(slab)
    r_i = slab.positions[i]
    r_j = slab.positions[j]

    # Bridge midpoint
    midpoint = 0.5 * (r_i + r_j)

    # Put carbon above the bridge midpoint, oxygen above carbon
    z_top = max(atom.position[2] for atom in slab if atom.symbol == "Pt")
    c_height = z_top + 1.85
    co_bond = 1.15

    carbon_pos = midpoint.copy()
    carbon_pos[2] = c_height

    oxygen_pos = carbon_pos + np.array([0.0, 0.0, co_bond])

    slab += read_co(carbon_pos, oxygen_pos)

    slab.wrap()
    slab.center(axis=2, vacuum=VACUUM_Z)

    # Sanity checks
    c_idx = len(slab) - 2
    o_idx = len(slab) - 1

    min_c_pt = min(
        slab.get_distance(c_idx, k, mic=True)
        for k, atom in enumerate(slab[:-2])
        if atom.symbol == "Pt"
    )

    print(f"[site] bridge pair: Pt indices ({i}, {j}), Pt-Pt distance = {d_ptpt:.3f} Å")
    print(f"[sanity] C height above top Pt layer = {carbon_pos[2] - z_top:.3f} Å")
    print(f"[sanity] initial CO bond = {slab.get_distance(c_idx, o_idx):.3f} Å")
    print(f"[sanity] min distance (C ↔ Pt) = {min_c_pt:.3f} Å")

    return slab


def read_co(carbon_pos, oxygen_pos):
    from ase import Atoms
    return Atoms(
        "CO",
        positions=[carbon_pos, oxygen_pos],
        pbc=(False, False, False),
    )


def main():
    atoms = build_co_on_bridge()

    atoms.calc = build_calc(OUT / "co_bridge.txt")

    opt = BFGS(
        atoms,
        trajectory=str(OUT / "co_bridge.traj"),
        logfile=str(OUT / "co_bridge_opt.log"),
    )
    opt.run(fmax=FMAX, steps=STEPS)

    energy = atoms.get_potential_energy()

    c_idx = len(atoms) - 2
    o_idx = len(atoms) - 1
    d_co = atoms.get_distance(c_idx, o_idx)

    write(OUT / "co_bridge_relaxed.xyz", atoms)
    atoms.calc.write(str(OUT / "co_bridge_relaxed.gpw"))

    print(f"[done] E(CO*/bridge) = {energy:.6f} eV")
    print(f"[done] CO bond length = {d_co:.3f} Å")


if __name__ == "__main__":
    main()