from ase.io import read
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.build import molecule, add_adsorbate

from gpaw import GPAW, PW, FermiDirac


def fix_bottom_two_layers(slab):
    z = slab.positions[:, 2]
    z_unique = sorted(set([round(v, 6) for v in z]))
    z_fix_max = z_unique[1]  # bottom 2 layers
    mask = z <= z_fix_max
    slab.set_constraint(FixAtoms(mask=mask))


def top_layer_indices(slab, tol=1e-3):
    z = slab.positions[:, 2]
    top_z = z.max()
    return [i for i, zi in enumerate(z) if abs(zi - top_z) < tol], top_z


def mic_midpoint(slab, i1, i2):
    """
    MIC-safe midpoint between atoms i1 and i2:
    r_mid = r1 + 0.5 * (r2 - r1)_MIC
    """
    r1 = slab.positions[i1].copy()
    # vector from i1 -> i2 with minimum-image convention
    v12 = slab.get_distance(i1, i2, mic=True, vector=True)
    rmid = r1 + 0.5 * v12
    return rmid


def min_c_pt_distance(atoms, mic=True):
    c_idx = [i for i, a in enumerate(atoms) if a.symbol == "C"][-1]
    pt_idx = [i for i, a in enumerate(atoms) if a.symbol == "Pt"]
    dmin = min(atoms.get_distance(c_idx, j, mic=mic) for j in pt_idx)
    return dmin


def c_height_above_top_layer(atoms, top_z):
    c_idx = [i for i, a in enumerate(atoms) if a.symbol == "C"][-1]
    return atoms.positions[c_idx, 2] - top_z


def main():
    slab = read("outputs/pt111_clean_initial.xyz")

    # For PW mode, keep full PBC (vacuum already in the z cell length)
    slab.pbc = (True, True, True)
    slab.wrap()

    fix_bottom_two_layers(slab)

    top_idx, top_z = top_layer_indices(slab)

    # choose the closest pair among top-layer Pt atoms
    best_pair = None
    best_d = 1e9
    for a_i in range(len(top_idx)):
        for b_i in range(a_i + 1, len(top_idx)):
            a = top_idx[a_i]
            b = top_idx[b_i]
            d = slab.get_distance(a, b, mic=True)
            if d < best_d:
                best_d = d
                best_pair = (a, b)

    i1, i2 = best_pair
    rmid = mic_midpoint(slab, i1, i2)

    # Add CO
    co = molecule("CO")
    add_adsorbate(slab, co, height=2.4, position=(rmid[0], rmid[1]))

    # Center slab in z, but DO NOT wrap in z after placing adsorbate
    slab.center(axis=2)

    # Recompute top_z AFTER centering (important!)
    z_pt = slab.positions[[i for i, a in enumerate(slab) if a.symbol == "Pt"], 2]
    top_z_new = z_pt.max()

    dmin = min_c_pt_distance(slab, mic=True)
    cheight = c_height_above_top_layer(slab, top_z_new)

    print(f"C height above top Pt layer (Å): {cheight:.3f}")
    print(f"Min C–Pt distance (Å): {dmin:.3f}")

    # If geometry is nonsense, stop early.
    if cheight < 1.2 or dmin < 1.6:
        raise RuntimeError(
            "Geometry looks too close to the surface. "
            "Increase height (e.g., 2.2 Å) or check bridge site selection."
        )

    calc = GPAW(
        mode=PW(400),
        xc="PBE",
        kpts=(4, 4, 1),
        occupations=FermiDirac(0.1),
        symmetry="off",
        txt="outputs/co_on_bridge.txt",
    )
    slab.calc = calc

    # Recommended: relax for meaningful on-top vs bridge comparison
    opt = BFGS(slab, trajectory="outputs/co_on_bridge.traj",
               logfile="outputs/co_on_bridge_opt.log")
    opt.run(fmax=0.05)

    E = slab.get_potential_energy()
    print("Relaxed CO* bridge energy (eV):", E)

    calc.write("outputs/co_on_bridge.gpw")
    slab.write("outputs/co_on_bridge.xyz")


if __name__ == "__main__":
    main()