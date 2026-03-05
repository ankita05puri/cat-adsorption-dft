"""
03_co_on_top.py
Pillar 3A: CO adsorption geometry on Pt(111) on-top site (relax)

Inputs:
  outputs/pt111_relaxed.xyz   (from slab relaxation)

Outputs:
  outputs/co_top.txt          (GPAW SCF log)
  outputs/co_top_opt.log      (BFGS log)
  outputs/co_top.traj         (trajectory)
  outputs/co_top_relaxed.xyz  (final geometry)
  outputs/co_top_relaxed.gpw  (restart)
"""

from pathlib import Path
import numpy as np

from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import BFGS

from gpaw import GPAW, PW, FermiDirac, Mixer


# ---------- paths ----------
BASE = Path(__file__).resolve().parents[1]   # project root
OUT = BASE / "outputs"
OUT.mkdir(exist_ok=True)


# ---------- helpers ----------
def fix_bottom_n_layers(atoms, n_layers=2, tol=1e-3):
    """Fix bottom n_layers by clustering z-coordinates (robust to float noise)."""
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

    if len(layers) < n_layers:
        raise RuntimeError(f"Found only {len(layers)} layers; increase tol or check slab.")

    z_fix_max = layers[n_layers - 1] + tol
    mask = z <= z_fix_max
    atoms.set_constraint(FixAtoms(mask=mask))
    return mask


def build_calc(txt_path):
    """Metallic settings consistent with slab relaxation."""
    return GPAW(
        mode=PW(400),
        xc="PBE",
        kpts={"size": (4, 4, 1), "gamma": True},
        occupations=FermiDirac(0.1),
        mixer=Mixer(0.05, 5, weight=100.0),
        symmetry={"point_group": False, "time_reversal": False},
        convergence={"energy": 1e-5, "density": 1e-4, "maximum iterations": 200},
        txt=str(txt_path),
    )


def main():
    slab = read(OUT / "pt111_relaxed.xyz")
    slab.pbc = (True, True, False)
    slab.wrap()
    slab.center(axis=2, vacuum=15.0)

    mask = fix_bottom_n_layers(slab, n_layers=2, tol=1e-3)
    print(f"[constraint] Fixed atoms: {int(mask.sum())} / {len(slab)}")

    # pick a top-layer Pt atom (highest z)
    z = slab.positions[:, 2]
    pt_idx = [i for i, a in enumerate(slab) if a.symbol == "Pt"]
    top_z = max(z[i] for i in pt_idx)
    top_sites = [i for i in pt_idx if abs(z[i] - top_z) < 1e-3]
    i_top = top_sites[0]
    x_top, y_top, z_top = slab.positions[i_top]

    # place CO (C down, O up)
    h_C = 1.85   # initial C height above top Pt (Å)
    d_CO = 1.15  # initial CO bond length (Å)

    slab.append("C")
    slab.positions[-1] = [x_top, y_top, z_top + h_C]

    slab.append("O")
    slab.positions[-1] = [x_top, y_top, z_top + h_C + d_CO]

    slab.center(axis=2)  # keep vacuum symmetric

    slab.calc = build_calc(OUT / "co_top.txt")

    opt = BFGS(
        slab,
        trajectory=str(OUT / "co_top.traj"),
        logfile=str(OUT / "co_top_opt.log"),
    )
    opt.run(fmax=0.05, steps=200)

    E = slab.get_potential_energy()
    print("[done] E(slab+CO on-top) (eV):", E)

    slab.calc.write(str(OUT / "co_top_relaxed.gpw"))
    write(OUT / "co_top_relaxed.xyz", slab)

    # quick sanity: C height above top Pt layer
    z_pt_top = max(slab.positions[i, 2] for i in pt_idx)
    c_idx = [i for i, a in enumerate(slab) if a.symbol == "C"][-1]
    print("[sanity] C height above top Pt (Å):", slab.positions[c_idx, 2] - z_pt_top)


if __name__ == "__main__":
    main()