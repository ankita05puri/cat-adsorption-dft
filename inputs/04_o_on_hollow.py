# 04_o_on_hollow.py

from pathlib import Path
import numpy as np

from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from gpaw import GPAW, PW, FermiDirac, Mixer

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "outputs"
OUT.mkdir(exist_ok=True)

VACUUM = 15.0

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

def build_calc(txt_path):
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

def top_and_second_layer_indices(slab, tol=1e-3):
    z = slab.positions[:, 2]
    pt_idx = [i for i, a in enumerate(slab) if a.symbol == "Pt"]
    z_pt = np.array([z[i] for i in pt_idx])

    z_top = z_pt.max()
    top = [pt_idx[i] for i in np.where(np.abs(z_pt - z_top) < tol)[0]]

    # second layer: next distinct z below top
    z_unique = np.unique(np.round(z_pt, 3))
    z_unique = np.sort(z_unique)
    z_second = z_unique[-2]
    second = [pt_idx[i] for i in np.where(np.abs(z_pt - z_second) < tol)[0]]

    return top, second

def nearest_three_top_atoms(slab, top_idx):
    """Pick 3 top atoms that form a local triangle (closest to first top atom)."""
    r = slab.positions[top_idx][:, :2]
    r0 = r[0]
    d = np.linalg.norm(r - r0, axis=1)
    order = np.argsort(d)
    # include itself + two nearest neighbors
    pick = [top_idx[order[0]], top_idx[order[1]], top_idx[order[2]]]
    return pick

def hollow_xy(slab, top_idx):
    tri = nearest_three_top_atoms(slab, top_idx)
    xy = slab.positions[tri, :2].mean(axis=0)
    return xy[0], xy[1]

def classify_hollow_site(slab, x, y, second_idx, cutoff=1.2):
    """
    If there's a 2nd-layer Pt nearly directly underneath (in xy), it's hcp.
    If not, it's fcc.
    """
    r2 = slab.positions[second_idx][:, :2]
    d = np.linalg.norm(r2 - np.array([x, y]), axis=1)
    if d.min() < cutoff:
        return "hcp"
    return "fcc"

def place_O_on_site(slab, site, height=1.2):
    top_idx, second_idx = top_and_second_layer_indices(slab)

    # get a hollow center in the cell
    x, y = hollow_xy(slab, top_idx)

    a1 = slab.cell[0, :2]
    a2 = slab.cell[1, :2]

    candidates = []
    shifts = [(0,0), (1/3,1/3), (2/3,1/3), (1/3,2/3), (2/3,2/3)]
    for s1, s2 in shifts:
        xx, yy = np.array([x, y]) + s1 * a1 + s2 * a2
        typ = classify_hollow_site(slab, xx, yy, second_idx)
        candidates.append((typ, xx, yy))

    # pick first candidate that matches requested site
    for typ, xx, yy in candidates:
        if typ == site:
            z_top = slab.positions[:, 2].max()
            slab.append("O")
            slab.positions[-1] = [xx, yy, z_top + height]
            return

    # fallback: place at first hollow (should rarely happen)
    z_top = slab.positions[:, 2].max()
    slab.append("O")
    slab.positions[-1] = [x, y, z_top + height]

def run_site(site):
    slab = read(OUT / "pt111_relaxed.xyz")
    slab.pbc = (True, True, False)
    slab.wrap()
    slab.center(axis=2, vacuum=VACUUM)

    mask = fix_bottom_n_layers(slab, n_layers=2)
    print(f"[{site}] [constraint] Fixed atoms: {int(mask.sum())} / {len(slab)}")

    place_O_on_site(slab, site=site, height=1.2)
    slab.center(axis=2, vacuum=VACUUM)

    slab.calc = build_calc(OUT / f"o_{site}.txt")

    opt = BFGS(slab, trajectory=str(OUT / f"o_{site}.traj"), logfile=str(OUT / f"o_{site}_opt.log"))
    opt.run(fmax=0.05, steps=200)

    E = slab.get_potential_energy()
    print(f"[{site}] [done] E(slab+O) (eV):", E)

    slab.calc.write(str(OUT / f"o_{site}_relaxed.gpw"))
    write(OUT / f"o_{site}_relaxed.xyz", slab)

def main():
    for site in ["fcc", "hcp"]:
        run_site(site)

if __name__ == "__main__":
    main()
