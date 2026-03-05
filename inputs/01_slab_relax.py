# inputs/01_slab_relax.py

import os
import numpy as np

from ase.build import fcc111
from ase.io import write, read
from ase.constraints import FixAtoms
from ase.optimize import BFGS

from gpaw import GPAW, PW, FermiDirac, Mixer


VACUUM = 15.0
SLAB_SIZE = (2, 2, 4)
A_Pt = 3.92
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

    if len(layers) < n_layers:
        raise RuntimeError(
            f"Could not identify {n_layers} layers (found {len(layers)})."
        )

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
        txt=txt_path,
    )


def main():
    os.makedirs("outputs", exist_ok=True)

    traj_path = "outputs/pt111_relaxed.traj"
    optlog_path = "outputs/pt111_relaxed_opt.log"
    gpawtxt_path = "outputs/pt111_relaxed.txt"

    # Restart from final relaxed structure if exists
    if os.path.exists(traj_path) and os.path.getsize(traj_path) > 0:
        slab = read(traj_path)
        print("[restart] Loaded existing relaxed slab:", traj_path)
    else:
        slab = fcc111("Pt", size=SLAB_SIZE, vacuum=VACUUM, a=A_Pt)
        print(f"[start] Built Pt(111): size={SLAB_SIZE}, vacuum={VACUUM} Å")

    slab.pbc = (True, True, False)
    slab.wrap()
    slab.center(axis=2, vacuum=VACUUM)

    print("[check] pbc:", slab.pbc)
    print("[check] cell lengths:", slab.cell.lengths())

    mask = fix_bottom_n_layers(slab, n_layers=2)
    print(f"[constraint] Fixed atoms: {int(mask.sum())} / {len(slab)}")

    slab.calc = build_calc(gpawtxt_path)

    opt = BFGS(slab, trajectory=traj_path, logfile=optlog_path)
    opt.run(fmax=FMAX, steps=STEPS)

    slab.calc.write("outputs/pt111_relaxed.gpw")
    slab.write("outputs/pt111_relaxed.xyz")

    E = slab.get_potential_energy()
    print("[done] Relaxed slab energy (eV):", E)


if __name__ == "__main__":
    main()