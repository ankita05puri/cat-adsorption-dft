from ase.io import read
from ase.constraints import FixAtoms
from ase.optimize import BFGS

from gpaw import GPAW, PW, FermiDirac


def fix_bottom_two_layers(slab):
    z = slab.positions[:, 2]
    z_unique = sorted(set([round(v, 6) for v in z]))
    z_fix_max = z_unique[1]  # bottom 2 layers
    slab.set_constraint(FixAtoms(mask=z <= z_fix_max))


def main():
    slab = read("outputs/pt111_clean_initial.xyz")

    slab.pbc = (True, True, True)
    slab.center(axis=2)
    fix_bottom_two_layers(slab)

    # Pick a top-layer Pt atom
    z = slab.positions[:, 2]
    pt_idx = [i for i, a in enumerate(slab) if a.symbol == "Pt"]
    top_z = max(z[i] for i in pt_idx)
    top_indices = [i for i in pt_idx if abs(z[i] - top_z) < 1e-3]
    i_top = top_indices[0]
    x_top, y_top, z_top = slab.positions[i_top]

    # Build CO manually (vertical)
    h_C = 1.85   # initial C height above top Pt (Å)
    d_CO = 1.15  # CO bond length (Å)

    slab.append("C")
    slab.positions[-1] = [x_top, y_top, z_top + h_C]

    slab.append("O")
    slab.positions[-1] = [x_top, y_top, z_top + h_C + d_CO]

    slab.center(axis=2)

    # Calculator (symmetry off)
    calc = GPAW(
        mode=PW(400),
        xc="PBE",
        kpts=(4, 4, 1),
        occupations=FermiDirac(0.1),
        symmetry="off",
        txt="outputs/co_on_top_relax.txt",
    )
    slab.calc = calc

    # Relax geometry
    opt = BFGS(
        slab,
        trajectory="outputs/co_on_top.traj",
        logfile="outputs/co_on_top_opt.log",
    )
    opt.run(fmax=0.05)

    # Final energy + save
    E = slab.get_potential_energy()
    print("Relaxed CO* on-top energy (eV):", E)

    calc.write("outputs/co_on_top_relaxed.gpw")
    slab.write("outputs/co_on_top_relaxed.xyz")

    # Sanity: C height above top Pt layer
    z_pt_top = max(slab.positions[i, 2] for i in pt_idx)
    c_idx = [i for i, a in enumerate(slab) if a.symbol == "C"][-1]
    print("C height above top Pt layer (Å):", slab.positions[c_idx, 2] - z_pt_top)


if __name__ == "__main__":
    main()