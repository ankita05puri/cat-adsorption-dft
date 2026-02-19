from ase.build import fcc111
from ase.constraints import FixAtoms
from ase.io import write

from gpaw import GPAW, PW, FermiDirac


def main():
    slab = fcc111("Pt", size=(2, 2, 4), vacuum=15.0, a=3.92)

    # Fix bottom 2 layers
    z = slab.positions[:, 2]
    z_sorted = sorted(set(z))
    z_fix_max = z_sorted[1]
    mask = z <= z_fix_max
    slab.set_constraint(FixAtoms(mask=mask))

    write("outputs/pt111_clean_initial.xyz", slab)

    calc = GPAW(
        mode=PW(400),
        xc="PBE",
        kpts=(3, 3, 1),
        occupations=FermiDirac(0.1),
        txt="outputs/pt111_clean.txt",
    )

    slab.calc = calc
    E = slab.get_potential_energy()

    calc.write("outputs/pt111_clean.gpw")
    print("Clean slab energy (eV):", E)


if __name__ == "__main__":
    main()