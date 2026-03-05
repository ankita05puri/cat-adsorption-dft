# 05_o2_gas.py

from ase import Atoms
from ase.optimize import BFGS
from ase.io import write
from gpaw import GPAW, PW, FermiDirac


def main():
    o2 = Atoms("O2", positions=[(0, 0, 0), (0, 0, 1.23)], pbc=False)
    o2.center(vacuum=8.0)

    calc = GPAW(
        mode=PW(400),
        xc="PBE",
        kpts=(1, 1, 1),
        occupations=FermiDirac(0.01),
        spinpol=True,
        hund=True,
        txt="outputs/o2_gas.txt",
    )
    o2.calc = calc

    opt = BFGS(o2, trajectory="outputs/o2_gas.traj", logfile="outputs/o2_gas_opt.log")
    opt.run(fmax=0.01, steps=80)

    E = o2.get_potential_energy()

    print("[done] O2(g) bond length (Å):", o2.get_distance(0, 1))
    print("[done] O2(g) energy (eV):", E)
    print("[done] magnetic moment (μB):", o2.get_magnetic_moment())

    calc.write("outputs/o2_gas.gpw")
    write("outputs/o2_gas.xyz", o2)


if __name__ == "__main__":
    main()
