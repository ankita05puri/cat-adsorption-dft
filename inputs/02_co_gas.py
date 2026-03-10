# 02_co_gas.py
import os
from ase import Atoms
from ase.io import write
from ase.optimize import BFGS
from gpaw import GPAW, PW, FermiDirac


def main():
    os.makedirs("outputs", exist_ok=True)

    # Initial CO 
    co = Atoms(
        "CO",
        positions=[(0.0, 0.0, 0.0),
                   (0.0, 0.0, 1.13)],
        cell=(15.0, 15.0, 15.0),
        pbc=False,
    )

    co.center(vacuum=6.0)      
    co.positions += 1e-3      

    calc = GPAW(
        mode=PW(400),
        xc="PBE",
        kpts=(1, 1, 1),                 
        occupations=FermiDirac(0.01),
        symmetry={"point_group": False, "time_reversal": False},
        txt="outputs/co_gas.txt",
    )
    co.calc = calc

    # Relaxation
    opt = BFGS(co, trajectory="outputs/co_gas.traj", logfile="outputs/co_gas_opt.log")
    opt.run(fmax=0.01, steps=50)

    E = co.get_potential_energy()
    d_co = co.get_distance(0, 1)

    calc.write("outputs/co_gas.gpw")
    write("outputs/co_gas.xyz", co)

    print("[done] CO(g) optimized bond length (Å):", round(d_co, 4))
    print("[done] CO(g) energy (eV):", E)


if __name__ == "__main__":
    main()