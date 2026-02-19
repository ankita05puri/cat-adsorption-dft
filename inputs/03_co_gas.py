from ase import Atoms
from gpaw import GPAW, PW, FermiDirac

def main():
    # CO in a big box (gas-phase)
    co = Atoms("CO",
               positions=[(0.0, 0.0, 0.0),
                          (0.0, 0.0, 1.13)],
               cell=(15.0, 15.0, 15.0),
               pbc=(False, False, False))

    calc = GPAW(
        mode=PW(400),
        xc="PBE",
        kpts=(1, 1, 1),
        occupations=FermiDirac(0.1),
        symmetry="off",
        txt="outputs/co_gas.txt",
    )
    co.calc = calc

    E = co.get_potential_energy()
    calc.write("outputs/co_gas.gpw")
    print("CO(g) energy (eV):", E)

if __name__ == "__main__":
    main()