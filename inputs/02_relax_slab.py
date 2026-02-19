from ase.io import read
from ase.constraints import FixAtoms
from ase.optimize import BFGS

from gpaw import GPAW

def main():
    slab = read("outputs/pt111_clean_initial.xyz")

    # Fix bottom 2 layers (same logic you used before)
    z = slab.positions[:, 2]
    z_sorted = sorted(set(z))
    z_fix_max = z_sorted[1]
    mask = z <= z_fix_max
    slab.set_constraint(FixAtoms(mask=mask))

    calc = GPAW("outputs/pt111_clean.gpw")  # restart from converged SCF
    slab.calc = calc

    opt = BFGS(slab, trajectory="outputs/pt111_relax.traj", logfile="outputs/pt111_relax_opt.log")
    opt.run(fmax=0.05)  # eV/Å (good first target)

    calc.write("outputs/pt111_relaxed.gpw")
    slab.write("outputs/pt111_relaxed.xyz")

    E = slab.get_potential_energy()
    print("Relaxed slab energy (E0 from restart calc):", E)

if __name__ == "__main__":
    main()