import re
from pathlib import Path

def read_energy_from_txt(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")

    E = None
    with open(path, "r") as f:
        for line in f:
            if line.strip().startswith("Extrapolated:"):
                parts = line.split()
                E = float(parts[1])

    if E is None:
        raise RuntimeError(f"Could not find 'Extrapolated' energy in {path}")
    return E


def main():
    E_slab   = read_energy_from_txt("outputs/pt111_clean.txt")
    E_COgas  = read_energy_from_txt("outputs/co_gas.txt")
    E_slabCO = read_energy_from_txt("outputs/co_on_top.txt")

    E_ads = E_slabCO - E_slab - E_COgas

    print("E_slab (eV):     ", E_slab)
    print("E_CO(g) (eV):    ", E_COgas)
    print("E_slab+CO (eV):  ", E_slabCO)
    print("E_ads (eV):      ", E_ads)
    print("\nConvention: negative E_ads = exothermic adsorption (binds).")

    # Save summary
    with open("outputs/energies_summary.txt", "w") as f:
        f.write(f"E_slab   {E_slab}\n")
        f.write(f"E_COgas  {E_COgas}\n")
        f.write(f"E_slabCO {E_slabCO}\n")
        f.write(f"E_ads    {E_ads}\n")

if __name__ == "__main__":
    main()