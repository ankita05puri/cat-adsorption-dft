# 06_extract_energies.py
<<<<<<< HEAD

from __future__ import annotations
=======
>>>>>>> aa75805 (Add adsorption site comparisons and update DFT workflow scripts)

from pathlib import Path
import json
from gpaw import GPAW

OUT = Path(__file__).resolve().parents[1] / "outputs"

def energy(gpw_path: Path) -> float:
    calc = GPAW(str(gpw_path))
    return float(calc.get_potential_energy())

def main():
    E_slab = energy(OUT / "pt111_relaxed.gpw")
    E_co_gas = energy(OUT / "co_gas.gpw")
    E_o2_gas = energy(OUT / "o2_gas.gpw")


    coads_file = OUT / "co_o_init_relaxed.gpw" 
    E_coads_tot = energy(coads_file)

    E_coads = E_coads_tot - E_slab - E_co_gas - 0.5 * E_o2_gas

    results = {
        "E_slab": E_slab,
        "E_co_gas": E_co_gas,
        "E_o2_gas": E_o2_gas,
        "E_slab_co_o": E_coads_tot,
        "E_coads": E_coads,
        "formula": "E(slab+CO+O) - E(slab) - E(CO_g) - 0.5*E(O2_g)",
        "units": "eV",
    }

    print(json.dumps(results, indent=2))

    out_json = OUT / "coadsorption_energy.json"
    out_json.write_text(json.dumps(results, indent=2))
    print(f"\n[saved] {out_json}")

if __name__ == "__main__":
    main()
