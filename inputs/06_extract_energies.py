# inputs/06_extract_energies.py
from __future__ import annotations

from pathlib import Path
import json
import csv
from datetime import datetime

from gpaw import GPAW


BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "outputs"
OUT.mkdir(exist_ok=True)

FILES = {
    "slab": OUT / "pt111_relaxed.gpw",
    "co_gas": OUT / "co_gas.gpw",
    "o2_gas": OUT / "o2_gas.gpw",
    "co_top": OUT / "co_top_relaxed.gpw",
    "o_fcc": OUT / "o_fcc_relaxed.gpw",
    "o_hcp": OUT / "o_hcp_relaxed.gpw",
}

JSON_OUT = OUT / "adsorption_energies.json"
CSV_OUT = OUT / "adsorption_energies.csv"


def energy_from_gpw(path: Path) -> float:
    """Load a GPAW .gpw file and return the potential energy in eV."""
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    calc = GPAW(str(path))
    return float(calc.get_potential_energy())


def main() -> None:
    # --- Load energies ---
    E_slab = energy_from_gpw(FILES["slab"])
    E_CO = energy_from_gpw(FILES["co_gas"])
    E_O2 = energy_from_gpw(FILES["o2_gas"])

    E_slab_CO_top = energy_from_gpw(FILES["co_top"])
    E_slab_O_fcc = energy_from_gpw(FILES["o_fcc"])
    E_slab_O_hcp = energy_from_gpw(FILES["o_hcp"])

    # --- Adsorption energies ---
    # E_ads(CO) = E(slab+CO) - E(slab) - E(CO)
    E_ads_CO_top = E_slab_CO_top - E_slab - E_CO

    # E_ads(O) = E(slab+O) - E(slab) - 1/2 E(O2)
    half_E_O2 = 0.5 * E_O2
    E_ads_O_fcc = E_slab_O_fcc - E_slab - half_E_O2
    E_ads_O_hcp = E_slab_O_hcp - E_slab - half_E_O2

    # Site preference
    dE_hcp_minus_fcc = E_ads_O_hcp - E_ads_O_fcc

    results = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "files": {k: str(v) for k, v in FILES.items()},
        "energies_eV": {
            "E_slab": E_slab,
            "E_CO_gas": E_CO,
            "E_O2_gas": E_O2,
            "E_slab_CO_top": E_slab_CO_top,
            "E_slab_O_fcc": E_slab_O_fcc,
            "E_slab_O_hcp": E_slab_O_hcp,
        },
        "adsorption_energies_eV": {
            "E_ads_CO_top": E_ads_CO_top,
            "E_ads_O_fcc": E_ads_O_fcc,
            "E_ads_O_hcp": E_ads_O_hcp,
            "delta_E_hcp_minus_fcc": dE_hcp_minus_fcc,
        },
        "notes": {
            "E_ads_CO_definition": "E(slab+CO) - E(slab) - E(CO)",
            "E_ads_O_definition": "E(slab+O) - E(slab) - 0.5*E(O2)",
            "delta_definition": "E_ads(O_hcp) - E_ads(O_fcc); positive => fcc more stable",
        },
    }

    # --- Save JSON ---
    with open(JSON_OUT, "w") as f:
        json.dump(results, f, indent=2)

    # --- Save CSV (flat, easy for spreadsheet) ---
    rows = [
        ("E_slab", E_slab),
        ("E_CO_gas", E_CO),
        ("E_O2_gas", E_O2),
        ("E_slab_CO_top", E_slab_CO_top),
        ("E_slab_O_fcc", E_slab_O_fcc),
        ("E_slab_O_hcp", E_slab_O_hcp),
        ("E_ads_CO_top", E_ads_CO_top),
        ("E_ads_O_fcc", E_ads_O_fcc),
        ("E_ads_O_hcp", E_ads_O_hcp),
        ("delta_E_hcp_minus_fcc", dE_hcp_minus_fcc),
    ]
    with open(CSV_OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "value_eV"])
        w.writerows(rows)

    # --- Print summary ---
    print("\n[OK] Energies extracted from .gpw files:")
    print(f"  E_slab        = {E_slab:.12f} eV")
    print(f"  E_CO(g)       = {E_CO:.12f} eV")
    print(f"  E_O2(g)       = {E_O2:.12f} eV")
    print(f"  E_slab+CO_top = {E_slab_CO_top:.12f} eV")
    print(f"  E_slab+O_fcc  = {E_slab_O_fcc:.12f} eV")
    print(f"  E_slab+O_hcp  = {E_slab_O_hcp:.12f} eV")

    print("\n[OK] Adsorption energies:")
    print(f"  E_ads(CO_top) = {E_ads_CO_top:.6f} eV")
    print(f"  E_ads(O_fcc)  = {E_ads_O_fcc:.6f} eV")
    print(f"  E_ads(O_hcp)  = {E_ads_O_hcp:.6f} eV")
    print(f"  ΔE (hcp-fcc)  = {dE_hcp_minus_fcc:.6f} eV  (positive => fcc preferred)")

    print("\n[SAVED]")
    print(f"  {JSON_OUT}")
    print(f"  {CSV_OUT}\n")


if __name__ == "__main__":
    main()