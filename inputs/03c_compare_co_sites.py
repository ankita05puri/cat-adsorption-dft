# 03c_compare_co_sites.py
#
# Compare total energies and adsorption energies for CO on top vs bridge sites
# using the same clean slab and CO(g) reference, and generate a bar chart.

from pathlib import Path
import json

import matplotlib.pyplot as plt
from gpaw import GPAW


BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "outputs"
OUT.mkdir(exist_ok=True)


def read_energy_from_gpw(path: Path) -> float:
    calc = GPAW(str(path))
    return calc.get_potential_energy()


def read_last_extrapolated(txt_path: Path) -> float:
    """
    Fallback parser for GPAW txt output if gpw file is missing.
    """
    lines = txt_path.read_text().splitlines()
    vals = []

    for line in lines:
        if "Extrapolated:" in line:
            try:
                vals.append(float(line.split()[-1]))
            except Exception:
                pass

    if not vals:
        raise ValueError(f"No 'Extrapolated:' energies found in {txt_path}")
    return vals[-1]


def get_energy(preferred_gpw: str, fallback_txt: str) -> float:
    gpw_path = OUT / preferred_gpw
    txt_path = OUT / fallback_txt

    if gpw_path.exists():
        return read_energy_from_gpw(gpw_path)
    if txt_path.exists():
        return read_last_extrapolated(txt_path)

    raise FileNotFoundError(
        f"Could not find either {gpw_path.name} or {txt_path.name}"
    )


def main():
    # References
    e_slab = get_energy("pt111_relaxed.gpw", "pt111_relaxed.txt")
    e_co_gas = get_energy("co_gas.gpw", "co_gas.txt")

    # Adsorbed states
    e_co_top = get_energy("co_top_relaxed.gpw", "co_top.txt")
    e_co_bridge = get_energy("co_bridge_relaxed.gpw", "co_bridge.txt")

    # Adsorption energies
    eads_top = e_co_top - e_slab - e_co_gas
    eads_bridge = e_co_bridge - e_slab - e_co_gas

    # Stability differences
    # Negative means top is lower in energy than bridge
    delta_total_top_minus_bridge = e_co_top - e_co_bridge
    delta_ads_top_minus_bridge = eads_top - eads_bridge

    result = {
        "E_slab_eV": e_slab,
        "E_co_gas_eV": e_co_gas,
        "E_co_top_total_eV": e_co_top,
        "E_co_bridge_total_eV": e_co_bridge,
        "Eads_co_top_eV": eads_top,
        "Eads_co_bridge_eV": eads_bridge,
        "delta_total_top_minus_bridge_eV": delta_total_top_minus_bridge,
        "delta_ads_top_minus_bridge_eV": delta_ads_top_minus_bridge,
        "interpretation": {
            "more_negative_Eads_is_more_stable": True,
            "negative_delta_top_minus_bridge_means": "top site is more stable than bridge",
            "positive_delta_top_minus_bridge_means": "bridge site is more stable than top",
        },
    }

    out_json = OUT / "co_site_comparison.json"
    out_txt = OUT / "co_site_comparison.txt"

    out_json.write_text(json.dumps(result, indent=2))

    lines = [
        "CO adsorption site comparison on Pt(111)",
        "",
        f"E(slab)              = {e_slab:.6f} eV",
        f"E(CO gas)            = {e_co_gas:.6f} eV",
        "",
        f"E(slab+CO_top)       = {e_co_top:.6f} eV",
        f"E(slab+CO_bridge)    = {e_co_bridge:.6f} eV",
        "",
        f"Eads(CO_top)         = {eads_top:.6f} eV",
        f"Eads(CO_bridge)      = {eads_bridge:.6f} eV",
        "",
        f"ΔE_total(top-bridge) = {delta_total_top_minus_bridge:.6f} eV",
        f"ΔE_ads(top-bridge)   = {delta_ads_top_minus_bridge:.6f} eV",
        "",
    ]

    if delta_ads_top_minus_bridge < 0:
        lines.append("Conclusion: CO on TOP is more stable than CO on BRIDGE.")
    elif delta_ads_top_minus_bridge > 0:
        lines.append("Conclusion: CO on BRIDGE is more stable than CO on TOP.")
    else:
        lines.append("Conclusion: TOP and BRIDGE are equal within this comparison.")

    out_txt.write_text("\n".join(lines))

    print("\n".join(lines))
    print(f"\n[done] wrote: {out_txt}")
    print(f"[done] wrote: {out_json}")

    # -----------------------------
    # Plot CO adsorption comparison
    # -----------------------------
    labels = ["CO (top)", "CO (bridge)"]
    energies = [eads_top, eads_bridge]

    plt.figure(figsize=(5,4))

    bars = plt.bar(
        labels,
        energies,
        width=0.55
    )

    plt.ylabel("Adsorption Energy (eV)")
    plt.title("CO Adsorption on Pt(111)")

    plt.axhline(0, linewidth=1)
    plt.grid(axis="y", linestyle="--", alpha=0.3)

    plt.ylim(min(energies)-0.2, 0.05)

    for bar, val in zip(bars, energies):
        plt.text(
            bar.get_x() + bar.get_width()/2,
            val - 0.06,
            f"{val:.2f}",
            ha="center",
            va="top",
            fontsize=10
        )

    plt.tight_layout()

    fig_path = OUT / "co_adsorption_site_comparison.png"
    plt.savefig(fig_path, dpi=300)
    plt.close()
    
    print(f"[done] wrote: {fig_path}")


if __name__ == "__main__":
    main()