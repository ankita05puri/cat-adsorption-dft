# 04b_compare_o_sites.py
#
# Compare total energies and adsorption energies for O on fcc vs hcp hollow sites
# using the same clean slab and O2(g) reference.
#
# Expected files:
#   outputs/pt111_relaxed.gpw
#   outputs/o2_gas.gpw
#   outputs/o_fcc_relaxed.gpw   (or outputs/o_fcc.txt fallback)
#   outputs/o_hcp_relaxed.gpw   (or outputs/o_hcp.txt fallback)

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
    e_o2_gas = get_energy("o2_gas.gpw", "o2_gas.txt")

    # Adsorbed states
    e_o_fcc = get_energy("o_fcc_relaxed.gpw", "o_fcc.txt")
    e_o_hcp = get_energy("o_hcp_relaxed.gpw", "o_hcp.txt")

    # Adsorption energies for atomic O using 1/2 E(O2)
    eads_fcc = e_o_fcc - e_slab - 0.5 * e_o2_gas
    eads_hcp = e_o_hcp - e_slab - 0.5 * e_o2_gas

    # Stability difference
    # negative means fcc is lower in energy than hcp
    delta_total_fcc_minus_hcp = e_o_fcc - e_o_hcp
    delta_ads_fcc_minus_hcp = eads_fcc - eads_hcp

    result = {
        "E_slab_eV": e_slab,
        "E_o2_gas_eV": e_o2_gas,
        "E_o_fcc_total_eV": e_o_fcc,
        "E_o_hcp_total_eV": e_o_hcp,
        "Eads_o_fcc_eV": eads_fcc,
        "Eads_o_hcp_eV": eads_hcp,
        "delta_total_fcc_minus_hcp_eV": delta_total_fcc_minus_hcp,
        "delta_ads_fcc_minus_hcp_eV": delta_ads_fcc_minus_hcp,
        "formula": "Eads(O*) = E(slab+O) - E(slab) - 0.5*E(O2_g)",
        "interpretation": {
            "more_negative_Eads_is_more_stable": True,
            "negative_delta_fcc_minus_hcp_means": "fcc hollow is more stable than hcp hollow",
            "positive_delta_fcc_minus_hcp_means": "hcp hollow is more stable than fcc hollow",
        },
    }

    out_json = OUT / "o_site_comparison.json"
    out_txt = OUT / "o_site_comparison.txt"

    out_json.write_text(json.dumps(result, indent=2))

    lines = [
        "O adsorption site comparison on Pt(111)",
        "",
        f"E(slab)              = {e_slab:.6f} eV",
        f"E(O2 gas)            = {e_o2_gas:.6f} eV",
        "",
        f"E(slab+O_fcc)        = {e_o_fcc:.6f} eV",
        f"E(slab+O_hcp)        = {e_o_hcp:.6f} eV",
        "",
        f"Eads(O_fcc)          = {eads_fcc:.6f} eV",
        f"Eads(O_hcp)          = {eads_hcp:.6f} eV",
        "",
        f"ΔE_total(fcc-hcp)    = {delta_total_fcc_minus_hcp:.6f} eV",
        f"ΔE_ads(fcc-hcp)      = {delta_ads_fcc_minus_hcp:.6f} eV",
        "",
    ]

    if delta_ads_fcc_minus_hcp < 0:
        lines.append("Conclusion: O on FCC hollow is more stable than O on HCP hollow.")
    elif delta_ads_fcc_minus_hcp > 0:
        lines.append("Conclusion: O on HCP hollow is more stable than O on FCC hollow.")
    else:
        lines.append("Conclusion: FCC and HCP are equal within this comparison.")

    out_txt.write_text("\n".join(lines))

    print("\n".join(lines))
    print(f"\n[done] wrote: {out_txt}")
    print(f"[done] wrote: {out_json}")

    # -----------------------------
    # Plot adsorption energy comparison
    # -----------------------------
    labels = ["O@FCC", "O@HCP"]
    energies = [eads_fcc, eads_hcp]

    plt.figure(figsize=(5,4))

    bars = plt.bar(
        labels,
        energies,
        width=0.55
    )

    plt.ylabel("Adsorption Energy (eV)")
    plt.title("O Adsorption on Pt(111)")

    # zero reference line
    plt.axhline(0, linewidth=1)

    # light grid improves readability
    plt.grid(axis="y", linestyle="--", alpha=0.3)

    plt.ylim(-1.25, 0.05)

    # annotate bars cleanly
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

    fig_path = OUT / "o_adsorption_site_comparison.png"
    plt.savefig(fig_path, dpi=300)
    plt.close()
    
    print(f"[done] wrote: {fig_path}")


if __name__ == "__main__":
    main()