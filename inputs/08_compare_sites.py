# inputs/08_compare_sites.py
# Compare CO adsorption energies: on-top vs bridge
# Convention: E_ads = E(slab+CO) - E(slab) - E(CO_gas)
# Negative E_ads => exothermic adsorption (strong binding)

def read_extrapolated_energy(txt_path: str) -> float:
    E = None
    with open(txt_path, "r") as f:
        for line in f:
            s = line.strip()
            if s.startswith("Extrapolated:"):
                # Format: "Extrapolated:  -112.280189"
                parts = s.split()
                E = float(parts[1])
    if E is None:
        raise RuntimeError(f"Could not find 'Extrapolated:' in {txt_path}")
    return E


def main():
    # Update these paths only if your filenames differ
    slab_txt = "outputs/pt111_clean.txt"
    co_gas_txt = "outputs/co_gas.txt"
    on_top_txt = "outputs/co_on_top_relax.txt"
    bridge_txt = "outputs/co_on_bridge.txt"

    E_slab = read_extrapolated_energy(slab_txt)
    E_co_gas = read_extrapolated_energy(co_gas_txt)
    E_top = read_extrapolated_energy(on_top_txt)
    E_bridge = read_extrapolated_energy(bridge_txt)

    Eads_top = E_top - E_slab - E_co_gas
    Eads_bridge = E_bridge - E_slab - E_co_gas

    print("Energies (eV) [from 'Extrapolated']:")
    print(f"  E_slab        = {E_slab: .6f}   ({slab_txt})")
    print(f"  E_CO(g)       = {E_co_gas: .6f}   ({co_gas_txt})")
    print(f"  E_slab+CO top  = {E_top: .6f}   ({on_top_txt})")
    print(f"  E_slab+CO br   = {E_bridge: .6f}   ({bridge_txt})")
    print()
    print("Adsorption energies (eV):  E_ads = E(slab+CO) - E(slab) - E(CO_g)")
    print(f"  E_ads(top)    = {Eads_top: .6f}")
    print(f"  E_ads(bridge) = {Eads_bridge: .6f}")
    print()

    if Eads_top < Eads_bridge:
        stronger = "on-top"
    else:
        stronger = "bridge"

    print(f"Stronger binding (more negative E_ads): {stronger}")

    # Optional: quick sanity check for weird positives
    if Eads_top > 0 or Eads_bridge > 0:
        print("\nNote: Positive E_ads means endothermic adsorption.")
        print("If this is unexpected, ensure BOTH site calculations are relaxed and consistent settings are used.")


if __name__ == "__main__":
    main()