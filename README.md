CO Adsorption on Pt(111) using DFT (GPAW + ASE)

Overview

This project investigates CO adsorption on a Pt(111) surface using plane-wave DFT with GPAW and ASE.

We compute:
	•	Clean slab energy
	•	CO(g) energy
	•	CO adsorbed on-top
	•	CO adsorbed at bridge site
	•	Adsorption energies for both geometries

The goal is to compare site preference and quantify binding strength.

⸻

Scientific Question

Which site binds CO more strongly on Pt(111)?

We compare:

E_{ads} = E_{slab+CO} - E_{slab} - E_{CO(g)}

More negative adsorption energy → stronger binding.

⸻

Computational Details
	•	Code: GPAW (plane-wave mode)
	•	XC functional: PBE
	•	Cutoff: 400 eV
	•	k-points: (4, 4, 1)
	•	Fermi smearing: 0.1 eV
	•	Slab: Pt(111), 3 layers
	•	Bottom two layers fixed during relaxation
	•	Geometry optimization: BFGS, fmax = 0.05 eV/Å

  
⸻

Results

Configuration  Total Energy (eV)  Adsorption Energy (eV)
Clean slab  -98.367235  —
CO(g)  -13.192740  —
CO on-top  -112.591177  -1.031202
CO bridge  -112.712397  -1.152422

⸻

Conclusion

Bridge site binds CO more strongly than on-top site.

E_{ads}^{bridge} < E_{ads}^{top}

This aligns with expected stronger multi-center metal–adsorbate interaction at bridge coordination.

⸻

Project Structure

cat-adsorption-dft/
│
├── inputs/
│   ├── 01_clean_slab.py
│   ├── 03_co_gas.py
│   ├── 04_co_on_top.py
│   ├── 07_co_on_bridge.py
│   ├── 08_compare_sites.py
│
├── outputs/           # (large .gpw files excluded via .gitignore)
├── energies_summary.txt
├── requirements.txt
└── README.md


⸻

How to Run

Create conda environment:

conda create -n gpaw-x64 python=3.10
conda activate gpaw-x64

Run in order:

python inputs/01_clean_slab.py
python inputs/03_co_gas.py
python inputs/04_co_on_top.py
python inputs/07_co_on_bridge.py
python inputs/08_compare_sites.py


⸻

Skills Demonstrated
	•	Surface slab construction
	•	Adsorbate placement (on-top and bridge)
	•	Geometry relaxation with constraints
	•	Extraction of DFT energies from output
	•	Adsorption energy comparison
	•	Debugging GPAW symmetry and geometry issues
	•	Reproducible workflow design

⸻

Next Extensions
	•	Add hollow site comparison
	•	Perform full relaxation for all sites
	•	Compute vibrational frequencies
	•	Add coverage dependence
	•	Build ML model to predict adsorption energies

⸻

This project represents a complete, reproducible DFT adsorption workflow suitable for heterogeneous catalysis research.
conda install -c conda-forge gpaw ase numpy

