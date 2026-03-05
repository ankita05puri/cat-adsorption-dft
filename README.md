# Periodic DFT Surface Modeling

## CO Adsorption and Reaction Setup on Pt(111)

This repository implements a reproducible periodic Density Functional Theory (DFT) workflow for studying adsorption and reaction configurations on catalytic metal surfaces.

The system studied is CO oxidation on Pt(111), a classical model reaction in heterogeneous catalysis. The workflow demonstrates how electronic structure calculations can be used to compute adsorption energetics and construct reaction geometries for further mechanistic modeling.

Calculations are performed using ASE (Atomic Simulation Environment) and GPAW.

## Modeling Workflow

The repository follows a typical surface catalysis modeling pipeline:

```bash
Pt(111) Surface Construction
        ↓
Slab Relaxation
        ↓
Gas-Phase Reference Calculations
        ↓
Adsorption Geometry Optimization
        ↓
Adsorption Energy Analysis
        ↓
Co-Adsorption Configuration
        ↓
Reaction Pathway Setup (NEB)
``` 

This workflow forms the electronic-structure foundation for microkinetic modeling and catalytic performance analysis.

## System Studied

Surface: Pt(111)
Reaction: CO Oxidation

Elementary processes explored:
	•	CO adsorption
	•	O adsorption
	•	CO + O co-adsorption
	•	CO oxidation reaction pathway setup

## Computational Details

Density Functional Theory calculations were performed using:

### DFT Method
	•	Exchange–correlation functional: PBE
	•	Basis: Plane-wave (PW400)
	•	Brillouin zone sampling: 4 × 4 × 1 k-point mesh

### Surface Model
	•	Periodic Pt(111) slab
	•	Four atomic layers
	•	Bottom layers fixed to represent bulk lattice
	•	Vacuum region (~15 Å) to avoid interactions between periodic images

### Geometry Optimization
	•	Optimizer: BFGS
	•	Convergence criterion: fmax < 0.05 eV/Å

## Adsorption Energy Defination

Adsorption energies are computed as: E_ads = E(slab + adsorbate) − E(slab) − E(gas)

For atomic oxygen adsorption: E_ads(O) = E(slab + O) − E(slab) − ½ E(O₂)

## Repository Structure

```bash
cat-adsorption-dft/
│
├── inputs/
│   01_slab_relax.py        # Build and relax Pt(111) slab
│   02_co_gas.py            # Gas-phase CO reference
│   03_co_on_top.py         # CO adsorption on top site
│   04_o_on_hollow.py       # O adsorption on hollow sites
│   05_o2_gas.py            # Gas-phase O2 reference
│   06_extract_energies.py  # Adsorption energy calculation
│   07_co_o_coadsorb.py     # CO + O co-adsorption configuration
│   08_neb_co_oxidation.py  # Reaction pathway setup (NEB)
│
├── outputs/
│   Optimized structures and calculation outputs
│
├── figures/
│   Surface structures and visualization images
│
├── data/
│   Extracted adsorption energies and analysis results
│
└── README.md
```
## Running the Workflow

Run calculations sequentially.

Surface and adsorption calculations:
```bash
python inputs/01_slab_relax.py
python inputs/02_co_gas.py
python inputs/03_co_on_top.py
python inputs/04_o_on_hollow.py
python inputs/05_o2_gas.py
python inputs/06_extract_energies.py
```

Advanced configurations:
```bash
python inputs/07_co_o_coadsorb.py
python inputs/08_neb_co_oxidation.py
```

## Purpose of the Project

This repository demonstrates how periodic DFT calculations can be used to quantify surface adsorption energetics and construct reaction configurations for catalytic systems.

The results generated here serve as inputs for microkinetic modeling of catalytic performance, forming the first stage of a multiscale modeling pipeline.

## Tools Used

Electronic Structure
	•	GPAW
	•	ASE

Scientific Programming
	•	Python
	•	NumPy
	•	Matplotlib

# Author

## Ankita Puri
Computational Materials Scientist
Heterogeneous Catalysis | Atomistic Modeling | Physics-Informed ML

📍 Portland, OR
📧 ankit05puri@gmail.com
🔗 https://linkedin.com/in/ankita-puri-phd

