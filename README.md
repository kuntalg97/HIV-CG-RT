# HIV-CG-RT

This repository contains all codes and input files required to simulate HIV-1 capsid uncoating driven by reverse transcription.

We use an integrative coarse-grained (CG) approach to dynamically simulate reverse transcription, which converts flexible ssRNA into rigid dsDNA, inside the HIV-1 capsid. Beyond a critical threshold, dsDNA growth generates sufficient mechanical stress to drive capsid rupture, leading to uncoating.

Corresponding paper: [*Proc. Natl. Acad. Sci. U.S.A.* 123 (43) e2617714123 (2026)](https://doi.org/10.1073/pnas.2617714123)  
Datasets and models: [Zenodo](https://doi.org/10.5281/zenodo.22169970)

Model and simulation details:
- Capsid model: [Hudait and Voth, PNAS, 121 (4) e2313737121, 2024](https://www.pnas.org/doi/10.1073/pnas.2313737121)
- MD engine: [LAMMPS](https://www.lammps.org/#gsc.tab=0)
- Coarse-Grained Kinetic Monte Carlo (CG-KMC) algorithm: inputs and codes provided
