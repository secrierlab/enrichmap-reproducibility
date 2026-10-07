# EnrichMap: Spatially-aware gene set enrichment

Notebooks for reproducing the results of the EnrichMap [preprint](https://www.biorxiv.org/content/10.1101/2025.05.30.656960v1). Figure numbers below follow the preprint. See [Setup](#setup), [Data](#data) and [Running the notebooks](#running-the-notebooks) before you start.

1. Simulated dataset
- [Figure 1](notebooks/05_00_simulated_dataset.ipynb)
- [Supplementary Figure 1a–d](notebooks/05_01_simulated_benchmarking.ipynb)
- [Supplementary Figure 1e–f](notebooks/05_00_simulation_classification.ipynb)

2. Visium mouse brain dataset
- [Figure 2](notebooks/01_visium_platform.ipynb)
- [Supplementary Figure 2](notebooks/06_performance.ipynb)
- [Figure 3](notebooks/18_visium_weights.ipynb)
- [Supplementary Figure 3](notebooks/01_visium_platform.ipynb) and [here](notebooks/14_unit_tests.ipynb)
- [Supplementary Figure 4](notebooks/17_benchmarking.ipynb)
- [Supplementary Figures 5–6](notebooks/13_gene_set_sizes.ipynb) and [here](notebooks/01_visium_platform.ipynb)
- [Supplementary Figure 7](notebooks/06_performance.ipynb)
- [Supplementary Figure 8](notebooks/07_normalisation_effect.ipynb)
- [Supplementary Figure 9](notebooks/12_visium_platform_n_neighbours.ipynb)

3. Visium human breast cancer dataset
- [Figure 4](notebooks/08_cancer_hallmarks.ipynb), [here](notebooks/09_00_G0_arrest_in_breast_cancer.ipynb) and [here](notebooks/09_01_HRD_in_breast_cancer.ipynb)
- [Supplementary Figure 10](notebooks/08_cancer_hallmarks.ipynb), [here](notebooks/09_00_G0_arrest_in_breast_cancer.ipynb) and [here](notebooks/09_01_HRD_in_breast_cancer.ipynb)

4. Other platforms
Figure 5:
- [Visium HD](notebooks/02_visiumhd_platform.ipynb)
- [Xenium](notebooks/03_xenium_platform.ipynb)
- [MERFISH](notebooks/04_merfish_platform.ipynb)
- [Imaging Mass Cytometry](notebooks/10_masscytometry_platform.ipynb)

[Supplementary Figure 11](notebooks/02_visiumhd_platform.ipynb), [here](notebooks/03_xenium_platform.ipynb), [here](notebooks/04_merfish_platform.ipynb) and [here](notebooks/10_masscytometry_platform.ipynb)

## Setup

Python ≥ 3.10 is required (tested with Python 3.11). Create an environment, install the dependencies and start Jupyter:

```bash
git clone https://github.com/secrierlab/enrichmap-reproducibility.git
cd enrichmap-reproducibility
python -m venv .venv && source .venv/bin/activate   # or use a conda environment
pip install -r requirements.txt
jupyter lab
```

`07_normalisation_effect.ipynb` additionally uses R (scran normalisation through `rpy2`). It needs R with the Bioconductor packages `scran` and `BiocParallel`, plus `pip install -r requirements-r.txt`. All other notebooks run without R.

Start Jupyter from the repository root or from `notebooks/`; the notebooks locate the repository root automatically, so no paths need to be edited.

## Data

| Data | Source | Used in |
|---|---|---|
| Visium human breast cancer (`adata_slides_0_3.h5ad`) | [Zenodo, 10.5281/zenodo.15438169](https://doi.org/10.5281/zenodo.15438169) (CC BY 4.0) | `09_00`, `09_01` |
| Xenium (`xenium_aligned_ductal_invasive.h5ad`) | [Zenodo, 10.5281/zenodo.15438169](https://doi.org/10.5281/zenodo.15438169) | `03` |
| Visium HD mouse brain (`visium_hd_mouse_brain.h5ad`) | [10X Genomics](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he-v4) | processed Visium HD data, see `02` |
| Visium mouse brain | `squidpy.datasets.visium_hne_adata()`, downloaded automatically | `00`, `01`, `06`, `07`, `12`, `14`, `15`, `18` |
| MERFISH | `squidpy.datasets.merfish()`, downloaded automatically | `04` |
| Imaging mass cytometry | `squidpy.datasets.imc()`, downloaded automatically | `10` |
| PathwayCommons v12, MSigDB (through OmniPath), HGNC gene symbols | public downloads made inside the notebooks | `08`, `13`, `11` |

Download the Zenodo files into a `processed_data/` folder in the repository root (the notebooks create this folder). `02_visiumhd_platform.ipynb` reads the Visium HD Space Ranger output from a `Visium_HD_Mouse_Brain/` folder in the repository root.

## Running the notebooks

Figures are written to `figures/`, tables to `results/`, and intermediate files to `processed_data/` and the repository root. These are generated and git-ignored. Some notebooks use files written by others, so run them in this order:

1. `00_gene_set.ipynb` first. It writes `Pyramidal_layer_signatures.pkl` and `Striatum_signatures.pkl`, which `01`, `02`, `06`, `07`, `12`, `15`, `17` and `18` read.
2. `08_cancer_hallmarks.ipynb` before `09_00_G0_arrest_in_breast_cancer.ipynb` (it writes `hallmarks_genelists.pkl`).
3. `05_01_simulated_benchmarking.ipynb` before `05_00_simulation_classification.ipynb` (it writes `results/sim_classification.csv`).

The stored notebook outputs were produced on the authors' machine; re-running with newer package versions can change numerical results slightly.

## License

The code in this repository is released under the GNU General Public License v3.0 (GPL-3.0-only); see [LICENSE](LICENSE). The Zenodo datasets are released under CC BY 4.0.

For more details of the EnrichMap package and learn how to use `EnrichMap` package please go to documentation at https://enrichmap.readthedocs.io/en/stable.