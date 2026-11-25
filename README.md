# EnrichMap: Spatially-aware gene set enrichment

Notebooks that reproduce the figures of the EnrichMap article. Each figure notebook is named after the panels it draws. See [Setup](#setup), [Data](#data) and [Running the notebooks](#running-the-notebooks) before you start.

## Figures and notebooks

| Figure | Panels | Notebook |
|---|---|---|
| Figure 1 | b–h | [`Fig1`](notebooks/Fig1.ipynb) |
| Figure 2 | b–f | [`Fig2`](notebooks/Fig2.ipynb) |
| Figure 3 | a–f | [`Fig3`](notebooks/Fig3.ipynb) |
| Figure 4 | a–h | [`Fig4_ExtendedDataFig8_1_scores`](notebooks/Fig4_ExtendedDataFig8_1_scores.ipynb) (scores), then [`Fig4_ExtendedDataFig8_2_panels`](notebooks/Fig4_ExtendedDataFig8_2_panels.ipynb) |
| Figure 5 | a–c | [`Fig5a-c_ExtendedDataFig9`](notebooks/Fig5a-c_ExtendedDataFig9.ipynb) |
| Figure 5 | d–n | [`Fig5d-n`](notebooks/Fig5d-n.ipynb) |
| Extended Data Fig. 1 | a–e | [`ExtendedDataFig1`](notebooks/ExtendedDataFig1.ipynb) |
| Extended Data Fig. 2 | a–f | [`ExtendedDataFig2`](notebooks/ExtendedDataFig2.ipynb) |
| Extended Data Fig. 3 | a–d | [`ExtendedDataFig3_ExtendedDataFig4a-b`](notebooks/ExtendedDataFig3_ExtendedDataFig4a-b.ipynb); gene sets of c from [`ExtendedDataFig3c_signatures`](notebooks/ExtendedDataFig3c_signatures.ipynb) |
| Extended Data Fig. 4 | a, b | [`ExtendedDataFig3_ExtendedDataFig4a-b`](notebooks/ExtendedDataFig3_ExtendedDataFig4a-b.ipynb) |
| Extended Data Fig. 4 | c–h | [`ExtendedDataFig4c-h_ExtendedDataFig5a-c`](notebooks/ExtendedDataFig4c-h_ExtendedDataFig5a-c.ipynb) |
| Extended Data Fig. 5 | a–c | [`ExtendedDataFig4c-h_ExtendedDataFig5a-c`](notebooks/ExtendedDataFig4c-h_ExtendedDataFig5a-c.ipynb) |
| Extended Data Fig. 5 | d | [`ExtendedDataFig5d`](notebooks/ExtendedDataFig5d.ipynb) |
| Extended Data Fig. 5 | e, f | [`ExtendedDataFig5e-f`](notebooks/ExtendedDataFig5e-f.ipynb) |
| Extended Data Fig. 6 | a–e | [`ExtendedDataFig6`](notebooks/ExtendedDataFig6.ipynb) |
| Extended Data Fig. 7 | a–d | [`ExtendedDataFig7`](notebooks/ExtendedDataFig7.ipynb) |
| Extended Data Fig. 8 | a, c–e | [`Fig4_ExtendedDataFig8_2_panels`](notebooks/Fig4_ExtendedDataFig8_2_panels.ipynb) |
| Extended Data Fig. 8 | b | [`ExtendedDataFig8b`](notebooks/ExtendedDataFig8b.ipynb) |
| Extended Data Fig. 9 | a–c | [`Fig5a-c_ExtendedDataFig9`](notebooks/Fig5a-c_ExtendedDataFig9.ipynb) |

Panels were reordered after the notebooks were run, so some output file names use the earlier letters:

| Output name in the notebook | Panel in the article |
|---|---|
| `ExtendedDataFig7a_knn_distance` | Extended Data Fig. 7b |
| `ExtendedDataFig7b_variograms` | Extended Data Fig. 7c |
| `ExtendedDataFig7c_maps` | Extended Data Fig. 7d |
| `ExtendedDataFig7d_ring_recovery` | Extended Data Fig. 7a |
| `ExtendedDataFig8b_bootstrap_rho` | Extended Data Fig. 8c |
| `ExtendedDataFig8c_comparator_maps` | Extended Data Fig. 8d |
| `ExtendedDataFig8d_comparators_vs_ground_truth` | Extended Data Fig. 8e |

## Setup

Python ≥ 3.10 is required (the published figures were made with Python 3.11). Create an environment, install the dependencies and start Jupyter:

```bash
git clone https://github.com/secrierlab/enrichmap-reproducibility.git
cd enrichmap-reproducibility
python -m venv .venv && source .venv/bin/activate   # or use a conda environment
pip install -r requirements.txt
jupyter lab
```

Start Jupyter from the repository root or from `notebooks/`. The notebooks locate the repository root themselves, so no paths need to be edited.

### R

Two parts of the analysis use R.

- **GSDensity** (a comparator method; the `REV_GSD_base_*` notebooks and `scripts/F2b_compute.py --part gsd`) runs in its own conda environment, `gsdensity_env`, built with `bash scripts/build_gsdensity_env.sh`. If the archived CRAN packages fail to install, run `bash scripts/build_gsdensity_env_step2.sh` afterwards. The notebooks call `scripts/run_gsdensity.R` after `conda activate gsdensity_env`. To activate R some other way, set `GSDENSITY_ACTIVATE` to the shell prefix to use, for example `export GSDENSITY_ACTIVATE="source ~/miniforge3/etc/profile.d/conda.sh && conda activate gsdensity_env && "`.
- **scran** size factors (`ExtendedDataFig6`) need an R with the Bioconductor packages `scran` and `BiocParallel`. The notebook runs `Rscript` from your `PATH`; set `RSCRIPT` to use another one.

## Data

Data files are not part of this repository.

| Data | Source | Where to put it | Used in |
|---|---|---|---|
| Visium human breast cancer (`adata_slides_0_3.h5ad`) | [Zenodo, 10.5281/zenodo.15438169](https://doi.org/10.5281/zenodo.15438169) (CC BY 4.0) | `processed_data/` | `Fig4_*`, `ExtendedDataFig5e-f`, `ExtendedDataFig8b`, `REV_GSD_base_breast` |
| Xenium (`xenium_aligned_ductal_invasive.h5ad`) | [Zenodo, 10.5281/zenodo.15438169](https://doi.org/10.5281/zenodo.15438169) | repository root | `Fig5d-n`, `REV_GSD_base_xenium` |
| Visium HD mouse brain (Space Ranger output) | [10X Genomics](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he-v4) | `Visium_HD_Mouse_Brain/` | `00_visiumhd_preprocessing` |
| Visium mouse brain | `squidpy.datasets.visium_hne_adata()`, downloaded automatically | `data/anndata/` | `00_visium_brain_preprocessing` |
| MERFISH, imaging mass cytometry | `squidpy.datasets.merfish()` and `imc()`, downloaded automatically | `data/anndata/` | `Fig5d-n`, `REV_GSD_base_merfish_imc` |
| PathwayCommons v12 | downloaded inside `00_hallmark_gene_sets` | repository root | `00_hallmark_gene_sets` |
| HRD and G0 arrest signatures | [`signature_all_phenotypes.csv`](signature_all_phenotypes.csv), in this repository | repository root | `Fig4_*`, `ExtendedDataFig5e-f`, `REV_GSD_base_breast` |

Simulated data (Fig. 1 and Extended Data Fig. 1) are generated inside the notebooks and scripts with fixed seeds.

## Running the notebooks

Figures are written to `figures/` and the values they plot to `source_data/`. Intermediate files go to `processed_data/`, `data/`, `logs/` and the repository root. All of these are generated and git-ignored. Run the notebooks in this order:

1. **Inputs.**
   - `00_gene_set` writes `Pyramidal_layer_signatures.pkl` and `Striatum_signatures.pkl`, the mouse brain region signatures.
   - `00_visium_brain_preprocessing` writes `processed_data/visium_brain_adata.h5ad`.
   - `00_visiumhd_preprocessing` (after `00_gene_set`) writes `processed_data/visiumhd_pyramidal_layer_scores.h5ad`.
   - `00_hallmark_gene_sets` writes `hallmarks_genelists.pkl`.
   - `ExtendedDataFig3c_signatures` writes `data/F1/ED3c_gene_sets.json`.
2. **GSDensity.** The `REV_GSD_base_*` notebooks write the GSDensity scores that the figure notebooks read from `data/gsd_cache/`. Run `REV_GSD_base_brain` after `ExtendedDataFig3c_signatures`.
3. **Simulations for Extended Data Fig. 1.** Run each part of `scripts/F2b_compute.py` over the 100 simulated slides; results go to `data/F2b/`. The parts are resumable and can run in parallel:
   ```bash
   python scripts/F2b_compute.py --part fast --seeds 0-99 --workers 2
   python scripts/F2b_compute.py --part ssgsea --seeds 0-99 --workers 2
   python scripts/F2b_compute.py --part gsva --seeds 0-99 --workers 3
   python scripts/F2b_compute.py --part gsva_noise --seeds 0-99 --workers 2   # after gsva
   python scripts/F2b_compute.py --part gsd --seeds 0-99 --workers 4
   python scripts/F2b_compute.py --part noiseclip --seeds 0-4 --workers 1     # after ssgsea and gsd
   ```
4. **Figures.** `Fig4_ExtendedDataFig8_1_scores` before `Fig4_ExtendedDataFig8_2_panels`; all other figure notebooks in any order.

The stored notebook outputs were produced on the authors' machine with the package versions noted in `requirements.txt`; re-running with newer versions can change numerical results slightly.

## License

The code in this repository is released under the GNU General Public License v3.0 (GPL-3.0-only); see [LICENSE](LICENSE). The Zenodo datasets are released under CC BY 4.0.

For details of the EnrichMap package and how to use it, see the documentation at https://enrichmap.readthedocs.io/en/stable.
