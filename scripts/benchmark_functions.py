"""
Benchmark functions of the original 17_benchmarking notebook (cell 7), copied verbatim. Used by
gsdensity_figure_caches.py to replay the Fig. 2d contamination draws and Fig. 2e noise matrices (seed 0),
so that GSDensity is run on exactly the perturbations the other methods saw.
"""

from typing import Callable, Optional, Sequence, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import issparse
from scipy.stats import spearmanr


def _resolve_labels(adata, label_source, positive_class):
    if isinstance(label_source, np.ndarray):
        labels = label_source.astype(int)
    elif isinstance(label_source, str):
        if label_source not in adata.obs.columns:
            raise KeyError(f"'{label_source}' not found in adata.obs")
        col = adata.obs[label_source]
        if positive_class is not None:
            labels = (col == positive_class).astype(int).values
        else:
            vals = col.unique()
            if set(vals).issubset({0, 1, True, False}):
                labels = col.astype(int).values
            else:
                raise ValueError(
                    f"Column '{label_source}' is not binary. Set positive_class."
                )
    else:
        raise TypeError("label_source must be str or np.ndarray.")
    if len(labels) != adata.n_obs:
        raise ValueError(
            f"Label length ({len(labels)}) != adata.n_obs ({adata.n_obs})."
        )
    n_pos = labels.sum()
    if n_pos == 0 or n_pos == len(labels):
        raise ValueError(f"Need both classes. Found {n_pos}/{len(labels)} positives.")
    return labels


def benchmark_noise_robustness(
    adata: AnnData,
    gene_set: list[str],
    scoring_functions: dict[str, callable],
    noise_levels: Sequence[float] = (0.0, 0.25, 0.5, 1.0, 1.5, 2.0),
    n_repeats: int = 5,
    batch_key: str | None = None,
    seed: int = 0,
) -> pd.DataFrame:
    """
    Noise robustness benchmark (17_benchmarking cell 7, verbatim).

    Parameters
    ----------
    adata : AnnData
        Log-normalised data.
    gene_set : list of str
        Signature.
    scoring_functions : dict
        Method name -> scorer(adata, genes).
    noise_levels : sequence of float
        Gaussian noise SD as a multiple of each signature gene's SD.
    n_repeats : int
        Repeats per level.
    batch_key : str or None
        Not used here.
    seed : int
        Seed of the random stream.

    Returns
    -------
    pd.DataFrame
        method, noise_level, repeat, correlation (Spearman with the clean score).
    """
    rng = np.random.default_rng(seed)
    common_genes = [g for g in gene_set if g in adata.var_names]

    clean_scores: dict[str, np.ndarray] = {}
    for method_name, score_fn in scoring_functions.items():
        adata_clean = adata.copy()
        clean_scores[method_name] = np.asarray(
            score_fn(adata_clean, common_genes),
            dtype=np.float64,
        ).flatten()

    gene_stds = np.zeros(len(common_genes))
    for i, gene in enumerate(common_genes):
        x = adata[:, gene].X
        x = x.toarray().flatten() if issparse(x) else x.flatten()
        gene_stds[i] = x.std()

    records: list[dict] = []
    for noise_level in noise_levels:
        for rep in range(n_repeats):
            adata_noisy = adata.copy()
            X_dense = (
                adata_noisy.X.toarray()
                if issparse(adata_noisy.X)
                else adata_noisy.X.copy()
            )
            gene_idx = [list(adata.var_names).index(g) for g in common_genes]
            for j, gi in enumerate(gene_idx):
                X_dense[:, gi] += rng.normal(
                    0, gene_stds[j] * noise_level, X_dense.shape[0]
                )
            adata_noisy.X = X_dense

            # Also update layers that scorers might read from
            for layer_name in adata_noisy.layers:
                adata_noisy.layers[layer_name] = X_dense.copy()

            for method_name, score_fn in scoring_functions.items():
                noisy_scores = np.asarray(
                    score_fn(adata_noisy.copy(), common_genes),
                    dtype=np.float64,
                ).flatten()
                if batch_key is not None:
                    corrs = []
                    for patient in adata.obs[batch_key].unique():
                        mask = (adata.obs[batch_key] == patient).values
                        rho, _ = spearmanr(
                            clean_scores[method_name][mask], noisy_scores[mask]
                        )
                        if not np.isnan(rho):
                            corrs.append(rho)
                    rho_avg = np.mean(corrs) if corrs else 0.0
                else:
                    rho_avg, _ = spearmanr(clean_scores[method_name], noisy_scores)
                records.append(
                    {
                        "method": method_name,
                        "noise_level": noise_level,
                        "repeat": rep,
                        "correlation": rho_avg,
                    }
                )

    result = pd.DataFrame(records)
    return result


def benchmark_contamination(
    adata: AnnData,
    gene_set: list[str],
    scoring_functions: dict[str, callable],
    contamination_fractions: Sequence[float] = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0),
    n_repeats: int = 10,
    seed: int = 0,
) -> pd.DataFrame:
    """
    Signature contamination benchmark (17_benchmarking cell 7, verbatim).

    Parameters
    ----------
    adata : AnnData
        Log-normalised data.
    gene_set : list of str
        Clean signature.
    scoring_functions : dict
        Method name -> scorer(adata, genes).
    contamination_fractions : sequence of float
        Random non-signature genes added, as a multiple of the signature size.
    n_repeats : int
        Repeats per fraction.
    seed : int
        Seed of the random stream.

    Returns
    -------
    pd.DataFrame
        method, contamination, repeat, correlation (Spearman with the clean score).
    """
    rng = np.random.default_rng(seed)
    common_genes = [g for g in gene_set if g in adata.var_names]
    n_sig = len(common_genes)

    non_sig = [g for g in adata.var_names if g not in common_genes]

    if len(non_sig) < n_sig * max(contamination_fractions):
        import warnings

        warnings.warn(
            "Not enough non-signature genes for highest contamination. "
            "Sampling with replacement.",
            UserWarning,
        )

    # Clean reference scores
    clean_scores: dict[str, np.ndarray] = {}
    for method_name, score_fn in scoring_functions.items():
        clean_scores[method_name] = np.asarray(
            score_fn(adata.copy(), common_genes),
            dtype=np.float64,
        ).flatten()

    records: list[dict] = []

    for frac in contamination_fractions:
        n_contam = int(round(n_sig * frac))

        for rep in range(n_repeats):
            replace = len(non_sig) < n_contam
            noise_genes = (
                list(rng.choice(non_sig, n_contam, replace=replace))
                if n_contam > 0
                else []
            )
            contaminated_set = common_genes + noise_genes

            for method_name, score_fn in scoring_functions.items():
                try:
                    contam_scores = np.asarray(
                        score_fn(adata.copy(), contaminated_set),
                        dtype=np.float64,
                    ).flatten()
                except Exception:
                    contam_scores = np.zeros(adata.n_obs)

                rho, _ = spearmanr(clean_scores[method_name], contam_scores)
                if np.isnan(rho):
                    rho = 0.0

                records.append(
                    {
                        "method": method_name,
                        "contamination": frac,
                        "repeat": rep,
                        "correlation": rho,
                    }
                )

    result = pd.DataFrame(records)
    return result
