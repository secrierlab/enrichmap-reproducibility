from __future__ import annotations
import os
from typing import Sequence
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import squidpy as sq
import decoupler as dc
import scanpy as sc
import gseapy as gp
import enrichmap as em
from anndata import AnnData
from scipy.sparse import issparse
from scipy.stats import spearmanr, ks_2samp
from sklearn.cluster import DBSCAN
from sklearn.metrics import (
    average_precision_score,
    f1_score,
    precision_recall_curve,
    roc_auc_score,
)
from matplotlib.collections import QuadMesh
import anndata
from scipy.sparse import issparse, csr_matrix
from scipy.spatial.distance import cdist
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans
from sklearn.metrics import (
    average_precision_score,
    f1_score,
    precision_recall_curve,
    roc_auc_score,
)
import colorsys
import json
import time
import anndata as ad
import matplotlib.transforms as mtransforms
from scipy.sparse import csr_matrix, issparse
from scipy.stats import norm, t as t_dist
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
)
import matplotlib.gridspec as gridspec
from pygam import LinearGAM, te
from scipy import stats
from matplotlib.colors import SymLogNorm
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score
from pygam import LinearGAM, s, te
import pygam
from scipy.sparse import csr_matrix, diags, identity, issparse
from scipy.stats import t as t_dist
from sklearn.metrics import roc_auc_score
from scipy.sparse import csr_matrix, csc_matrix, issparse
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr, spearmanr, rankdata, t as t_dist
from scipy.stats import t as student_t
from esda.geary import Geary
from esda.moran import Moran
from libpysal.weights import KNN
from scipy.spatial import KDTree
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import f1_score, precision_recall_curve, roc_auc_score
from enrichmap.tools._score import _signed_weights
from libpysal.weights.spatial_lag import lag_spatial
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist, squareform
from scipy.stats import linregress, pearsonr, rankdata, spearmanr
import zlib
from matplotlib.colors import to_rgba
from scipy.stats import mannwhitneyu, rankdata
from sklearn.metrics import roc_auc_score, roc_curve
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from enrichmap.tools import _score
import subprocess
from scipy.io import mmwrite
import re
from sklearn.metrics import precision_recall_curve
import warnings
from adjustText import adjust_text
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, average_precision_score
import urllib.parse
import urllib.request
import mpmath


# Shell prefix that puts GSDensity's R on PATH before Rscript runs (see README, "GSDensity").
GSDENSITY_ACTIVATE = os.environ.get(
    "GSDENSITY_ACTIVATE", 'eval "$(conda shell.bash hook)" && conda activate gsdensity_env && '
)

def benchmark_noise_robustness(
    adata: AnnData,
    gene_set: list[str],
    scoring_functions: dict[str, callable],
    noise_levels: Sequence[float] = (0.0, 0.25, 0.5, 1.0, 1.5, 2.0),
    n_repeats: int = 5,
    batch_key: str | None = None,
    seed: int = 0,
    plot: bool = True,
    figsize: tuple[float, float] = (6, 4),
    palette: str = "muted",
    save: str | None = None,
) -> pd.DataFrame:
    """Noise robustness benchmark. See previous version for full docstring."""
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
    if plot and not result.empty:
        _plot_noise_robustness(result, figsize=figsize, palette=palette, save=save)
    return result


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


def _plot_noise_robustness(df, figsize=(6, 4), palette="muted", save=None):
    fig, ax = plt.subplots(figsize=figsize)
    sns.lineplot(
        data=df,
        x="noise_level",
        y="correlation",
        hue="method",
        style="method",
        markers=True,
        dashes=False,
        err_style="band",
        errorbar="sd",
        palette=palette,
        ax=ax,
    )
    ax.set_xlabel("Noise level (× gene σ)", fontsize=9)
    ax.set_ylabel("Spearman ρ (noisy vs clean)", fontsize=9)
    ax.set_title("Noise robustness", fontsize=10)
    ax.set_ylim(0, 1.05)
    ax.tick_params(labelsize=8)
    ax.legend(title="", fontsize=8, frameon=False)
    ax.grid(False)
    sns.despine(ax=ax)
    plt.tight_layout()
    if save:
        os.makedirs(os.path.dirname(save) or "figures", exist_ok=True)
        path = os.path.join("figures", save) if not os.path.dirname(save) else save
        plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.show()
    return ax


def benchmark_contamination(
    adata: AnnData,
    gene_set: list[str],
    scoring_functions: dict[str, callable],
    contamination_fractions: Sequence[float] = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0),
    n_repeats: int = 10,
    seed: int = 0,
    plot: bool = True,
    figsize: tuple[float, float] = (6, 4),
    palette: str = "muted",
    save: str | None = None,
) -> pd.DataFrame:
    """
    Benchmark scoring methods under signature contamination.

    Progressively adds random irrelevant genes to the signature and
    measures how much each method's score degrades. A method with gene
    weighting (EnrichMap) should be robust because it downweights the
    irrelevant genes. Methods that treat all genes equally degrade
    proportionally because each random gene dilutes the signal.

    The metric is the Spearman correlation between each method's
    contaminated-signature score and its own clean-signature score.
    This is method-internal consistency — it does not depend on ground
    truth labels or classification performance.

    Parameters
    ----------
    adata : AnnData
    gene_set : list of str
        Clean gene signature.
    scoring_functions : dict of str -> callable
        Each takes ``(adata, gene_list)`` and returns a 1D score array.
    contamination_fractions : sequence of float
        Random genes to add as a fraction of signature size. 1.0 means
        add as many random genes as real genes (50% contamination).
    n_repeats : int, default 10
    seed : int, default 0
    plot, figsize, palette, save
        Plotting options.

    Returns
    -------
    pd.DataFrame
        Columns: ``method``, ``contamination``, ``repeat``, ``correlation``.
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

    if plot and not result.empty:
        _plot_contamination(result, n_sig, figsize=figsize, palette=palette, save=save)

    return result


def _plot_contamination(df, n_sig, figsize=(6, 4), palette="muted", save=None):
    """Degradation curves: score consistency vs contamination."""
    fig, ax = plt.subplots(figsize=figsize)
    sns.lineplot(
        data=df,
        x="contamination",
        y="correlation",
        hue="method",
        style="method",
        markers=True,
        dashes=False,
        err_style="band",
        errorbar="sd",
        palette=palette,
        ax=ax,
    )
    ax.set_xlabel(f"Random genes added (\u00d7 {n_sig} signature genes)", fontsize=9)
    ax.set_ylabel("Spearman \u03c1 (contaminated vs clean)", fontsize=9)
    ax.set_title("Signature contamination robustness", fontsize=10)
    ax.set_ylim(0, 1.05)
    ax.tick_params(labelsize=8)
    ax.legend(title="", fontsize=8, frameon=False)
    ax.grid(False)
    sns.despine(ax=ax)
    plt.tight_layout()

    if save:
        os.makedirs(os.path.dirname(save) or "figures", exist_ok=True)
        path = os.path.join("figures", save) if not os.path.dirname(save) else save
        plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.show()
    return ax


def benchmark_spatial_concordance(
    adata: AnnData,
    score_keys: Sequence[str],
    label_source: str | np.ndarray,
    positive_class: str | None = None,
    n_neighbors: int = 6,
    smoothing_iterations: Sequence[int] = (1, 2, 3),
    plot: bool = True,
    figsize: tuple[float, float] = (5, 4),
    palette: str = "muted",
    save: str | None = None,
) -> pd.DataFrame:
    """
    Evaluate how well each method's continuous scores match the true
    spatial domain structure.

    Cell-type domains in tissue do not have razor-sharp boundaries.
    Spots deep inside a domain are unambiguously positive; spots at the
    boundary are in a transition zone; spots far away are clearly
    negative. Hard binary labels (0/1) are a discretisation of this
    continuous spatial reality, and metrics that binarise predictions
    (F1, AUPRC) penalise methods that produce soft, biologically
    realistic boundaries.

    This function creates a continuous "spatial ground truth" by
    iteratively averaging the binary labels over the spatial neighbours
    graph (equivalent to spatial smoothing of the labels). The result
    is a continuous field with values near 1.0 deep inside the domain,
    intermediate values at boundaries and near 0.0 far from the domain.
    The Spearman correlation between each method's scores and this
    continuous spatial ground truth is then computed.

    Multiple smoothing iterations are used to test concordance at
    different spatial scales: iteration=1 captures the immediate
    boundary zone, iteration=3 captures broader transition regions.

    A method that produces spatially coherent, domain-shaped scores
    (like EnrichMap with spatial smoothing) will correlate highly with
    the continuous spatial ground truth. A method that produces patchy,
    per-spot scores will correlate less well because its noise does not
    match the smooth domain structure.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with spatial coordinates in ``adata.obsm``
        and pre-computed scores in ``adata.obs``.

    score_keys : sequence of str
        Column names in ``adata.obs`` for each scoring method's output.

    label_source : str or np.ndarray
        Ground truth labels (same format as other benchmark functions).

    positive_class : str or None
        Category to treat as positive.

    n_neighbors : int, default 6
        Number of spatial neighbours for constructing the graph.

    smoothing_iterations : sequence of int, default (1, 2, 3)
        Number of neighbour-averaging iterations applied to the binary
        labels to create the continuous spatial ground truth. Each
        iteration spreads the domain boundary further. Results for each
        iteration are reported separately.

    plot : bool, default True
    figsize, palette, save
        Plotting options.

    Returns
    -------
    pd.DataFrame
        Tidy dataframe with columns ``method``, ``smoothing_iter``,
        ``concordance``.

    Examples
    --------
    >>> results = benchmark_spatial_concordance(
    ...     adata,
    ...     score_keys=["enrichmap_score", "aucell_score", "scanpy_score"],
    ...     label_source="cell_type",
    ...     positive_class="Pyramidal_layer",
    ... )
    """
    labels = _resolve_labels(adata, label_source, positive_class)

    # Build spatial graph
    adata_tmp = adata.copy()
    sq.gr.spatial_neighbors(
        adata_tmp,
        n_neighs=n_neighbors,
        coord_type="generic",
        key_added="spatial",
    )
    W = adata_tmp.obsp["spatial_connectivities"]

    # Row-normalise: each row sums to 1
    row_sums = np.array(W.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1.0
    W_norm = W.multiply(1.0 / row_sums[:, None])

    records: list[dict] = []

    for n_iter in smoothing_iterations:
        # Create continuous spatial ground truth by iteratively smoothing labels
        spatial_gt = labels.astype(np.float64).copy()
        for _ in range(n_iter):
            spatial_gt = W_norm.dot(spatial_gt)

        for method in score_keys:
            if method not in adata.obs.columns:
                continue

            scores = adata.obs[method].values.astype(np.float64)
            valid = ~np.isnan(scores)
            rho, _ = spearmanr(scores[valid], spatial_gt[valid])
            if np.isnan(rho):
                rho = 0.0

            records.append(
                {
                    "method": method,
                    "smoothing_iter": n_iter,
                    "concordance": rho,
                }
            )

    result = pd.DataFrame(records)

    if plot and not result.empty:
        _plot_spatial_concordance(
            result,
            figsize=figsize,
            palette=palette,
            save=save,
        )

    return result


def _plot_spatial_concordance(df, figsize=(5, 4), palette="muted", save=None):
    """
    Bar chart (single iteration) or line plot (multiple iterations)
    of spatial concordance by method.
    """
    df = df.copy()
    df["method_label"] = df["method"].str.replace("_score", "", regex=False)

    n_iters = df["smoothing_iter"].nunique()

    if n_iters == 1:
        fig, ax = plt.subplots(figsize=figsize)
        sns.barplot(
            data=df,
            x="method_label",
            y="concordance",
            hue="method_label",
            palette=palette,
            edgecolor="white",
            linewidth=0.5,
            legend=False,
            ax=ax,
        )
        ax.set_ylabel("Spatial concordance (Spearman ρ)", fontsize=9)
        ax.set_xlabel("")
        ax.set_title("Neighbourhood enrichment concordance", fontsize=10)
        ax.tick_params(axis="x", rotation=45, labelsize=8)
    else:
        fig, ax = plt.subplots(figsize=figsize)
        sns.lineplot(
            data=df,
            x="smoothing_iter",
            y="concordance",
            hue="method_label",
            style="method_label",
            markers=True,
            dashes=False,
            palette=palette,
            ax=ax,
        )
        ax.set_xlabel("Smoothing iterations", fontsize=9)
        ax.set_ylabel("Spatial concordance (Spearman ρ)", fontsize=9)
        ax.set_title("Neighbourhood enrichment concordance", fontsize=10)
        ax.legend(title="", fontsize=8, frameon=False)

    ax.set_ylim(bottom=0)
    ax.tick_params(labelsize=8)
    ax.grid(False)
    sns.despine(ax=ax)
    plt.tight_layout()
    ax.tick_params(labelsize=8)
    ax.grid(False)
    sns.despine(ax=ax)
    plt.tight_layout()

    if save:
        os.makedirs(os.path.dirname(save) or "figures", exist_ok=True)
        path = os.path.join("figures", save) if not os.path.dirname(save) else save
        plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.show()
    return ax


# 7. Scoring functions for the robustness benchmarks
def make_enrichmap_null_scorer():
    def scorer(ad, gene_list):
        em.tl.score(
            ad,
            gene_set=gene_list,
            score_key="_bench",
            smoothing=True,
            null_model=True,
            correct_spatial_covariates=False,
        )
        return ad.obs["_bench_score"].values

    return scorer


# EnrichMap: smoothing-only
def make_enrichmap_scorer():
    def scorer(ad, gene_list):
        em.tl.score(
            ad,
            gene_set=gene_list,
            score_key="_bench",
            smoothing=True,
            null_model=False,
            correct_spatial_covariates=False,
        )
        return ad.obs["_bench_score"].values

    return scorer


def make_scanpy_scorer():
    def scorer(ad, gene_list):
        sc.tl.score_genes(ad, gene_list=gene_list, score_name="_bench", use_raw=False)
        return ad.obs["_bench"].values

    return scorer


def make_aucell_scorer():
    def scorer(ad, gene_list):
        net = pd.DataFrame({"source": "_bench", "target": gene_list})
        dc.mt.aucell(ad, net=net)
        return (
            dc.pp.get_obsm(adata=ad, key="score_aucell")
            .obsm["score_aucell"]["_bench"]
            .values
        )

    return scorer


def make_zscore_scorer():
    def scorer(ad, gene_list):
        X = ad.X.toarray() if issparse(ad.X) else ad.X
        expr_df = pd.DataFrame(X.T, index=ad.var_names, columns=ad.obs_names)
        common = [g for g in gene_list if g in expr_df.index]
        return (expr_df.loc[common].sum(axis=0) / np.sqrt(len(common))).values

    return scorer


def make_gsva_scorer():
    def scorer(ad, gene_list):
        X = ad.X.toarray() if issparse(ad.X) else ad.X
        expr_df = pd.DataFrame(X.T, index=ad.var_names, columns=ad.obs_names)
        result = gp.gsva(
            data=expr_df, gene_sets={"_bench": gene_list}, outdir=None, min_size=3
        )
        es = result.res2d.pivot(index="Term", columns="Name", values="ES")
        return es.loc["_bench"].T.astype(float).reindex(ad.obs_names).values

    return scorer


def make_ssgsea_scorer():
    def scorer(ad, gene_list):
        X = ad.X.toarray() if issparse(ad.X) else ad.X
        expr_df = pd.DataFrame(X.T, index=ad.var_names, columns=ad.obs_names)
        expr_df.index = expr_df.index.str.upper()
        gl_upper = [g.upper() for g in gene_list]
        result = gp.ssgsea(
            data=expr_df,
            gene_sets={"_bench": gl_upper},
            outdir=None,
            sample_norm_method="rank",
            no_plot=True,
            min_size=3,
        )
        nes = result.res2d.pivot(index="Term", columns="Name", values="NES")
        return nes.loc["_bench"].T.astype(float).reindex(ad.obs_names).values

    return scorer


# --------------------------------------------------------------------------------------------------
# _common
# --------------------------------------------------------------------------------------------------


def set_box_ratio(fig=None, ratio=1.25, axes=None):
    """
    Set the width / height ratio of data axes boxes: 1.25 (x longer than y) for benchmark plots,
    1 (square) for correlation plots (one quantity against another, e.g. Spearman scatter plots).

    Spatial maps (equal aspect), heatmaps (QuadMesh), images and colour bars keep their shape.

    Parameters
    ----------
    fig : matplotlib.figure.Figure or None
        Figure to adjust (default: the current figure).
    ratio : float
        Width / height of each axes box (1.25 benchmark, 1 correlation).
    axes : list of matplotlib.axes.Axes or None
        Only these axes (default: every axes of `fig`), for figures that mix both kinds of panel.
    """
    fig = plt.gcf() if fig is None else fig
    for ax in fig.axes if axes is None else axes:
        if (
            ax.get_label() == "<colorbar>"
            or ax.images
            or any(isinstance(c, QuadMesh) for c in ax.collections)
            or ax.get_aspect() == 1.0
        ):
            continue
        ax.set_box_aspect(1 / ratio)


# --------------------------------------------------------------------------------------------------
# AUD_cv_vs_logfc
# --------------------------------------------------------------------------------------------------


def simulate_visium(
    seed=0,
    n_genes=3000,
    n_signature_genes=20,
    n_clusters=4,
    target_cluster=1,
    n_spots_range=(1500, 4000),
    spot_spacing_um=100,
    spot_diameter_um=55,
    mean_counts_per_spot=5000,
    dispersion_range=(0.1, 2.0),
    fold_change_range=(1.5, 5.0),
    dropout_rate=0.1,
):
    """
    Simulate a 10x Visium spatial transcriptomics dataset.

    The grid layout matches the Visium platform: a square capture area
    with hexagonally arranged spots at 100 um centre-to-centre spacing.
    The number of spots is drawn randomly from ``n_spots_range`` to
    mimic variable tissue coverage across real slides. This is achieved
    by varying the capture area dimensions per simulation.

    Returns (adata, sig_genes, fold_changes, target_label).
    """
    rng = np.random.default_rng(seed)

    # ── Variable spot count via variable capture area ────────────
    # Draw a target spot count, then compute the grid side that
    # produces approximately that many spots in a hexagonal grid.
    # For a hex grid: n_spots ≈ n_rows × n_cols ≈ side^2 / 0.866
    n_spots_target = rng.integers(
        n_spots_range[0],
        n_spots_range[1] + 1,
    )
    grid_side = int(np.ceil(np.sqrt(n_spots_target * 0.866)))
    capture_area_um = grid_side * spot_spacing_um

    # ── Hexagonal grid within square capture area ────────────────
    n_cols = int(capture_area_um / spot_spacing_um)
    n_rows = int(capture_area_um / (spot_spacing_um * np.sqrt(3) / 2))

    x_coords, y_coords = [], []
    for row in range(n_rows):
        for col in range(n_cols):
            x = col * spot_spacing_um + (row % 2) * (spot_spacing_um / 2)
            y = row * spot_spacing_um * np.sqrt(3) / 2
            if x <= capture_area_um and y <= capture_area_um:
                x_coords.append(x)
                y_coords.append(y)

    x_coords = np.array(x_coords)
    y_coords = np.array(y_coords)
    n_spots = len(x_coords)
    coords = np.column_stack([x_coords, y_coords])

    # ── Spatially coherent clusters via Voronoi tessellation ─────
    margin = capture_area_um * 0.15
    seed_coords = np.column_stack(
        [
            rng.uniform(margin, capture_area_um - margin, n_clusters),
            rng.uniform(margin, capture_area_um - margin, n_clusters),
        ]
    )

    dists = cdist(coords, seed_coords)
    dists += rng.normal(0, spot_spacing_um * 0.8, dists.shape)
    cluster_ids = np.argmin(dists, axis=1)
    cluster_labels = np.array([f"cluster_{i + 1}" for i in cluster_ids])

    target_label = f"cluster_{target_cluster + 1}"

    # ── Gene names ──────────────────────────────────────────────
    gene_names = [f"GENE_{i + 1:05d}" for i in range(n_genes)]

    # ── Signature genes with graded fold changes ────────────────
    sig_indices = rng.choice(n_genes, n_signature_genes, replace=False)
    fold_changes = rng.uniform(
        fold_change_range[0],
        fold_change_range[1],
        n_signature_genes,
    )
    order = np.argsort(-fold_changes)
    sig_indices = sig_indices[order]
    fold_changes = fold_changes[order]
    sig_genes = [gene_names[i] for i in sig_indices]

    # ── Expression matrix (negative binomial) ───────────────────
    gene_means = np.exp(rng.normal(1.5, 1.0, n_genes))
    gene_means = np.clip(gene_means, 0.1, 100)

    dispersions = rng.uniform(
        dispersion_range[0],
        dispersion_range[1],
        n_genes,
    )

    lib_size_factors = np.exp(rng.normal(0, 0.3, n_spots))
    lib_size_factors *= mean_counts_per_spot / (
        gene_means.sum() * lib_size_factors.mean()
    )

    X = np.zeros((n_spots, n_genes), dtype=np.float32)
    target_mask = cluster_labels == target_label

    for g in range(n_genes):
        mu = gene_means[g] * lib_size_factors

        if g in sig_indices:
            fc_idx = np.where(sig_indices == g)[0][0]
            mu[target_mask] *= fold_changes[fc_idx]

        r = dispersions[g]
        p = r / (r + mu)
        counts = rng.negative_binomial(r, p)

        dropout_mask = rng.random(n_spots) < dropout_rate
        counts[dropout_mask] = 0

        X[:, g] = counts

    # ── Assemble AnnData ────────────────────────────────────────
    adata = anndata.AnnData(X=csr_matrix(X))
    adata.var_names = gene_names
    adata.obs_names = [f"spot_{i}" for i in range(n_spots)]
    adata.obs["cluster"] = pd.Categorical(cluster_labels)
    adata.obsm["spatial"] = coords

    adata.raw = adata.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    return adata, sig_genes, fold_changes, target_label


def cv_weight_gene_table(seed):
    """
    Per-gene CV weight, SD and baseline expression for one simulated slide.

    Parameters
    ----------
    seed : int
        Simulation seed.

    Returns
    -------
    pd.DataFrame
        One row per signature gene.
    """
    adata, sig_genes, fold_changes, target_label = simulate_visium(seed=seed)
    weights = em.tl.infer_gene_weights(adata, sig_genes)

    X = adata[:, sig_genes].X
    X = X.toarray() if issparse(X) else np.asarray(X)
    is_target = (adata.obs["cluster"] == target_label).values

    return pd.DataFrame(
        {
            "seed": seed,
            "gene": sig_genes,
            "fold_change": fold_changes,
            "log2_fold_change": np.log2(fold_changes),
            "cv_weight": [weights[g] for g in sig_genes],
            "sd": X.std(axis=0),
            "baseline_mean": X[~is_target].mean(axis=0),
        }
    )


# --------------------------------------------------------------------------------------------------
# R4_D1_F1_audit
# --------------------------------------------------------------------------------------------------


def median_f1_ceiling(prevalence):
    """
    Maximum F1 when half of the spots are called positive.

    Parameters
    ----------
    prevalence : float or array-like
        Fraction of positive spots p.

    Returns
    -------
    float or np.ndarray
        2 min(p, 1/2) / (1/2 + p), i.e. 4p / (1 + 2p) for p <= 1/2.
    """
    p = np.asarray(prevalence, dtype=np.float64)
    return 2 * np.minimum(p, 0.5) / (0.5 + p)


def all_positive_f1(prevalence):
    """
    F1 of calling every spot positive.

    Parameters
    ----------
    prevalence : float or array-like
        Fraction of positive spots p.

    Returns
    -------
    float or np.ndarray
        2p / (1 + p).
    """
    p = np.asarray(prevalence, dtype=np.float64)
    return 2 * p / (1 + p)


def _benchmark_pooled(adata, score_keys, labels, metrics, n_neighbors=6):
    records: list[dict] = []
    spatial_info = None
    if any(m in metrics for m in ("spatial_coherence", "boundary_f1")):
        spatial_info = _build_spatial_info(adata, labels, n_neighbors)
    for method in score_keys:
        if method not in adata.obs.columns:
            continue
        scores = adata.obs[method].values.astype(np.float64)
        valid = ~np.isnan(scores)
        s, y = scores[valid], labels[valid]
        if y.sum() == 0 or y.sum() == len(y):
            continue
        if "auroc" in metrics:
            try:
                val = roc_auc_score(y, s)
            except ValueError:
                val = np.nan
            records.append(
                {"method": method, "metric": "AUROC", "value": val, "patient": "pooled"}
            )
        if "auprc" in metrics:
            try:
                val = average_precision_score(y, s)
            except ValueError:
                val = np.nan
            records.append(
                {"method": method, "metric": "AUPRC", "value": val, "patient": "pooled"}
            )
        if "f1" in metrics:
            f1_opt, _ = _optimal_f1(y, s)
            records.append(
                {
                    "method": method,
                    "metric": "F1 (optimal)",
                    "value": f1_opt,
                    "patient": "pooled",
                }
            )
        if any(m in metrics for m in ("spatial_coherence", "boundary_f1")):
            _, best_thr = _optimal_f1(y, s)
            q = float(np.clip(np.mean(s <= best_thr), 0.01, 0.99))
            preds = (scores > np.quantile(scores, q)).astype(int)
            if "spatial_coherence" in metrics and spatial_info is not None:
                val = _spatial_coherence(
                    preds, spatial_info["W"], spatial_info["gt_morans_i"]
                )
                records.append(
                    {
                        "method": method,
                        "metric": "Spatial coherence",
                        "value": val,
                        "patient": "pooled",
                    }
                )
            if "boundary_f1" in metrics and spatial_info is not None:
                val = _boundary_f1(labels, preds, spatial_info["boundary_mask"])
                records.append(
                    {
                        "method": method,
                        "metric": "Boundary F1",
                        "value": val,
                        "patient": "pooled",
                    }
                )
    return pd.DataFrame(records)


def _build_spatial_info(adata, labels, n_neighbors):
    adata_tmp = adata.copy()
    sq.gr.spatial_neighbors(
        adata_tmp, n_neighs=n_neighbors, coord_type="generic", key_added="spatial"
    )
    W = adata_tmp.obsp["spatial_connectivities"]
    gt_morans_i = _morans_i_continuous(labels.astype(float), W)
    boundary_mask = np.zeros(len(labels), dtype=bool)
    for i in range(len(labels)):
        neighbours = W[i].nonzero()[1]
        if len(neighbours) > 0 and np.any(labels[neighbours] != labels[i]):
            boundary_mask[i] = True
    return {"W": W, "gt_morans_i": gt_morans_i, "boundary_mask": boundary_mask}


def _morans_i_continuous(values, W):
    values = np.asarray(values, dtype=np.float64)
    n = len(values)
    z = values - values.mean()
    z2 = z @ z
    if z2 == 0:
        return 0.0
    return float((n / W.sum()) * (z @ W.dot(z) / z2))


def _spatial_coherence(predictions, W, gt_morans_i):
    pred_morans_i = _morans_i_continuous(predictions.astype(float), W)
    if gt_morans_i <= 0:
        return max(pred_morans_i, 0.0)
    return float(np.clip(pred_morans_i / gt_morans_i, 0.0, 1.0))


def _boundary_f1(y_true, y_pred, boundary_mask):
    if boundary_mask.sum() == 0:
        return np.nan
    y_b, p_b = y_true[boundary_mask], y_pred[boundary_mask]
    if y_b.sum() == 0 or y_b.sum() == len(y_b):
        return np.nan
    return float(f1_score(y_b, p_b, zero_division=0))


def _optimal_f1(y_true, scores):
    precision, recall, thresholds = precision_recall_curve(y_true, scores)
    precision, recall = precision[:-1], recall[:-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        f1 = 2 * (precision * recall) / (precision + recall)
    f1 = np.nan_to_num(f1, nan=0.0)
    best_idx = np.argmax(f1)
    return float(f1[best_idx]), float(thresholds[best_idx])


def d1_f1_calls(y, pred):
    """
    F1 of binary calls.

    Parameters
    ----------
    y : array-like
        Binary labels.
    pred : array-like
        Binary calls.

    Returns
    -------
    float
        2 TP / (n called + n positive), 0 if both are empty.
    """
    y = np.asarray(y).astype(bool)
    pred = np.asarray(pred).astype(bool)
    denom = pred.sum() + y.sum()
    return float(2 * (y & pred).sum() / denom) if denom else 0.0


def d1_f1_median(y, s):
    """F1 when the score is thresholded at its median."""
    return d1_f1_calls(y, s > np.median(s))


def d1_f1_topk(y, s):
    """F1 when the k highest-scoring spots are called positive, k = number of positives."""
    k = int(np.asarray(y).sum())
    pred = np.zeros(len(y), bool)
    pred[np.argsort(-np.asarray(s), kind="stable")[:k]] = True
    return d1_f1_calls(y, pred)


def d1_otsu_threshold(s, nbins=256):
    """
    Otsu threshold of a score (label-free).

    Parameters
    ----------
    s : array-like
        Scores.
    nbins : int
        Histogram bins.

    Returns
    -------
    float
        Bin centre that maximises the between-class variance.
    """
    hist, edges = np.histogram(np.asarray(s, dtype=np.float64), bins=nbins)
    mids = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * mids) / np.maximum(w0, 1)
    m1 = (np.sum(hist * mids) - np.cumsum(hist * mids)) / np.maximum(w1, 1)
    return float(mids[int(np.argmax(w0 * w1 * (m0 - m1) ** 2))])


def d1_f1_otsu(y, s):
    """F1 at the Otsu threshold."""
    return d1_f1_calls(y, np.asarray(s) > d1_otsu_threshold(s))


def d1_spatial_folds(coords, n_blocks=25, n_folds=5, seed=0):
    """
    Spatial folds: k-means blocks of spots, blocks assigned to folds at random.

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates.
    n_blocks, n_folds : int
        Number of k-means blocks and folds (n_folds = n_blocks gives one id per block).
    seed : int
        Seed of k-means and of the block-to-fold assignment.

    Returns
    -------
    np.ndarray
        Fold of every spot.
    """
    km = KMeans(n_clusters=n_blocks, n_init=10, random_state=seed).fit(
        np.asarray(coords)
    )
    rng = np.random.default_rng(seed)
    return rng.permutation(np.arange(n_blocks) % n_folds)[km.labels_]


def d1_f1_spatial_cv(y, s, folds):
    """
    F1 with the F1-optimal threshold learnt on the other folds and applied to each held-out fold (calls pooled).

    Parameters
    ----------
    y, s : np.ndarray
        Labels and scores.
    folds : np.ndarray
        Fold of every spot.

    Returns
    -------
    float
        F1 of the pooled held-out calls.
    """
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=np.float64)
    pred = np.zeros(len(y), bool)
    for f in np.unique(folds):
        test = folds == f
        if y[~test].sum() == 0:
            continue
        _, thr = _optimal_f1(y[~test], s[~test])
        pred[test] = s[test] >= thr
    return d1_f1_calls(y, pred)


def d1_threshold_rules(y, s, partitions):
    """
    AUROC, AUPRC and F1 under the five shared threshold rules.

    Parameters
    ----------
    y, s : np.ndarray
        Labels and scores.
    partitions : list of np.ndarray
        Spatial fold partitions; the spatial-block CV F1 is the mean over them.

    Returns
    -------
    dict
        Metric name -> value.
    """
    return {
        "AUROC": roc_auc_score(y, s),
        "AUPRC": average_precision_score(y, s),
        "F1 (oracle threshold)": _optimal_f1(y, s)[0],
        "F1 (spatial-block CV threshold)": float(
            np.mean([d1_f1_spatial_cv(y, s, f) for f in partitions])
        ),
        "F1 (top-k, prevalence-matched)": d1_f1_topk(y, s),
        "F1 (Otsu, label-free)": d1_f1_otsu(y, s),
        "F1 (median threshold)": d1_f1_median(y, s),
    }


def d1_spatial_graph(coords, n_neighbors=6):
    """
    Binary spatial kNN graph as built inside em.tl.score.

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates.
    n_neighbors : int
        Number of neighbours.

    Returns
    -------
    scipy.sparse.csr_matrix
        Spatial connectivities.
    """
    tmp = anndata.AnnData(X=np.zeros((len(coords), 1), dtype=np.float32))
    tmp.obsm["spatial"] = np.asarray(coords)
    sq.gr.spatial_neighbors(
        tmp, n_neighs=n_neighbors, coord_type="generic", key_added="spatial"
    )
    return csr_matrix(tmp.obsp["spatial_connectivities"])


def d1_fixed_rule_calls(null_z, p_empirical, gsd_positive, gsd_rule):
    """
    Binary calls of the fixed rules (the methods' own thresholds).

    Parameters
    ----------
    null_z : np.ndarray
        Null-model z score of every spot.
    p_empirical : np.ndarray
        Null-model empirical p of every spot.
    gsd_positive : np.ndarray
        GSDensity antimode label of every spot (True = positive).
    gsd_rule : str
        Name of the GSDensity rule.

    Returns
    -------
    dict
        (method, rule) -> boolean calls of every spot (index them for a bootstrap resample).
    """
    return {
        ("EnrichMap (null)", "z > 0"): null_z > 0,
        ("EnrichMap (null)", "z > 1.96"): null_z > 1.96,
        ("EnrichMap (null)", "p_empirical <= 0.05"): p_empirical <= 0.05,
        ("GSDensity", gsd_rule): np.asarray(gsd_positive),
    }


def d1_paired_block_bootstrap(
    y, method_scores, fixed_calls, blocks, metric_fns, n_boot=1000, seed=7, min_pos=3
):
    """
    Paired spatial block bootstrap of every method and fixed rule.

    Parameters
    ----------
    y : array-like
        Binary labels of the task.
    method_scores : dict
        Method -> score of every spot.
    fixed_calls : dict
        (method, rule) -> boolean calls of every spot (d1_fixed_rule_calls).
    blocks : np.ndarray
        Spatial block of every spot (the resampling unit).
    metric_fns : dict
        Metric name -> function(y, s).
    n_boot : int
        Number of replicates.
    seed : int
        Seed of the block resampling.
    min_pos : int
        Replicates with fewer region spots are skipped.

    Returns
    -------
    pd.DataFrame
        replicate, method, metric, value.
    """
    rng = np.random.default_rng(seed)
    y = np.asarray(y).astype(int)
    ub = np.unique(blocks)
    idx_by_block = {b: np.flatnonzero(blocks == b) for b in ub}
    records = []
    for r in range(n_boot):
        pick = rng.choice(ub, size=len(ub), replace=True)
        idx = np.concatenate([idx_by_block[b] for b in pick])
        yy = y[idx]
        if yy.sum() < min_pos or yy.sum() == len(yy):
            continue
        for m, s in method_scores.items():
            ss = s[idx]
            for name, fn in metric_fns.items():
                records.append((r, m, name, float(fn(yy, ss))))
        for (m, rule), calls in fixed_calls.items():
            records.append((r, m, f"F1 ({rule})", d1_f1_calls(yy, calls[idx])))
    return pd.DataFrame(records, columns=["replicate", "method", "metric", "value"])


def d1_paired_differences(piv, wide, metrics, methods, reference):
    """
    Paired bootstrap differences of every other method to a reference method.

    Parameters
    ----------
    piv : pd.DataFrame
        Bootstrap values, index = (task, metric, replicate), columns = methods.
    wide : pd.DataFrame
        Full-data values, index = (task, method), columns = metrics.
    metrics : list of str
        Metrics to compare.
    methods : list of str
        Methods.
    reference : str
        Reference method.

    Returns
    -------
    pd.DataFrame
        One row per task, metric and method: point difference, 95% percentile interval, share of replicates with the
        same sign.
    """
    out = []
    for (ct, met), d in piv.groupby(level=[0, 1]):
        if met not in metrics:
            continue
        for m in methods:
            if m == reference:
                continue
            dd = d[m] - d[reference]
            point = wide.loc[(ct, m), met] - wide.loc[(ct, reference), met]
            lo, hi = np.quantile(dd, [0.025, 0.975])
            out.append(
                {
                    "task": ct,
                    "metric": met,
                    "method": m,
                    "difference": point,
                    "lo": lo,
                    "hi": hi,
                    "CI excludes 0": (lo > 0) or (hi < 0),
                    "share same sign": float(np.mean(np.sign(dd) == np.sign(point))),
                }
            )
    return pd.DataFrame(out)


def d1_plot_metric_pair(table, x, y, ax, colours, ref_x=None, ref_y=None):
    """
    One metric against another, one filled circle per method, with whiskers from f"{x}_lo"/f"{x}_hi" columns.

    Parameters
    ----------
    table : pd.DataFrame
        One row per method with a "method" column and the metric columns.
    x, y : str
        Metric columns (e.g. "AUROC" and "AUPRC").
    ax : matplotlib.axes.Axes
        Target axes.
    colours : dict
        Method -> colour.
    ref_x, ref_y : float or None
        Dotted reference lines.
    """
    for _, r in table.iterrows():
        xerr = (
            [[r[x] - r[f"{x}_lo"]], [r[f"{x}_hi"] - r[x]]]
            if f"{x}_lo" in table
            else None
        )
        yerr = (
            [[r[y] - r[f"{y}_lo"]], [r[f"{y}_hi"] - r[y]]]
            if f"{y}_lo" in table
            else None
        )
        ax.errorbar(
            r[x],
            r[y],
            xerr=xerr,
            yerr=yerr,
            fmt="o",
            color=colours[r["method"]],
            ms=5,
            mec="white",
            mew=0.4,
            elinewidth=0.6,
            capsize=0,
        )
    if ref_x is not None:
        ax.axvline(ref_x, ls=":", c="0.6", lw=0.6)
    if ref_y is not None:
        ax.axhline(ref_y, ls=":", c="0.6", lw=0.6)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.grid(False)


def d1_plot_benchmark_summary(table, ax, colours):
    """
    Methods x metrics summary: filled circles in the method colour, sized by each metric scaled 0-1 across methods,
    values printed, mean rank on the right.

    Parameters
    ----------
    table : pd.DataFrame
        Index = methods, columns = metrics (higher is better).
    ax : matplotlib.axes.Axes
        Target axes.
    colours : dict
        Method -> colour.
    """
    scaled = table.astype(float)
    scaled = (scaled - scaled.min()) / (scaled.max() - scaled.min()).replace(0, np.nan)
    scaled = scaled.fillna(1.0)
    ranks = scaled.rank(ascending=False).mean(axis=1)
    for i, m in enumerate(table.index):
        face = colours[m]
        dark = np.dot(to_rgb(face), [0.299, 0.587, 0.114]) < 0.55
        for j, c in enumerate(table.columns):
            s = scaled.loc[m, c]
            ax.scatter(j, i, s=30 + 220 * s, marker="o", color=face, edgecolors="none")
            ax.text(
                j,
                i,
                f"{table.loc[m, c]:.2g}",
                ha="center",
                va="center",
                fontsize=5,
                color="white" if (s > 0.6 and dark) else "black",
            )
        ax.text(
            len(table.columns) + 0.2,
            i,
            f"{ranks[m]:.1f}",
            ha="left",
            va="center",
            fontsize=6,
        )
    ax.set_xticks(range(len(table.columns)))
    ax.set_xticklabels(table.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(table.index)))
    ax.set_yticklabels(table.index)
    for lab in ax.get_yticklabels():
        lab.set_color(colours[lab.get_text()])
    ax.text(len(table.columns) + 0.2, -0.9, "mean rank", fontsize=6, ha="left")
    ax.set_xlim(-0.6, len(table.columns) + 0.9)
    ax.set_ylim(len(table.index) - 0.4, -1.2)
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


def d1_dot_ci(ax, k, v, lo, hi, c, ms=3.2, orient="h"):
    """
    One estimate (filled circle) with its interval.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    k : float
        Position on the categorical axis.
    v, lo, hi : float
        Estimate and interval limits.
    c : str
        Colour.
    ms : float
        Marker size.
    orient : str
        "h" (value on x) or "v" (value on y).
    """
    if orient == "h":
        ax.plot([lo, hi], [k, k], color=c, lw=0.9, solid_capstyle="butt")
        ax.plot(v, k, "o", ms=ms, color=c, mec="white", mew=0.3, zorder=3)
    else:
        ax.plot([k, k], [lo, hi], color=c, lw=0.9, solid_capstyle="butt")
        ax.plot(k, v, "o", ms=ms, color=c, mec="white", mew=0.3, zorder=3)


def d1_format_axis(ax, met, n_rows):
    """
    x range, ticks and grid of one metric axes (rows = methods or rules).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    met : str
        Metric name.
    n_rows : int
        Number of rows.
    """
    if met == "AUROC":
        ax.set_xlim(0.5, 1.02)
        ax.set_xticks([0.5, 0.75, 1.0])
        ax.set_xticklabels(["0.5", "", "1"])
    else:
        ax.set_xlim(0, 1.04)
        ax.set_xticks([0, 0.5, 1.0])
        ax.set_xticklabels(["0", "0.5", "1"])
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.tick_params(labelsize=6)
    ax.grid(False)
    ax.grid(axis="x", lw=0.3, color="0.88")
    ax.set_axisbelow(True)
    sns.despine(ax=ax)


def d1_fixed_rules_axis(ax, table, colours, fixed_labels):
    """
    F1 of the fixed rules with intervals; a vertical tick marks the oracle F1 of the same method.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    table : pd.DataFrame
        One row per fixed rule with method, rule, estimate, lo, hi and oracle.
    colours : dict
        Method -> colour.
    fixed_labels : dict
        Rule -> axis label.
    """
    for k, r in enumerate(table.itertuples()):
        c = colours[r.method]
        d1_dot_ci(ax, k, r.estimate, r.lo, r.hi, c)
        ax.plot([r.oracle, r.oracle], [k - 0.22, k + 0.22], color=c, lw=0.8)
    d1_format_axis(ax, "F1", len(table))
    ax.set_yticks(range(len(table)))
    ax.set_yticklabels([fixed_labels[r] for r in table["rule"]], fontsize=6)
    for lab, m in zip(ax.get_yticklabels(), table["method"]):
        lab.set_color(colours[m])


def d1_legend_handles(names, colours):
    """Legend handles, one filled circle per method."""
    return [
        Line2D([], [], marker="o", ls="", color=colours[m], ms=4, label=m)
        for m in names
    ]


def d1_draw_brain_panel(
    ax, summary, wide, prevalence, methods, colours, ct, met, show_y
):
    """
    One task x metric panel of the full Sup. Fig. 2e replacement.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    summary : pd.DataFrame
        task, method, metric, value, lo, hi, n_rep (estimates with bootstrap intervals).
    wide : pd.DataFrame
        Full-data values, index = (task, method); holds "F1 (spatial-block CV) s.d. over partitions".
    prevalence : float
        Fraction of region spots of the task.
    methods : list of str
        Methods, in plotting order.
    colours : dict
        Method -> colour.
    ct : str
        Task (region).
    met : str
        Metric name.
    show_y : bool
        Whether to label the method axis.

    Returns
    -------
    pd.DataFrame
        The plotted values, one row per method.
    """
    d = (
        summary[(summary.task == ct) & (summary.metric == met)]
        .set_index("method")
        .reindex(methods)
    )
    p = prevalence
    plotted = []
    for k, m in enumerate(methods):
        v = d.loc[m, "value"]
        if met == "F1 (spatial-block CV threshold)":
            spread = wide.loc[(ct, m), "F1 (spatial-block CV) s.d. over partitions"]
            lo, hi = v - spread, v + spread
            bar, n = "estimate ± 1 s.d. over 20 block partitions", 20
        else:
            lo, hi = d.loc[m, "lo"], d.loc[m, "hi"]
            bar, n = (
                "95% percentile interval, paired block bootstrap",
                int(d.loc[m, "n_rep"]),
            )
        d1_dot_ci(ax, k, v, lo, hi, colours[m])
        plotted.append(
            {
                "task": ct,
                "metric": met,
                "method": m,
                "estimate": v,
                "bar_lo": lo,
                "bar_hi": hi,
                "error bar": bar,
                "n": n,
            }
        )
    if met == "F1 (median threshold)":
        ax.axvline(median_f1_ceiling(p), ls="--", c="0.35", lw=0.6)
    if met == "AUPRC":
        ax.axvline(p, ls=":", c="0.35", lw=0.6)
    d1_format_axis(ax, met, len(methods))
    ax.set_yticks(range(len(methods)))
    ax.set_yticklabels(methods if show_y else [])
    for lab in ax.get_yticklabels() if show_y else []:
        lab.set_color(colours[lab.get_text()])
    if not show_y:
        ax.tick_params(axis="y", length=0)
    return pd.DataFrame(plotted)


def d1_dodge_plot(
    ax, sim_all, methods, colours, rng, metrics, xlabels, ylim=(0, 1.02), centre="mean"
):
    """
    Per-slide dots with a centre and interval per method, dodged within each metric.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    sim_all : pd.DataFrame
        seed, method, metric, value of every simulated slide.
    methods : list of str
        Methods, in plotting order.
    colours : dict
        Method -> colour.
    rng : np.random.Generator
        Generator of the horizontal jitter.
    metrics : list of str
        Metrics on the x axis (values of sim_all.metric).
    xlabels : list of str
        Tick labels.
    ylim : tuple of float
        y-axis limits.
    centre : str
        "mean" (mean ± 1.96 s.e.m.) or "median" (median and IQR).

    Returns
    -------
    pd.DataFrame
        One row per slide and method, with the plotted centre and interval of its group.
    """
    w = 0.8 / len(methods)
    plotted = []
    for i, met in enumerate(metrics):
        for k, m in enumerate(methods):
            x0 = i - 0.4 + w * (k + 0.5)
            g = sim_all[(sim_all.method == m) & (sim_all.metric == met)].sort_values(
                "seed"
            )
            v = g["value"].values
            ax.scatter(
                x0 + rng.uniform(-w * 0.3, w * 0.3, len(v)),
                v,
                s=1.0,
                color=colours[m],
                alpha=0.18,
                lw=0,
                zorder=1,
            )
            if centre == "mean":
                c = v.mean()
                half = 1.96 * v.std(ddof=1) / np.sqrt(len(v))
                lo, hi = c - half, c + half
            else:
                c = np.median(v)
                lo, hi = np.quantile(v, [0.25, 0.75])
            d1_dot_ci(ax, x0, c, lo, hi, colours[m], ms=2.6, orient="v")
            plotted.append(
                pd.DataFrame(
                    {
                        "metric": met,
                        "method": m,
                        "seed": g["seed"].values,
                        "slide value": v,
                        "centre": c,
                        "bar_lo": lo,
                        "bar_hi": hi,
                    }
                )
            )
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(xlabels)
    ax.set_ylim(*ylim)
    ax.grid(False)
    ax.grid(axis="y", lw=0.3, color="0.88")
    ax.set_axisbelow(True)
    sns.despine(ax=ax)
    return pd.concat(plotted, ignore_index=True)


def d1_smooth(values, W):
    """
    Mean of the spatial neighbours (the smoothing step of EnrichMap).

    Parameters
    ----------
    values : np.ndarray
        Score per spot.
    W : scipy.sparse matrix
        Binary spatial connectivities.

    Returns
    -------
    np.ndarray
        Neighbourhood mean per spot.
    """
    d = np.maximum(np.asarray(W.sum(axis=1)).ravel(), 1e-10)
    return W.dot(np.asarray(values, dtype=np.float64)) / d


def d1_derive_signature(adata_train, ct, n=10):
    """
    Top up-regulated genes of a region vs the rest of the training spots.

    Parameters
    ----------
    adata_train : AnnData
        Training spots.
    ct : str
        Region.
    n : int
        Signature size.

    Returns
    -------
    list of str
        Genes with padj < 0.01 and log2FC > 1, ranked by log2FC.
    """
    a = anndata.AnnData(
        X=adata_train.X,
        obs=pd.DataFrame(
            {"grp": np.where(adata_train.obs["cluster"] == ct, "in", "rest")},
            index=adata_train.obs_names,
        ),
        var=pd.DataFrame(index=adata_train.var_names),
    )
    a.obs["grp"] = a.obs["grp"].astype("category")
    sc.tl.rank_genes_groups(
        a,
        groupby="grp",
        groups=["in"],
        reference="rest",
        method="wilcoxon",
        use_raw=False,
    )
    df = sc.get.rank_genes_groups_df(a, group="in")
    df = df[(df.pvals_adj < 0.01) & (df.logfoldchanges > 1)].sort_values(
        "logfoldchanges", ascending=False
    )
    return df.names.head(n).tolist()


# --------------------------------------------------------------------------------------------------
# R4_F1_GSDensity
# --------------------------------------------------------------------------------------------------


def read_gsdensity_task(task_dir, spots=None):
    """
    Read the outputs of one GSDensity run (scripts/run_gsdensity.R).

    Parameters
    ----------
    task_dir : str
        Folder with gsd_scores.csv, gsd_labels.csv, gsd_timing.json and, if the spatial test was run,
        gsd_spatial_kld.csv.
    spots : list or None
        Spot order to reindex the per-spot tables to.

    Returns
    -------
    dict
        "scores" (spots x sets relevance), "labels" (spots x sets antimode calls), "spatial_kld"
        (sets x statistics or None) and "timing" (seconds per step).
    """
    scores = pd.read_csv(os.path.join(task_dir, "gsd_scores.csv"), index_col=0)
    labels = pd.read_csv(os.path.join(task_dir, "gsd_labels.csv"), index_col=0)
    kld_path = os.path.join(task_dir, "gsd_spatial_kld.csv")
    kld = pd.read_csv(kld_path, index_col=0) if os.path.exists(kld_path) else None
    timing = json.load(open(os.path.join(task_dir, "gsd_timing.json")))
    if spots is not None:
        scores, labels = scores.reindex(spots), labels.reindex(spots)
    return {"scores": scores, "labels": labels, "spatial_kld": kld, "timing": timing}


def gsd_peak_rss_gb(log):
    """
    Peak resident memory reported by /usr/bin/time -v in a GSDensity log.

    Parameters
    ----------
    log : str
        Path of the log file.

    Returns
    -------
    float
        Peak RSS in GB (NaN if not reported).
    """
    for line in open(log):
        if "Maximum resident set size" in line:
            return int(line.split(":")[1]) / 1024**2
    return np.nan


def dense(X):
    """
    Dense numpy version of a (possibly sparse) matrix.

    Parameters
    ----------
    X : array-like or sparse matrix

    Returns
    -------
    np.ndarray
    """
    return X.toarray() if issparse(X) else np.asarray(X)


def raw_enrichmap(adata, genes, weights=None):
    """
    EnrichMap Eq. 1 (CV-weighted mean of log-expression) before z-scoring, smoothing and the null model.

    Parameters
    ----------
    adata : AnnData
        Log-normalised expression.
    genes : list
        Signature genes (those missing from adata are dropped).
    weights : dict or None
        Gene weights; default em.tl.infer_gene_weights.

    Returns
    -------
    np.ndarray
        Raw score per spot.
    """
    genes = [g for g in dict.fromkeys(genes) if g in adata.var_names]
    w = weights if weights is not None else em.tl.infer_gene_weights(adata, genes)
    wv = np.array([w.get(g, 1.0) for g in genes], dtype=np.float64)
    X = dense(adata[:, genes].X)
    return X @ wv / np.abs(wv).sum()


def raw_zscore(adata, genes):
    """
    Z score: sum of log-expression of the signature genes divided by sqrt(number of genes).

    Parameters
    ----------
    adata : AnnData
        Log-normalised expression.
    genes : list
        Signature genes.

    Returns
    -------
    np.ndarray
        Score per spot.
    """
    genes = [g for g in dict.fromkeys(genes) if g in adata.var_names]
    return dense(adata[:, genes].X).sum(axis=1) / np.sqrt(len(genes))


def raw_scanpy(adata, genes):
    """
    sc.tl.score_genes of one signature on a light copy of `adata`.

    Parameters
    ----------
    adata : AnnData
        Log-normalised expression.
    genes : list
        Signature genes.

    Returns
    -------
    np.ndarray
        Score per spot.
    """
    a = ad.AnnData(
        X=adata.X,
        obs=pd.DataFrame(index=adata.obs_names),
        var=pd.DataFrame(index=adata.var_names),
    )
    sc.tl.score_genes(a, gene_list=list(genes), score_name="_bench", use_raw=False)
    return a.obs["_bench"].values.astype(np.float64)


def aucell_multi(adata, sets):
    """
    decoupler AUCell for many gene sets in one call (AUCell ranks genes within each spot, so a set's
    score does not depend on the other sets).

    Parameters
    ----------
    adata : AnnData
        Log-normalised expression.
    sets : dict
        Set name -> genes.

    Returns
    -------
    pd.DataFrame
        Spots x sets.
    """
    a = ad.AnnData(
        X=adata.X,
        obs=pd.DataFrame(index=adata.obs_names),
        var=pd.DataFrame(index=adata.var_names),
    )
    net = pd.DataFrame(
        [(k, g) for k, v in sets.items() for g in v], columns=["source", "target"]
    )
    dc.mt.aucell(a, net=net, tmin=1)
    scores = dc.pp.get_obsm(adata=a, key="score_aucell").obsm["score_aucell"]
    return pd.DataFrame(
        {k: scores[k].values.astype(np.float64) for k in sets}, index=adata.obs_names
    )


def score_sets(adata, sets, methods=("EnrichMap", "z-score", "scanpy", "AUCell")):
    """
    Raw (unsmoothed) scores of many gene sets with the fast scorers.

    Parameters
    ----------
    adata : AnnData
        Log-normalised expression.
    sets : dict
        Set name -> genes.
    methods : tuple of str
        Any of "EnrichMap", "z-score", "scanpy", "AUCell".

    Returns
    -------
    dict
        Method -> DataFrame (spots x sets).
    """
    scorers = {"EnrichMap": raw_enrichmap, "z-score": raw_zscore, "scanpy": raw_scanpy}
    out = {}
    for m in methods:
        if m == "AUCell":
            out[m] = aucell_multi(adata, sets)
        else:
            out[m] = pd.DataFrame(
                {
                    k: np.asarray(scorers[m](adata, v), dtype=np.float64).ravel()
                    for k, v in sets.items()
                },
                index=adata.obs_names,
            )
    return out


def expression_frame(adata, upper=False):
    """
    Genes x spots expression table for gseapy.

    Parameters
    ----------
    adata : AnnData
        Log-normalised expression.
    upper : bool
        Upper-case gene names (ssGSEA).

    Returns
    -------
    pd.DataFrame
    """
    df = pd.DataFrame(dense(adata.X).T, index=adata.var_names, columns=adata.obs_names)
    if upper:
        df.index = df.index.str.upper()
    return df


def gseapy_multi(adata, sets, which, threads=1):
    """
    GSVA (ES) or ssGSEA (NES) of several gene sets in one gseapy call.

    Parameters
    ----------
    adata : AnnData
        Log-normalised expression.
    sets : dict
        Set name -> genes.
    which : str
        "GSVA" or "ssGSEA".
    threads : int
        gseapy threads.

    Returns
    -------
    pd.DataFrame
        Spots x sets.
    """
    if which == "GSVA":
        res = gp.gsva(
            data=expression_frame(adata),
            gene_sets={k: list(v) for k, v in sets.items()},
            outdir=None,
            min_size=3,
            threads=threads,
        )
        tab = res.res2d.pivot(index="Term", columns="Name", values="ES")
    else:
        res = gp.ssgsea(
            data=expression_frame(adata, upper=True),
            gene_sets={k: [g.upper() for g in v] for k, v in sets.items()},
            outdir=None,
            sample_norm_method="rank",
            no_plot=True,
            min_size=3,
            threads=threads,
        )
        tab = res.res2d.pivot(index="Term", columns="Name", values="NES")
    return pd.DataFrame(
        {k: tab.loc[k].T.astype(float).reindex(adata.obs_names).values for k in sets},
        index=adata.obs_names,
    )


def noisy_adata(adata, genes, level, rep, region_idx):
    """
    Gaussian noise with SD = level x gene SD added to the log-expression of the signature genes only
    (as benchmark_noise_robustness in 17_benchmarking, with one seed per region, level and repeat).

    Parameters
    ----------
    adata : AnnData
        Log-normalised expression with obs["cluster"] and obsm["spatial"].
    genes : list
        Signature genes to perturb.
    level : float
        Noise SD in units of each gene's SD.
    rep : int
        Repeat index.
    region_idx : int
        Index of the region signature (seed offset).

    Returns
    -------
    AnnData
        Perturbed copy (sparse float32 X, obs["cluster"], obsm["spatial"]).
    """
    rng = np.random.default_rng(
        20_000 + 1000 * region_idx + int(round(100 * level)) + 7 * rep
    )
    X = dense(adata.X).astype(np.float64)
    idx = [adata.var_names.get_loc(g) for g in genes]
    sd = X[:, idx].std(axis=0)
    X[:, idx] += rng.normal(0, 1, (X.shape[0], len(idx))) * (sd * level)[None, :]
    a = ad.AnnData(
        X=csr_matrix(X.astype(np.float32)),
        obs=adata.obs[["cluster"]].copy(),
        var=pd.DataFrame(index=adata.var_names),
    )
    a.obsm["spatial"] = np.asarray(adata.obsm["spatial"]).copy()
    return a


def spatial_conn(coords, n_neighbors=6):
    """
    Binary spatial kNN graph as built inside em.tl.score.

    Parameters
    ----------
    coords : array-like
        Spot coordinates (n x 2).
    n_neighbors : int
        Number of neighbours.

    Returns
    -------
    scipy.sparse.csr_matrix
        Binary connectivities.
    """
    tmp = ad.AnnData(X=np.zeros((len(coords), 1), dtype=np.float32))
    tmp.obsm["spatial"] = np.asarray(coords)
    sq.gr.spatial_neighbors(
        tmp, n_neighs=n_neighbors, coord_type="generic", key_added="spatial"
    )
    return csr_matrix(tmp.obsp["spatial_connectivities"])


def smooth(Y, conn, n_iter=1):
    """
    EnrichMap Eq. 3: neighbourhood average with the row-normalised graph.

    Parameters
    ----------
    Y : array-like
        Score vector or spots x k matrix.
    conn : sparse matrix
        Binary spatial graph.
    n_iter : int
        Number of smoothing passes.

    Returns
    -------
    np.ndarray
        Smoothed values.
    """
    Y = np.asarray(Y, dtype=np.float64)
    d = np.maximum(np.asarray(conn.sum(axis=1)).ravel(), 1e-10)
    out = Y
    for _ in range(n_iter):
        out = conn.dot(out) / (d[:, None] if out.ndim == 2 else d)
    return out


def spatial_folds(coords, n_blocks=25, n_folds=5, seed=0):
    """
    Spatial cross-validation folds: k-means blocks of spots assigned to folds.

    Parameters
    ----------
    coords : array-like
        Spot coordinates.
    n_blocks, n_folds : int
        Number of k-means blocks and folds.
    seed : int
        Random seed.

    Returns
    -------
    np.ndarray
        Fold per spot.
    """
    km = KMeans(n_clusters=n_blocks, n_init=10, random_state=seed).fit(
        np.asarray(coords)
    )
    rng = np.random.default_rng(seed)
    block_fold = rng.permutation(np.arange(n_blocks) % n_folds)
    return block_fold[km.labels_]


def hop_distance(conn, mask, max_hops=12):
    """
    Graph distance (neighbour rings) of every spot to the nearest spot in `mask`.

    Parameters
    ----------
    conn : sparse matrix
        Spatial graph.
    mask : array-like of bool
        Region spots.
    max_hops : int
        Largest ring.

    Returns
    -------
    np.ndarray
        Hops per spot (inf beyond max_hops).
    """
    mask = np.asarray(mask, bool)
    A = csr_matrix((conn + conn.T) > 0).astype(np.int8)
    dist = np.full(len(mask), np.inf)
    dist[mask] = 0
    frontier = mask.copy()
    seen = mask.copy()
    for h in range(1, max_hops + 1):
        nxt = (A.dot(frontier.astype(np.int8)) > 0) & ~seen
        dist[nxt] = h
        seen |= nxt
        frontier = nxt
        if not nxt.any():
            break
    return dist


def leakage_profile(score, mask, dist, max_ring=8, far_from=6):
    """
    Standardised score by ring outside a region, L(d) = (ring mean - far mean) / (inside mean - far mean).

    Parameters
    ----------
    score : array-like
        Score per spot.
    mask : array-like of bool
        Region spots.
    dist : np.ndarray
        Hop distance to the region.
    max_ring, far_from : int
        Last ring and the distance from which spots count as far background.

    Returns
    -------
    dict
        Ring -> L(d).
    """
    score = np.asarray(score, dtype=np.float64)
    inside = score[np.asarray(mask, bool)].mean()
    far = score[dist >= far_from].mean()
    span = inside - far
    prof = {}
    for d in range(1, max_ring + 1):
        sel = dist == d
        prof[d] = (
            float((score[sel].mean() - far) / span)
            if sel.any() and span != 0
            else np.nan
        )
    return prof


def leakage_summary(score, mask, dist, halo=(1, 2, 3)):
    """
    Scalar leakage metrics of a score around a region.

    Parameters
    ----------
    score : array-like
        Score per spot.
    mask : array-like of bool
        Region spots.
    dist : np.ndarray
        Hop distance to the region.
    halo : tuple of int
        Rings forming the halo.

    Returns
    -------
    dict
        Halo leakage index, region-vs-halo AUROC, top-k share in the halo and inside the region.
    """
    score = np.asarray(score, dtype=np.float64)
    mask = np.asarray(mask, bool)
    prof = leakage_profile(score, mask, dist)
    halo_sel = np.isin(dist, halo)
    y = np.r_[np.ones(mask.sum()), np.zeros(halo_sel.sum())]
    sv = np.r_[score[mask], score[halo_sel]]
    k = int(mask.sum())
    top = np.argsort(-score, kind="stable")[:k]
    return {
        "halo leakage index": float(np.nanmean([prof[d] for d in halo])),
        "region-vs-halo AUROC": float(roc_auc_score(y, sv)),
        "top-k share in halo": float(np.isin(dist[top], halo).mean()),
        "top-k share inside": float(mask[top].mean()),
    }


def morans_i(values, W):
    """
    Moran's I with binary spatial weights.

    Parameters
    ----------
    values : array-like
        Score per spot.
    W : sparse matrix
        Binary spatial weights.

    Returns
    -------
    float
    """
    v = np.asarray(values, dtype=np.float64)
    z = v - v.mean()
    z2 = z @ z
    if z2 == 0:
        return 0.0
    return float((len(v) / W.sum()) * (z @ W.dot(z)) / z2)


def morans_i_cols(V, W):
    """
    Moran's I of every column of a spots x k matrix.

    Parameters
    ----------
    V : array-like
        Spots x k values.
    W : sparse matrix
        Binary spatial weights.

    Returns
    -------
    np.ndarray
        One Moran's I per column.
    """
    V = np.asarray(V, dtype=np.float64)
    Z = V - V.mean(axis=0)
    num = (Z * W.dot(Z)).sum(axis=0)
    den = (Z * Z).sum(axis=0)
    return (V.shape[0] / W.sum()) * num / den


def gearys_c(values, W):
    """
    Geary's C with binary spatial weights (< 1 = positive autocorrelation).

    Parameters
    ----------
    values : array-like
        Score per spot.
    W : sparse matrix
        Binary spatial weights.

    Returns
    -------
    float
        (n - 1) sum_ij w_ij (x_i - x_j)^2 / (2 W sum_i (x_i - mean)^2).
    """
    z = np.asarray(values, dtype=np.float64)
    z = z - z.mean()
    z2 = z @ z
    if z2 == 0:
        return np.nan
    W = W.tocsr()
    deg = np.asarray(W.sum(axis=1)).ravel() + np.asarray(W.sum(axis=0)).ravel()
    num = (z**2 * deg).sum() - 2 * z @ W.dot(z)
    return float((len(z) - 1) * num / (2 * W.sum() * z2))


def perm_matrix(n, B, seed=0):
    """
    Spot permutations.

    Parameters
    ----------
    n : int
        Number of spots.
    B : int
        Number of permutations.
    seed : int
        Random seed.

    Returns
    -------
    np.ndarray
        n x B index matrix.
    """
    rng = np.random.default_rng(seed)
    return np.stack([rng.permutation(n) for _ in range(B)], axis=1)


def emp_p(obs, null):
    """
    Empirical upper-tail p-value, (1 + #null >= obs) / (1 + #null).

    Parameters
    ----------
    obs : float
        Observed statistic.
    null : array-like
        Null statistics.

    Returns
    -------
    float
    """
    null = np.asarray(null)
    return float((1 + (null >= obs).sum()) / (1 + len(null)))


def competitive_p(obs, null):
    """
    Competitive p-value of a statistic against matched control sets.

    Parameters
    ----------
    obs : float
        Statistic of the named signature.
    null : array-like
        Statistic of the matched control sets (non-finite values dropped).

    Returns
    -------
    float
        (1 + #controls >= obs) / (1 + #controls), NaN if obs is undefined.
    """
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    return (
        np.nan if not np.isfinite(obs) else (1 + (null >= obs).sum()) / (1 + len(null))
    )


def to_weights(s):
    """
    Shift a signed score to be non-negative and normalise it to sum 1 (lowest spot gets weight 0).

    Parameters
    ----------
    s : array-like
        Score vector or spots x k matrix.

    Returns
    -------
    np.ndarray
    """
    s = np.asarray(s, dtype=np.float64)
    w = s - s.min(axis=0)
    return w / w.sum(axis=0)


def _bandwidth_nrd(x):
    """
    MASS::bandwidth.nrd, 4 * 1.06 * min(sd, IQR / 1.34) * n^(-1/5).

    Parameters
    ----------
    x : array-like
        Coordinate values.

    Returns
    -------
    float
    """
    x = np.asarray(x, dtype=np.float64)
    q = np.quantile(x, [0.25, 0.75])
    return 4 * 1.06 * min(np.std(x, ddof=1), (q[1] - q[0]) / 1.34) * len(x) ** (-0.2)


class SpatialKLD:
    """
    Python port of gsdensity::compute.spatial.kld: KL divergence between the weighted and the unweighted
    2-D KDE of spot positions on an n x n grid (bandwidth.nrd / 4 per axis).

    Parameters
    ----------
    coords : array-like
        Spot coordinates (N x 2).
    n : int
        Grid size per axis.
    """

    def __init__(self, coords, n=10):
        x, y = (
            np.asarray(coords, dtype=np.float64)[:, 0],
            np.asarray(coords, dtype=np.float64)[:, 1],
        )
        gx, gy = np.linspace(x.min(), x.max(), n), np.linspace(y.min(), y.max(), n)
        h = np.array([_bandwidth_nrd(x), _bandwidth_nrd(y)]) / 4
        self.ax = norm.pdf((gx[:, None] - x[None, :]) / h[0])
        self.ay = norm.pdf((gy[:, None] - y[None, :]) / h[1])
        self.h = h
        self.N = len(x)
        self.Q = self.dens(np.full(self.N, 1.0 / self.N)) + 1e-300

    def dens(self, W):
        """
        Weighted KDE on the grid.

        Parameters
        ----------
        W : array-like
            Weights per spot, or spots x k.

        Returns
        -------
        np.ndarray
            Grid densities (n^2, or n^2 x k).
        """
        W = np.asarray(W, dtype=np.float64)
        one = W.ndim == 1
        W = W[:, None] if one else W
        Z = np.einsum("ik,jk,kc->ijc", self.ax, self.ay, W, optimize=True) / (
            W.sum(axis=0) * self.h[0] * self.h[1]
        )
        Z = Z.reshape(-1, W.shape[1])
        return Z[:, 0] if one else Z

    def kld(self, W):
        """
        Unnormalised KL divergence sum P log(P / Q) of the weighted vs the uniform KDE.

        Parameters
        ----------
        W : array-like
            Weights per spot, or spots x k.

        Returns
        -------
        float or np.ndarray
        """
        P = self.dens(W)
        Q = self.Q if P.ndim == 1 else self.Q[:, None]
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(P > 0, P * np.log(P / Q), 0.0)
        return t.sum(axis=0)

    def test(self, w, null_W):
        """
        GSDensity's statistic (log KLD), its Gaussian p-value against a null, and an empirical p-value.

        Parameters
        ----------
        w : array-like
            Observed weights per spot.
        null_W : array-like
            Spots x B null weight vectors (undefined null draws are dropped and counted).

        Returns
        -------
        dict
            spatial_kld, rspatial_kld (null mean), rspatial_kld_sd, delta_kld, p_gauss, p_emp, n_null_nan.
        """
        with np.errstate(invalid="ignore", divide="ignore"):
            obs = np.log(self.kld(w))
            nl = np.log(self.kld(null_W))
        n_nan = int((~np.isfinite(nl)).sum())
        nl = nl[np.isfinite(nl)]
        mu, sd = nl.mean(), nl.std(ddof=1)
        ok = np.isfinite(obs)
        return dict(
            spatial_kld=obs,
            rspatial_kld=mu,
            rspatial_kld_sd=sd,
            delta_kld=obs - mu,
            p_gauss=float(norm.sf(obs, mu, sd)) if ok else np.nan,
            p_emp=float((1 + (nl >= obs).sum()) / (1 + len(nl))) if ok else np.nan,
            n_null_nan=n_nan,
        )


def f1_from_calls(y, pred):
    """
    F1 of binary calls.

    Parameters
    ----------
    y : array-like
        True labels (0/1).
    pred : array-like
        Calls (bool or 0/1).

    Returns
    -------
    float
    """
    y = np.asarray(y).astype(bool)
    pred = np.asarray(pred).astype(bool)
    tp = (y & pred).sum()
    denom = pred.sum() + y.sum()
    return float(2 * tp / denom) if denom else 0.0


def f1_oracle(y, s):
    """
    F1 at the F1-maximising threshold chosen on the evaluated spots.

    Parameters
    ----------
    y : array-like
        True labels.
    s : array-like
        Score.

    Returns
    -------
    tuple
        (F1, threshold).
    """
    p, r, t = precision_recall_curve(y, s)
    p, r = p[:-1], r[:-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.nan_to_num(2 * p * r / (p + r))
    i = int(np.argmax(f))
    return float(f[i]), float(t[i])


def f1_median(y, s):
    """
    F1 with the score thresholded at its median.

    Parameters
    ----------
    y : array-like
        True labels.
    s : array-like
        Score.

    Returns
    -------
    float
    """
    return f1_from_calls(y, s > np.median(s))


def f1_topk(y, s):
    """
    F1 when the k highest-scoring spots are called positive, k = number of positives.

    Parameters
    ----------
    y : array-like
        True labels.
    s : array-like
        Score.

    Returns
    -------
    float
    """
    y = np.asarray(y).astype(bool)
    k = int(y.sum())
    order = np.argsort(-np.asarray(s), kind="stable")
    pred = np.zeros(len(y), bool)
    pred[order[:k]] = True
    return f1_from_calls(y, pred)


def otsu_threshold(s, nbins=256):
    """
    Otsu threshold of a score.

    Parameters
    ----------
    s : array-like
        Score.
    nbins : int
        Histogram bins.

    Returns
    -------
    float
    """
    s = np.asarray(s, dtype=np.float64)
    hist, edges = np.histogram(s, bins=nbins)
    mids = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * mids) / np.maximum(w0, 1)
    m1 = (np.sum(hist * mids) - np.cumsum(hist * mids)) / np.maximum(w1, 1)
    between = w0 * w1 * (m0 - m1) ** 2
    return float(mids[int(np.argmax(between))])


def f1_otsu(y, s):
    """
    F1 at the label-free Otsu threshold.

    Parameters
    ----------
    y : array-like
        True labels.
    s : array-like
        Score.

    Returns
    -------
    float
    """
    return f1_from_calls(y, np.asarray(s) > otsu_threshold(s))


def f1_spatial_cv(y, s, folds):
    """
    F1 with the F1-optimal threshold learnt on training spatial blocks and applied to held-out blocks.

    Parameters
    ----------
    y : array-like
        True labels.
    s : array-like
        Score.
    folds : np.ndarray
        Fold per spot.

    Returns
    -------
    float
    """
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=np.float64)
    pred = np.zeros(len(y), bool)
    for f in np.unique(folds):
        test = folds == f
        train = ~test
        if y[train].sum() == 0:
            continue
        _, thr = f1_oracle(y[train], s[train])
        pred[test] = s[test] >= thr
    return f1_from_calls(y, pred)


def spatial_coherence(pred, W, gt_morans):
    """
    Moran's I of a binary prediction relative to that of the truth (as in 17_benchmarking).

    Parameters
    ----------
    pred : array-like
        Binary calls.
    W : sparse matrix
        Spatial weights.
    gt_morans : float
        Moran's I of the true labels.

    Returns
    -------
    float
    """
    mi = morans_i(np.asarray(pred, dtype=float), W)
    return float(np.clip(mi / gt_morans, 0, 1)) if gt_morans > 0 else max(mi, 0.0)


def classification_metrics(y, s, W=None, folds=None):
    """
    Per-spot classification metrics of one score for one binary task.

    Parameters
    ----------
    y : array-like
        True labels (0/1).
    s : array-like
        Score.
    W : sparse matrix or None
        Spatial weights (adds spatial coherence and Moran's I).
    folds : np.ndarray or None
        Spatial folds (adds the cross-validated threshold F1).

    Returns
    -------
    dict
        AUROC, AUPRC and F1 under the oracle, median, top-k, Otsu (and spatial-CV) thresholds.
    """
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=np.float64)
    if np.isnan(s).any() or np.nanstd(s) == 0:
        return {}
    f1o, thr = f1_oracle(y, s)
    out = {
        "AUROC": roc_auc_score(y, s),
        "AUPRC": average_precision_score(y, s),
        "F1 (oracle threshold)": f1o,
        "F1 (median threshold)": f1_median(y, s),
        "F1 (top-k, prevalence-matched)": f1_topk(y, s),
        "F1 (Otsu, label-free)": f1_otsu(y, s),
    }
    if folds is not None:
        out["F1 (spatial-block CV threshold)"] = f1_spatial_cv(y, s, folds)
    if W is not None:
        gt = morans_i(y.astype(float), W)
        out["Spatial coherence (oracle thr)"] = spatial_coherence(s >= thr, W, gt)
        out["Moran's I"] = morans_i(s, W)
    return out


def median_f1_upper_bound(prevalence):
    """
    Largest F1 reachable when half of the spots are called positive, 4p / (1 + 2p).

    Parameters
    ----------
    prevalence : float
        Fraction of positive spots.

    Returns
    -------
    float
    """
    p = np.asarray(prevalence, dtype=np.float64)
    return 4 * p / (1 + 2 * p)


def block_bootstrap(y, s, blocks, metric_fns, B=500, seed=0):
    """
    Percentile 95% CIs from resampling spatial blocks with replacement.

    Parameters
    ----------
    y : array-like
        True labels.
    s : array-like
        Score (or 0/1 calls).
    blocks : np.ndarray
        Block per spot.
    metric_fns : dict
        Name -> function(y, s).
    B : int
        Number of resamples.
    seed : int
        Random seed (the same seed gives the same resamples for every method).

    Returns
    -------
    dict
        Name -> (2.5th, 97.5th percentile).
    """
    rng = np.random.default_rng(seed)
    ub = np.unique(blocks)
    members = [np.where(blocks == b)[0] for b in ub]
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=np.float64)
    out = {k: [] for k in metric_fns}
    for _ in range(B):
        idx = np.concatenate([members[i] for i in rng.integers(0, len(ub), len(ub))])
        yy, ss = y[idx], s[idx]
        if yy.sum() == 0 or yy.sum() == len(yy):
            continue
        for k, fn in metric_fns.items():
            out[k].append(fn(yy, ss))
    return {
        k: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
        for k, v in out.items()
    }


def robustness_arm_scores(raw_scores, z_null, arms, conn):
    """
    Final score of every robustness arm: a raw score passed through one spatial operator, or the
    null-model score.

    Parameters
    ----------
    raw_scores : dict
        Raw-score source ("EnrichMap", "z-score", "scanpy", "AUCell", "GSDensity") -> raw vector.
    z_null : array-like
        Null-model score of the same (perturbed) signature and matrix.
    arms : dict
        Arm name -> (raw-score source, operator), operator "raw", "smooth" or "null model".
    conn : sparse matrix
        Spatial graph for smoothing.

    Returns
    -------
    dict
        Arm name -> final score (GSDensity relevance is log(W + 1e-6) before smoothing).
    """
    out = {}
    for name, (source, op) in arms.items():
        if op == "null model":
            out[name] = np.asarray(z_null, float)
            continue
        v = np.asarray(raw_scores[source], float)
        if op == "raw":
            out[name] = v
            continue
        if source == "GSDensity":
            v = np.log(v + 1e-6)
        out[name] = smooth(v, conn)
    return out


def robustness_area_table(rob, metric):
    """
    Normalised area of every robustness curve: per repeat the mean over the perturbation levels
    (level 0 excluded), then the mean and 95% t-interval over repeats.

    Parameters
    ----------
    rob : pd.DataFrame
        One row per region, experiment, level, repeat and method.
    metric : str
        "spearman_vs_clean" or "AUROC".

    Returns
    -------
    pd.DataFrame
        One row per region, experiment and method.
    """
    per_rep = (
        rob[rob.level > 0]
        .groupby(["region", "experiment", "method", "rep"])[metric]
        .mean()
        .rename("area")
        .reset_index()
    )
    out = (
        per_rep.groupby(["region", "experiment", "method"])["area"]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    half = t_dist.ppf(0.975, out["count"] - 1) * out["std"] / np.sqrt(out["count"])
    out["lo"], out["hi"] = out["mean"] - half, out["mean"] + half
    out["metric"] = metric
    return out


def time_python_scorers(adata, gene_sets, pools, enrichmap_variants, n_repeats=3):
    """
    Wall-clock seconds to score N gene sets with the EnrichMap variants (one em.tl.score call with all sets,
    6-NN graph built once) and the fast comparators, on a fresh copy of the expression each time.

    Parameters
    ----------
    adata : AnnData
        Log-normalised expression with obsm["spatial"].
    gene_sets : dict
        Set name -> genes.
    pools : dict
        N -> list of set names to score.
    enrichmap_variants : dict
        Method label -> em.tl.score keyword arguments.
    n_repeats : int
        Timing repeats.

    Returns
    -------
    pd.DataFrame
        method, n_sets, repeat, seconds, load_avg_1min.
    """
    rows = []
    for rep in range(n_repeats):
        for n, names in pools.items():
            sets = {k: gene_sets[k] for k in names}
            for method in list(enrichmap_variants) + ["Z score", "Scanpy", "AUCell"]:
                a = ad.AnnData(
                    X=adata.X,
                    obs=pd.DataFrame(index=adata.obs_names),
                    var=pd.DataFrame(index=adata.var_names),
                )
                a.obsm["spatial"] = np.asarray(adata.obsm["spatial"])
                t0 = time.time()
                if method in enrichmap_variants:
                    em.tl.score(a, gene_set=sets, **enrichmap_variants[method])
                elif method == "Z score":
                    [raw_zscore(a, v) for v in sets.values()]
                elif method == "Scanpy":
                    [raw_scanpy(a, v) for v in sets.values()]
                else:
                    aucell_multi(a, sets)
                rows.append(
                    dict(
                        method=method,
                        n_sets=n,
                        repeat=rep,
                        seconds=time.time() - t0,
                        load_avg_1min=os.getloadavg()[0],
                    )
                )
    return pd.DataFrame(rows)


def gsdensity_shades(base_colour, n):
    """
    `n` shades of a colour's hue (dark to light), for the GSDensity parameter settings.

    Parameters
    ----------
    base_colour : str
        Hex colour.
    n : int
        Number of shades.

    Returns
    -------
    list of str
        Hex colours.
    """
    hue, _, sat = colorsys.rgb_to_hls(
        *(int(base_colour[i : i + 2], 16) / 255 for i in (1, 3, 5))
    )
    return [
        "#%02x%02x%02x"
        % tuple(round(v * 255) for v in colorsys.hls_to_rgb(hue, lt, sat))
        for lt in np.linspace(0.30, 0.78, n)
    ]


def plot_metric_pair(table, x, y, ax, colours, ref_x=None, ref_y=None, invert_y=False):
    """
    One metric against another, one filled circle per method, with 95% CI whiskers when f"{x}_lo"/f"{x}_hi" and
    f"{y}_lo"/f"{y}_hi" columns exist.

    Parameters
    ----------
    table : pd.DataFrame
        One row per method (and task) with a "method" column and the metric columns.
    x, y : str
        Metric columns.
    ax : plt.Axes
        Target axes.
    colours : dict
        Method -> colour.
    ref_x, ref_y : float or None
        Dotted reference lines.
    invert_y : bool
        Invert the y axis (lower-is-better metrics).

    Returns
    -------
    plt.Axes
    """
    for _, r in table.iterrows():
        xerr = yerr = None
        if f"{x}_lo" in table:
            xerr = [[r[x] - r[f"{x}_lo"]], [r[f"{x}_hi"] - r[x]]]
        if f"{y}_lo" in table:
            yerr = [[r[y] - r[f"{y}_lo"]], [r[f"{y}_hi"] - r[y]]]
        ax.errorbar(
            r[x],
            r[y],
            xerr=xerr,
            yerr=yerr,
            fmt="o",
            color=colours[r["method"]],
            ms=5,
            mec="white",
            mew=0.4,
            elinewidth=0.6,
            capsize=0,
            label=r["method"],
        )
    if ref_x is not None:
        ax.axvline(ref_x, ls=":", c="0.6", lw=0.6)
    if ref_y is not None:
        ax.axhline(ref_y, ls=":", c="0.6", lw=0.6)
    if invert_y:
        ax.invert_yaxis()
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.grid(False)
    return ax


def plot_benchmark_summary(table, ax, colours, lower_is_better=(), annotate=True):
    """
    Methods x metrics summary: filled circles in the method colour, sized by each metric scaled 0-1 across
    methods (best = 1), raw values printed, and a mean-rank column.

    Parameters
    ----------
    table : pd.DataFrame
        Index = methods, columns = metrics.
    ax : plt.Axes
        Target axes.
    colours : dict
        Method -> colour of the circles and the row label (grey circles and black label if missing).
    lower_is_better : tuple of str
        Metrics where smaller is better.
    annotate : bool
        Print the raw value in each cell.

    Returns
    -------
    plt.Axes
    """
    scaled = table.copy().astype(float)
    for c in table.columns:
        col = table[c].astype(float)
        span = col.max() - col.min()
        s = (col - col.min()) / span if span > 0 else col * 0 + 1
        scaled[c] = 1 - s if c in lower_is_better else s
    ranks = scaled.rank(ascending=False).mean(axis=1)
    for i, m in enumerate(table.index):
        face = colours.get(m, "0.6")
        dark = np.dot(to_rgb(face), [0.299, 0.587, 0.114]) < 0.55
        for j, c in enumerate(table.columns):
            s = scaled.loc[m, c]
            ax.scatter(j, i, s=30 + 220 * s, marker="o", color=face, edgecolors="none")
            if annotate:
                ax.text(
                    j,
                    i,
                    f"{table.loc[m, c]:.2g}",
                    ha="center",
                    va="center",
                    fontsize=5,
                    color="white" if (s > 0.6 and dark) else "black",
                )
        ax.text(
            len(table.columns) + 0.2,
            i,
            f"{ranks[m]:.1f}",
            ha="left",
            va="center",
            fontsize=6,
        )
    ax.set_xticks(range(len(table.columns)))
    ax.set_xticklabels(table.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(table.index)))
    ax.set_yticklabels(table.index)
    for lab in ax.get_yticklabels():
        lab.set_color(colours.get(lab.get_text(), "black"))
    ax.text(len(table.columns) + 0.2, -0.9, "mean rank", fontsize=6, ha="left")
    ax.set_xlim(-0.6, len(table.columns) + 0.9)
    ax.set_ylim(len(table.index) - 0.4, -1.2)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    return ax


def label_colour(name, colours):
    """
    Method colour of a test or statistic label ("GSDensity ...", "EnrichMap (null) ...", ...).

    Parameters
    ----------
    name : str
        Label starting with the method.
    colours : dict
        Method -> colour.

    Returns
    -------
    str
        Colour (black if no method matches).
    """
    if name.startswith("GSDensity"):
        return colours["GSDensity"]
    for tag, m in (
        ("(null)", "EnrichMap (null)"),
        ("(smoothing-only)", "EnrichMap (smoothing-only)"),
        ("(raw)", "EnrichMap (raw, no smoothing)"),
    ):
        if tag in name:
            return colours[m]
    return "black"


def annotated_value_heatmap(ax, table, cmap, vmin, vmax, fmt, colours, neg_log10=False):
    """
    Annotated heat map with row labels in method colours.

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    table : pd.DataFrame
        Rows = tests or statistics, columns = groups or sets.
    cmap : str
        Colour map.
    vmin, vmax : float
        Colour limits (of the transformed values).
    fmt : str
        Format of the printed values.
    colours : dict
        Method -> colour (row labels).
    neg_log10 : bool
        Colour by -log10 of the values (p-values); vmin and vmax are then on that scale.
    """
    v = table.values.astype(float)
    shade = -np.log10(v) if neg_log10 else v
    ax.imshow(shade, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    for i in range(v.shape[0]):
        for j in range(v.shape[1]):
            txt = "–" if not np.isfinite(v[i, j]) else format(v[i, j], fmt)
            dark = (
                np.isfinite(shade[i, j]) and (shade[i, j] - vmin) / (vmax - vmin) > 0.6
            )
            ax.text(
                j,
                i,
                txt,
                ha="center",
                va="center",
                fontsize=5,
                color="white" if dark else "black",
            )
    ax.set_xticks(range(v.shape[1]))
    ax.set_xticklabels(table.columns, rotation=45, ha="right", fontsize=6)
    ax.set_yticks(range(v.shape[0]))
    ax.set_yticklabels(table.index, fontsize=5.5)
    for lab in ax.get_yticklabels():
        lab.set_color(label_colour(lab.get_text(), colours))
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)


def annotate_points_right(
    ax, x, values, labels, text_colours=None, dx=0.12, min_gap=0.07, fontsize=4.5
):
    """
    Short text labels to the right of points that share one x position, pushed apart vertically so that they do
    not overlap (call it after everything else is drawn, so that the y limits are final).

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    x : float
        x position of the points (data units).
    values : array-like
        y value of every point (NaN: no label).
    labels : list of str
        Label of every point.
    text_colours : list of str or None
        Colour of every label (default black).
    dx : float
        Horizontal offset of the labels (data units).
    min_gap : float
        Minimum vertical distance between labels (fraction of the axes height).
    fontsize : float
        Font size.
    """
    values = np.asarray(values, dtype=float)
    y0, y1 = ax.get_ylim()
    frac = (values - y0) / (y1 - y0)
    order = [i for i in np.argsort(frac) if np.isfinite(frac[i])]
    pos, last = {}, -np.inf
    for i in order:
        pos[i] = max(frac[i], last + min_gap)
        last = pos[i]
    shift = max(last - 1.0, 0.0) if pos else 0.0
    trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    for i, p in pos.items():
        ax.text(
            x + dx,
            p - shift,
            labels[i],
            transform=trans,
            ha="left",
            va="center",
            fontsize=fontsize,
            color="black" if text_colours is None else text_colours[i],
        )


# --------------------------------------------------------------------------------------------------
# R4_F2b_ablation_simulations
# --------------------------------------------------------------------------------------------------


def f2b_spatial_conn(coords, n_neighbors=6):
    """
    Binary spatial kNN connectivity, built as inside em.tl.score.

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates (n x 2).
    n_neighbors : int
        Number of spatial neighbours.

    Returns
    -------
    scipy.sparse.csr_matrix
        Spatial connectivities.
    """
    tmp = anndata.AnnData(X=np.zeros((len(coords), 1), dtype=np.float32))
    tmp.obsm["spatial"] = np.asarray(coords)
    sq.gr.spatial_neighbors(
        tmp, n_neighs=n_neighbors, coord_type="generic", key_added="spatial"
    )
    return csr_matrix(tmp.obsp["spatial_connectivities"])


def f2b_smooth(Y, conn, n_iter=1):
    """
    Neighbourhood mean with the row-normalised adjacency (EnrichMap Eq. 3), optionally iterated.

    Parameters
    ----------
    Y : np.ndarray
        Score per spot (or spots x columns).
    conn : scipy.sparse matrix
        Spatial connectivity.
    n_iter : int
        Number of smoothing iterations.

    Returns
    -------
    np.ndarray
        Smoothed score.
    """
    Y = np.asarray(Y, dtype=np.float64)
    d = np.maximum(np.asarray(conn.sum(axis=1)).ravel(), 1e-10)
    out = Y
    for _ in range(n_iter):
        out = conn.dot(out) / (d[:, None] if out.ndim == 2 else d)
    return out


def f2b_spatial_gam(coords):
    """
    The default spatial GAM of em.tl.score, LinearGAM(te(0, 1)), written as a linear smoother H = B (B'B + P)^-1 B'.

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates (n x 2).

    Returns
    -------
    dict
        Spline basis "B" and penalised normal matrix "A" = B'B + P (with pyGAM's sqrt(eps) ridge).
    """
    coords = np.asarray(coords, dtype=np.float64)
    gam = LinearGAM(te(0, 1)).fit(coords, coords[:, 0])
    B = gam._modelmat(coords)
    B = B.toarray() if hasattr(B, "toarray") else np.asarray(B)
    P = gam._P()
    P = P.toarray() if hasattr(P, "toarray") else np.asarray(P)
    return {
        "B": B,
        "A": B.T @ B + P + np.sqrt(np.finfo(float).eps) * np.eye(B.shape[1]),
    }


def f2b_gam_fitted(gam, Y):
    """
    GAM fit H y.

    Parameters
    ----------
    gam : dict
        Output of `f2b_spatial_gam`.
    Y : np.ndarray
        Score per spot.

    Returns
    -------
    np.ndarray
        Fitted values.
    """
    return gam["B"] @ np.linalg.solve(
        gam["A"], gam["B"].T @ np.asarray(Y, dtype=np.float64)
    )


def f2b_gam_residual(gam, Y):
    """
    GAM residual (I - H) y, the spatial-covariate correction of em.tl.score.

    Parameters
    ----------
    gam : dict
        Output of `f2b_spatial_gam`.
    Y : np.ndarray
        Score per spot.

    Returns
    -------
    np.ndarray
        Residual.
    """
    return np.asarray(Y, dtype=np.float64) - f2b_gam_fitted(gam, Y)


def f2b_configurations(raw, conn, gam):
    """
    The 2 x 2 factorial {smoothing} x {GAM} of one raw score.

    Parameters
    ----------
    raw : np.ndarray
        Raw score per spot.
    conn : scipy.sparse matrix
        6-NN connectivity.
    gam : dict
        Output of `f2b_spatial_gam`.

    Returns
    -------
    dict
        raw (r), smooth (S r), GAM ((I - H) r) and smooth+GAM ((I - H) S r).
    """
    raw = np.asarray(raw, dtype=np.float64)
    sm = f2b_smooth(raw, conn)
    return {
        "raw": raw,
        "smooth": sm,
        "GAM": f2b_gam_residual(gam, raw),
        "smooth+GAM": f2b_gam_residual(gam, sm),
    }


def f2b_spatial_folds(coords, n_blocks=25, n_folds=5, seed=0):
    """
    Fold of every spot: k-means spatial blocks assigned at random to folds.

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates.
    n_blocks, n_folds : int
        Number of spatial blocks and folds.
    seed : int
        Seed of the k-means and of the block-to-fold assignment.

    Returns
    -------
    np.ndarray
        Fold index per spot.
    """
    km = KMeans(n_clusters=n_blocks, n_init=10, random_state=seed).fit(
        np.asarray(coords)
    )
    return np.random.default_rng(seed).permutation(np.arange(n_blocks) % n_folds)[
        km.labels_
    ]


def f2b_boundary_mask(conn, y):
    """
    Spots with at least one neighbour of the other label (as run_classification in 05_01_simulated_benchmarking).

    Parameters
    ----------
    conn : scipy.sparse matrix
        6-NN connectivity.
    y : np.ndarray
        Binary labels.

    Returns
    -------
    np.ndarray
        Boolean mask of boundary spots.
    """
    A = csr_matrix(conn)
    pos_nb = np.asarray(A.dot(y.astype(float))).ravel()
    deg = np.asarray(A.sum(axis=1)).ravel()
    return ((y == 1) & (pos_nb < deg)) | ((y == 0) & (pos_nb > 0))


def f2b_perturbations(adata, sig_genes, seed, fractions, noise_levels, n_repeats):
    """
    Contaminated gene sets and noise draws of 05_01_simulated_benchmarking (stream default_rng(seed + 1000),
    contamination draws first, then noise draws).

    Parameters
    ----------
    adata : anndata.AnnData
        Simulated slide (log-normalised expression).
    sig_genes : list
        Signature genes.
    seed : int
        Simulation seed.
    fractions : list
        Random genes added, as multiples of the signature size.
    noise_levels : list
        Gaussian noise SD, as multiples of each signature gene's SD.
    n_repeats : int
        Repeats per fraction and per noise level.

    Returns
    -------
    tuple
        {(fraction, repeat): gene list}, {(level, repeat): spots x signature-gene noise matrix} and the column
        indices of the signature genes.
    """
    rng = np.random.default_rng(seed + 1000)
    var = list(adata.var_names)
    common = [g for g in sig_genes if g in adata.var_names]
    common_set = set(common)
    non_sig = [g for g in var if g not in common_set]
    contam = {}
    for frac in fractions:
        n_contam = int(round(len(common) * frac))
        for rep in range(n_repeats):
            noise_genes = (
                list(rng.choice(non_sig, n_contam, replace=len(non_sig) < n_contam))
                if n_contam > 0
                else []
            )
            contam[(frac, rep)] = common + noise_genes
    X = adata.X.toarray() if issparse(adata.X) else np.asarray(adata.X)
    gene_idx = [var.index(g) for g in common]
    sds = np.array([X[:, gi].std() for gi in gene_idx])
    noise = {}
    for level in noise_levels:
        for rep in range(n_repeats):
            noise[(level, rep)] = np.column_stack(
                [
                    rng.normal(0, sds[j] * level, X.shape[0])
                    for j in range(len(gene_idx))
                ]
            )
    return contam, noise, gene_idx


def f2b_noisy_adata(adata, gene_idx, noise):
    """
    Slide with Gaussian noise added to the signature genes (X + N(0, (sd_g * level)^2), column by column).

    Parameters
    ----------
    adata : anndata.AnnData
        Simulated slide (log-normalised expression, obsm["spatial"]).
    gene_idx : list
        Column indices of the signature genes.
    noise : np.ndarray
        Spots x signature-gene noise matrix (from `f2b_perturbations`).

    Returns
    -------
    anndata.AnnData
        Dense noisy expression with the same spots, genes and coordinates.
    """
    Xn = adata.X.toarray() if issparse(adata.X) else np.asarray(adata.X).copy()
    for j, gi in enumerate(gene_idx):
        Xn[:, gi] += noise[:, j]
    out = anndata.AnnData(
        X=Xn,
        obs=pd.DataFrame(index=adata.obs_names),
        var=pd.DataFrame(index=adata.var_names),
    )
    out.obsm["spatial"] = np.asarray(adata.obsm["spatial"])
    return out


def f2b_f1(y, pred):
    """
    F1 of binary calls.

    Parameters
    ----------
    y, pred : np.ndarray
        Binary labels and calls.

    Returns
    -------
    float
        F1 (0 when there are no labels and no calls).
    """
    y, pred = np.asarray(y).astype(bool), np.asarray(pred).astype(bool)
    denom = pred.sum() + y.sum()
    return float(2 * (y & pred).sum() / denom) if denom else 0.0


def f2b_f1_oracle(y, s):
    """
    F1 at the F1-maximising threshold chosen on the evaluated spots (an optimistic upper bound).

    Parameters
    ----------
    y : np.ndarray
        Binary labels.
    s : np.ndarray
        Score per spot.

    Returns
    -------
    tuple of float
        Best F1 and its threshold.
    """
    p, r, t = precision_recall_curve(y, s)
    p, r = p[:-1], r[:-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.nan_to_num(2 * p * r / (p + r))
    i = int(np.argmax(f))
    return float(f[i]), float(t[i])


def f2b_otsu_threshold(s, nbins=256):
    """
    Label-free Otsu threshold of a score.

    Parameters
    ----------
    s : np.ndarray
        Score per spot.
    nbins : int
        Histogram bins.

    Returns
    -------
    float
        Threshold.
    """
    hist, edges = np.histogram(np.asarray(s, dtype=np.float64), bins=nbins)
    mids = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * mids) / np.maximum(w0, 1)
    m1 = (np.sum(hist * mids) - np.cumsum(hist * mids)) / np.maximum(w1, 1)
    return float(mids[int(np.argmax(w0 * w1 * (m0 - m1) ** 2))])


def f2b_f1_spatial_cv(y, s, folds):
    """
    F1 with the F1-optimal threshold learnt on the training blocks and applied to the held-out blocks.

    Parameters
    ----------
    y : np.ndarray
        Binary labels.
    s : np.ndarray
        Score per spot.
    folds : np.ndarray
        Fold per spot (`f2b_spatial_folds`).

    Returns
    -------
    float
        F1 of the pooled held-out calls.
    """
    y, s = np.asarray(y).astype(int), np.asarray(s, dtype=np.float64)
    pred = np.zeros(len(y), bool)
    for k in np.unique(folds):
        test = folds == k
        if y[~test].sum() == 0:
            continue
        pred[test] = s[test] >= f2b_f1_oracle(y[~test], s[~test])[1]
    return f2b_f1(y, pred)


def f2b_morans_i(values, W):
    """
    Moran's I with a binary spatial weight matrix (as _morans_i_continuous in 17_benchmarking).

    Parameters
    ----------
    values : np.ndarray
        Value per spot.
    W : scipy.sparse matrix
        Spatial weights.

    Returns
    -------
    float
        Moran's I (0 for a constant vector).
    """
    v = np.asarray(values, dtype=np.float64)
    z = v - v.mean()
    z2 = z @ z
    if z2 == 0:
        return 0.0
    return float((len(v) / W.sum()) * (z @ W.dot(z)) / z2)


def f2b_spatial_coherence(pred, W, gt_morans):
    """
    Moran's I of binary calls divided by Moran's I of the truth (17_benchmarking), clipped to [0, 1].

    Parameters
    ----------
    pred : np.ndarray
        Binary calls.
    W : scipy.sparse matrix
        Spatial weights.
    gt_morans : float
        Moran's I of the ground truth.

    Returns
    -------
    float
        Spatial coherence.
    """
    mi = f2b_morans_i(np.asarray(pred, dtype=float), W)
    return float(np.clip(mi / gt_morans, 0, 1)) if gt_morans > 0 else max(mi, 0.0)


def f2b_classification_metrics(y, s, conn, folds, bmask):
    """
    Classification metrics of one score against the simulated target domain (the factorial metric set).

    Parameters
    ----------
    y : np.ndarray
        Binary labels.
    s : np.ndarray
        Score per spot.
    conn : scipy.sparse matrix
        6-NN connectivity (Moran's I and spatial coherence).
    folds : np.ndarray
        Spatial-block folds (`f2b_spatial_folds`).
    bmask : np.ndarray
        Boundary spots (`f2b_boundary_mask`).

    Returns
    -------
    dict
        Metric -> value (empty for a constant or NaN score).
    """
    y, s = np.asarray(y).astype(int), np.asarray(s, dtype=np.float64)
    if np.isnan(s).any() or np.nanstd(s) == 0:
        return {}
    f1o, _ = f2b_f1_oracle(y, s)
    gt = f2b_morans_i(y.astype(float), conn)
    top = np.zeros(len(y), int)
    top[np.argsort(-s, kind="stable")[: int(y.sum())]] = 1
    return {
        "AUROC": roc_auc_score(y, s),
        "AUPRC": average_precision_score(y, s),
        "F1 (oracle threshold)": f1o,
        "F1 (median threshold)": f2b_f1(y, s > np.median(s)),
        "F1 (top-k, prevalence-matched)": f2b_f1(y, top),
        "F1 (Otsu, label-free)": f2b_f1(y, s > f2b_otsu_threshold(s)),
        "F1 (spatial-block CV threshold)": f2b_f1_spatial_cv(y, s, folds),
        "Moran's I": f2b_morans_i(s, conn),
        "Boundary F1 (top-k)": float(f1_score(y[bmask], top[bmask], zero_division=0)),
        "Spatial coherence (top-k)": f2b_spatial_coherence(top, conn, gt),
        "KS": float(stats.ks_2samp(s[y == 1], s[y == 0]).statistic),
    }


def f2b_concordance_labels_smoothed(scores, labels, conn, iters):
    """
    Sup. Fig. 1a as implemented: Spearman rho between the score and the labels smoothed k times.

    Parameters
    ----------
    scores : np.ndarray
        Score per spot.
    labels : np.ndarray
        Binary labels.
    conn : scipy.sparse matrix
        6-NN connectivity.
    iters : list of int
        Smoothing iterations reported.

    Returns
    -------
    dict
        k -> rho.
    """
    cur = np.asarray(labels, dtype=np.float64)
    d = np.maximum(np.asarray(conn.sum(axis=1)).ravel(), 1e-10)
    out = {}
    for k in range(1, max(iters) + 1):
        cur = conn.dot(cur) / d
        if k in iters:
            rho = stats.spearmanr(scores, cur).correlation
            out[k] = 0.0 if np.isnan(rho) else float(rho)
    return out


def f2b_concordance_scores_smoothed(scores, labels, conn, iters):
    """
    Methods-text variant of Sup. Fig. 1a: the score is smoothed k times and correlated with the binary labels.

    Parameters
    ----------
    scores : np.ndarray
        Score per spot.
    labels : np.ndarray
        Binary labels.
    conn : scipy.sparse matrix
        6-NN connectivity.
    iters : list of int
        Smoothing iterations reported.

    Returns
    -------
    dict
        k -> rho.
    """
    cur = np.asarray(scores, dtype=np.float64)
    d = np.maximum(np.asarray(conn.sum(axis=1)).ravel(), 1e-10)
    out = {}
    for k in range(1, max(iters) + 1):
        cur = conn.dot(cur) / d
        if k in iters:
            rho = stats.spearmanr(cur, labels).correlation
            out[k] = 0.0 if np.isnan(rho) else float(rho)
    return out


def f2b_gearys_c(values, W):
    """
    Geary's C with a binary spatial weight matrix (1 = no autocorrelation, < 1 = positive autocorrelation).

    Parameters
    ----------
    values : np.ndarray
        Score per spot.
    W : scipy.sparse matrix
        Spatial connectivities (6-NN; the graph of Moran's I).

    Returns
    -------
    float
        Geary's C.
    """
    v = np.asarray(values, dtype=np.float64)
    z = v - v.mean()
    z2 = z @ z
    if z2 == 0:
        return np.nan
    W = W.tocsr()
    rows = np.asarray(W.sum(axis=1)).ravel()
    cols = np.asarray(W.sum(axis=0)).ravel()
    sq_diff = rows @ (v**2) + cols @ (v**2) - 2 * v @ W.dot(v)
    return float((len(v) - 1) * sq_diff / (2 * W.sum() * z2))


def f2b_hop_distance(conn, mask, max_hops=12):
    """
    Graph distance (rings of 6-NN neighbours) of every spot to the nearest spot in `mask`.

    Parameters
    ----------
    conn : scipy.sparse matrix
        6-NN connectivity.
    mask : np.ndarray
        Boolean mask of the target domain.
    max_hops : int
        Largest ring computed (farther spots stay at inf).

    Returns
    -------
    np.ndarray
        Ring index per spot (0 inside the domain).
    """
    mask = np.asarray(mask, bool)
    A = csr_matrix((conn + conn.T) > 0).astype(np.int8)
    dist = np.full(len(mask), np.inf)
    dist[mask] = 0
    frontier, seen = mask.copy(), mask.copy()
    for h in range(1, max_hops + 1):
        nxt = (A.dot(frontier.astype(np.int8)) > 0) & ~seen
        dist[nxt] = h
        seen |= nxt
        frontier = nxt
        if not nxt.any():
            break
    return dist


def f2b_leakage_profile(score, mask, dist, max_ring=8, far_from=6):
    """
    Standardised score by ring outside the domain: L(d) = (ring-d mean - far mean) / (inside mean - far mean).

    Parameters
    ----------
    score : np.ndarray
        Score per spot.
    mask : np.ndarray
        Boolean mask of the target domain.
    dist : np.ndarray
        Ring index per spot (`f2b_hop_distance`).
    max_ring : int
        Last ring reported.
    far_from : int
        Far background = spots at least this many rings away.

    Returns
    -------
    dict
        Ring -> L(d) (0 = no leakage, 1 = as high as inside).
    """
    score = np.asarray(score, dtype=np.float64)
    inside = score[np.asarray(mask, bool)].mean()
    far = score[dist >= far_from].mean()
    span = inside - far
    prof = {}
    for d in range(1, max_ring + 1):
        sel = dist == d
        prof[d] = (
            float((score[sel].mean() - far) / span)
            if sel.any() and span != 0
            else np.nan
        )
    return prof


def f2b_leakage_summary(score, mask, dist, halo=(1, 2, 3)):
    """
    Scalar leakage metrics: halo leakage index, region-vs-halo AUROC and the top-k shares in the halo and inside.

    Parameters
    ----------
    score : np.ndarray
        Score per spot.
    mask : np.ndarray
        Boolean mask of the target domain.
    dist : np.ndarray
        Ring index per spot (`f2b_hop_distance`).
    halo : tuple of int
        Rings that form the halo.

    Returns
    -------
    dict
        Metric -> value.
    """
    score = np.asarray(score, dtype=np.float64)
    mask = np.asarray(mask, bool)
    prof = f2b_leakage_profile(score, mask, dist)
    halo_sel = np.isin(dist, halo)
    y = np.r_[np.ones(mask.sum()), np.zeros(halo_sel.sum())]
    sv = np.r_[score[mask], score[halo_sel]]
    top = np.argsort(-score, kind="stable")[: int(mask.sum())]
    return {
        "halo leakage index": float(np.nanmean([prof[d] for d in halo])),
        "region-vs-halo AUROC": float(roc_auc_score(y, sv)),
        "top-k share in halo": float(np.isin(dist[top], halo).mean()),
        "top-k share inside": float(mask[top].mean()),
    }


def f2b_summarise(df, metrics, methods, configs):
    """
    Mean, SD, SEM and median of per-slide metrics.

    Parameters
    ----------
    df : pd.DataFrame
        Long table with `method`, `config`, `metric`, `value` (one row per slide).
    metrics, methods, configs : list
        Levels to summarise.

    Returns
    -------
    pd.DataFrame
        One row per method x configuration x metric.
    """
    out = []
    for m in methods:
        for c in configs:
            for met in metrics:
                v = df[(df.method == m) & (df.config == c) & (df.metric == met)][
                    "value"
                ].dropna()
                if len(v):
                    out.append(
                        {
                            "method": m,
                            "config": c,
                            "metric": met,
                            "n": len(v),
                            "mean": v.mean(),
                            "sd": v.std(ddof=1),
                            "sem": v.std(ddof=1) / np.sqrt(len(v)),
                            "median": v.median(),
                        }
                    )
    return pd.DataFrame(out)


def f2b_boot_ci(values, n_boot=2000, seed=0):
    """
    95% percentile bootstrap CI of the mean over slides.

    Parameters
    ----------
    values : np.ndarray
        One value per slide.
    n_boot : int
        Bootstrap resamples.
    seed : int
        Seed of the resampling.

    Returns
    -------
    tuple
        (lower, upper).
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    rng = np.random.default_rng(seed)
    means = v[rng.integers(0, len(v), (n_boot, len(v)))].mean(axis=1)
    return tuple(np.percentile(means, [2.5, 97.5]))


def f2b_mean_ci(per_slide, metrics, method_order, config_order):
    """
    Mean and 95% bootstrap CI over slides per method x configuration.

    Parameters
    ----------
    per_slide : pd.DataFrame
        Columns `method`, `config`, `seed` and the metrics.
    metrics : list
        Metric columns.
    method_order, config_order : list
        Row order of the output.

    Returns
    -------
    pd.DataFrame
        method, config, and for each metric: mean, `_lo`, `_hi` and the number of slides.
    """
    rows = []
    for (m, c), g in per_slide.groupby(["method", "config"], sort=False):
        row = {"method": m, "config": c}
        for k in metrics:
            v = g[k].dropna().values if k in g else np.array([])
            mu = v.mean() if len(v) else np.nan
            lo, hi = f2b_boot_ci(v) if len(v) > 1 else (mu, mu)
            row[k], row[f"{k}_lo"], row[f"{k}_hi"], row[f"{k} (n slides)"] = (
                mu,
                min(lo, mu),
                max(hi, mu),
                len(v),
            )
        rows.append(row)
    out = pd.DataFrame(rows)
    order = {m: i for i, m in enumerate(method_order)}
    corder = {c: i for i, c in enumerate(config_order)}
    out["_o"] = [
        order.get(m, 99) * 10 + corder.get(c, 9) for m, c in zip(out.method, out.config)
    ]
    return out.sort_values("_o").drop(columns="_o").reset_index(drop=True)


def f2b_paired_table(df, metrics, methods, contrasts):
    """
    Paired-by-slide differences between configurations (95% t-interval, Wilcoxon signed-rank p).

    Parameters
    ----------
    df : pd.DataFrame
        Long table with `seed`, `method`, `config`, `metric`, `value`.
    metrics, methods : list
        Levels to test.
    contrasts : dict
        Label -> (configuration without, configuration with).

    Returns
    -------
    pd.DataFrame
        One row per metric x method x contrast.
    """
    out = []
    for met in metrics:
        piv = df[df.metric == met].pivot_table(
            index=["method", "seed"], columns="config", values="value"
        )
        for m in methods:
            if m not in piv.index.get_level_values(0):
                continue
            pm = piv.loc[m]
            for lab_, (a, b) in contrasts.items():
                if a not in pm or b not in pm:
                    continue
                d = (pm[b] - pm[a]).dropna()
                if len(d) < 3:
                    continue
                ci = stats.t.interval(
                    0.95,
                    len(d) - 1,
                    loc=d.mean(),
                    scale=d.std(ddof=1) / np.sqrt(len(d)),
                )
                pw = stats.wilcoxon(d).pvalue if (d != 0).any() else np.nan
                out.append(
                    {
                        "metric": met,
                        "method": m,
                        "contrast": lab_,
                        "n seeds": len(d),
                        "mean without": pm[a][d.index].mean(),
                        "mean with": pm[b][d.index].mean(),
                        "mean diff": d.mean(),
                        "sd diff": d.std(ddof=1),
                        "95% CI low": ci[0],
                        "95% CI high": ci[1],
                        "median diff": d.median(),
                        "frac seeds improved": (d > 0).mean(),
                        "wilcoxon p": pw,
                    }
                )
    return pd.DataFrame(out)


def f2b_vs_enrichmap(df, metric, cfg_em, cfg_other, methods):
    """
    Paired-by-slide difference EnrichMap[cfg_em] - competitor[cfg_other].

    Parameters
    ----------
    df : pd.DataFrame
        Long classification table.
    metric : str
        Metric compared.
    cfg_em, cfg_other : str
        Configurations of EnrichMap and of the competitor.
    methods : list
        Competitors.

    Returns
    -------
    list of dict
        One row per competitor.
    """
    piv = df[df.metric == metric].pivot_table(
        index="seed", columns=["method", "config"], values="value"
    )
    out = []
    for m in methods:
        if (m, cfg_other) not in piv:
            continue
        d = (piv[("EnrichMap", cfg_em)] - piv[(m, cfg_other)]).dropna()
        ci = stats.t.interval(
            0.95, len(d) - 1, loc=d.mean(), scale=d.std(ddof=1) / np.sqrt(len(d))
        )
        out.append(
            {
                "metric": metric,
                "comparison": f"EnrichMap[{cfg_em}] − {m}[{cfg_other}]",
                "method": m,
                "n": len(d),
                "mean diff": d.mean(),
                "95% CI low": ci[0],
                "95% CI high": ci[1],
                "frac EnrichMap better": (d > 0).mean(),
                "wilcoxon p": stats.wilcoxon(d).pvalue if (d != 0).any() else np.nan,
            }
        )
    return out


def f2b_robustness_summary(df):
    """
    Average the repeats per slide, then summarise contamination or noise records over slides.

    Parameters
    ----------
    df : pd.DataFrame
        Records with `method`, `config`, `level`, `seed`, `rep`, `rho_self`, `auroc`, `auprc`.

    Returns
    -------
    tuple of pd.DataFrame
        Per-slide means and the mean / SD / count over slides.
    """
    per_seed = (
        df.groupby(["method", "config", "level", "seed"])[
            ["rho_self", "auroc", "auprc"]
        ]
        .mean()
        .reset_index()
    )
    out = per_seed.groupby(["method", "config", "level"])[
        ["rho_self", "auroc", "auprc"]
    ].agg(["mean", "std", "count"])
    out.columns = [f"{a}_{b}" for a, b in out.columns]
    return per_seed, out.reset_index()


def f2b_show_robustness(rob, val, name, levels, methods, configs):
    """
    Print mean ± SD over slides of one robustness metric per scorer x configuration x level.

    Parameters
    ----------
    rob : pd.DataFrame
        Output of `f2b_robustness_summary` (second element).
    val : str
        Metric (`rho_self` or `auroc`).
    name : str
        Table title.
    levels, methods, configs : list
        Levels, scorers and configurations shown.
    """
    t = rob[
        rob.method.isin(methods) & rob.config.isin(configs) & rob.level.isin(levels)
    ]
    t = t.assign(
        v=[f"{a:.3f}±{b:.3f}" for a, b in zip(t[val + "_mean"], t[val + "_std"])]
    )
    t = t.pivot_table(
        index=["method", "config"], columns="level", values="v", aggfunc="first"
    )
    print(f"\n=== {name} ===")
    print(
        t.reindex(
            [(m, c) for m in methods for c in configs if (m, c) in t.index]
        ).to_string()
    )


def f2b_robustness_area(ps, val):
    """
    Normalised area of one robustness curve per slide: mean over the perturbed levels (> 0), complete slides only.

    Parameters
    ----------
    ps : pd.DataFrame
        Per-slide means (`method`, `config`, `level`, `seed`, metrics).
    val : str
        Metric (`rho_self` or `auroc`).

    Returns
    -------
    pd.Series
        Area indexed by method, config and seed.
    """
    p = ps[ps.level > 0]
    g = p.groupby(["method", "config", "seed"])[val]
    return g.mean()[g.count() == p.level.nunique()]


def f2b_headline_table(table, cols, headline):
    """
    Headline configuration of each method (rows in the published order) from a method x configuration table.

    Parameters
    ----------
    table : pd.DataFrame
        Columns `method`, `config` and the metrics (output of `f2b_mean_ci`).
    cols : dict
        Metric column -> short column label.
    headline : list of tuple
        (method, configuration, published label) per row.

    Returns
    -------
    pd.DataFrame
        Index = published method labels, columns = short labels.
    """
    t = table.merge(
        pd.DataFrame(headline, columns=["method", "config", "label"]),
        on=["method", "config"],
    )
    return (
        t.set_index("label")[list(cols)]
        .rename(columns=cols)
        .reindex([label for *_, label in headline])
    )


def f2b_categorical_to_str(df, labels):
    """
    Plain strings for categorical columns and published scorer labels (source data).

    Parameters
    ----------
    df : pd.DataFrame
        Table to export.
    labels : dict
        Scorer name -> published label.

    Returns
    -------
    pd.DataFrame
        Copy with a default index.
    """
    df = df.reset_index(drop=True).copy()
    for col in df.columns:
        if isinstance(df[col].dtype, pd.CategoricalDtype):
            df[col] = df[col].astype(str)
    if "method" in df:
        df["method"] = df["method"].replace(labels)
    return df


def f2b_plot_metric_pair(
    table, x, y, ax, colours, ref_x=None, ref_y=None, invert_y=False, marker_col=None
):
    """
    One metric against another, one point per row (palette colours), with 95% CI whiskers when f"{x}_lo"/f"{x}_hi"
    and f"{y}_lo"/f"{y}_hi" exist.

    Parameters
    ----------
    table : pd.DataFrame
        One row per point, with a "method" column holding the palette key.
    x, y : str
        Metric columns.
    ax : matplotlib.axes.Axes
        Target axes.
    colours : dict
        Palette.
    ref_x, ref_y : float or None
        Dotted reference lines.
    invert_y : bool
        Invert the y axis (lower-is-better metrics).
    marker_col : str or None
        Column holding a marker per row.

    Returns
    -------
    matplotlib.axes.Axes
    """
    for _, r in table.iterrows():
        xerr = (
            [[r[x] - r[f"{x}_lo"]], [r[f"{x}_hi"] - r[x]]]
            if f"{x}_lo" in table
            else None
        )
        yerr = (
            [[r[y] - r[f"{y}_lo"]], [r[f"{y}_hi"] - r[y]]]
            if f"{y}_lo" in table
            else None
        )
        ax.errorbar(
            r[x],
            r[y],
            xerr=xerr,
            yerr=yerr,
            fmt=r[marker_col] if marker_col else "o",
            color=colours[r["method"]],
            ms=5,
            mec="white",
            mew=0.4,
            elinewidth=0.6,
            capsize=0,
            label=r["method"],
        )
    if ref_x is not None:
        ax.axvline(ref_x, ls=":", c="0.6", lw=0.6)
    if ref_y is not None:
        ax.axhline(ref_y, ls=":", c="0.6", lw=0.6)
    if invert_y:
        ax.invert_yaxis()
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.grid(False)
    return ax


def f2b_plot_benchmark_summary(table, ax, colours, lower_is_better=(), annotate=True):
    """
    Methods x metrics summary: circles sized and shaded by each metric scaled to 0-1 across methods (best = 1), raw
    values annotated, plus a mean-rank column.

    Parameters
    ----------
    table : pd.DataFrame
        Index = methods (published labels), columns = metrics.
    ax : matplotlib.axes.Axes
        Target axes.
    colours : dict
        Palette (colours the method labels).
    lower_is_better : tuple of str
        Metrics where smaller is better.
    annotate : bool
        Print the raw value in each cell.

    Returns
    -------
    matplotlib.axes.Axes
    """
    scaled = table.copy().astype(float)
    for c in table.columns:
        col = table[c].astype(float)
        rng_ = col.max() - col.min()
        s = (col - col.min()) / rng_ if rng_ > 0 else col * 0 + 1
        scaled[c] = 1 - s if c in lower_is_better else s
    ranks = scaled.rank(ascending=False).mean(axis=1)
    cmap = plt.get_cmap("Purples")
    for i, m in enumerate(table.index):
        for j, c in enumerate(table.columns):
            s = scaled.loc[m, c]
            ax.scatter(
                j, i, s=30 + 220 * s, color=cmap(0.25 + 0.75 * s), edgecolors="none"
            )
            if annotate:
                ax.text(
                    j,
                    i,
                    f"{table.loc[m, c]:.2g}",
                    ha="center",
                    va="center",
                    fontsize=5,
                    color="white" if s > 0.6 else "black",
                )
        ax.text(
            len(table.columns) + 0.2,
            i,
            f"{ranks[m]:.1f}",
            ha="left",
            va="center",
            fontsize=6,
        )
    ax.set_xticks(range(len(table.columns)))
    ax.set_xticklabels(table.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(table.index)))
    ax.set_yticklabels(table.index)
    for lab in ax.get_yticklabels():
        lab.set_color(colours.get(lab.get_text(), "black"))
    ax.text(len(table.columns) + 0.2, -0.9, "mean rank", fontsize=6, ha="left")
    ax.set_xlim(-0.6, len(table.columns) + 0.9)
    ax.set_ylim(len(table.index) - 0.4, -1.2)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    return ax


def f2b_pair_rows(table, x, y, methods, configs, style):
    """
    Rows of a factorial metric-pair panel, with palette key and marker (null model last, drawn on top).

    Parameters
    ----------
    table : pd.DataFrame
        Output of `f2b_mean_ci`.
    x, y : str
        Metrics on the axes.
    methods, configs : list
        Points drawn.
    style : dict
        Plot style (see the section header).

    Returns
    -------
    pd.DataFrame
        Plotted rows with `colour key` and `marker`.
    """
    t = (
        table[table.method.isin(methods) & table.config.isin(configs)]
        .dropna(subset=[x, y])
        .copy()
    )
    t = pd.concat(
        [t[t.config != style["null_config"]], t[t.config == style["null_config"]]]
    )
    t["colour key"] = [style["colour_key"][(m, c)] for m, c in zip(t.method, t.config)]
    t["marker"] = t.config.map(style["markers"])
    return t


def f2b_pair_source(table, x, y, methods, configs, style):
    """
    Source data of a factorial metric-pair panel (published scorer labels).

    Parameters
    ----------
    table : pd.DataFrame
        Output of `f2b_mean_ci`.
    x, y : str
        Metrics on the axes.
    methods, configs : list
        Points drawn.
    style : dict
        Plot style.

    Returns
    -------
    pd.DataFrame
        method, config, colour key and each metric with its CI and number of slides.
    """
    t = f2b_pair_rows(table, x, y, methods, configs, style)
    cols = ["method", "config", "colour key"] + [
        c for k in (x, y) for c in (k, f"{k}_lo", f"{k}_hi", f"{k} (n slides)")
    ]
    return (
        t[cols].assign(method=t.method.replace(style["labels"])).reset_index(drop=True)
    )


def f2b_factorial_pair(
    ax, table, x, y, methods, configs, style, ref_x=None, ref_y=None, invert_y=False
):
    """
    Compact factorial panel: colour = scorer (EnrichMap: the matching variant), marker = configuration, whiskers =
    95% bootstrap CI over slides.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    table : pd.DataFrame
        Output of `f2b_mean_ci`.
    x, y : str
        Metrics on the axes.
    methods, configs : list
        Points drawn.
    style : dict
        Plot style.
    ref_x, ref_y : float or None
        Dotted reference lines.
    invert_y : bool
        Invert the y axis.

    Returns
    -------
    pd.DataFrame
        The plotted rows (source data).
    """
    t = f2b_pair_rows(table, x, y, methods, configs, style)
    f2b_plot_metric_pair(
        t.assign(method=t["colour key"]),
        x,
        y,
        ax,
        style["colours"],
        ref_x=ref_x,
        ref_y=ref_y,
        invert_y=invert_y,
        marker_col="marker",
    )
    sns.despine(ax=ax)
    return f2b_pair_source(table, x, y, methods, configs, style)


def f2b_method_handles(methods, configs, style, kind="dot"):
    """
    Legend handles of the scorers (EnrichMap: one entry per variant colour).

    Parameters
    ----------
    methods : list
        Scorers in legend order.
    configs : list
        Configurations shown (decides which EnrichMap entries appear).
    style : dict
        Plot style (optional "enrichmap_labels": configuration -> legend text of EnrichMap).
    kind : str
        "dot" (marker only) or "line" (line + marker; oracle truth dashed).

    Returns
    -------
    list of matplotlib.lines.Line2D
    """
    colour = {k: style["colours"][v] for k, v in style["colour_key"].items()}
    hs = []
    for m in methods:
        if m == "EnrichMap":
            entries = [
                (style["null_config"], "EnrichMap (null)"),
                ("smooth", "EnrichMap (smoothing-only) = smooth"),
                ("smooth+GAM", "EnrichMap (GAM) = smooth+GAM"),
                ("raw", "EnrichMap unsmoothed (raw, GAM)"),
            ]
            for cfg, text in entries:
                if cfg in configs or (cfg == "raw" and "GAM" in configs):
                    hs.append(
                        Line2D(
                            [],
                            [],
                            color=colour[(m, cfg)],
                            marker="o",
                            ls="-" if kind == "line" else "",
                            lw=0.9,
                            ms=3,
                            label=style.get("enrichmap_labels", {}).get(cfg, text),
                        )
                    )
        else:
            ls = ("--" if m == "Truth (oracle)" else "-") if kind == "line" else ""
            hs.append(
                Line2D(
                    [],
                    [],
                    color=colour[(m, "raw")],
                    marker="o",
                    ls=ls,
                    lw=0.9,
                    ms=3,
                    label=style["labels"].get(m, m),
                )
            )
    return hs


def f2b_config_handles(configs, style):
    """
    Legend handles of the configurations (marker shapes).

    Parameters
    ----------
    configs : list
        Configurations.
    style : dict
        Plot style.

    Returns
    -------
    list of matplotlib.lines.Line2D
    """
    return [
        Line2D(
            [],
            [],
            color="0.35",
            marker=style["markers"][c],
            ls="",
            ms=3.5,
            label=style["config_labels"].get(c, c),
        )
        for c in configs
    ]


def f2b_legend(
    ax,
    methods,
    configs,
    style,
    kind="dot",
    loc="upper left",
    bbox=(1.02, 1.0),
    fontsize=5.5,
):
    """
    One legend to the right of an axis: scorers (colour), then configurations (marker shape).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes holding the legend.
    methods, configs : list
        Scorers and configurations shown.
    style : dict
        Plot style.
    kind : str
        "dot" or "line" handles.
    loc : str
        Legend location.
    bbox : tuple
        bbox_to_anchor.
    fontsize : float
        Font size.
    """
    handles = f2b_method_handles(methods, configs, style, kind)
    if kind == "dot":
        handles += [
            Line2D([], [], ls="", label=" "),
            Line2D([], [], ls="", label="Configuration (marker)"),
        ]
        handles += f2b_config_handles(configs, style)
    ax.legend(
        handles=handles,
        title="Scorer (colour)",
        frameon=False,
        fontsize=fontsize,
        title_fontsize=6,
        loc=loc,
        bbox_to_anchor=bbox,
        handlelength=1.2,
        alignment="left",
    )


def f2b_dotplot(
    ax, summary, metric, methods, configs, style, err="sd", ylim=None, title=None
):
    """
    Mean ± error per scorer (x) and configuration (dodged; colour = scorer, marker = configuration).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    summary : pd.DataFrame
        Table with `method`, `config`, `metric`, `mean` and the `err` column.
    metric : str
        Metric to plot.
    methods, configs : list
        Order of the x groups and of the dodged configurations.
    style : dict
        Plot style.
    err : str
        Column used for the error bars (`sd` or `sem`).
    ylim : tuple, optional
        y limits.
    title : str, optional
        Axes title (default: the metric).
    """
    w = 0.8 / len(configs)
    for i, m in enumerate(methods):
        for j, c in enumerate(configs):
            r = summary[
                (summary.method == m)
                & (summary.config == c)
                & (summary.metric == metric)
            ]
            if not len(r):
                continue
            ax.errorbar(
                i - 0.4 + w * (j + 0.5),
                r["mean"].iloc[0],
                yerr=r[err].iloc[0],
                fmt=style["markers"][c],
                ms=3,
                color=style["colours"][style["colour_key"][(m, c)]],
                mec="white",
                mew=0.3,
                elinewidth=0.8,
                capsize=1.2,
            )
    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(
        [
            style["labels"].get(m, m) if m != "Truth (oracle)" else "Truth"
            for m in methods
        ],
        rotation=40,
        ha="right",
    )
    ax.set_title(metric if title is None else title)
    if ylim:
        ax.set_ylim(*ylim)


def f2b_lineplot_by_config(
    axrow,
    df,
    xcol,
    methods,
    configs,
    style,
    ylabel,
    xlabel,
    err="std",
    ylim=None,
    titles=True,
    xlabels=True,
    bands=True,
):
    """
    One axes per configuration; mean line ± `err` band per scorer.

    Parameters
    ----------
    axrow : sequence of matplotlib.axes.Axes
        One axes per configuration.
    df : pd.DataFrame
        Table with `method`, `config`, `xcol`, `mean` and `err`.
    xcol : str
        x variable (`k` or `level`).
    methods, configs : list
        Scorers drawn on every axes and configurations (one per axes).
    style : dict
        Plot style.
    ylabel, xlabel : str
        Axis labels.
    err : str
        Column for the band.
    ylim : tuple, optional
        y limits.
    titles, xlabels, bands : bool
        Draw configuration titles / x labels / ± err bands.
    """
    for i, (ax, cfg) in enumerate(zip(axrow, configs)):
        for m in methods:
            s = df[(df.method == m) & (df.config == cfg)].sort_values(xcol)
            if not len(s):
                continue
            col = style["colours"][style["colour_key"][(m, cfg)]]
            ax.plot(
                s[xcol],
                s["mean"],
                "--" if m == "Truth (oracle)" else "-",
                marker="o",
                ms=2,
                lw=0.9,
                color=col,
                label=style["labels"].get(m, m),
            )
            if bands:
                ax.fill_between(
                    s[xcol],
                    s["mean"] - s[err],
                    s["mean"] + s[err],
                    color=col,
                    alpha=0.12,
                    lw=0,
                )
        if titles:
            ax.set_title(style["config_labels"].get(cfg, cfg), fontsize=6.5)
        if xlabels:
            ax.set_xlabel(xlabel)
        if ylim:
            ax.set_ylim(*ylim)
        if i > 0:
            ax.tick_params(labelleft=False)
    axrow[0].set_ylabel(ylabel)


def f2b_dot_table(summary, metric, methods, configs, labels, err="sd"):
    """
    The numbers `f2b_dotplot` draws: mean and error-bar bounds per scorer and configuration.

    Parameters
    ----------
    summary : pd.DataFrame
        Table with `method`, `config`, `metric`, `n`, `mean` and the `err` column.
    metric : str
        Plotted metric.
    methods, configs : list
        Plotted scorers and configurations.
    labels : dict
        Scorer name -> published label.
    err : str
        Error-bar column (`sd` or `sem`).

    Returns
    -------
    pd.DataFrame
        One row per point.
    """
    rows = []
    for m in methods:
        for c in configs:
            r = summary[
                (summary.method == m)
                & (summary.config == c)
                & (summary.metric == metric)
            ]
            if not len(r):
                continue
            mu, e = r["mean"].iloc[0], r[err].iloc[0]
            rows.append(
                {
                    "method": labels.get(m, m),
                    "configuration": c,
                    "n slides": r["n"].iloc[0],
                    "mean": mu,
                    err: e,
                    "error bar low": mu - e,
                    "error bar high": mu + e,
                }
            )
    return pd.DataFrame(rows)


def f2b_line_table(df, xcol, methods, cfg, n_col, labels, err="std"):
    """
    The numbers `f2b_lineplot_by_config` draws on one axes: mean and band bounds per scorer and x.

    Parameters
    ----------
    df : pd.DataFrame
        Table passed to `f2b_lineplot_by_config`.
    xcol : str
        x variable.
    methods : list
        Plotted scorers.
    cfg : str
        Configuration of this axes.
    n_col : str
        Column holding the number of slides.
    labels : dict
        Scorer name -> published label.
    err : str
        Band column.

    Returns
    -------
    pd.DataFrame
        One row per point.
    """
    out = []
    for m in methods:
        s = df[(df.method == m) & (df.config == cfg)].sort_values(xcol)
        for _, r in s.iterrows():
            out.append(
                {
                    "method": labels.get(m, m),
                    "config": cfg,
                    xcol: r[xcol],
                    "n slides": r[n_col],
                    "mean": r["mean"],
                    "sd": r[err],
                    "band low": r["mean"] - r[err],
                    "band high": r["mean"] + r[err],
                }
            )
    return pd.DataFrame(out)


def f2b_fixed_rule_strip(ax, table):
    """
    Fixed, label-free calls next to the rule summary: circles sized and shaded by the F1 (0-1), with values; empty
    where a rule does not exist for the method (rows aligned with `f2b_plot_benchmark_summary`).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    table : pd.DataFrame
        Index = methods (same order as the summary), columns = fixed rules, NaN where not applicable.
    """
    cmap = plt.get_cmap("Purples")
    for i, m in enumerate(table.index):
        for j, c in enumerate(table.columns):
            v = table.loc[m, c]
            if not np.isfinite(v):
                continue
            ax.scatter(
                j, i, s=30 + 220 * v, color=cmap(0.25 + 0.75 * v), edgecolors="none"
            )
            ax.text(
                j,
                i,
                f"{v:.2g}",
                ha="center",
                va="center",
                fontsize=5,
                color="white" if v > 0.6 else "black",
            )
    ax.set_xticks(range(len(table.columns)))
    ax.set_xticklabels(table.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(table.index)))
    ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-0.6, len(table.columns) - 0.4)
    ax.set_ylim(len(table.index) - 0.4, -1.2)
    ax.text(
        (len(table.columns) - 1) / 2, -0.9, "fixed, label-free", fontsize=6, ha="center"
    )
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)


def f2b_headline_curves(
    ax, df, xcol, ylabel, xlabel, headline, colours, null_config, ylim=None, truth=True
):
    """
    Curve panel of the nine methods at their headline configuration (and the oracle truth, raw).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    df : pd.DataFrame
        Table with `method`, `config`, `xcol`, `mean`, `std`, `count`.
    xcol : str
        x variable.
    ylabel, xlabel : str
        Axis labels.
    headline : list of tuple
        (method, configuration, published label) per curve.
    colours : dict
        Palette (keys = published labels).
    null_config : str
        Name of the null-model configuration (drawn thicker, underneath).
    ylim : tuple, optional
        y limits.
    truth : bool
        Draw the oracle truth (raw) as a dashed black line.

    Returns
    -------
    pd.DataFrame
        The plotted points (source data).
    """
    rows = headline + ([("Truth (oracle)", "raw", "Truth (oracle)")] if truth else [])
    out = []
    for m, c, label in rows:
        s = df[(df.method == m) & (df.config == c)].sort_values(xcol)
        if not len(s):
            continue
        ax.plot(
            s[xcol],
            s["mean"],
            "--" if m in ("GSDensity", "Truth (oracle)") else "-",
            marker="o",
            ms=3 if c == null_config else 2,
            lw=2.0 if c == null_config else 0.9,
            color=colours[label],
            zorder=2 if c == null_config else 3,
        )
        out.append(
            s.assign(method=label)[["method", "config", xcol, "mean", "std", "count"]]
        )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if ylim:
        ax.set_ylim(*ylim)
    ax.grid(False)
    sns.despine(ax=ax)
    return pd.concat(out, ignore_index=True)


def f2b_compact_layout():
    """
    Compact benchmark layout: a, b (metric panels), c (summary) and the legend to the right of c.

    Returns
    -------
    tuple
        Figure, [axes a, axes b], axes c and the (empty) legend axes.
    """
    fig = plt.figure(figsize=(10.0, 2.5))
    outer = gridspec.GridSpec(1, 3, figure=fig, width_ratios=[1, 1, 2.0], wspace=0.75)
    inner = gridspec.GridSpecFromSubplotSpec(
        1, 2, subplot_spec=outer[2], width_ratios=[1, 0.9], wspace=0.75
    )
    lax = fig.add_subplot(inner[1])
    lax.axis("off")
    return (
        fig,
        [fig.add_subplot(outer[0]), fig.add_subplot(outer[1])],
        fig.add_subplot(inner[0]),
        lax,
    )


# --------------------------------------------------------------------------------------------------
# R4_F2c_when_GAM_helps
# --------------------------------------------------------------------------------------------------


def f2c_cond_label(arm, field, strength):
    """
    Short label of one artefact condition.

    Parameters
    ----------
    arm : str
        Artefact arm (clean, random, worst or global).
    field : str
        Field shape.
    strength : float
        Natural-log range of the field.

    Returns
    -------
    str
        "clean" or "arm|field|strength".
    """
    return "clean" if arm == "clean" else f"{arm}|{field}|{strength:g}"


def f2c_simulate_visium_artefact(
    seed=0,
    n_genes=3000,
    n_signature_genes=20,
    n_clusters=4,
    target_cluster=1,
    n_spots_range=(1500, 4000),
    spot_spacing_um=100,
    spot_diameter_um=55,
    mean_counts_per_spot=5000,
    dispersion_range=(0.1, 2.0),
    fold_change_range=(1.5, 5.0),
    dropout_rate=0.1,
    domain="voronoi",
    domain_width_spots=1.5,
    n_patches=12,
    patch_radius_spots=3.0,
    tech_field=None,
    tech_strength=0.0,
    tech_gene_frac=0.5,
    tech_scale=0.5,
    tech_on_signature=None,
    return_truth=False,
):
    """
    Simulated Visium slide of 05_01_simulated_benchmarking, with domain geometries and technical fields.

    With the default arguments the random stream and the output are identical to simulate_visium of
    05_01_simulated_benchmarking. The extensions are drawn from an independent random stream.

    Parameters
    ----------
    seed : int
        Random seed.
    n_genes, n_signature_genes, n_clusters, target_cluster : int
        Genes, signature genes, Voronoi clusters and the target cluster (0-based).
    n_spots_range : tuple of int
        Range of the number of spots.
    spot_spacing_um, spot_diameter_um : float
        Visium spot spacing and diameter (um).
    mean_counts_per_spot : float
        Mean library size.
    dispersion_range, fold_change_range : tuple of float
        Negative-binomial dispersion and signature fold-change ranges.
    dropout_rate : float
        Fraction of counts set to zero.
    domain : str
        "voronoi" (target = cluster_2), "stripe" (thin sinusoidal band, layer-like) or "patches"
        (n_patches discs of radius patch_radius_spots, niche-like).
    domain_width_spots : float
        Half-width of the stripe in spot spacings.
    n_patches : int
        Number of patches.
    patch_radius_spots : float
        Patch radius in spot spacings.
    tech_field : str or None
        Smooth multiplicative technical field F(x, y) on a fraction `tech_gene_frac` of genes before count
        sampling: "gradient" (linear trend in a random direction), "edge" (capture loss towards the border) or
        "grf" (Gaussian random field with length scale `tech_scale` x side).
    tech_strength : float
        Range of log F (natural log).
    tech_gene_frac : float
        Fraction of genes affected.
    tech_scale : float
        Length scale of the Gaussian random field.
    tech_on_signature : bool or None
        None (signature genes affected at the same rate as the others), True (all) or False (none).
    return_truth : bool
        Also return log F and the affected genes.

    Returns
    -------
    tuple
        (adata, sig_genes, fold_changes, target_label), plus {"log_field", "affected"} if return_truth.
    """
    rng = np.random.default_rng(seed)
    rng_ext = np.random.default_rng(10_000_019 + seed)

    n_spots_target = rng.integers(n_spots_range[0], n_spots_range[1] + 1)
    grid_side = int(np.ceil(np.sqrt(n_spots_target * 0.866)))
    capture_area_um = grid_side * spot_spacing_um
    n_cols = int(capture_area_um / spot_spacing_um)
    n_rows = int(capture_area_um / (spot_spacing_um * np.sqrt(3) / 2))
    x_coords, y_coords = [], []
    for row in range(n_rows):
        for col in range(n_cols):
            x = col * spot_spacing_um + (row % 2) * (spot_spacing_um / 2)
            y = row * spot_spacing_um * np.sqrt(3) / 2
            if x <= capture_area_um and y <= capture_area_um:
                x_coords.append(x)
                y_coords.append(y)
    x_coords = np.array(x_coords)
    y_coords = np.array(y_coords)
    n_spots = len(x_coords)
    coords = np.column_stack([x_coords, y_coords])

    margin = capture_area_um * 0.15
    seed_coords = np.column_stack(
        [
            rng.uniform(margin, capture_area_um - margin, n_clusters),
            rng.uniform(margin, capture_area_um - margin, n_clusters),
        ]
    )
    dists = cdist(coords, seed_coords)
    dists += rng.normal(0, spot_spacing_um * 0.8, dists.shape)
    cluster_ids = np.argmin(dists, axis=1)
    cluster_labels = np.array([f"cluster_{i + 1}" for i in cluster_ids])
    target_label = f"cluster_{target_cluster + 1}"

    if domain == "stripe":
        L = capture_area_um
        phase = rng_ext.uniform(0, 2 * np.pi)
        amp = rng_ext.uniform(0.08, 0.15) * L
        y0 = rng_ext.uniform(0.35, 0.65) * L
        curve = y0 + amp * np.sin(2 * np.pi * coords[:, 0] / (0.8 * L) + phase)
        in_target = np.abs(coords[:, 1] - curve) <= domain_width_spots * spot_spacing_um
        cluster_labels = np.where(
            in_target,
            target_label,
            np.where(cluster_labels == target_label, "cluster_bg", cluster_labels),
        )
    elif domain == "patches":
        L = capture_area_um
        centres = rng_ext.uniform(0.08 * L, 0.92 * L, size=(n_patches, 2))
        dd = cdist(coords, centres).min(axis=1)
        in_target = dd <= patch_radius_spots * spot_spacing_um
        cluster_labels = np.where(
            in_target,
            target_label,
            np.where(cluster_labels == target_label, "cluster_bg", cluster_labels),
        )
    elif domain != "voronoi":
        raise ValueError(domain)

    gene_names = [f"GENE_{i + 1:05d}" for i in range(n_genes)]
    sig_indices = rng.choice(n_genes, n_signature_genes, replace=False)
    fold_changes = rng.uniform(
        fold_change_range[0], fold_change_range[1], n_signature_genes
    )
    order = np.argsort(-fold_changes)
    sig_indices = sig_indices[order]
    fold_changes = fold_changes[order]
    sig_genes = [gene_names[i] for i in sig_indices]

    gene_means = np.exp(rng.normal(1.5, 1.0, n_genes))
    gene_means = np.clip(gene_means, 0.1, 100)
    dispersions = rng.uniform(dispersion_range[0], dispersion_range[1], n_genes)
    lib_size_factors = np.exp(rng.normal(0, 0.3, n_spots))
    lib_size_factors *= mean_counts_per_spot / (
        gene_means.sum() * lib_size_factors.mean()
    )

    # technical field (independent random stream)
    log_field = np.zeros(n_spots)
    affected = np.zeros(n_genes, bool)
    if tech_field is not None and tech_strength > 0:
        u = coords / capture_area_um
        if tech_field == "gradient":
            ang = rng_ext.uniform(0, 2 * np.pi)
            proj = u @ np.array([np.cos(ang), np.sin(ang)])
            base = (proj - proj.min()) / (proj.max() - proj.min())
        elif tech_field == "edge":
            d_edge = np.minimum.reduce([u[:, 0], 1 - u[:, 0], u[:, 1], 1 - u[:, 1]])
            base = 1 - np.exp(-d_edge / 0.12)
        elif tech_field == "grf":
            # random Fourier features of an RBF-kernel Gaussian field, length scale tech_scale x side
            K = 200
            freqs = rng_ext.normal(0, 1.0 / tech_scale, size=(K, 2))
            phases = rng_ext.uniform(0, 2 * np.pi, K)
            f = np.cos(u @ freqs.T + phases).sum(axis=1)
            base = (f - f.min()) / (f.max() - f.min())
        else:
            raise ValueError(tech_field)
        log_field = tech_strength * (base - base.mean())
        affected = rng_ext.random(n_genes) < tech_gene_frac
        if tech_on_signature is True:
            affected[sig_indices] = True
        elif tech_on_signature is False:
            affected[sig_indices] = False
    field = np.exp(log_field)

    X = np.zeros((n_spots, n_genes), dtype=np.float32)
    target_mask = cluster_labels == target_label
    for g in range(n_genes):
        mu = gene_means[g] * lib_size_factors
        if g in sig_indices:
            fc_idx = np.where(sig_indices == g)[0][0]
            mu[target_mask] *= fold_changes[fc_idx]
        if affected[g]:
            mu = mu * field
        r = dispersions[g]
        p = r / (r + mu)
        counts = rng.negative_binomial(r, p)
        dropout_mask = rng.random(n_spots) < dropout_rate
        counts[dropout_mask] = 0
        X[:, g] = counts

    adata = ad.AnnData(X=csr_matrix(X))
    adata.var_names = gene_names
    adata.obs_names = [f"spot_{i}" for i in range(n_spots)]
    adata.obs["cluster"] = pd.Categorical(cluster_labels)
    adata.obsm["spatial"] = coords
    adata.layers["counts"] = adata.X.copy()
    adata.obs["total_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["tech_log_field"] = log_field
    adata.var["tech_affected"] = affected
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    if return_truth:
        return (
            adata,
            sig_genes,
            fold_changes,
            target_label,
            {"log_field": log_field, "affected": affected},
        )
    return adata, sig_genes, fold_changes, target_label


def f2c_spatial_conn(coords, n_neighbors=6):
    """
    Binary spatial kNN connectivity, built as inside em.tl.score.

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates (n x 2).
    n_neighbors : int
        Number of neighbours.

    Returns
    -------
    scipy.sparse.csr_matrix
        Connectivity matrix.
    """
    tmp = ad.AnnData(X=np.zeros((len(coords), 1), dtype=np.float32))
    tmp.obsm["spatial"] = np.asarray(coords)
    sq.gr.spatial_neighbors(
        tmp, n_neighs=n_neighbors, coord_type="generic", key_added="spatial"
    )
    return csr_matrix(tmp.obsp["spatial_connectivities"])


def f2c_smooth(Y, conn, n_iter=1):
    """
    Neighbourhood average with the row-normalised adjacency (EnrichMap Eq. 3).

    Parameters
    ----------
    Y : np.ndarray
        Values per spot (vector or spots x columns).
    conn : scipy.sparse matrix
        Spatial connectivity.
    n_iter : int
        Number of smoothing iterations.

    Returns
    -------
    np.ndarray
        Smoothed values.
    """
    Y = np.asarray(Y, dtype=np.float64)
    d = np.maximum(np.asarray(conn.sum(axis=1)).ravel(), 1e-10)
    out = Y
    for _ in range(n_iter):
        out = conn.dot(out) / (d[:, None] if out.ndim == 2 else d)
    return out


def f2c_morans_i(values, W):
    """
    Moran's I with a binary spatial weight matrix.

    Parameters
    ----------
    values : np.ndarray
        Values per spot.
    W : scipy.sparse matrix
        Spatial weights.

    Returns
    -------
    float
        Moran's I (0 for a constant vector).
    """
    v = np.asarray(values, dtype=np.float64)
    z = v - v.mean()
    z2 = z @ z
    if z2 == 0:
        return 0.0
    return float((len(v) / W.sum()) * (z @ W.dot(z)) / z2)


def f2c_raw_enrichmap(adata, genes, weights=None):
    """
    Raw EnrichMap score (Eq. 1, CV-weighted mean) before smoothing and GAM.

    Parameters
    ----------
    adata : AnnData
        Log-normalised data.
    genes : list of str
        Signature genes.
    weights : dict, optional
        Gene weights (default: em.tl.infer_gene_weights).

    Returns
    -------
    np.ndarray
        Score per spot.
    """
    genes = [g for g in dict.fromkeys(genes) if g in adata.var_names]
    w = weights if weights is not None else em.tl.infer_gene_weights(adata, genes)
    wv = np.array([w.get(g, 1.0) for g in genes], dtype=np.float64)
    X = adata[:, genes].X
    X = X.toarray() if issparse(X) else np.asarray(X)
    return X @ wv / np.abs(wv).sum()


def f2c_raw_zscore(adata, genes):
    """
    Raw z-score scorer: sum of the signature genes / sqrt(number of genes).

    Parameters
    ----------
    adata : AnnData
        Log-normalised data.
    genes : list of str
        Signature genes.

    Returns
    -------
    np.ndarray
        Score per spot.
    """
    genes = [g for g in dict.fromkeys(genes) if g in adata.var_names]
    X = adata[:, genes].X
    X = X.toarray() if issparse(X) else np.asarray(X)
    return X.sum(axis=1) / np.sqrt(len(genes))


def f2c_f1(y, pred):
    """
    F1 score of binary calls.

    Parameters
    ----------
    y : np.ndarray
        True labels.
    pred : np.ndarray
        Binary calls.

    Returns
    -------
    float
        F1 (0 when there are no positives and no calls).
    """
    y = np.asarray(y).astype(bool)
    pred = np.asarray(pred).astype(bool)
    tp = (y & pred).sum()
    denom = pred.sum() + y.sum()
    return float(2 * tp / denom) if denom else 0.0


def f2c_spearman(a, b):
    """
    Spearman correlation, NaN when either vector is constant.

    Parameters
    ----------
    a, b : np.ndarray
        Values.

    Returns
    -------
    float
        Spearman rho.
    """
    if np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(spearmanr(a, b).correlation)


def f2c_target_metrics(s, y, log_field, clean_ref=None, clean_contrast=None):
    """
    F2c metrics of one score against the true domain and the technical field.

    Parameters
    ----------
    s : np.ndarray
        Score per spot.
    y : np.ndarray
        Target-domain indicator.
    log_field : np.ndarray
        log F per spot.
    clean_ref : np.ndarray, optional
        Uncorrected score of the artefact-free twin (rho_clean).
    clean_contrast : float, optional
        In-vs-out contrast of the artefact-free twin (retention).

    Returns
    -------
    dict
        AUROC, AUPRC, F1 top-k, contrast, rho_bg_field, rho_all_field (and retention, rho_clean).
    """
    y = np.asarray(y, bool)
    bg = ~y
    contrast = float(s[y].mean() - s[bg].mean())
    k = int(y.sum())
    top = np.argsort(-s, kind="stable")[:k]
    m = {
        "AUROC": roc_auc_score(y, s),
        "AUPRC": average_precision_score(y, s),
        "F1 top-k": float(y[top].mean()),
        "contrast": contrast,
        "rho_bg_field": f2c_spearman(s[bg], log_field[bg]),
        "rho_all_field": f2c_spearman(s, log_field),
    }
    if clean_contrast is not None:
        m["retention"] = contrast / clean_contrast
    if clean_ref is not None:
        m["rho_clean"] = f2c_spearman(s, clean_ref)
    return m


def f2c_null_metrics(s, log_field, W):
    """
    Spatial structure of a null-signature score and its correlation with the technical field.

    Parameters
    ----------
    s : np.ndarray
        Null-signature score per spot.
    log_field : np.ndarray
        log F per spot.
    W : scipy.sparse matrix
        6-NN connectivity (Moran's I weights).

    Returns
    -------
    dict
        null_moran, null_abs_rho_field, null_rho_field.
    """
    return {
        "null_moran": f2c_morans_i(s, W),
        "null_abs_rho_field": abs(f2c_spearman(s, log_field))
        if np.std(log_field) > 0
        else np.nan,
        "null_rho_field": f2c_spearman(s, log_field),
    }


def f2c_null_genes_sim(adata, sig, seed):
    """
    The 20-random-gene null signature of a simulated slide (random stream of F2c_compute.py).

    Parameters
    ----------
    adata : AnnData
        Simulated slide.
    sig : list of str
        Signature genes (excluded).
    seed : int
        Simulation seed.

    Returns
    -------
    list of str
        Null-signature genes.
    """
    gidx = {g: i for i, g in enumerate(adata.var_names)}
    rng_null = np.random.default_rng(777 + seed)
    pool = np.setdiff1d(np.arange(adata.n_vars), [gidx[g] for g in sig])
    return list(adata.var_names[rng_null.choice(pool, len(sig), replace=False)])


def f2c_null_model_gene_args(gene_sets, scorer):
    """
    em.tl.score gene-set arguments of the two F2c scorers.

    Parameters
    ----------
    gene_sets : dict
        Score key -> genes.
    scorer : str
        "EnrichMap" (CV gene weights, the package default) or "z-score" (uniform gene weights).

    Returns
    -------
    dict
        gene_set=... or gene_weights=...
    """
    if scorer == "z-score":
        return dict(gene_weights={k: {g: 1.0 for g in v} for k, v in gene_sets.items()})
    return dict(gene_set={k: list(v) for k, v in gene_sets.items()})


def f2c_null_model_output(adata, key, genes, scorer):
    """
    Null-model z of one gene set, its numerator in score units and its empirical p.

    Parameters
    ----------
    adata : AnnData
        Data scored with em.tl.score(null_model=True).
    key : str
        Score key.
    genes : list of str
        Genes of the set.
    scorer : str
        "EnrichMap" or "z-score" (the z-score unit of the F2c raw scorer is sum / sqrt(n), the uniform-weight
        mean is sum / n).

    Returns
    -------
    dict
        z, num (z x null SD in score units), p (p_empirical), overlap, offset_ratio, unit.
    """
    s = adata.obs[f"{key}_score"].values.astype(np.float64)
    unit = (
        np.sqrt(len([g for g in dict.fromkeys(genes) if g in adata.var_names]))
        if scorer == "z-score"
        else 1.0
    )
    info = adata.uns["null_model"][key]
    return {
        "z": s,
        "num": s * info["null_sd"].astype(np.float64) * unit,
        "unit": unit,
        "p": info["p_empirical"].astype(np.float64),
        "overlap": float(info["control_overlap"]),
        "offset_ratio": float(info["level_offset_ratio"]),
    }


def f2c_null_model_target_metrics(
    o, y, log_field, clean_none, clean_contrast, clean_same
):
    """
    F2c target metrics of a null-model score (z for ranks, numerator for score-unit contrasts) plus z calls.

    Parameters
    ----------
    o : dict
        Output of f2c_null_model_output for the target signature.
    y : np.ndarray
        Target-domain indicator.
    log_field : np.ndarray
        log F per spot.
    clean_none : np.ndarray
        Uncorrected (smoothed raw) score of the artefact-free twin (rho_clean reference).
    clean_contrast : float
        In-vs-out contrast of `clean_none` (retention reference).
    clean_same : dict or None
        Same option on the artefact-free twin (rho_clean_same, retention_z); None for the twin itself.

    Returns
    -------
    dict
        Metric -> value.
    """
    y = np.asarray(y, bool)
    z, num, p = o["z"], o["num"], o["p"]
    m = f2c_target_metrics(z, y, log_field, clean_ref=clean_none)
    m["contrast_z"] = m.pop("contrast")
    m["contrast"] = float(num[y].mean() - num[~y].mean())
    m["retention"] = m["contrast"] / clean_contrast
    ref = clean_same if clean_same is not None else o
    m["rho_clean_same"] = f2c_spearman(z, ref["z"])
    m["retention_z"] = m["contrast_z"] / float(ref["z"][y].mean() - ref["z"][~y].mean())
    m["F1 z>1.96"] = f2c_f1(y, z > 1.96)
    m["frac_tgt_z196"] = float((z[y] > 1.96).mean())
    m["frac_bg_z196"] = float((z[~y] > 1.96).mean())
    m["frac_bg_zm196"] = float((z[~y] < -1.96).mean())
    m["frac_bg_p05"] = float((p[~y] <= 0.05).mean())
    m["control_overlap"] = o["overlap"]
    m["level_offset_ratio"] = o["offset_ratio"]
    return m


def f2c_null_model_signature_metrics(o, log_field, conn):
    """
    F2c null-signature metrics of a null-model score plus its call rates.

    Parameters
    ----------
    o : dict
        Output of f2c_null_model_output for the null signature.
    log_field : np.ndarray
        log F per spot.
    conn : scipy.sparse matrix
        6-NN connectivity.

    Returns
    -------
    dict
        null_moran, null_abs_rho_field, null_rho_field, null_frac_z196, null_frac_absz196, null_frac_p05.
    """
    m = f2c_null_metrics(o["z"], log_field, conn)
    m["null_frac_z196"] = float((o["z"] > 1.96).mean())
    m["null_frac_absz196"] = float((np.abs(o["z"]) > 1.96).mean())
    m["null_frac_p05"] = float((o["p"] <= 0.05).mean())
    return m


def f2c_read_gsdensity(path, index):
    """
    GSDensity relevance and native positive labels of one gene set, in spot order.

    Parameters
    ----------
    path : str
        CSV with columns spot_id, relevance and label.
    index : array-like
        Spot names in the order wanted.

    Returns
    -------
    tuple of np.ndarray
        Relevance per spot and boolean positive label per spot.
    """
    d = pd.read_csv(path).set_index("spot_id").reindex(index)
    return d["relevance"].values.astype(np.float64), d["label"].astype(
        str
    ).values == "positive"


def f2c_pass_through_slope(delta, log_field, ref_sd):
    """
    Artefact pass-through slope: d(score change) / d(log F), in units of a reference SD.

    Parameters
    ----------
    delta : np.ndarray
        Spiked minus clean score per spot (same correction).
    log_field : np.ndarray
        Spiked-in log F per spot.
    ref_sd : float
        SD of the clean reference score.

    Returns
    -------
    float
        Least-squares slope of delta on log F divided by `ref_sd`.
    """
    lfc = log_field - log_field.mean()
    return float((lfc @ (delta - delta.mean())) / (lfc @ lfc) / ref_sd)


def f2c_lognorm(counts, var_names, obs_names):
    """
    AnnData of counts, normalised to 1e4 per spot and log1p-transformed (counts kept in layers).

    Parameters
    ----------
    counts : scipy.sparse matrix
        Counts (spots x genes).
    var_names, obs_names : array-like
        Gene and spot names.

    Returns
    -------
    AnnData
        Log-normalised data.
    """
    a = ad.AnnData(
        X=counts.copy(),
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=var_names),
    )
    a.layers["counts"] = counts.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a


def f2c_brain_base_field(kind, coords, rng):
    """
    Smooth field shape in [0, 1] on the tissue ("edge" uses the distance to the tissue border).

    Parameters
    ----------
    kind : str
        "gradient", "edge" or "GRF-0.5".
    coords : np.ndarray
        Spot coordinates.
    rng : np.random.Generator
        Random stream (gradient direction, GRF features).

    Returns
    -------
    np.ndarray
        Field shape per spot.
    """
    L = np.ptp(coords, axis=0).max()
    u = (coords - coords.min(axis=0)) / L
    if kind == "gradient":
        ang = rng.uniform(0, 2 * np.pi)
        proj = u @ np.array([np.cos(ang), np.sin(ang)])
        base = (proj - proj.min()) / np.ptp(proj)
    elif kind == "edge":
        tree = cKDTree(coords)
        dnn, _ = tree.query(coords, k=7)
        spacing = np.median(dnn[:, 1])
        border = (dnn[:, 1:] < 1.2 * spacing).sum(axis=1) < 6
        d_edge = cKDTree(coords[border]).query(coords)[0] / L
        base = 1 - np.exp(-d_edge / 0.12)
        base = (base - base.min()) / np.ptp(base)
    elif kind == "GRF-0.5":
        K = 200
        freqs = rng.normal(0, 1.0 / 0.5, size=(K, 2))
        phases = rng.uniform(0, 2 * np.pi, K)
        f = np.cos(u @ freqs.T + phases).sum(axis=1)
        base = (f - f.min()) / np.ptp(f)
    else:
        raise ValueError(kind)
    return base


def f2c_spike(counts, base, strength, affected, rng):
    """
    Multiplicative capture field exp(strength x (base - mean)) spiked in by binomial thinning.

    Counts of affected genes are thinned with retention p = exp(strength x (base - 1)); the global constant is
    removed by normalize_total.

    Parameters
    ----------
    counts : scipy.sparse matrix
        Observed counts (spots x genes).
    base : np.ndarray
        Field shape per spot in [0, 1].
    strength : float
        Range of log F (natural log).
    affected : np.ndarray
        Boolean per gene.
    rng : np.random.Generator
        Random stream of the thinning.

    Returns
    -------
    scipy.sparse.csr_matrix
        Spiked counts.
    """
    p = np.exp(strength * (base - 1.0))
    C = counts.tocsr(copy=True)
    rows = np.repeat(np.arange(C.shape[0]), np.diff(C.indptr))
    sel = affected[C.indices]
    C.data[sel] = rng.binomial(C.data[sel].astype(np.int64), p[rows[sel]]).astype(
        np.float32
    )
    C.eliminate_zeros()
    return C


def f2c_brain_null_genes(adata, genes, rep):
    """
    The 10-gene brain null signature of one replicate (genes detected in >= 10% of spots, no signature genes).

    Parameters
    ----------
    adata : AnnData
        Mouse brain with layers["counts"].
    genes : dict
        Region -> signature genes.
    rep : int
        Replicate.

    Returns
    -------
    list of str
        Null-signature genes.
    """
    det = np.asarray((adata.layers["counts"] > 0).mean(axis=0)).ravel()
    sig = set(sum(genes.values(), []))
    pool = [g for g, d in zip(adata.var_names, det) if d >= 0.10 and g not in sig]
    return list(np.random.default_rng(5000 + rep).choice(pool, 10, replace=False))


def f2c_ci95(x):
    """
    Half-width of the normal-approximation 95% CI of a mean.

    Parameters
    ----------
    x : array-like
        Values (NaN ignored).

    Returns
    -------
    float
        1.96 x SD / sqrt(n), NaN when fewer than two values.
    """
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    return 1.96 * x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan


def f2c_plot_phase(
    summ,
    metric,
    scorer,
    title,
    corr,
    labels,
    domains,
    cond_labels,
    cond_parts,
    scorer_labels,
    missing_colour,
    ci_col=None,
    cmap="RdBu",
    vmax=0.5,
    lin=0.01,
    diverging=True,
    panel_titles=None,
):
    """
    Domain x condition heat map of one summary metric, one panel per correction.

    Parameters
    ----------
    summ : pd.DataFrame
        Summary per scorer, domain, condition and correction.
    metric : str
        Column of `summ` to plot.
    scorer : str
        Scorer (EnrichMap or z-score).
    title : str
        Colour-bar label.
    corr : list of str
        Corrections, one panel each.
    labels : dict
        Correction -> panel title.
    domains : list of str
        Domain geometries (rows).
    cond_labels : list of str
        Condition labels (columns).
    cond_parts : dict
        Condition label -> (arm, field, strength).
    scorer_labels : dict
        Scorer -> label.
    missing_colour : str
        Colour of cells that were not run.
    ci_col : str, optional
        Column with the paired 95% CI half-width; cells are marked when it excludes 0.
    cmap : str
        Colour map.
    vmax : float
        Upper colour limit (and -vmax for diverging scales).
    lin : float
        Linear range of the symmetric-log scale.
    diverging : bool
        Symmetric-log diverging scale (True) or linear 0-vmax scale (False).
    panel_titles : dict, optional
        Panel titles replacing `labels`.

    Returns
    -------
    tuple
        (figure, tidy table of the plotted cells of every panel).
    """
    titles = {**labels, **(panel_titles or {})}
    n_rows = int(np.ceil(len(corr) / 2))
    fig, axes = plt.subplots(
        n_rows,
        2,
        figsize=(7.2, 1.25 * n_rows + 1.2),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    norm = (
        SymLogNorm(linthresh=lin, linscale=0.6, vmin=-vmax, vmax=vmax, base=10)
        if diverging
        else plt.Normalize(0, vmax)
    )
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad(missing_colour)
    cells = pd.MultiIndex.from_product([domains, cond_labels], names=["domain", "cond"])
    cols = (
        ["domain", "cond", "arm", "field", "strength", metric]
        + ([ci_col] if ci_col else [])
        + ["n_seeds"]
    )
    tables = []
    for ax, c in zip(axes.flat, corr):
        d = summ[(summ.scorer == scorer) & (summ.correction == c)]
        M = d.pivot(index="domain", columns="cond", values=metric).reindex(
            index=domains, columns=cond_labels
        )
        im = ax.imshow(
            np.ma.masked_invalid(M.values), cmap=cm, norm=norm, aspect="auto"
        )
        t = d.set_index(["domain", "cond"]).reindex(cells).reset_index()
        for j, part in enumerate(("arm", "field", "strength")):
            t[part] = [cond_parts[x][j] for x in t.cond]
        t = t[cols]
        if ci_col is not None:
            C = d.pivot(index="domain", columns="cond", values=ci_col).reindex(
                index=domains, columns=cond_labels
            )
            for i in range(M.shape[0]):
                for j in range(M.shape[1]):
                    m, ci = M.values[i, j], C.values[i, j]
                    if abs(m) >= 0.01 and abs(m) > ci:
                        ax.text(
                            j,
                            i,
                            "+" if m > 0 else "−",
                            ha="center",
                            va="center",
                            fontsize=5,
                            color="white" if abs(m) > 0.1 else "black",
                        )
            # same rule as the markers
            sig = (t[metric].abs() >= 0.01) & (t[metric].abs() > t[ci_col])
            t["ci95_low"] = t[metric] - t[ci_col]
            t["ci95_high"] = t[metric] + t[ci_col]
            t["marker"] = np.where(sig, np.where(t[metric] > 0, "+", "−"), "")
        t.insert(0, "correction", c)
        tables.append(t)
        ax.set_title(titles[c], fontsize=7, pad=2)
        ax.set_yticks(range(len(domains)))
        ax.set_yticklabels(domains, fontsize=6)
        for x in (0.5, 12.5, 15.5):
            ax.axvline(x, color="k", lw=0.5)
        ax.grid(False)
    for ax in axes.flat[len(corr) :]:
        ax.set_visible(False)
    # condition labels under the last panel of each column
    for col in range(2):
        ax = axes.flat[max(range(col, len(corr), 2))]
        ax.tick_params(labelbottom=True)
        ax.set_xticks(range(len(cond_labels)))
        ax.set_xticklabels(
            [l.replace("|", " ") for l in cond_labels], rotation=90, fontsize=5
        )
    cbar = fig.colorbar(im, ax=axes, shrink=0.5, pad=0.01)
    cbar.set_label(title, fontsize=6)
    cbar.ax.tick_params(labelsize=5)
    fig.suptitle(f"{scorer_labels[scorer]}: {title}", fontsize=8)
    return fig, pd.concat(tables, ignore_index=True)


def f2c_point_table(m, correction, x_name):
    """
    Tidy table of one plotted error-bar series.

    Parameters
    ----------
    m : pd.DataFrame
        Indexed by the x variable, with columns mean, ci95 and n.
    correction : str
        Correction of the series.
    x_name : str
        Name of the x variable.

    Returns
    -------
    pd.DataFrame
        One row per point with mean, 95% CI bounds and n.
    """
    t = m.rename_axis(x_name).reset_index()
    t.insert(0, "correction", correction)
    t["ci95_low"] = t["mean"] - t["ci95"]
    t["ci95_high"] = t["mean"] + t["ci95"]
    return t


def f2c_errorbar(ax, m, c, style, x=None, ms=2.5, lw=0.9):
    """
    One correction's mean ± 95% CI series in its palette colour, marker and line style.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    m : pd.DataFrame
        Indexed by x, with columns mean and ci95.
    c : str
        Correction.
    style : dict
        "colours", "markers", "lines" and "labels" per correction.
    x : array-like, optional
        x values (default: the index of `m`).
    ms, lw : float
        Marker size and line width.
    """
    ax.errorbar(
        m.index if x is None else x,
        m["mean"],
        yerr=m["ci95"],
        color=style["colours"][c],
        marker=style["markers"][c],
        ls=style["lines"][c],
        ms=ms * (1.6 if c == "GSDensity" else 1.0),
        mew=0,
        lw=lw * (1.4 if c in ("null model", "GSDensity") else 1.0),
        capsize=1,
        label=style["labels"][c],
        zorder=4
        if c in ("null model", "GSDensity", "none", "GAM default (lam 0.6)")
        else 3,
    )


def f2c_legend_handles(corr, style):
    """
    Legend entries (colour, marker, line style) of the corrections.

    Parameters
    ----------
    corr : list of str
        Corrections in legend order.
    style : dict
        "colours", "markers", "lines" and "labels" per correction.

    Returns
    -------
    list of Line2D
    """
    return [
        Line2D(
            [],
            [],
            color=style["colours"][c],
            marker=style["markers"][c],
            ls=style["lines"][c],
            lw=0.9,
            ms=4 if c != "GSDensity" else 5.5,
            markeredgecolor="k",
            markeredgewidth=0.2,
            label=style["labels"][c],
        )
        for c in corr
    ]


def f2c_wide_table(
    summ, metric, corr_all, short, domains, cond_labels, scorer="EnrichMap", arms=None
):
    """
    Domain x condition rows of one summary metric, one column per correction.

    Parameters
    ----------
    summ : pd.DataFrame
        Summary per scorer, domain, condition and correction.
    metric : str
        Column of `summ`.
    corr_all : list of str
        Corrections (column order).
    short : dict
        Correction -> short name.
    domains : list of str
        Domain geometries.
    cond_labels : list of str
        Condition labels.
    scorer : str
        Scorer.
    arms : list of str, optional
        Restrict to these artefact arms.

    Returns
    -------
    pd.DataFrame
        Wide table with short correction names.
    """
    d = summ[(summ.scorer == scorer) & (summ.correction.isin(corr_all))]
    if arms is not None:
        d = d[d.arm.isin(arms)]
    t = d.pivot_table(
        index=["domain", "cond"], columns="correction", values=metric
    ).reindex(columns=corr_all)
    t = t.reindex(
        [(dm, c) for dm in domains for c in cond_labels if (dm, c) in t.index]
    )
    return t.rename(columns=short)


def f2c_tradeoff_table(d):
    """
    Confounding removed and signal retained per scorer, domain and correction.

    Parameters
    ----------
    d : pd.DataFrame
        Per-seed rows (artefact conditions only), including the uncorrected ("none") rows.

    Returns
    -------
    pd.DataFrame
        One row per scorer, domain and correction.
    """
    rows = []
    for (scn, dm), dd in d.groupby(["scorer", "domain"]):
        ref = dd[dd.correction == "none"].abs_rho_bg.mean()
        for c, x in dd.groupby("correction"):
            rows.append(
                {
                    "scorer": scn,
                    "domain": dm,
                    "correction": c,
                    "confounding_removed": 1 - x.abs_rho_bg.mean() / ref,
                    "residual_abs_rho_bg": x.abs_rho_bg.mean(),
                    "retention": x.retention.mean(),
                    "retention_sd_over_seeds": x.groupby("seed").retention.mean().std(),
                    "dAUROC": x.d_AUROC.mean(),
                    "null_moran": x.null_moran.mean(),
                    "n_simulations": len(x),
                }
            )
    return pd.DataFrame(rows)


def f2c_domain_values_text(summ, correction, cond, domains, metric="AUROC_mean"):
    """
    One correction's value per domain in one condition, as text.

    Parameters
    ----------
    summ : pd.DataFrame
        Summary of one scorer.
    correction : str
        Correction.
    cond : str
        Condition label.
    domains : list of str
        Domain geometries.
    metric : str
        Column of `summ`.

    Returns
    -------
    str
    """
    d = (
        summ[(summ.correction == correction) & (summ.cond == cond)]
        .set_index("domain")[metric]
        .reindex(domains)
    )
    return ", ".join(f"{dm} {v:.3f}" for dm, v in d.items())


def f2c_plot_metric_pair(
    table, x, y, ax, colours, ref_x=None, ref_y=None, invert_y=False, marker_col=None
):
    """
    One metric against another, one point per method, with 95% CI whiskers from the _lo/_hi columns.

    Parameters
    ----------
    table : pd.DataFrame
        One row per method (and condition), with a "method" column and the metric columns.
    x, y : str
        Metric columns.
    ax : matplotlib.axes.Axes
        Target axes.
    colours : dict
        Method -> colour.
    ref_x, ref_y : float, optional
        Dotted reference lines.
    invert_y : bool
        Invert the y axis (lower-is-better metrics).
    marker_col : str, optional
        Column holding a marker per row (default circles).

    Returns
    -------
    matplotlib.axes.Axes
    """
    for _, r in table.iterrows():
        xerr = yerr = None
        if f"{x}_lo" in table:
            xerr = [[r[x] - r[f"{x}_lo"]], [r[f"{x}_hi"] - r[x]]]
        if f"{y}_lo" in table:
            yerr = [[r[y] - r[f"{y}_lo"]], [r[f"{y}_hi"] - r[y]]]
        ax.errorbar(
            r[x],
            r[y],
            xerr=xerr,
            yerr=yerr,
            fmt=r[marker_col] if marker_col else "o",
            color=colours[r["method"]],
            ms=5,
            mec="white",
            mew=0.4,
            elinewidth=0.6,
            capsize=0,
            label=r["method"],
        )
    if ref_x is not None:
        ax.axvline(ref_x, ls=":", c="0.6", lw=0.6)
    if ref_y is not None:
        ax.axhline(ref_y, ls=":", c="0.6", lw=0.6)
    if invert_y:
        ax.invert_yaxis()
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.grid(False)
    return ax


def f2c_plot_benchmark_summary(table, ax, colours, lower_is_better=(), annotate=True):
    """
    Methods x metrics summary: circles sized and shaded by each metric scaled to 0-1 across methods (best = 1),
    raw values annotated, plus a mean-rank column.

    Parameters
    ----------
    table : pd.DataFrame
        Index = methods, columns = metrics.
    ax : matplotlib.axes.Axes
        Target axes.
    colours : dict
        Method -> colour (tick labels).
    lower_is_better : tuple of str
        Metrics where smaller is better.
    annotate : bool
        Print the raw value in each cell.

    Returns
    -------
    matplotlib.axes.Axes
    """
    scaled = table.copy().astype(float)
    for c in table.columns:
        col = table[c].astype(float)
        rng_ = col.max() - col.min()
        s = (col - col.min()) / rng_ if rng_ > 0 else col * 0 + 1
        scaled[c] = 1 - s if c in lower_is_better else s
    ranks = scaled.rank(ascending=False).mean(axis=1)
    cmap = plt.get_cmap("Purples")
    for i, m in enumerate(table.index):
        for j, c in enumerate(table.columns):
            s = scaled.loc[m, c]
            ax.scatter(
                j, i, s=30 + 220 * s, color=cmap(0.25 + 0.75 * s), edgecolors="none"
            )
            if annotate:
                ax.text(
                    j,
                    i,
                    f"{table.loc[m, c]:.2g}",
                    ha="center",
                    va="center",
                    fontsize=5,
                    color="white" if s > 0.6 else "black",
                )
        ax.text(
            len(table.columns) + 0.2,
            i,
            f"{ranks[m]:.1f}",
            ha="left",
            va="center",
            fontsize=6,
        )
    ax.set_xticks(range(len(table.columns)))
    ax.set_xticklabels(table.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(table.index)))
    ax.set_yticklabels(table.index)
    for lab in ax.get_yticklabels():
        lab.set_color(colours[lab.get_text()])
    ax.text(len(table.columns) + 0.2, -0.9, "mean rank", fontsize=6, ha="left")
    ax.set_xlim(-0.6, len(table.columns) + 0.9)
    ax.set_ylim(len(table.index) - 0.4, -1.2)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    return ax


# --------------------------------------------------------------------------------------------------
# R4_F2d_GAM_mechanism
# --------------------------------------------------------------------------------------------------


def f2d_spatial_conn(coords, n_neighbors=6):
    """
    Binary spatial kNN connectivity, built as inside em.tl.score.

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates (n_spots x 2).
    n_neighbors : int
        Number of spatial neighbours.

    Returns
    -------
    scipy.sparse.csr_matrix
        Binary connectivity (n_spots x n_spots).
    """
    tmp = anndata.AnnData(X=np.zeros((len(coords), 1), dtype=np.float32))
    tmp.obsm["spatial"] = np.asarray(coords)
    sq.gr.spatial_neighbors(
        tmp, n_neighs=n_neighbors, coord_type="generic", key_added="spatial"
    )
    return csr_matrix(tmp.obsp["spatial_connectivities"])


def f2d_smooth(Y, conn):
    """
    Neighbourhood average with the row-normalised adjacency (EnrichMap Eq. 3).

    Parameters
    ----------
    Y : np.ndarray
        Values (n_spots,) or (n_spots x n_columns).
    conn : scipy.sparse matrix
        Binary spatial connectivity.

    Returns
    -------
    np.ndarray
        Smoothed values, same shape as `Y`.
    """
    Y = np.asarray(Y, dtype=np.float64)
    d = np.maximum(np.asarray(conn.sum(axis=1)).ravel(), 1e-10)
    return conn.dot(Y) / (d[:, None] if Y.ndim == 2 else d)


def f2d_gam_smoother(coords, basis="te", n_splines=None, lam=None):
    """
    The spatial-covariate GAM of EnrichMap written as a linear smoother, H = B (B'B + P)^-1 B'.

    With basis="te" and n_splines = lam = None this is LinearGAM(te(0, 1)) as called by em.tl.score
    (pyGAM defaults: 10 x 10 cubic B-splines, lam = 0.6 per margin, intercept); basis="additive" gives
    s(0) + s(1) (pyGAM default 20 splines per margin). The sqrt(machine-eps) ridge that pyGAM adds before
    solving is added to B'B + P, so that the edof matches pyGAM's statistics_["edof"].

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates (n_spots x 2).
    basis : str
        "te" (tensor product) or "additive".
    n_splines : int or None
        Splines per margin (None: pyGAM default).
    lam : float or None
        Penalty per margin (None: pyGAM default).

    Returns
    -------
    dict
        Model matrix "B", penalised normal matrix "A", effective degrees of freedom "edof" (trace of H)
        and the pyGAM "terms".
    """
    coords = np.asarray(coords, dtype=np.float64)
    kw = {}
    if n_splines is not None:
        kw["n_splines"] = n_splines
    if lam is not None:
        kw["lam"] = lam
    term = te(0, 1, **kw) if basis == "te" else s(0, **kw) + s(1, **kw)
    gam = LinearGAM(term).fit(coords, coords[:, 0])
    B = gam._modelmat(coords)
    B = B.toarray() if hasattr(B, "toarray") else np.asarray(B)
    P = gam._P()
    P = P.toarray() if hasattr(P, "toarray") else np.asarray(P)
    A = B.T @ B + P + np.sqrt(np.finfo(float).eps) * np.eye(B.shape[1])
    edof = float(np.trace(np.linalg.solve(A, B.T @ B)))
    return {"B": B, "A": A, "edof": edof, "terms": str(gam.terms)}


def f2d_gam_fitted(gam, Y):
    """
    GAM fit H y.

    Parameters
    ----------
    gam : dict
        Output of `f2d_gam_smoother`.
    Y : np.ndarray
        Values (n_spots,) or (n_spots x n_columns).

    Returns
    -------
    np.ndarray
        Fitted values, same shape as `Y`.
    """
    Y = np.asarray(Y, dtype=np.float64)
    return gam["B"] @ np.linalg.solve(gam["A"], gam["B"].T @ Y)


def f2d_gam_residual(gam, Y):
    """
    GAM residual (I - H) y, the spatial-covariate correction of em.tl.score.

    Parameters
    ----------
    gam : dict
        Output of `f2d_gam_smoother`.
    Y : np.ndarray
        Values (n_spots,) or (n_spots x n_columns).

    Returns
    -------
    np.ndarray
        Residuals, same shape as `Y`.
    """
    return np.asarray(Y, dtype=np.float64) - f2d_gam_fitted(gam, Y)


def f2d_gam_r2(gam, y):
    """
    Share of the centred variance of y absorbed by the GAM, R2 = 1 - ||y - H y||^2 / ||y - mean(y)||^2.

    Parameters
    ----------
    gam : dict
        Output of `f2d_gam_smoother`.
    y : np.ndarray
        Values (n_spots,).

    Returns
    -------
    float
        Absorbed share.
    """
    y = np.asarray(y, dtype=np.float64)
    yc = y - y.mean()
    return float(1 - ((y - f2d_gam_fitted(gam, y)) ** 2).sum() / (yc**2).sum())


def f2d_raw_enrichmap(adata, genes):
    """
    EnrichMap Eq. 1: CV-weighted mean expression, before smoothing and the GAM.

    Parameters
    ----------
    adata : AnnData
        Log-normalised expression.
    genes : list
        Signature genes (genes missing from the data are skipped).

    Returns
    -------
    np.ndarray
        Raw score (n_spots,).
    """
    genes = [g for g in dict.fromkeys(genes) if g in adata.var_names]
    w = em.tl.infer_gene_weights(adata, genes)
    wv = np.array([w.get(g, 1.0) for g in genes], dtype=np.float64)
    X = adata[:, genes].X
    X = X.toarray() if issparse(X) else np.asarray(X)
    return X @ wv / np.abs(wv).sum()


def f2d_mode_energy(y, evecs):
    """
    Energy of a centred, unit-norm indicator in each graph-Laplacian mode.

    Parameters
    ----------
    y : np.ndarray
        Region indicator (n_spots,).
    evecs : np.ndarray
        Laplacian eigenvectors (n_spots x n_modes), ascending eigenvalue.

    Returns
    -------
    np.ndarray
        Squared projection on every mode (sums to 1).
    """
    yc = (y - y.mean()) / np.linalg.norm(y - y.mean())
    return (evecs.T @ yc) ** 2


def f2d_region_depth(conn, mask, max_hops=60):
    """
    Mean graph distance (neighbour rings) of in-region spots to the nearest spot outside the region.

    1 = every in-region spot touches the boundary (a one-spot-thick layer).

    Parameters
    ----------
    conn : scipy.sparse matrix
        Spatial connectivity (symmetrised inside).
    mask : np.ndarray
        Region membership (n_spots,), bool.
    max_hops : int
        Largest number of rings explored.

    Returns
    -------
    float
        Mean number of rings from an in-region spot to the region boundary.
    """
    mask = np.asarray(mask, bool)
    A = csr_matrix((conn + conn.T) > 0).astype(np.int8)
    outside = ~mask
    dist = np.full(len(mask), np.inf)
    dist[outside] = 0
    frontier = outside.copy()
    seen = outside.copy()
    for h in range(1, max_hops + 1):
        nxt = (A.dot(frontier.astype(np.int8)) > 0) & ~seen
        dist[nxt] = h
        seen |= nxt
        frontier = nxt
        if not nxt.any():
            break
    return float(np.mean(dist[mask]))


def f2d_spot_pitch(coords):
    """
    Median distance between a spot and its nearest neighbour.

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates (n_spots x 2).

    Returns
    -------
    float
        Spot pitch in coordinate units.
    """
    d, _ = cKDTree(coords).query(coords, k=2)
    return float(np.median(d[:, 1]))


def f2d_label_concordance(scores, labels, conn, iters=range(1, 11)):
    """
    Neighbourhood concordance of Fig. 2c / Sup. Fig. 1a (05_01 run_concordance).

    The ground-truth labels are smoothed k times with the row-normalised graph and Spearman-correlated
    with the unchanged scores.

    Parameters
    ----------
    scores : np.ndarray
        Score (n_spots,).
    labels : np.ndarray
        Region label (n_spots,), 0/1.
    conn : scipy.sparse matrix
        Binary spatial connectivity.
    iters : iterable of int
        Smoothing iterations k reported.

    Returns
    -------
    dict
        k -> Spearman rho (0 when undefined).
    """
    cur = np.asarray(labels, dtype=np.float64).copy()
    d = np.maximum(np.asarray(conn.sum(axis=1)).ravel(), 1e-10)
    out = {}
    for k in range(1, max(iters) + 1):
        cur = conn.dot(cur) / d
        if k in iters:
            rho = spearmanr(scores, cur).correlation
            out[k] = 0.0 if np.isnan(rho) else float(rho)
    return out


def f2d_classification_metrics(y, score):
    """
    AUROC, AUPRC and prevalence-matched (top-k) F1 of a score against a region label.

    Parameters
    ----------
    y : np.ndarray
        Region label (n_spots,), 0/1.
    score : np.ndarray
        Score (n_spots,).

    Returns
    -------
    dict
        "AUROC", "AUPRC" and "F1 top-k" (the k highest-scoring spots called, k = number of region spots).
    """
    y = np.asarray(y).astype(int)
    score = np.asarray(score, dtype=np.float64)
    pred = np.zeros(len(y), bool)
    pred[np.argsort(-score, kind="stable")[: int(y.sum())]] = True
    tp = (y.astype(bool) & pred).sum()
    return {
        "AUROC": roc_auc_score(y, score),
        "AUPRC": average_precision_score(y, score),
        "F1 top-k": float(2 * tp / (pred.sum() + y.sum())),
    }


def f2d_topk_call(score, k):
    """
    The k highest-scoring spots (stable order).

    Parameters
    ----------
    score : np.ndarray
        Score (n_spots,).
    k : int
        Number of spots to call (the number of region spots).

    Returns
    -------
    np.ndarray
        Boolean call (n_spots,).
    """
    call = np.zeros(len(score), bool)
    call[np.argsort(-np.asarray(score), kind="stable")[:k]] = True
    return call


def f2d_call_metrics(y, call):
    """
    Precision, recall and F1 of a binary call against the region label.

    Parameters
    ----------
    y : np.ndarray
        Region label (n_spots,), 0/1.
    call : np.ndarray
        Spots called positive (n_spots,), bool.

    Returns
    -------
    dict
        Number of spots called, precision, recall and F1 (0 when nothing is called).
    """
    y = np.asarray(y).astype(bool)
    call = np.asarray(call).astype(bool)
    tp = int((y & call).sum())
    precision = tp / call.sum() if call.sum() else 0.0
    recall = tp / y.sum()
    f1 = 2 * precision * recall / (precision + recall) if tp else 0.0
    return {
        "n called": int(call.sum()),
        "precision": precision,
        "recall": recall,
        "F1": f1,
    }


def f2d_hue_shades(base_colour, n, lightest=0.86, darkest=0.30):
    """
    `n` shades of the hue of `base_colour`, from dark to light (for sweeps of GAM settings).

    Parameters
    ----------
    base_colour : str
        Hex colour whose hue and saturation are kept.
    n : int
        Number of shades.
    lightest, darkest : float
        HLS lightness range.

    Returns
    -------
    list of str
        Hex colours.
    """
    h, _, sat = colorsys.rgb_to_hls(
        *(int(base_colour[i : i + 2], 16) / 255 for i in (1, 3, 5))
    )
    return [
        "#%02x%02x%02x"
        % tuple(round(v * 255) for v in colorsys.hls_to_rgb(h, lightness, sat))
        for lightness in np.linspace(darkest, lightest, n)
    ]


def f2d_plot_metric_pair(
    table, x, y, ax, colours, ref_x=None, ref_y=None, marker_col=None
):
    """
    One metric against another, one point per method, with 95% CI whiskers from f"{x}_lo"/f"{x}_hi" and
    f"{y}_lo"/f"{y}_hi".

    Parameters
    ----------
    table : pd.DataFrame
        One row per method, with a "method" column, the metric columns and their CI columns.
    x, y : str
        Metric columns, e.g. "AUROC" and "AUPRC".
    ax : matplotlib Axes
        Axes to draw on.
    colours : dict
        Method -> colour.
    ref_x, ref_y : float or None
        Dotted reference lines (e.g. AUROC 0.5; AUPRC = prevalence).
    marker_col : str or None
        Column holding a marker per row (default circles).

    Returns
    -------
    matplotlib Axes
        The axes.
    """
    for _, r in table.iterrows():
        xerr = [[r[x] - r[f"{x}_lo"]], [r[f"{x}_hi"] - r[x]]]
        yerr = [[r[y] - r[f"{y}_lo"]], [r[f"{y}_hi"] - r[y]]]
        ax.errorbar(
            r[x],
            r[y],
            xerr=xerr,
            yerr=yerr,
            fmt=r[marker_col] if marker_col else "o",
            color=colours[r["method"]],
            ms=5,
            mec="white",
            mew=0.4,
            elinewidth=0.6,
            capsize=0,
            label=r["method"],
        )
    if ref_x is not None:
        ax.axvline(ref_x, ls=":", c="0.6", lw=0.6)
    if ref_y is not None:
        ax.axhline(ref_y, ls=":", c="0.6", lw=0.6)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.grid(False)
    return ax


def f2d_plot_benchmark_summary(table, ax, colours):
    """
    Methods x metrics summary: circles sized and shaded by each metric scaled 0-1 across methods (best = 1),
    raw values printed with three decimals ("n/a" for a missing metric, drawn empty), plus a mean-rank column.

    Parameters
    ----------
    table : pd.DataFrame
        Index = methods, columns = metrics (higher = better).
    ax : matplotlib Axes
        Axes to draw on.
    colours : dict
        Method -> colour of its tick label.

    Returns
    -------
    matplotlib Axes
        The axes.
    """
    scaled = table.copy().astype(float)
    for c in table.columns:
        col = table[c].astype(float)
        span = col.max() - col.min()
        scaled[c] = (col - col.min()) / span if span > 0 else col * 0 + 1
    ranks = scaled.rank(ascending=False).mean(axis=1)
    cmap = plt.get_cmap("Purples")
    for i, m in enumerate(table.index):
        for j, c in enumerate(table.columns):
            sv = scaled.loc[m, c]
            v = table.loc[m, c]
            ax.scatter(
                j, i, s=30 + 220 * sv, color=cmap(0.25 + 0.75 * sv), edgecolors="none"
            )
            ax.text(
                j,
                i,
                "n/a" if np.isnan(v) else f"{v:.3f}",
                ha="center",
                va="center",
                fontsize=5,
                color="white" if sv > 0.6 else "black",
            )
        ax.text(
            len(table.columns) + 0.2,
            i,
            f"{ranks[m]:.1f}",
            ha="left",
            va="center",
            fontsize=6,
        )
    ax.set_xticks(range(len(table.columns)))
    ax.set_xticklabels(table.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(table.index)))
    ax.set_yticklabels(table.index)
    for lab in ax.get_yticklabels():
        lab.set_color(colours[lab.get_text()])
    ax.text(len(table.columns) + 0.2, -0.9, "mean rank", fontsize=6, ha="left")
    ax.set_xlim(-0.6, len(table.columns) + 0.9)
    ax.set_ylim(len(table.index) - 0.4, -1.2)
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    return ax


# --------------------------------------------------------------------------------------------------
# R4_F2e_leakage
# --------------------------------------------------------------------------------------------------


def f2e_dense(X):
    """
    Dense copy of a (possibly sparse) matrix.

    Parameters
    ----------
    X : array-like or scipy.sparse matrix
        Matrix.

    Returns
    -------
    np.ndarray
        Dense array.
    """
    return X.toarray() if issparse(X) else np.asarray(X)


def f2e_spatial_conn(coords, n_neighbors=6):
    """
    Binary spatial kNN connectivity, built as inside em.tl.score.

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates (n x 2).
    n_neighbors : int
        Number of spatial neighbours.

    Returns
    -------
    scipy.sparse.csr_matrix
        n x n connectivity.
    """
    tmp = anndata.AnnData(X=np.zeros((len(coords), 1), dtype=np.float32))
    tmp.obsm["spatial"] = np.asarray(coords)
    sq.gr.spatial_neighbors(
        tmp, n_neighs=n_neighbors, coord_type="generic", key_added="spatial"
    )
    return csr_matrix(tmp.obsp["spatial_connectivities"])


def f2e_smooth(Y, conn, n_iter=1):
    """
    EnrichMap Eq. 3: neighbourhood average with the row-normalised adjacency (the spot itself excluded).

    Parameters
    ----------
    Y : np.ndarray
        Per-spot values (n or n x k).
    conn : scipy.sparse matrix
        Spatial connectivity.
    n_iter : int
        Number of smoothing iterations.

    Returns
    -------
    np.ndarray
        Smoothed values.
    """
    Y = np.asarray(Y, dtype=np.float64)
    d = np.maximum(np.asarray(conn.sum(axis=1)).ravel(), 1e-10)
    out = Y
    for _ in range(n_iter):
        out = conn.dot(out) / (d[:, None] if out.ndim == 2 else d)
    return out


def f2e_smooth_self_inclusive(Y, conn):
    """
    Variant of EnrichMap Eq. 3 that includes the spot itself (mean of the spot and its neighbours).

    Parameters
    ----------
    Y : np.ndarray
        Per-spot values.
    conn : scipy.sparse matrix
        Spatial connectivity.

    Returns
    -------
    np.ndarray
        Smoothed values.
    """
    return f2e_smooth(Y, csr_matrix(conn + identity(conn.shape[0], format="csr")))


class F2eSpatialGAM:
    """
    The spatial-covariate GAM of EnrichMap written as a linear smoother, fitted = H y with H = B (B'B + P)^-1 B'.

    With basis "te", n_splines None and lam None this is LinearGAM(te(0, 1)) as called by em.tl.score (10 x 10 cubic
    B-splines, lam 0.6 per margin, intercept); basis "additive" gives s(0) + s(1) (20 splines per margin). te(0, 1)
    is a partition of unity, so B'B + P is singular along one direction; the sqrt(machine-eps) ridge that pyGAM adds
    before solving is applied twice, as in the original F2e analysis. `edof` = trace(H).

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates (n x 2).
    basis : str
        "te" (tensor product) or "additive".
    n_splines : int, optional
        Splines per margin (pyGAM default when None).
    lam : float, optional
        Penalty per margin (pyGAM default 0.6 when None).
    """

    def __init__(self, coords, basis="te", n_splines=None, lam=None):
        coords = np.asarray(coords, dtype=np.float64)
        kw = {}
        if n_splines is not None:
            kw["n_splines"] = n_splines
        if lam is not None:
            kw["lam"] = lam
        if basis == "te":
            term = pygam.te(0, 1, **kw)
        else:
            term = pygam.s(0, **kw) + pygam.s(1, **kw)
        gam = pygam.LinearGAM(term).fit(coords, coords[:, 0])
        B = gam._modelmat(coords)
        B = B.toarray() if hasattr(B, "toarray") else np.asarray(B)
        P = gam._P()
        P = P.toarray() if hasattr(P, "toarray") else np.asarray(P)
        ridge = np.sqrt(np.finfo(float).eps) * np.eye(B.shape[1])
        self.B, self.P = B, P
        self._A = B.T @ B + P + ridge
        self._A = self._A + ridge
        self.edof = float(np.trace(np.linalg.solve(self._A, B.T @ B)))

    def fitted(self, Y):
        """
        GAM fit H y.

        Parameters
        ----------
        Y : np.ndarray
            Per-spot values (n or n x k).

        Returns
        -------
        np.ndarray
            Fitted values.
        """
        Y = np.asarray(Y, dtype=np.float64)
        return self.B @ np.linalg.solve(self._A, self.B.T @ Y)

    def residual(self, Y):
        """
        GAM residual (I - H) y.

        Parameters
        ----------
        Y : np.ndarray
            Per-spot values (n or n x k).

        Returns
        -------
        np.ndarray
            Residuals.
        """
        return np.asarray(Y, dtype=np.float64) - self.fitted(Y)

    def r2(self, y):
        """
        Fraction of the centred variance of y that the GAM absorbs.

        Parameters
        ----------
        y : np.ndarray
            Per-spot values.

        Returns
        -------
        float
            R2 of the GAM fit.
        """
        y = np.asarray(y, dtype=np.float64)
        yc = y - y.mean()
        return float(1 - ((y - self.fitted(y)) ** 2).sum() / (yc**2).sum())


def f2e_make_gams(coords):
    """
    The three GAMs compared in F2e: default te, additive s(x)+s(y) and low-flexibility te (lam 60).

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates (n x 2).

    Returns
    -------
    dict
        GAM name -> F2eSpatialGAM.
    """
    return {
        "te (default)": F2eSpatialGAM(coords),
        "additive s(x)+s(y)": F2eSpatialGAM(coords, basis="additive"),
        "low-flex te": F2eSpatialGAM(coords, basis="te", n_splines=10, lam=60.0),
    }


def f2e_statsmodels_additive_gam_residual(coords, Y, df=10, degree=3):
    """
    Residual of the GAM of the first submission (statsmodels GLMGam, additive cubic B-splines, df 10 + 10, unpenalised).

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates (n x 2).
    Y : np.ndarray
        Per-spot values (n or n x k).
    df, degree : int
        B-spline degrees of freedom and degree per coordinate.

    Returns
    -------
    np.ndarray
        Y minus the fitted values.
    """
    from statsmodels.gam.api import BSplines, GLMGam

    coords = np.asarray(coords, dtype=np.float64)
    bs = BSplines(coords, df=[df, df], degree=[degree, degree])
    Y = np.atleast_2d(np.asarray(Y, dtype=np.float64).T).T
    out = np.empty_like(Y)
    for j in range(Y.shape[1]):
        res = GLMGam.from_formula("y ~ 1", data={"y": Y[:, j]}, smoother=bs).fit()
        out[:, j] = Y[:, j] - res.fittedvalues
    return out.ravel() if out.shape[1] == 1 else out


def f2e_hop_distance(conn, mask, max_hops=12):
    """
    Graph distance (in neighbour rings) of every spot to the nearest spot in `mask`.

    Parameters
    ----------
    conn : scipy.sparse matrix
        Spatial connectivity.
    mask : np.ndarray
        Boolean region label per spot.
    max_hops : int
        Largest ring computed (farther spots stay at inf).

    Returns
    -------
    np.ndarray
        Rings from the region (0 = inside).
    """
    mask = np.asarray(mask, bool)
    A = csr_matrix((conn + conn.T) > 0).astype(np.int8)
    dist = np.full(len(mask), np.inf)
    dist[mask] = 0
    frontier = mask.copy()
    seen = mask.copy()
    for h in range(1, max_hops + 1):
        nxt = (A.dot(frontier.astype(np.int8)) > 0) & ~seen
        dist[nxt] = h
        seen |= nxt
        frontier = nxt
        if not nxt.any():
            break
    return dist


def f2e_region_depth(conn, mask):
    """
    Mean graph distance (rings) of in-region spots to the nearest spot outside (1 = a one-spot-thick layer).

    Parameters
    ----------
    conn : scipy.sparse matrix
        Spatial connectivity.
    mask : np.ndarray
        Boolean region label per spot.

    Returns
    -------
    float
        Region depth in rings.
    """
    mask = np.asarray(mask, bool)
    d_out = f2e_hop_distance(conn, ~mask, max_hops=60)
    return float(np.mean(d_out[mask]))


def f2e_raw_enrichmap(adata, genes, weights=None):
    """
    EnrichMap Eq. 1 (CV-weighted mean of log-expression) before smoothing, the null model and the GAM.

    Parameters
    ----------
    adata : AnnData
        Log-normalised expression.
    genes : list
        Signature genes (absent genes are dropped).
    weights : dict, optional
        Gene weights; em.tl.infer_gene_weights when None.

    Returns
    -------
    np.ndarray
        Raw score per spot.
    """
    genes = [g for g in dict.fromkeys(genes) if g in adata.var_names]
    w = weights if weights is not None else em.tl.infer_gene_weights(adata, genes)
    wv = np.array([w.get(g, 1.0) for g in genes], dtype=np.float64)
    X = f2e_dense(adata[:, genes].X)
    return X @ wv / np.abs(wv).sum()


def f2e_raw_zscore(adata, genes):
    """
    Z score: sum of the signature genes' log-expression divided by sqrt(number of genes).

    Parameters
    ----------
    adata : AnnData
        Log-normalised expression.
    genes : list
        Signature genes (absent genes are dropped).

    Returns
    -------
    np.ndarray
        Score per spot.
    """
    genes = [g for g in dict.fromkeys(genes) if g in adata.var_names]
    return f2e_dense(adata[:, genes].X).sum(axis=1) / np.sqrt(len(genes))


def f2e_matched_control_sets(adata, genes, n_sets=20, n_bins=25, seed=0):
    """
    Expression-matched random control gene sets (same size, same mean-expression bins, signature genes excluded).

    Parameters
    ----------
    adata : AnnData
        Log-normalised expression.
    genes : list
        Signature genes.
    n_sets : int
        Number of control sets.
    n_bins : int
        Number of mean-expression bins.
    seed : int
        Random seed.

    Returns
    -------
    list of list
        Control gene sets.
    """
    rng = np.random.default_rng(seed)
    genes = [g for g in dict.fromkeys(genes) if g in adata.var_names]
    means = pd.Series(
        np.asarray(f2e_dense(adata.X).mean(axis=0)).ravel(), index=adata.var_names
    )
    bins = pd.qcut(means.rank(method="first"), n_bins, labels=False)
    sig = set(genes)
    pool = {
        b: [g for g in bins.index[bins == b] if g not in sig] for b in range(n_bins)
    }
    need = pd.Series([bins[g] for g in genes]).value_counts()
    sets = []
    for _ in range(n_sets):
        cs = []
        for b, k in need.items():
            cs += list(rng.choice(pool[b], size=min(k, len(pool[b])), replace=False))
        sets.append(cs)
    return sets


def f2e_reference_corrected(
    adata,
    genes,
    conn,
    gam,
    n_sets=20,
    seed=0,
    smooth_first=True,
    scorer=f2e_raw_enrichmap,
):
    """
    Reference correction: S r_sig - H (mean of S r_c over expression-matched control sets).

    The GAM is fitted to the control score, not to the signature score, so it removes only the spatial trend shared
    with random genes of matched expression.

    Parameters
    ----------
    adata : AnnData
        Log-normalised expression.
    genes : list
        Signature genes.
    conn : scipy.sparse matrix
        Spatial connectivity.
    gam : F2eSpatialGAM
        GAM operator.
    n_sets : int
        Number of control sets.
    seed : int
        Random seed of the control sets.
    smooth_first : bool
        Smooth signature and control scores before the GAM.
    scorer : callable
        Raw scorer scorer(adata, genes).

    Returns
    -------
    np.ndarray
        Reference-corrected score per spot.
    """
    r = scorer(adata, genes)
    ctrl = np.mean(
        [
            scorer(adata, cs)
            for cs in f2e_matched_control_sets(adata, genes, n_sets=n_sets, seed=seed)
        ],
        axis=0,
    )
    if smooth_first:
        r, ctrl = f2e_smooth(r, conn), f2e_smooth(ctrl, conn)
    return r - gam.fitted(ctrl)


def f2e_configs(raw, conn, co, gams, ref=None):
    """
    Score maps derived from one raw score: raw, smoothing (package and self-inclusive) and the GAM variants.

    Parameters
    ----------
    raw : np.ndarray
        Raw per-spot score (or the region indicator).
    conn : scipy.sparse matrix
        6-NN spatial connectivity.
    co : np.ndarray
        Spot coordinates.
    gams : dict
        GAM operators from f2e_make_gams.
    ref : np.ndarray, optional
        Reference-corrected score, added as "smooth + reference correction".

    Returns
    -------
    dict
        Configuration name -> per-spot score.
    """
    gte = gams["te (default)"]
    sm = f2e_smooth(raw, conn)
    out = {
        "raw": raw,
        "raw + te GAM": gte.residual(raw),
        "smooth": sm,
        "smooth (self-inclusive)": f2e_smooth_self_inclusive(raw, conn),
        "smooth + te GAM (default)": gte.residual(sm),
        "smooth + additive GAM (pyGAM)": gams["additive s(x)+s(y)"].residual(sm),
        "smooth + additive GAM (statsmodels, 1st submission)": f2e_statsmodels_additive_gam_residual(
            co, sm
        ),
        "smooth + low-flex te GAM": gams["low-flex te"].residual(sm),
    }
    if ref is not None:
        out["smooth + reference correction"] = ref
    return out


def f2e_distance_zones():
    """
    Distance zones (neighbour rings from the region) used for the top-k shares and the binary calls.

    Returns
    -------
    dict
        Zone name -> (first ring, last ring).
    """
    return {
        "inside": (0, 0),
        "halo (1-3)": (1, 3),
        "near (4-8)": (4, 8),
        "far (>=9)": (9, np.inf),
    }


def f2e_leakage_profile(score, mask, dist, max_ring=8, far_from=6):
    """
    Standardised ring profile L(d) = (ring-d mean - far mean) / (inside mean - far mean).

    Parameters
    ----------
    score : np.ndarray
        Per-spot score.
    mask : np.ndarray
        Boolean region label per spot.
    dist : np.ndarray
        Rings from the region.
    max_ring : int
        Last ring.
    far_from : int
        Far background = spots at least this many rings away.

    Returns
    -------
    dict
        Ring d -> L(d) (0 = no leakage, 1 = as high as inside).
    """
    score = np.asarray(score, dtype=np.float64)
    inside = score[np.asarray(mask, bool)].mean()
    far = score[dist >= far_from].mean()
    span = inside - far
    prof = {}
    for d in range(1, max_ring + 1):
        sel = dist == d
        prof[d] = (
            float((score[sel].mean() - far) / span)
            if sel.any() and span != 0
            else np.nan
        )
    return prof


def f2e_leakage_summary(score, mask, dist, halo=(1, 2, 3)):
    """
    Halo leakage index (mean L over the halo rings), region-vs-halo AUROC and the top-k shares in halo and region.

    Parameters
    ----------
    score : np.ndarray
        Per-spot score.
    mask : np.ndarray
        Boolean region label per spot.
    dist : np.ndarray
        Rings from the region.
    halo : tuple
        Halo rings.

    Returns
    -------
    dict
        Metric name -> value.
    """
    score = np.asarray(score, dtype=np.float64)
    mask = np.asarray(mask, bool)
    prof = f2e_leakage_profile(score, mask, dist)
    halo_sel = np.isin(dist, halo)
    y = np.r_[np.ones(mask.sum()), np.zeros(halo_sel.sum())]
    sv = np.r_[score[mask], score[halo_sel]]
    k = int(mask.sum())
    top = np.argsort(-score, kind="stable")[:k]
    return {
        "halo leakage index": float(np.nanmean([prof[d] for d in halo])),
        "region-vs-halo AUROC": float(roc_auc_score(y, sv)),
        "top-k share in halo": float(np.isin(dist[top], halo).mean()),
        "top-k share inside": float(mask[top].mean()),
    }


def f2e_leakage_metrics(score, mask, dist, max_ring=8, far_from=6):
    """
    All leakage metrics of one score map.

    L(d) and the halo index are not invariant to a smooth additive surface (a GAM can lower them by pushing the
    neighbourhood of a region below the far background, a negative moat), so offset-robust metrics are added:
    AUROC inside vs ring 1 (boundary sharpness), AUROC halo vs far (< 0.5 = moat), Z(d) in far-background SD, the
    moat (min L(d), min Z(d) and their rings), the zones of the top-k spots (k = region size) and AUROC region vs rest.

    Parameters
    ----------
    score : np.ndarray
        Per-spot score.
    mask : np.ndarray
        Boolean region label per spot.
    dist : np.ndarray
        Rings from the region.
    max_ring : int
        Last ring of the profiles.
    far_from : int
        Far background = spots at least this many rings away.

    Returns
    -------
    dict
        Metric name -> value.
    """
    s = np.asarray(score, dtype=np.float64)
    mask = np.asarray(mask, bool)
    out = {}
    prof = f2e_leakage_profile(s, mask, dist, max_ring=max_ring, far_from=far_from)
    for d, v in prof.items():
        out[f"L{d}"] = v
    out.update(f2e_leakage_summary(s, mask, dist))
    far = dist >= far_from
    mu_far, sd_far = s[far].mean(), s[far].std()
    for d in range(1, max_ring + 1):
        sel = dist == d
        out[f"Z{d}"] = (
            float((s[sel].mean() - mu_far) / sd_far)
            if (sel.any() and sd_far > 0)
            else np.nan
        )
    Ls = np.array([prof[d] for d in range(1, max_ring + 1)])
    Zs = np.array([out[f"Z{d}"] for d in range(1, max_ring + 1)])
    out["moat L min"] = float(np.nanmin(Ls)) if np.isfinite(Ls).any() else np.nan
    out["moat L ring"] = int(np.nanargmin(Ls) + 1) if np.isfinite(Ls).any() else -1
    out["moat Z min"] = float(np.nanmin(Zs)) if np.isfinite(Zs).any() else np.nan
    out["moat Z ring"] = int(np.nanargmin(Zs) + 1) if np.isfinite(Zs).any() else -1
    out["inside Z"] = (
        float((s[mask].mean() - mu_far) / sd_far) if sd_far > 0 else np.nan
    )
    r1 = dist == 1
    out["AUROC inside vs ring1"] = float(
        roc_auc_score(
            np.r_[np.ones(mask.sum()), np.zeros(r1.sum())], np.r_[s[mask], s[r1]]
        )
    )
    halo = np.isin(dist, (1, 2, 3))
    out["AUROC halo vs far"] = float(
        roc_auc_score(
            np.r_[np.ones(halo.sum()), np.zeros(far.sum())], np.r_[s[halo], s[far]]
        )
    )
    k = int(mask.sum())
    top = np.argsort(-s, kind="stable")[:k]
    for z, (lo, hi) in f2e_distance_zones().items():
        out[f"top-k share {z}"] = float(((dist[top] >= lo) & (dist[top] <= hi)).mean())
    out["AUROC"] = float(roc_auc_score(mask.astype(int), s))
    return out


def f2e_far_units(s, far):
    """
    Score in far-background SD units.

    Parameters
    ----------
    s : np.ndarray
        Per-spot score.
    far : np.ndarray
        Boolean mask of far-background spots (rings >= 6).

    Returns
    -------
    np.ndarray
        (s - mean of far background) / SD of far background.
    """
    return (s - s[far].mean()) / s[far].std()


def f2e_topk_call(s, k):
    """
    The k highest-scoring spots (stable order).

    Parameters
    ----------
    s : np.ndarray
        Per-spot score.
    k : int
        Number of spots to call (the number of region spots).

    Returns
    -------
    np.ndarray
        Boolean call per spot.
    """
    call = np.zeros(len(s), bool)
    call[np.argsort(-np.asarray(s, dtype=np.float64), kind="stable")[:k]] = True
    return call


def f2e_call_zones(call, mask, dist, negative=None):
    """
    Where the called spots lie, and the share of each zone that is called.

    Parameters
    ----------
    call : np.ndarray
        Boolean positive call per spot.
    mask : np.ndarray
        Boolean region label per spot.
    dist : np.ndarray
        Rings from the region (0 = inside).
    negative : np.ndarray, optional
        Boolean opposite call ("depleted", e.g. z < -1.96) per spot, for the halo and far rates.

    Returns
    -------
    dict
        Number of calls, precision, recall, F1, share of calls and call rate per zone, negative-call rates.
    """
    call = np.asarray(call, bool)
    mask = np.asarray(mask, bool)
    zones = f2e_distance_zones()
    n = int(call.sum())
    out = {
        "n region spots": int(mask.sum()),
        "n called": n,
        "called / region size": n / mask.sum(),
    }
    for z, (lo, hi) in zones.items():
        zone = (dist >= lo) & (dist <= hi)
        out[f"share of calls {z}"] = float((call & zone).sum() / n) if n else np.nan
        out[f"call rate {z}"] = float(call[zone].mean()) if zone.any() else np.nan
    tp = int((call & mask).sum())
    precision = tp / n if n else 0.0
    recall = tp / mask.sum()
    out["precision"] = precision
    out["recall"] = recall
    out["F1"] = 2 * precision * recall / (precision + recall) if tp else 0.0
    for z in ("halo (1-3)", "far (>=9)"):
        lo, hi = zones[z]
        zone = (dist >= lo) & (dist <= hi)
        out[f"negative-call rate {z}"] = (
            float(np.asarray(negative, bool)[zone].mean())
            if negative is not None
            else np.nan
        )
    return out


def f2e_binary_calls(scores, p_emp, gsd_label, mask, dist):
    """
    Binary calls of the three EnrichMap variants and GSDensity, by distance zone.

    Rules: prevalence-matched top-k for every method; the null model's z > 0, z > 1.96 and p_emp <= 0.05; the GAM
    residual's > 0; GSDensity's antimode labels. The smoothing-only score is a non-negative weighted mean of
    log-expression, so it has no sign rule.

    Parameters
    ----------
    scores : dict
        "EnrichMap (null)", "EnrichMap (smoothing-only)", "EnrichMap (GAM)" and "GSDensity" (relevance) ->
        per-spot score.
    p_emp : np.ndarray
        Empirical p-value of the null model per spot.
    gsd_label : np.ndarray
        GSDensity antimode label per spot ("positive" / "negative").
    mask : np.ndarray
        Boolean region label per spot.
    dist : np.ndarray
        Rings from the region.

    Returns
    -------
    list of dict
        One row per (method, rule) with the f2e_call_zones output.
    """
    methods = [
        "EnrichMap (null)",
        "EnrichMap (smoothing-only)",
        "EnrichMap (GAM)",
        "GSDensity",
    ]
    k = int(np.asarray(mask, bool).sum())
    z = np.asarray(scores["EnrichMap (null)"], dtype=np.float64)
    gam = np.asarray(scores["EnrichMap (GAM)"], dtype=np.float64)
    calls = [(m, "top-k", f2e_topk_call(scores[m], k), None) for m in methods]
    calls += [
        ("EnrichMap (null)", "score > 0", z > 0, z < 0),
        ("EnrichMap (null)", "z > 1.96", z > 1.96, z < -1.96),
        ("EnrichMap (null)", "p_emp ≤ 0.05", np.asarray(p_emp) <= 0.05, None),
        ("EnrichMap (GAM)", "score > 0", gam > 0, gam < 0),
        (
            "GSDensity",
            "antimode",
            np.asarray(gsd_label).astype(str) == "positive",
            None,
        ),
    ]
    return [
        {"method": m, "rule": r, **f2e_call_zones(c, mask, dist, negative=neg)}
        for m, r, c, neg in calls
    ]


def f2e_simulate_visium(
    seed=0,
    n_genes=3000,
    n_signature_genes=20,
    n_clusters=4,
    target_cluster=1,
    n_spots_range=(1500, 4000),
    spot_spacing_um=100,
    spot_diameter_um=55,
    mean_counts_per_spot=5000,
    dispersion_range=(0.1, 2.0),
    fold_change_range=(1.5, 5.0),
    dropout_rate=0.1,
    domain="voronoi",
    domain_width_spots=1.5,
    n_patches=12,
    patch_radius_spots=3.0,
    tech_field=None,
    tech_strength=0.0,
    tech_gene_frac=0.5,
    tech_scale=0.5,
    tech_on_signature=None,
    return_truth=False,
):
    """
    simulate_visium of 05_01_simulated_benchmarking (Sup. Fig. 1), extended with other target domains and smooth
    technical fields.

    With the default arguments the random stream and the output are those of the original function; the extensions
    draw from an independent random generator.

    Parameters
    ----------
    seed : int
        Random seed.
    n_genes, n_signature_genes, n_clusters, target_cluster : int
        As in 05_01_simulated_benchmarking.
    n_spots_range, spot_spacing_um, spot_diameter_um, mean_counts_per_spot, dispersion_range, fold_change_range,
    dropout_rate
        As in 05_01_simulated_benchmarking.
    domain : str
        "voronoi" (original: Voronoi cells, target = cluster_2), "stripe" (a thin sinusoidal band about
        2 x domain_width_spots spots wide) or "patches" (n_patches discs of radius patch_radius_spots).
    domain_width_spots, n_patches, patch_radius_spots : float
        Stripe half-width and patch layout.
    tech_field : str, optional
        None, "gradient", "edge" or "grf": smooth multiplicative technical field F(x, y) on a fraction
        `tech_gene_frac` of genes (log F has range `tech_strength`; "grf" has length scale `tech_scale` x side).
    tech_strength, tech_gene_frac, tech_scale : float
        Technical field strength, share of affected genes and length scale.
    tech_on_signature : bool, optional
        None (signature genes affected at the same rate), True (all affected) or False (none affected).
    return_truth : bool
        Also return the true log-field and the affected genes.

    Returns
    -------
    tuple
        (adata, signature genes, fold changes, target label) and, with `return_truth`, a dict with "log_field" and
        "affected".
    """
    rng = np.random.default_rng(seed)
    rng_ext = np.random.default_rng(10_000_019 + seed)

    n_spots_target = rng.integers(n_spots_range[0], n_spots_range[1] + 1)
    grid_side = int(np.ceil(np.sqrt(n_spots_target * 0.866)))
    capture_area_um = grid_side * spot_spacing_um
    n_cols = int(capture_area_um / spot_spacing_um)
    n_rows = int(capture_area_um / (spot_spacing_um * np.sqrt(3) / 2))
    x_coords, y_coords = [], []
    for row in range(n_rows):
        for col in range(n_cols):
            x = col * spot_spacing_um + (row % 2) * (spot_spacing_um / 2)
            y = row * spot_spacing_um * np.sqrt(3) / 2
            if x <= capture_area_um and y <= capture_area_um:
                x_coords.append(x)
                y_coords.append(y)
    x_coords = np.array(x_coords)
    y_coords = np.array(y_coords)
    n_spots = len(x_coords)
    coords = np.column_stack([x_coords, y_coords])

    margin = capture_area_um * 0.15
    seed_coords = np.column_stack(
        [
            rng.uniform(margin, capture_area_um - margin, n_clusters),
            rng.uniform(margin, capture_area_um - margin, n_clusters),
        ]
    )
    dists = cdist(coords, seed_coords)
    dists += rng.normal(0, spot_spacing_um * 0.8, dists.shape)
    cluster_ids = np.argmin(dists, axis=1)
    cluster_labels = np.array([f"cluster_{i + 1}" for i in cluster_ids])
    target_label = f"cluster_{target_cluster + 1}"

    if domain == "stripe":
        L = capture_area_um
        phase = rng_ext.uniform(0, 2 * np.pi)
        amp = rng_ext.uniform(0.08, 0.15) * L
        y0 = rng_ext.uniform(0.35, 0.65) * L
        curve = y0 + amp * np.sin(2 * np.pi * coords[:, 0] / (0.8 * L) + phase)
        in_target = np.abs(coords[:, 1] - curve) <= domain_width_spots * spot_spacing_um
        cluster_labels = np.where(
            in_target,
            target_label,
            np.where(cluster_labels == target_label, "cluster_bg", cluster_labels),
        )
    elif domain == "patches":
        L = capture_area_um
        centres = rng_ext.uniform(0.08 * L, 0.92 * L, size=(n_patches, 2))
        dd = cdist(coords, centres).min(axis=1)
        in_target = dd <= patch_radius_spots * spot_spacing_um
        cluster_labels = np.where(
            in_target,
            target_label,
            np.where(cluster_labels == target_label, "cluster_bg", cluster_labels),
        )

    gene_names = [f"GENE_{i + 1:05d}" for i in range(n_genes)]
    sig_indices = rng.choice(n_genes, n_signature_genes, replace=False)
    fold_changes = rng.uniform(
        fold_change_range[0], fold_change_range[1], n_signature_genes
    )
    order = np.argsort(-fold_changes)
    sig_indices = sig_indices[order]
    fold_changes = fold_changes[order]
    sig_genes = [gene_names[i] for i in sig_indices]

    gene_means = np.exp(rng.normal(1.5, 1.0, n_genes))
    gene_means = np.clip(gene_means, 0.1, 100)
    dispersions = rng.uniform(dispersion_range[0], dispersion_range[1], n_genes)
    lib_size_factors = np.exp(rng.normal(0, 0.3, n_spots))
    lib_size_factors *= mean_counts_per_spot / (
        gene_means.sum() * lib_size_factors.mean()
    )

    # technical field (independent random generator)
    log_field = np.zeros(n_spots)
    affected = np.zeros(n_genes, bool)
    if tech_field is not None and tech_strength > 0:
        u = coords / capture_area_um
        if tech_field == "gradient":
            ang = rng_ext.uniform(0, 2 * np.pi)
            proj = u @ np.array([np.cos(ang), np.sin(ang)])
            base = (proj - proj.min()) / (proj.max() - proj.min())
        elif tech_field == "edge":
            d_edge = np.minimum.reduce([u[:, 0], 1 - u[:, 0], u[:, 1], 1 - u[:, 1]])
            base = 1 - np.exp(-d_edge / 0.12)
        else:
            # random Fourier features of an RBF-kernel Gaussian field, length scale tech_scale x side
            K = 200
            freqs = rng_ext.normal(0, 1.0 / tech_scale, size=(K, 2))
            phases = rng_ext.uniform(0, 2 * np.pi, K)
            f = np.cos(u @ freqs.T + phases).sum(axis=1)
            base = (f - f.min()) / (f.max() - f.min())
        log_field = tech_strength * (base - base.mean())
        affected = rng_ext.random(n_genes) < tech_gene_frac
        if tech_on_signature is True:
            affected[sig_indices] = True
        elif tech_on_signature is False:
            affected[sig_indices] = False
    field = np.exp(log_field)

    X = np.zeros((n_spots, n_genes), dtype=np.float32)
    target_mask = cluster_labels == target_label
    for g in range(n_genes):
        mu = gene_means[g] * lib_size_factors
        if g in sig_indices:
            fc_idx = np.where(sig_indices == g)[0][0]
            mu[target_mask] *= fold_changes[fc_idx]
        if affected[g]:
            mu = mu * field
        r = dispersions[g]
        p = r / (r + mu)
        counts = rng.negative_binomial(r, p)
        dropout_mask = rng.random(n_spots) < dropout_rate
        counts[dropout_mask] = 0
        X[:, g] = counts

    adata = anndata.AnnData(X=csr_matrix(X))
    adata.var_names = gene_names
    adata.obs_names = [f"spot_{i}" for i in range(n_spots)]
    adata.obs["cluster"] = pd.Categorical(cluster_labels)
    adata.obsm["spatial"] = coords
    adata.layers["counts"] = adata.X.copy()
    adata.obs["total_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["tech_log_field"] = log_field
    adata.var["tech_affected"] = affected
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    if return_truth:
        return (
            adata,
            sig_genes,
            fold_changes,
            target_label,
            {"log_field": log_field, "affected": affected},
        )
    return adata, sig_genes, fold_changes, target_label


def f2e_spot_pitch(coords):
    """
    Median distance of each spot to its nearest neighbour.

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates (n x 2).

    Returns
    -------
    float
        Spot pitch.
    """
    d, _ = cKDTree(coords).query(coords, k=2)
    return float(np.median(d[:, 1]))


def f2e_diffusion_operator(coords, rho, sigma_pitch, pitch=None, truncate=3.0):
    """
    Lateral transcript diffusion M = (1 - rho) I + rho K (cf. SpotClean).

    K is a column-stochastic Gaussian kernel (self excluded, truncated at `truncate` sigma): K_ij = share of the
    transcripts leaving spot j that are captured at spot i, so the total is conserved.

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates (n x 2).
    rho : float
        Share of each spot's counts that diffuses to its neighbours.
    sigma_pitch : float
        Kernel SD in spot pitches.
    pitch : float, optional
        Spot pitch (median nearest-neighbour distance when None).
    truncate : float
        Kernel truncation in SDs.

    Returns
    -------
    tuple of scipy.sparse.csr_matrix
        (M, K).
    """
    coords = np.asarray(coords, dtype=np.float64)
    pitch = f2e_spot_pitch(coords) if pitch is None else pitch
    sig = sigma_pitch * pitch
    tree = cKDTree(coords)
    pairs = tree.query_pairs(truncate * sig, output_type="ndarray")
    i, j = pairs[:, 0], pairs[:, 1]
    d2 = ((coords[i] - coords[j]) ** 2).sum(1)
    w = np.exp(-d2 / (2 * sig**2))
    n = len(coords)
    W = csr_matrix((np.r_[w, w], (np.r_[i, j], np.r_[j, i])), shape=(n, n))
    colsum = np.asarray(W.sum(axis=0)).ravel()
    K = W @ diags(1 / np.maximum(colsum, 1e-300))
    M = (1 - rho) * diags(np.ones(n)) + rho * K
    return csr_matrix(M), csr_matrix(K)


def f2e_counts_to_adata(X, template):
    """
    AnnData from a count matrix, library-size normalised (1e4) and log1p-transformed as in f2e_simulate_visium.

    Parameters
    ----------
    X : np.ndarray
        Counts (spots x genes).
    template : AnnData
        Source of obs, var names and coordinates.

    Returns
    -------
    AnnData
        Log-normalised expression with counts in layers["counts"].
    """
    a = anndata.AnnData(
        X=csr_matrix(X.astype(np.float32)),
        obs=template.obs.copy(),
        var=pd.DataFrame(index=template.var_names),
    )
    a.obsm["spatial"] = template.obsm["spatial"].copy()
    a.layers["counts"] = a.X.copy()
    a.obs["total_counts"] = np.asarray(a.X.sum(axis=1)).ravel()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a


def f2e_apply_diffusion(adata, rho, sigma_pitch, seed=0):
    """
    Diffuse the counts of a simulated slide with M, round them stochastically and normalise again.

    Parameters
    ----------
    adata : AnnData
        Simulated slide with layers["counts"].
    rho : float
        Share of each spot's counts that diffuses.
    sigma_pitch : float
        Kernel SD in spot pitches.
    seed : int
        Seed of the stochastic rounding.

    Returns
    -------
    tuple
        (diffused AnnData, M, K).
    """
    rng = np.random.default_rng(20_000_003 + seed)
    X = f2e_dense(adata.layers["counts"]).astype(np.float64)
    M, K = f2e_diffusion_operator(adata.obsm["spatial"], rho, sigma_pitch)
    Xd = M @ X
    fl = np.floor(Xd)
    Xd = fl + (rng.random(Xd.shape) < (Xd - fl))
    return f2e_counts_to_adata(Xd, adata), M, K


def f2e_sim_gsdensity(arm, seed, spots, gsd_sim_dir):
    """
    GSDensity relevance and antimode label of one F2e slide (voronoi arm: the Sup. Fig. 1 base run of the same slide).

    Parameters
    ----------
    arm : str
        Simulation arm.
    seed : int
        Slide seed.
    spots : array-like
        Spot names of the slide.
    gsd_sim_dir : str
        GSDensity directory of the simulations (revision/data/gsd_cache/sims/).

    Returns
    -------
    pd.DataFrame
        Columns relevance, label and runtime_s (NaN for the voronoi arm), indexed by spot.
    """
    if arm == "voronoi":
        g = pd.read_csv(
            gsd_sim_dir + f"SupFig1_seed{seed:03d}_cluster_2.csv", index_col=0
        )
        g["runtime_s"] = np.nan
    else:
        name = arm.replace("(", "_").replace(")", "").replace(",", "-")
        g = pd.read_csv(gsd_sim_dir + f"R4_F2e/{name}__{seed}.csv", index_col=0)
    return g.reindex(spots)[["relevance", "label", "runtime_s"]]


def f2e_score_sim_slide(arm, seed, spec, variants, gsd_sim_dir):
    """
    Simulate one F2e slide, score the signature with the EnrichMap variants (6-NN smoothing) and read GSDensity.

    Parameters
    ----------
    arm : str
        Simulation arm.
    seed : int
        Slide seed.
    spec : dict
        Arm definition: "domain", optional "tech" (technical field arguments) and "diff" ((rho, sigma) of the lateral
        diffusion).
    variants : dict
        Variant name ("null", "smooth", "gam") -> em.tl.score arguments.
    gsd_sim_dir : str
        GSDensity directory of the simulations.

    Returns
    -------
    tuple
        (scored AnnData, dict of per-spot arrays: spot, coords, conn, mask, dist, raw, one array per variant, p_emp,
        gsd, gsd_label, plus control_overlap, level_offset_ratio and gsd_runtime).
    """
    adata, sig, fold_changes, target, truth = f2e_simulate_visium(
        seed=seed, domain=spec["domain"], return_truth=True, **spec.get("tech", {})
    )
    mask = (adata.obs["cluster"] == target).values
    coords = adata.obsm["spatial"].astype(float)
    conn = f2e_spatial_conn(coords)
    dist = f2e_hop_distance(conn, mask, max_hops=40)
    if spec.get("diff"):
        rho, sigma = spec["diff"]
        adata, _, _ = f2e_apply_diffusion(adata, rho, sigma, seed=seed)
    for variant, kwargs in variants.items():
        em.tl.score(
            adata,
            gene_set=list(sig),
            score_key=f"target_{variant}",
            smoothing=True,
            n_neighbors=6,
            **kwargs,
        )
    info = adata.uns["null_model"]["target_null"]
    gsd = f2e_sim_gsdensity(arm, seed, adata.obs_names, gsd_sim_dir)
    slide = {
        "spot": np.asarray(adata.obs_names),
        "coords": coords,
        "conn": conn,
        "mask": mask,
        "dist": dist,
        "raw": f2e_raw_enrichmap(adata, sig),
        **{
            v: adata.obs[f"target_{v}_score"].values.astype(np.float64)
            for v in variants
        },
        "p_emp": np.asarray(info["p_empirical"], dtype=np.float64),
        "control_overlap": float(info["control_overlap"]),
        "level_offset_ratio": float(info["level_offset_ratio"]),
        "gsd": gsd["relevance"].values.astype(np.float64),
        "gsd_label": gsd["label"].values.astype(str),
        "gsd_runtime": float(gsd["runtime_s"].iloc[0]),
    }
    return adata, slide


def f2e_sim_slide_tables(arm, seed, spec, variants, gsd_sim_dir):
    """
    Leakage metrics, null-model and GSDensity diagnostics and binary calls of one F2e slide.

    The package smoothing-only and GAM scores are also compared with the "smooth" and "smooth + te GAM (default)"
    configurations of the raw EnrichMap score (max absolute difference).

    Parameters
    ----------
    arm : str
        Simulation arm.
    seed : int
        Slide seed.
    spec : dict
        Arm definition (see f2e_score_sim_slide).
    variants : dict
        Variant name ("null", "smooth", "gam") -> em.tl.score arguments.
    gsd_sim_dir : str
        GSDensity directory of the simulations.

    Returns
    -------
    tuple of pd.DataFrame
        (leakage metrics and slide diagnostics: one row per method, binary calls: one row per method and rule).
    """
    adata, m = f2e_score_sim_slide(arm, seed, spec, variants, gsd_sim_dir)
    mask, dist, z, p, lab = m["mask"], m["dist"], m["null"], m["p_emp"], m["gsd_label"]
    halo, far = np.isin(dist, (1, 2, 3)), dist >= 6
    sm = f2e_smooth(m["raw"], m["conn"])
    meta = {
        "arm": arm,
        "seed": seed,
        "n_spots": adata.n_obs,
        "prevalence": mask.mean(),
        "null control overlap": m["control_overlap"],
        "null level offset ratio": m["level_offset_ratio"],
        "null empty spots": int(((z == 0) & (p == 1)).sum()),
        "null z>1.96 inside": float(np.mean(z[mask] > 1.96)),
        "null z>1.96 halo (1-3)": float(np.mean(z[halo] > 1.96)),
        "null z<-1.96 halo (1-3)": float(np.mean(z[halo] < -1.96)),
        "null z>1.96 far (>=6)": float(np.mean(z[far] > 1.96)),
        "null z<-1.96 far (>=6)": float(np.mean(z[far] < -1.96)),
        "null p_emp<=0.025 far (>=6)": float(np.mean(p[far] <= 0.025)),
        "max |smoothing-only - smooth|": float(np.abs(m["smooth"] - sm).max()),
        "max |GAM - smooth + te GAM|": float(
            np.abs(m["gam"] - F2eSpatialGAM(m["coords"]).residual(sm)).max()
        ),
        "GSDensity runtime_s": m["gsd_runtime"],
        "GSDensity zero relevance": int((m["gsd"] == 0).sum()),
        "GSDensity positive labels": int((lab == "positive").sum()),
        "GSDensity positive inside": float(np.mean(lab[mask] == "positive")),
        "GSDensity positive halo (1-3)": float(np.mean(lab[halo] == "positive")),
        "GSDensity positive far (>=6)": float(np.mean(lab[far] == "positive")),
    }
    scores = {
        "EnrichMap (null)": m["null"],
        "EnrichMap (smoothing-only)": m["smooth"],
        "EnrichMap (GAM)": m["gam"],
        "GSDensity": m["gsd"],
    }
    metrics = pd.DataFrame(
        [
            {**meta, "method": k, **f2e_leakage_metrics(s, mask, dist, max_ring=12)}
            for k, s in scores.items()
        ]
    )
    calls = pd.DataFrame(
        [
            {"arm": arm, "seed": seed, "n_spots": adata.n_obs, **row}
            for row in f2e_binary_calls(scores, p, lab, mask, dist)
        ]
    )
    return metrics, calls


def f2e_mean_ci(v):
    """
    Mean and 95% t-interval of the mean over seeds.

    Parameters
    ----------
    v : pd.Series
        Per-seed values.

    Returns
    -------
    tuple
        (mean, lower, upper).
    """
    v = v.dropna()
    m = v.mean()
    h = t_dist.ppf(0.975, len(v) - 1) * v.std() / np.sqrt(len(v)) if len(v) > 1 else 0.0
    return m, m - h, m + h


def f2e_mean_sd(df, metrics, arm_order, cfg_order):
    """
    "mean ± SD" over seeds per arm and configuration.

    Parameters
    ----------
    df : pd.DataFrame
        Per-seed leakage metrics with arm and config columns.
    metrics : list
        Metric columns.
    arm_order, cfg_order : list
        Row order.

    Returns
    -------
    pd.DataFrame
        "mean ± SD" strings indexed by (arm, config).
    """
    g = df.groupby(["arm", "config"])[metrics]
    m, s = g.mean(), g.std()
    out = pd.DataFrame(
        {k: [f"{a:.3f} ± {b:.3f}" for a, b in zip(m[k], s[k])] for k in metrics},
        index=m.index,
    )
    return out.reindex(
        [(a, c) for a in arm_order for c in cfg_order if (a, c) in out.index]
    )


def f2e_colour_shades(base_colour, n, lightest=0.86, darkest=0.30):
    """
    Shades of the hue of one colour, from dark to light.

    Parameters
    ----------
    base_colour : str
        Hex colour whose hue and saturation are kept.
    n : int
        Number of shades.
    lightest, darkest : float
        HLS lightness range.

    Returns
    -------
    list of str
        Hex colours.
    """
    h, _, s = colorsys.rgb_to_hls(
        *(int(base_colour[i : i + 2], 16) / 255 for i in (1, 3, 5))
    )
    return [
        "#%02x%02x%02x"
        % tuple(round(v * 255) for v in colorsys.hls_to_rgb(h, lightness, s))
        for lightness in np.linspace(darkest, lightest, n)
    ]


def f2e_plot_metric_pair(
    table, x, y, ax, colours, ref_x=None, ref_y=None, invert_y=False, marker_col=None
):
    """
    One metric against another, one point per method (and task), with 95% CI whiskers when f"{x}_lo"/f"{x}_hi" and
    f"{y}_lo"/f"{y}_hi" columns exist.

    Parameters
    ----------
    table : pd.DataFrame
        One row per method (and task) with a "method" column and the metric columns.
    x, y : str
        Metric columns.
    ax : matplotlib Axes
        Target axes.
    colours : dict
        Method -> colour.
    ref_x, ref_y : float, optional
        Dotted reference lines.
    invert_y : bool
        Invert the y axis (lower-is-better metrics).
    marker_col : str, optional
        Column with a marker per row.

    Returns
    -------
    matplotlib Axes
    """
    for _, r in table.iterrows():
        xerr = yerr = None
        if f"{x}_lo" in table:
            xerr = [[r[x] - r[f"{x}_lo"]], [r[f"{x}_hi"] - r[x]]]
        if f"{y}_lo" in table:
            yerr = [[r[y] - r[f"{y}_lo"]], [r[f"{y}_hi"] - r[y]]]
        ax.errorbar(
            r[x],
            r[y],
            xerr=xerr,
            yerr=yerr,
            fmt=r[marker_col] if marker_col else "o",
            color=colours[r["method"]],
            ms=5,
            mec="white",
            mew=0.4,
            elinewidth=0.6,
            capsize=0,
            label=r["method"],
        )
    if ref_x is not None:
        ax.axvline(ref_x, ls=":", c="0.6", lw=0.6)
    if ref_y is not None:
        ax.axhline(ref_y, ls=":", c="0.6", lw=0.6)
    if invert_y:
        ax.invert_yaxis()
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.grid(False)
    return ax


def f2e_plot_benchmark_summary(table, ax, colours, lower_is_better=(), annotate=True):
    """
    Methods x metrics summary: circles sized and shaded by each metric scaled to 0-1 across methods (best = 1), raw
    values annotated, and a mean-rank column.

    Parameters
    ----------
    table : pd.DataFrame
        Index = methods, columns = metrics.
    ax : matplotlib Axes
        Target axes.
    colours : dict
        Method -> colour of the tick label.
    lower_is_better : tuple of str
        Metrics where smaller is better.
    annotate : bool
        Print the raw value in each cell.

    Returns
    -------
    matplotlib Axes
    """
    scaled = table.copy().astype(float)
    for c in table.columns:
        col = table[c].astype(float)
        rng_ = col.max() - col.min()
        s = (col - col.min()) / rng_ if rng_ > 0 else col * 0 + 1
        scaled[c] = 1 - s if c in lower_is_better else s
    ranks = scaled.rank(ascending=False).mean(axis=1)
    cmap = plt.get_cmap("Purples")
    for i, m in enumerate(table.index):
        for j, c in enumerate(table.columns):
            s = scaled.loc[m, c]
            ax.scatter(
                j, i, s=30 + 220 * s, color=cmap(0.25 + 0.75 * s), edgecolors="none"
            )
            if annotate:
                ax.text(
                    j,
                    i,
                    f"{table.loc[m, c]:.2g}",
                    ha="center",
                    va="center",
                    fontsize=5,
                    color="white" if s > 0.6 else "black",
                )
        ax.text(
            len(table.columns) + 0.2,
            i,
            f"{ranks[m]:.1f}",
            ha="left",
            va="center",
            fontsize=6,
        )
    ax.set_xticks(range(len(table.columns)))
    ax.set_xticklabels(table.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(table.index)))
    ax.set_yticklabels(table.index)
    for lab in ax.get_yticklabels():
        lab.set_color(colours[lab.get_text()])
    ax.text(len(table.columns) + 0.2, -0.9, "mean rank", fontsize=6, ha="left")
    ax.set_xlim(-0.6, len(table.columns) + 0.9)
    ax.set_ylim(len(table.index) - 0.4, -1.2)
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    return ax


def f2e_plot_call_zones(
    ax, d, title, show_labels, method_order, rule_order, zone_colours, colours
):
    """
    Stacked bars of the share of called spots per distance zone, one bar per method and rule.

    Parameters
    ----------
    ax : matplotlib Axes
        Target axes.
    d : pd.DataFrame
        One row per (method, rule) with "share of calls <zone>" and "called / region size" columns.
    title : str
        Panel title.
    show_labels : bool
        Draw the method / rule tick labels.
    method_order, rule_order : list
        Bar order.
    zone_colours : dict
        Zone -> colour.
    colours : dict
        Method -> colour of the tick label.

    Returns
    -------
    list of dict
        Plotted rows for the source data.
    """
    idx = d.set_index(["method", "rule"]).index
    d = d.set_index(["method", "rule"]).reindex(
        [(m, r) for m in method_order for r in rule_order if (m, r) in idx]
    )
    rows = []
    for i, ((m, r), row) in enumerate(d.iterrows()):
        left = 0.0
        for z, col in zone_colours.items():
            v = row[f"share of calls {z}"]
            v = 0.0 if not np.isfinite(v) else v
            ax.barh(i, v, left=left, color=col, height=0.75, lw=0)
            left += v
        ax.text(
            1.02,
            i,
            f"{row['called / region size']:.2g}",
            va="center",
            ha="left",
            fontsize=4.5,
        )
        rows.append({"method": m, "rule": r, **row.to_dict()})
    ax.set_yticks(range(len(d)))
    if show_labels:
        ax.set_yticklabels(
            [
                f"{m.replace('EnrichMap (', '').replace(')', '')}: {r}"
                for m, r in d.index
            ],
            fontsize=5,
        )
        for lab, (m, r) in zip(ax.get_yticklabels(), d.index):
            lab.set_color(colours[m])
    else:
        ax.set_yticklabels([])
    ax.set_ylim(len(d) - 0.5, -0.5)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.5, 1])
    ax.set_title(title, fontsize=6)
    ax.grid(False)
    sns.despine(ax=ax)
    return rows


# --------------------------------------------------------------------------------------------------
# R4_F2f_downstream_breast
# --------------------------------------------------------------------------------------------------


def recover_counts_from_normalised(adata):
    """
    Integer UMI counts of the breast slides, recovered from the stored library-size normalised X.

    The stored X sums to 1e4 per spot, so counts = X * total_counts / 1e4 (rounded). Call before
    `sc.pp.normalize_total` and `sc.pp.log1p`.

    Parameters
    ----------
    adata : AnnData
        Slides as stored, with `obs["total_counts"]`.

    Returns
    -------
    scipy.sparse.csr_matrix
        Spots x genes counts (float32).
    """
    tc = adata.obs["total_counts"].values.astype(np.float64)
    C = csr_matrix(adata.X.multiply((tc / 1e4)[:, None]))
    C.data = np.round(C.data)
    C.eliminate_zeros()
    return C.astype(np.float32)


def lognorm_adata_from_counts(C, template):
    """
    Log-normalised AnnData (normalize_total 1e4 + log1p) from a count matrix, with the template's obs, var and spatial.

    Parameters
    ----------
    C : sparse matrix
        Spots x genes counts.
    template : AnnData
        Data whose obs, var_names and `obsm["spatial"]` are copied.

    Returns
    -------
    AnnData
    """
    b = ad.AnnData(
        X=csr_matrix(C, dtype=np.float32),
        obs=template.obs.copy(),
        var=pd.DataFrame(index=template.var_names),
    )
    b.obsm["spatial"] = np.asarray(template.obsm["spatial"]).copy()
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    b.X = csr_matrix(b.X)
    return b


def binomial_count_split(C, seed, eps=0.5):
    """
    Binomial count splitting (data thinning): C_A ~ Bin(C, eps), C_B = C - C_A.

    Parameters
    ----------
    C : sparse matrix
        Integer counts.
    seed : int
        Seed of the split.
    eps : float
        Fraction of counts sent to half A.

    Returns
    -------
    tuple of csr_matrix
        Counts of half A and half B.
    """
    rng = np.random.default_rng(seed)
    C = csr_matrix(C)
    A = C.copy()
    A.data = rng.binomial(C.data.astype(np.int64), eps).astype(np.float32)
    B = C.copy()
    B.data = (C.data - A.data).astype(np.float32)
    A.eliminate_zeros()
    B.eliminate_zeros()
    return A, B


def sparse_gene_mean_sd(X):
    """
    Per-gene mean and population SD (ddof=0) of a sparse matrix.

    Parameters
    ----------
    X : sparse matrix
        Spots x genes.

    Returns
    -------
    tuple of np.ndarray
        Gene means and SDs (float64).
    """
    X = csr_matrix(X)
    m = np.asarray(X.mean(axis=0)).ravel().astype(np.float64)
    m2 = np.asarray(X.multiply(X).mean(axis=0)).ravel().astype(np.float64)
    return m, np.sqrt(np.maximum(m2 - m**2, 0))


def present_genes_in_order(genes, var_names):
    """
    Signature genes present in the data, de-duplicated, in input order.

    Parameters
    ----------
    genes : list
        Gene symbols.
    var_names : pd.Index
        Genes of the data.

    Returns
    -------
    list
    """
    return [g for g in dict.fromkeys(genes) if g in var_names]


def cv_weights_in_order(adata, genes):
    """
    CV weights of `em.tl.infer_gene_weights`, keyed in input gene order (the package returns set order).

    Parameters
    ----------
    adata : AnnData
        Log-normalised data.
    genes : list
        Signature genes.

    Returns
    -------
    dict
        Gene -> CV weight.
    """
    w = em.tl.infer_gene_weights(adata, genes)
    return {g: w[g] for g in genes if g in w}


def signature_weight_matrix(var_names, weight_dicts):
    """
    Sparse genes x signatures matrix; column j = weight dict j normalised to sum |w| = 1 (EnrichMap Eq. 1).

    Parameters
    ----------
    var_names : pd.Index
        Genes of the data (absent genes are dropped).
    weight_dicts : list of dict
        Gene -> weight per signature.

    Returns
    -------
    scipy.sparse.csc_matrix
    """
    pos = pd.Index(var_names)
    rows, cols, vals = [], [], []
    for j, w in enumerate(weight_dicts):
        g = [k for k in w if k in pos]
        v = np.array([w[k] for k in g], dtype=np.float64)
        s = np.abs(v).sum()
        rows += list(pos.get_indexer(g))
        cols += [j] * len(g)
        vals += list(v / s if s > 0 else v)
    return csc_matrix((vals, (rows, cols)), shape=(len(pos), len(weight_dicts)))


def index_weight_matrix(n_genes, idx_lists, weight_lists):
    """
    Sparse genes x sets matrix from gene-index and weight lists, each column normalised to sum |w| = 1.

    Parameters
    ----------
    n_genes : int
        Number of genes.
    idx_lists : list of array-like
        Gene indices per set.
    weight_lists : list of array-like
        Weights per set.

    Returns
    -------
    scipy.sparse.csc_matrix
    """
    rows, cols, vals = [], [], []
    for j, (ix, w) in enumerate(zip(idx_lists, weight_lists)):
        w = np.asarray(w, dtype=np.float64)
        rows += list(ix)
        cols += [j] * len(ix)
        vals += list(w / np.abs(w).sum())
    return csc_matrix((vals, (rows, cols)), shape=(n_genes, len(idx_lists)))


def linear_scores(X, W):
    """
    Raw linear scores R = X W (EnrichMap Eq. 1 before standardisation).

    Parameters
    ----------
    X : sparse matrix
        Spots x genes.
    W : sparse matrix
        Genes x signatures.

    Returns
    -------
    np.ndarray
        Spots x signatures (float64).
    """
    R = csr_matrix(X).dot(W)
    R = R.toarray() if issparse(R) else np.asarray(R)
    return R.astype(np.float64)


def signed_weights_vectorised(X, var_names, up, down, cv):
    """
    Vectorised `em.tl.build_signed_weights`: |w| = CV, sign = Spearman(g, mean(up)) - Spearman(g, mean(down)),
    divided by the size of the gene's own group.

    Parameters
    ----------
    X : sparse matrix
        Spots x genes (log-normalised).
    var_names : pd.Index
        Genes of the data.
    up, down : list
        Up and down genes.
    cv : np.ndarray
        CV of every gene of the data.

    Returns
    -------
    dict
        Gene -> signed weight.
    """
    var = pd.Index(var_names)
    up_p = [g for g in up if g in var]
    down_p = [g for g in down if g in var]
    iu, idn = var.get_indexer(up_p), var.get_indexer(down_p)
    Xc = csc_matrix(X)
    Xu = Xc[:, iu].toarray().astype(np.float64)
    Xd = Xc[:, idn].toarray().astype(np.float64)
    ref_u, ref_d = rankdata(Xu.mean(1)), rankdata(Xd.mean(1))
    Rg = rankdata(np.hstack([Xu, Xd]), axis=0)
    Rg = Rg - Rg.mean(0)
    dr = []
    for v in (ref_u, ref_d):
        v = v - v.mean()
        with np.errstate(divide="ignore", invalid="ignore"):
            r = (Rg * v[:, None]).sum(0) / np.sqrt((Rg**2).sum(0) * (v**2).sum())
        dr.append(np.nan_to_num(r))
    dr = dr[0] - dr[1]
    cvv = cv[np.r_[iu, idn]]
    n_up, n_dn = len(up_p), len(down_p)
    return {
        g: cvv[k] * dr[k] / (n_up if k < n_up else n_dn)
        for k, g in enumerate(up_p + down_p)
    }


def mean_expression_bins(gene_means, n_bins=25):
    """
    Bins of mean-expression rank (as Scanpy `score_genes`).

    Parameters
    ----------
    gene_means : np.ndarray
        Mean expression per gene.
    n_bins : int
        Number of bins.

    Returns
    -------
    np.ndarray
        Bin per gene.
    """
    return pd.qcut(
        pd.Series(gene_means).rank(method="first"), n_bins, labels=False
    ).values


def matched_control_index(bins, sig_idx, n_sets, seed, exclude=()):
    """
    Expression-matched random gene sets: the gene matched to signature gene i is drawn from i's expression bin,
    without replacement within a set (unless the pool is too small), excluding the signature genes and `exclude`.

    Parameters
    ----------
    bins : np.ndarray
        Expression bin per gene.
    sig_idx : array-like
        Gene indices of the signature.
    n_sets : int
        Number of sets.
    seed : int
        Seed.
    exclude : array-like
        Further gene indices never drawn.

    Returns
    -------
    np.ndarray
        n_sets x len(sig_idx) gene indices.
    """
    rng = np.random.default_rng(seed)
    sig_idx = np.asarray(sig_idx)
    excl = np.zeros(len(bins), bool)
    excl[sig_idx] = True
    excl[np.asarray(list(exclude), dtype=int)] = True
    out = np.empty((n_sets, len(sig_idx)), dtype=np.int64)
    sig_bins = bins[sig_idx]
    pools = {b: np.where((bins == b) & ~excl)[0] for b in np.unique(sig_bins)}
    for s in range(n_sets):
        for b, pool in pools.items():
            pos = np.where(sig_bins == b)[0]
            out[s, pos] = rng.choice(pool, size=len(pos), replace=len(pos) > len(pool))
    return out


def weight_transfer_control_column(n_genes, sig_idx, sig_w, ctrl_idx):
    """
    Weight-transfer control: the control gene matched to signature gene i receives weight w_i; the mean control
    score over the sets is X @ (this column), on the scale of the signature's raw score.

    Parameters
    ----------
    n_genes : int
        Number of genes.
    sig_idx : array-like
        Gene indices of the signature.
    sig_w : array-like
        Weights of the signature genes.
    ctrl_idx : np.ndarray
        Control gene indices (n_sets x len(sig_idx)).

    Returns
    -------
    np.ndarray
        Weight per gene.
    """
    sig_w = np.asarray(sig_w, dtype=np.float64)
    col = np.zeros(n_genes)
    for s in range(ctrl_idx.shape[0]):
        np.add.at(col, ctrl_idx[s], sig_w)
    return col / (ctrl_idx.shape[0] * np.abs(sig_w).sum())


def reference_control_score(X, var_names, weights, bins, n_ctrl, seed):
    """
    Expression-matched, weight-transferred control score (mean of `n_ctrl` sets) for the reference-GAM arm.

    Parameters
    ----------
    X : sparse matrix
        Spots x genes.
    var_names : pd.Index
        Genes of the data.
    weights : dict
        Gene -> weight of the signature.
    bins : np.ndarray
        Expression bin per gene.
    n_ctrl : int
        Number of control sets.
    seed : int
        Seed for drawing the control genes.

    Returns
    -------
    np.ndarray
        Raw Eq. 1 control score per spot.
    """
    vpos = pd.Index(var_names)
    g = [x for x in weights if x in vpos]
    idx = vpos.get_indexer(g)
    ci = matched_control_index(bins, idx, n_ctrl, seed=seed)
    col = weight_transfer_control_column(len(vpos), idx, [weights[x] for x in g], ci)
    return linear_scores(X, csc_matrix(col[:, None]))[:, 0]


def spatial_knn_connectivity(coords, n_neighbors=6):
    """
    Binary spatial kNN connectivity, built as inside `em.tl.score`.

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates.
    n_neighbors : int
        Number of neighbours.

    Returns
    -------
    scipy.sparse.csr_matrix
    """
    tmp = ad.AnnData(X=np.zeros((len(coords), 1), dtype=np.float32))
    tmp.obsm["spatial"] = np.asarray(coords)
    sq.gr.spatial_neighbors(
        tmp, n_neighs=n_neighbors, coord_type="generic", key_added="spatial"
    )
    return csr_matrix(tmp.obsp["spatial_connectivities"])


def neighbour_average(Y, conn):
    """
    EnrichMap Eq. 3: neighbourhood average with the row-normalised adjacency.

    Parameters
    ----------
    Y : np.ndarray
        Values per spot (vector or spots x k).
    conn : sparse matrix
        Spatial connectivity.

    Returns
    -------
    np.ndarray
    """
    Y = np.asarray(Y, dtype=np.float64)
    d = np.maximum(np.asarray(conn.sum(axis=1)).ravel(), 1e-10)
    return conn.dot(Y) / (d[:, None] if Y.ndim == 2 else d)


class SpatialGAMSmoother:
    """
    The spatial-covariate GAM of EnrichMap, LinearGAM(te(0, 1)) with pyGAM defaults, written as a linear smoother.

    fitted = H y and residual = (I - H) y with H = B (B'B + P)^-1 B'; `edof` is trace(H).

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates of one slide.
    """

    def __init__(self, coords):
        coords = np.asarray(coords, dtype=np.float64)
        gam = LinearGAM(te(0, 1)).fit(coords, coords[:, 0])
        B = gam._modelmat(coords)
        B = B.toarray() if hasattr(B, "toarray") else np.asarray(B)
        P = gam._P()
        P = P.toarray() if hasattr(P, "toarray") else np.asarray(P)
        self.B = B
        self._A = B.T @ B + P + np.sqrt(np.finfo(float).eps) * np.eye(B.shape[1])
        self.edof = float(np.trace(np.linalg.solve(self._A, B.T @ B)))

    def fitted(self, Y):
        """GAM fit of each column of Y."""
        Y = np.asarray(Y, dtype=np.float64)
        return self.B @ np.linalg.solve(self._A, self.B.T @ Y)

    def residual(self, Y):
        """GAM residual of each column of Y."""
        return np.asarray(Y, dtype=np.float64) - self.fitted(Y)

    def r2(self, y):
        """Fraction of the centred variance of y absorbed by the GAM."""
        y = np.asarray(y, dtype=np.float64)
        yc = y - y.mean()
        return float(1 - ((y - self.fitted(y)) ** 2).sum() / (yc**2).sum())


class SlideOperators:
    """
    Per-slide EnrichMap post-processing of raw scores, as `em.tl.score` with `batch_key`.

    Arms: "raw" (per-slide z-score), "smooth" (z-score + 6-NN smoothing), "gam" (smoothing + GAM residual) and
    "ref" (smoothing minus the GAM fit of an expression-matched control score).

    Parameters
    ----------
    adata : AnnData
        Data with `obs[batch_key]` and `obsm["spatial"]`.
    batch_key : str
        Slide column.
    n_neighbors : int
        Neighbours for smoothing.
    """

    def __init__(self, adata, batch_key="batch", n_neighbors=6):
        b = adata.obs[batch_key].astype(str).values
        self.batches = list(pd.unique(b))
        self.idx = {k: np.where(b == k)[0] for k in self.batches}
        coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
        self.conn = {
            k: spatial_knn_connectivity(coords[i], n_neighbors)
            for k, i in self.idx.items()
        }
        self.gam = {k: SpatialGAMSmoother(coords[i]) for k, i in self.idx.items()}

    def _per(self, R, fn):
        R = np.asarray(R, dtype=np.float64)
        out = np.empty_like(R)
        for k, i in self.idx.items():
            out[i] = fn(k, R[i])
        return out

    def zscore(self, R):
        """Per-slide z-score of a score (or of each column)."""

        def f(k, x):
            sd = x.std(0)
            return (x - x.mean(0)) / np.where(sd > 0, sd, 1.0)

        return self._per(R, f)

    def smooth(self, R):
        """Per-slide neighbourhood average."""
        return self._per(R, lambda k, x: neighbour_average(x, self.conn[k]))

    def process(self, raw, arm, ctrl=None):
        """
        Apply one post-processing arm to raw score(s).

        Parameters
        ----------
        raw : np.ndarray
            Raw scores (spots,) or (spots, k).
        arm : str
            "raw", "smooth", "gam" or "ref".
        ctrl : np.ndarray or None
            Raw control scores (same shape), for "ref".

        Returns
        -------
        np.ndarray
        """
        R = np.asarray(raw, dtype=np.float64)
        Z = self.zscore(R)
        if arm == "raw":
            return Z
        SZ = self.smooth(Z)
        if arm == "smooth":
            return SZ
        if arm == "gam":
            return self._per(SZ, lambda k, x: self.gam[k].residual(x))
        C = np.asarray(ctrl, dtype=np.float64)
        out = np.empty_like(R)
        for k, i in self.idx.items():
            sd = R[i].std(0)
            cz = (C[i] - C[i].mean(0)) / np.where(sd > 0, sd, 1.0)
            v = SZ[i] - self.gam[k].fitted(neighbour_average(cz, self.conn[k]))
            out[i] = v - v.mean(0)
        return out

    def gam_r2(self, R):
        """Per slide: fraction of the variance of the smoothed z-score absorbed by the GAM."""
        SZ = self.smooth(self.zscore(R))
        return {k: self.gam[k].r2(SZ[i]) for k, i in self.idx.items()}


def null_model_score_matrix(adata, gene_weights=None, gene_set=None, batch_key="batch"):
    """
    EnrichMap (null) scores of many signatures, scored on a light copy of `adata`.

    Parameters
    ----------
    adata : AnnData
        Log-normalised data with `obs[batch_key]` and `obsm["spatial"]`.
    gene_weights : dict or None
        Signature -> {gene: weight} with fixed weights (inherited by the control sets).
    gene_set : dict or None
        Signature -> {"up": [...], "down": [...]} for signed signatures (weights refitted for every control set).
    batch_key : str
        Slide column.

    Returns
    -------
    pd.DataFrame
        Spots x signatures (fixed-weight signatures first, then signed ones).
    """
    light = ad.AnnData(
        X=adata.X,
        obs=adata.obs[[batch_key]].copy(),
        var=pd.DataFrame(index=adata.var_names),
    )
    light.obsm["spatial"] = np.asarray(adata.obsm["spatial"])
    names = []
    if gene_weights:
        em.tl.score(
            light,
            gene_weights=gene_weights,
            batch_key=batch_key,
            null_model=True,
            correct_spatial_covariates=False,
        )
        names += list(gene_weights)
    if gene_set:
        em.tl.score(
            light,
            gene_set=gene_set,
            batch_key=batch_key,
            null_model=True,
            correct_spatial_covariates=False,
        )
        names += list(gene_set)
    return pd.DataFrame(
        {k: light.obs[f"{k}_score"].values for k in names}, index=adata.obs_names
    )


def libpysal_knn_lag_matrix(coords, k=6):
    """
    Row-standardised libpysal KNN(k) weights, as `em.pl.cross_moran_scatter` uses for the spatial lag.

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates.
    k : int
        Number of neighbours.

    Returns
    -------
    scipy.sparse.csr_matrix
    """
    from libpysal.weights import KNN

    W = KNN.from_array(np.asarray(coords), k=k)
    W.transform = "r"
    return csr_matrix(W.sparse)


def pearson_with_columns(x, Y):
    """
    Pearson r between vector x and each column of Y.

    Parameters
    ----------
    x : np.ndarray
        Vector.
    Y : np.ndarray
        Vector or matrix.

    Returns
    -------
    np.ndarray
    """
    x = np.asarray(x, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64)
    if Y.ndim == 1:
        Y = Y[:, None]
    xc = x - x.mean()
    Yc = Y - Y.mean(0)
    return (xc @ Yc) / np.sqrt((xc @ xc) * (Yc**2).sum(0))


def spearman_corr_matrix(M):
    """
    Spearman correlation matrix of the columns of M.

    Parameters
    ----------
    M : np.ndarray
        Spots x scores.

    Returns
    -------
    np.ndarray
    """
    return np.corrcoef(rankdata(np.asarray(M), axis=0), rowvar=False)


class DutilleulModifiedTTest:
    """
    Modified t-test for the correlation of two spatially autocorrelated variables (Clifford, Richardson & Hemon
    1989; Dutilleul 1993; as SpatialPack::modified.ttest), with the distance classes of one set of spots
    pre-computed for repeated tests.

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates.
    nclass : int
        Equal-width distance classes.
    """

    def __init__(self, coords, nclass=13):
        D = pdist(np.asarray(coords, dtype=np.float64))
        edges = np.linspace(0, D.max() * (1 + 1e-9), nclass + 1)
        self.cls = np.digitize(D, edges[1:-1]).astype(np.int32)
        self.i, self.j = np.triu_indices(len(coords), 1)
        self.Nk = np.bincount(self.cls, minlength=nclass).astype(np.float64)
        self.n, self.nclass = len(coords), nclass

    def test(self, x, y):
        """
        Test the correlation of x and y.

        Parameters
        ----------
        x, y : np.ndarray
            Values on the spots.

        Returns
        -------
        dict
            r, n, effective sample size n_eff, t and two-sided p.
        """
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        n = self.n
        xc, yc = x - x.mean(), y - y.mean()
        Cx = np.bincount(
            self.cls, weights=xc[self.i] * xc[self.j], minlength=self.nclass
        ) / np.maximum(self.Nk, 1)
        Cy = np.bincount(
            self.cls, weights=yc[self.i] * yc[self.j], minlength=self.nclass
        ) / np.maximum(self.Nk, 1)
        c0x, c0y = (xc @ xc) / n, (yc @ yc) / n
        den = n * c0x * c0y + 2 * np.sum(self.Nk * Cx * Cy)
        n_eff = 1 + (n * c0x) * (n * c0y) / den if den > 0 else np.nan
        r = float((xc @ yc) / np.sqrt((xc @ xc) * (yc @ yc)))
        dfree = max(n_eff - 2, 1)
        t = r * np.sqrt(dfree / max(1 - r**2, 1e-12))
        return {
            "r": r,
            "n": n,
            "n_eff": float(n_eff),
            "t": float(t),
            "p": float(2 * t_dist.sf(abs(t), dfree)),
        }


def steiger_z(r12, r13, r23, n):
    r_bar = (r12 + r13) / 2
    f = (1 - r23) / (2 * (1 - r_bar**2))
    h = (1 - f * r_bar**2) / (1 - r_bar**2)
    z1, z2 = np.arctanh(r12), np.arctanh(r13)
    Z = (z1 - z2) * np.sqrt((n - 3) / (2 * (1 - r23) * h))
    p = 2 * (1 - stats.norm.cdf(abs(Z)))
    return Z, p


def kmeans_spatial_blocks(coords, n_blocks=30, seed=0):
    """
    Spatial blocks of one slide (k-means of the coordinates) for the block bootstrap.

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates.
    n_blocks : int
        Number of blocks.
    seed : int
        k-means seed.

    Returns
    -------
    np.ndarray
        Block label per spot.
    """
    return (
        KMeans(n_clusters=n_blocks, n_init=10, random_state=seed)
        .fit(np.asarray(coords))
        .labels_
    )


def block_bootstrap_pearson(gt, S, blocks, n_boot=1000, seed=0):
    """
    Spatial block bootstrap of Pearson r(gt, S[:, j]) for every column j (paired across columns).

    Parameters
    ----------
    gt : np.ndarray
        Ground truth per spot.
    S : np.ndarray
        Spots x test scores.
    blocks : np.ndarray
        Block label per spot.
    n_boot : int
        Resamples.
    seed : int
        Seed.

    Returns
    -------
    np.ndarray
        n_boot x k correlations.
    """
    rng = np.random.default_rng(seed)
    ub = np.unique(blocks)
    members = [np.where(blocks == b)[0] for b in ub]
    out = np.empty((n_boot, S.shape[1]))
    for bi in range(n_boot):
        pick = rng.integers(0, len(ub), len(ub))
        idx = np.concatenate([members[p] for p in pick])
        out[bi] = pearson_with_columns(gt[idx], S[idx])
    return out


def partial_spearman_rank_residuals(x, y, Z):
    """
    Partial Spearman correlation of x and y given the columns of Z (Pearson r of the rank residuals).

    Parameters
    ----------
    x, y : np.ndarray
        Variables to correlate.
    Z : np.ndarray
        Covariates (spots x k).

    Returns
    -------
    float
    """
    rx, ry = rankdata(x), rankdata(y)
    Zr = np.column_stack([np.ones(len(x))] + [rankdata(z) for z in Z.T])
    ex = rx - Zr @ np.linalg.lstsq(Zr, rx, rcond=None)[0]
    ey = ry - Zr @ np.linalg.lstsq(Zr, ry, rcond=None)[0]
    return pearsonr(ex, ey)[0]


def hrd_half_conditions(adata, hrd_gt, hrd_pos, hrd_neg):
    """
    HRD conditions of one count-split half (weights C and D inferred on that half).

    Parameters
    ----------
    adata : AnnData
        Log-normalised half.
    hrd_gt, hrd_pos, hrd_neg : dict
        Jacobson weights: all, positive and negative genes.

    Returns
    -------
    dict
        Condition -> {gene: weight}.
    """
    up, dn = list(hrd_pos), list(hrd_neg)
    return {
        "GT": hrd_gt,
        "up_part": hrd_pos,
        "down_part": {g: -v for g, v in hrd_neg.items()},
        "B_pos": hrd_pos,
        "U_unweighted": {g: 1.0 for g in up},
        "C_inferred": em.tl.infer_gene_weights(adata, up),
        "D_signed": em.tl.build_signed_weights(
            adata, up_genes=up, down_genes=dn, score_key="s"
        )["s"],
        "E_naive_signed": {
            **{g: 1.0 / len(up) for g in up},
            **{g: -1.0 / len(dn) for g in dn},
        },
    }


def evaluate_hrd_against_gt(
    gt_by_arm, tests_by_arm, label, batch, slides, patient, blocks, n_boot=1000
):
    """
    Per-slide correlation of every HRD test score with a ground truth, and the D-B, B-C, D-C and E-D contrasts.

    Parameters
    ----------
    gt_by_arm : dict or np.ndarray
        GT per arm (dict) or one GT used for every arm.
    tests_by_arm : dict
        Arm -> DataFrame of test scores (spots x methods).
    label : str
        GT name written to the output.
    batch : np.ndarray
        Slide per spot.
    slides : list
        Slides.
    patient : dict
        Slide -> patient label.
    blocks : dict
        Slide -> spatial block per spot of that slide.
    n_boot : int
        Spatial block-bootstrap resamples for the contrasts.

    Returns
    -------
    tuple of pd.DataFrame
        Per-slide Pearson r and Spearman rho; contrasts with 95% bootstrap CIs and Steiger Z.
    """
    rows, diffs = [], []
    for arm, T in tests_by_arm.items():
        g_all = gt_by_arm[arm] if isinstance(gt_by_arm, dict) else gt_by_arm
        for b in slides:
            m = batch == b
            g = np.asarray(g_all)[m]
            Tm = T.loc[m]
            rs = pearson_with_columns(g, Tm.values)
            for c, r in zip(Tm.columns, rs):
                rows.append(
                    {
                        "GT": label,
                        "arm": arm,
                        "slide": b,
                        "patient": patient[b],
                        "method": c,
                        "r": r,
                        "rho": spearmanr(g, Tm[c].values)[0],
                    }
                )
            cols = ["D_signed", "B_pos", "C_inferred", "E_naive_signed"]
            bb = block_bootstrap_pearson(
                g, Tm[cols].values, blocks[b], n_boot=n_boot, seed=1
            )
            ci = dict(zip(cols, range(4)))
            rr = dict(zip(Tm.columns, rs))
            for x, y in (
                ("D_signed", "B_pos"),
                ("B_pos", "C_inferred"),
                ("D_signed", "C_inferred"),
                ("E_naive_signed", "D_signed"),
            ):
                dd = bb[:, ci[x]] - bb[:, ci[y]]
                Z, p = steiger_z(
                    rr[x], rr[y], np.corrcoef(Tm[x], Tm[y])[0, 1], int(m.sum())
                )
                diffs.append(
                    {
                        "GT": label,
                        "arm": arm,
                        "slide": b,
                        "patient": patient[b],
                        "contrast": f"{x} - {y}",
                        "diff": rr[x] - rr[y],
                        "boot lo": np.quantile(dd, 0.025),
                        "boot hi": np.quantile(dd, 0.975),
                        "Steiger Z (n=spots)": Z,
                        "Steiger p": p,
                    }
                )
    return pd.DataFrame(rows), pd.DataFrame(diffs)


def hrd_r_table(res, gt, order, arms, value="r"):
    """
    Method x (patient, arm) table of correlations with one GT (mean over count-split replicates).

    Parameters
    ----------
    res : pd.DataFrame
        Output rows of `evaluate_hrd_against_gt`.
    gt : str
        GT label.
    order : list
        Row (method) order.
    arms : list
        Column (arm) order.
    value : str
        Column to tabulate ("r" or "rho").

    Returns
    -------
    pd.DataFrame
    """
    sub = (
        res[res.GT == gt]
        .groupby(["arm", "patient", "method"], as_index=False)[value]
        .mean()
    )
    t = sub.pivot_table(index="method", columns=["patient", "arm"], values=value)
    return t.reindex([c for c in order if c in t.index]).reindex(columns=arms, level=1)


def hrd_contrast_table(dif, gt, arms):
    """
    Contrast x patient table of "diff [lo, hi]" strings per arm (mean over count-split replicates).

    Parameters
    ----------
    dif : pd.DataFrame
        Contrast rows of `evaluate_hrd_against_gt`.
    gt : str
        GT label.
    arms : list
        Column (arm) order.

    Returns
    -------
    pd.DataFrame
    """
    agg = (
        dif[dif.GT == gt]
        .groupby(["contrast", "patient", "arm"])
        .agg(
            diff=("diff", "mean"),
            lo=("boot lo", "mean"),
            hi=("boot hi", "mean"),
        )
        .reset_index()
    )
    agg["cell"] = [
        f"{d:+.3f} [{l:+.3f},{h:+.3f}]"
        for d, l, h in zip(agg["diff"], agg["lo"], agg["hi"])
    ]
    return agg.pivot_table(
        index=["contrast", "patient"], columns="arm", values="cell", aggfunc="first"
    ).reindex(columns=arms)


def mean_r_with_gt(res, gt, method, arm, slide=None):
    """
    Pearson r of one test score with one GT (mean over count-split replicates; over both slides when slide is None).

    Parameters
    ----------
    res : pd.DataFrame
        Output rows of `evaluate_hrd_against_gt`.
    gt, method, arm : str
        GT label, test column and arm.
    slide : str or None
        Slide.

    Returns
    -------
    float
    """
    q = res[(res.GT == gt) & (res.method == method) & (res.arm == arm)]
    if slide is not None:
        q = q[q.slide == slide]
    return float(q.r.mean())


def plot_metric_pair_with_markers(table, x, y, ax, colours, marker_col=None):
    """
    One metric against another, one point per row, coloured by method.

    Parameters
    ----------
    table : pd.DataFrame
        One row per point, with a "method" column and the metric columns.
    x, y : str
        Metric columns.
    ax : matplotlib Axes
    colours : dict
        Method -> colour.
    marker_col : str or None
        Column holding a marker per row (default circles).

    Returns
    -------
    matplotlib Axes
    """
    for _, r in table.iterrows():
        ax.errorbar(
            r[x],
            r[y],
            fmt=r[marker_col] if marker_col else "o",
            color=colours[r["method"]],
            ms=5,
            mec="white",
            mew=0.4,
            elinewidth=0.6,
            capsize=0,
            label=r["method"],
        )
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.grid(False)
    return ax


def plot_benchmark_summary_circles(table, ax, colours):
    """
    Methods x metrics summary: circles sized and shaded by each metric scaled to 0-1 across methods (best = 1), raw
    values annotated, plus a mean-rank column.

    Parameters
    ----------
    table : pd.DataFrame
        Index = methods, columns = metrics (higher is better).
    ax : matplotlib Axes
    colours : dict
        Method -> colour of the row label.

    Returns
    -------
    matplotlib Axes
    """
    scaled = table.copy().astype(float)
    for c in table.columns:
        col = table[c].astype(float)
        rng_ = col.max() - col.min()
        scaled[c] = (col - col.min()) / rng_ if rng_ > 0 else col * 0 + 1
    ranks = scaled.rank(ascending=False).mean(axis=1)
    cmap = plt.get_cmap("Purples")
    for i, m in enumerate(table.index):
        for j, c in enumerate(table.columns):
            s = scaled.loc[m, c]
            ax.scatter(
                j, i, s=30 + 220 * s, color=cmap(0.25 + 0.75 * s), edgecolors="none"
            )
            ax.text(
                j,
                i,
                f"{table.loc[m, c]:.2g}",
                ha="center",
                va="center",
                fontsize=5,
                color="white" if s > 0.6 else "black",
            )
        ax.text(
            len(table.columns) + 0.2,
            i,
            f"{ranks[m]:.1f}",
            ha="left",
            va="center",
            fontsize=6,
        )
    ax.set_xticks(range(len(table.columns)))
    ax.set_xticklabels(table.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(table.index)))
    ax.set_yticklabels(table.index)
    for lab in ax.get_yticklabels():
        lab.set_color(colours[lab.get_text()])
    ax.text(len(table.columns) + 0.2, -0.9, "mean rank", fontsize=6, ha="left")
    ax.set_xlim(-0.6, len(table.columns) + 0.9)
    ax.set_ylim(len(table.index) - 0.4, -1.2)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    return ax


# --------------------------------------------------------------------------------------------------
# R4_F2g_interpretable_scores
# --------------------------------------------------------------------------------------------------


def simulate_visium_domains(
    seed=0,
    n_genes=3000,
    n_signature_genes=20,
    n_clusters=4,
    target_cluster=1,
    n_spots_range=(1500, 4000),
    spot_spacing_um=100,
    spot_diameter_um=55,
    mean_counts_per_spot=5000,
    dispersion_range=(0.1, 2.0),
    fold_change_range=(1.5, 5.0),
    dropout_rate=0.1,
    domain="voronoi",
    domain_width_spots=1.5,
):
    """
    Simulate a Visium slide with the simulate_visium generator of 05_01_simulated_benchmarking, with a choice of
    target domain.

    With domain="voronoi" the random stream and the output are those of the original function (Sup. Fig. 1). With
    domain="stripe" the target is a thin sinusoidal band about 2 x domain_width_spots spots wide (layer-like), drawn
    from an independent random stream so that the original stream is not perturbed.

    Parameters
    ----------
    seed : int
        Simulation seed.
    n_genes, n_signature_genes, n_clusters, target_cluster : int
        Genes, signature genes, Voronoi clusters and the target cluster (0-based).
    n_spots_range : tuple
        Range of the number of spots.
    spot_spacing_um, spot_diameter_um : float
        Visium spot spacing and diameter.
    mean_counts_per_spot : float
        Mean library size.
    dispersion_range, fold_change_range : tuple
        Negative binomial dispersions and signature fold changes.
    dropout_rate : float
        Fraction of counts set to zero.
    domain : str
        "voronoi" or "stripe".
    domain_width_spots : float
        Half-width of the stripe in spots.

    Returns
    -------
    tuple
        (adata, signature genes, fold changes, target label); adata is log-normalised.
    """
    rng = np.random.default_rng(seed)
    rng_ext = np.random.default_rng(10_000_019 + seed)

    n_spots_target = rng.integers(n_spots_range[0], n_spots_range[1] + 1)
    grid_side = int(np.ceil(np.sqrt(n_spots_target * 0.866)))
    capture_area_um = grid_side * spot_spacing_um
    n_cols = int(capture_area_um / spot_spacing_um)
    n_rows = int(capture_area_um / (spot_spacing_um * np.sqrt(3) / 2))
    x_coords, y_coords = [], []
    for row in range(n_rows):
        for col in range(n_cols):
            x = col * spot_spacing_um + (row % 2) * (spot_spacing_um / 2)
            y = row * spot_spacing_um * np.sqrt(3) / 2
            if x <= capture_area_um and y <= capture_area_um:
                x_coords.append(x)
                y_coords.append(y)
    x_coords = np.array(x_coords)
    y_coords = np.array(y_coords)
    n_spots = len(x_coords)
    coords = np.column_stack([x_coords, y_coords])

    margin = capture_area_um * 0.15
    seed_coords = np.column_stack(
        [
            rng.uniform(margin, capture_area_um - margin, n_clusters),
            rng.uniform(margin, capture_area_um - margin, n_clusters),
        ]
    )
    dists = cdist(coords, seed_coords)
    dists += rng.normal(0, spot_spacing_um * 0.8, dists.shape)
    cluster_ids = np.argmin(dists, axis=1)
    cluster_labels = np.array([f"cluster_{i + 1}" for i in cluster_ids])
    target_label = f"cluster_{target_cluster + 1}"

    if domain == "stripe":
        L = capture_area_um
        phase = rng_ext.uniform(0, 2 * np.pi)
        amp = rng_ext.uniform(0.08, 0.15) * L
        y0 = rng_ext.uniform(0.35, 0.65) * L
        curve = y0 + amp * np.sin(2 * np.pi * coords[:, 0] / (0.8 * L) + phase)
        in_target = np.abs(coords[:, 1] - curve) <= domain_width_spots * spot_spacing_um
        cluster_labels = np.where(
            in_target,
            target_label,
            np.where(cluster_labels == target_label, "cluster_bg", cluster_labels),
        )

    gene_names = [f"GENE_{i + 1:05d}" for i in range(n_genes)]
    sig_indices = rng.choice(n_genes, n_signature_genes, replace=False)
    fold_changes = rng.uniform(
        fold_change_range[0], fold_change_range[1], n_signature_genes
    )
    order = np.argsort(-fold_changes)
    sig_indices = sig_indices[order]
    fold_changes = fold_changes[order]
    sig_genes = [gene_names[i] for i in sig_indices]

    gene_means = np.exp(rng.normal(1.5, 1.0, n_genes))
    gene_means = np.clip(gene_means, 0.1, 100)
    dispersions = rng.uniform(dispersion_range[0], dispersion_range[1], n_genes)
    lib_size_factors = np.exp(rng.normal(0, 0.3, n_spots))
    lib_size_factors *= mean_counts_per_spot / (
        gene_means.sum() * lib_size_factors.mean()
    )

    X = np.zeros((n_spots, n_genes), dtype=np.float32)
    target_mask = cluster_labels == target_label
    for g in range(n_genes):
        mu = gene_means[g] * lib_size_factors
        if g in sig_indices:
            fc_idx = np.where(sig_indices == g)[0][0]
            mu[target_mask] *= fold_changes[fc_idx]
        r = dispersions[g]
        p = r / (r + mu)
        counts = rng.negative_binomial(r, p)
        dropout_mask = rng.random(n_spots) < dropout_rate
        counts[dropout_mask] = 0
        X[:, g] = counts

    adata = ad.AnnData(X=csr_matrix(X))
    adata.var_names = gene_names
    adata.obs_names = [f"spot_{i}" for i in range(n_spots)]
    adata.obs["cluster"] = pd.Categorical(cluster_labels)
    adata.obsm["spatial"] = coords
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata, sig_genes, fold_changes, target_label


def simulated_slide_with_null_signature(seed, domain):
    """
    One simulated slide with its 20-gene target signature and a 20-gene random (null) signature.

    Parameters
    ----------
    seed : int
        Simulation seed; the null signature is drawn with seed + 7.
    domain : str
        "voronoi" (Sup. Fig. 1) or "stripe" (thin sinusoidal band).

    Returns
    -------
    tuple
        (adata, {"target": genes, "null signature": genes}, target label).
    """
    sim, sig_genes, _, target_label = simulate_visium_domains(seed=seed, domain=domain)
    rng = np.random.default_rng(seed + 7)
    null_genes = list(
        rng.choice([g for g in sim.var_names if g not in sig_genes], 20, replace=False)
    )
    return sim, {"target": list(sig_genes), "null signature": null_genes}, target_label


def sign_calls(method, s, p_emp=None, label=None):
    """
    Binary calls of one score.

    Parameters
    ----------
    method : str
        Method label.
    s : np.ndarray
        Score.
    p_emp : np.ndarray or None
        Empirical p-value of the null model (used for the null model only).
    label : np.ndarray or None
        GSDensity antimode label, True = positive (used for GSDensity only; None when compute.cell.label found no
        antimode, i.e. no call).

    Returns
    -------
    dict
        Rule -> boolean call per spot.
    """
    calls = {"> 0": s > 0}
    if method == "EnrichMap (null)":
        calls["> 1.96"] = s > 1.96
        calls["p_emp <= 0.05"] = p_emp <= 0.05
    if method == "GSDensity" and label is not None:
        calls["antimode label"] = np.asarray(label, dtype=bool)
    return calls


def f1_of_calls(y, pred):
    """
    F1 score of binary calls.

    Parameters
    ----------
    y : np.ndarray
        Binary ground truth.
    pred : np.ndarray
        Binary calls.

    Returns
    -------
    float
        2 TP / (called + positives), 0 when both are empty.
    """
    y = np.asarray(y).astype(bool)
    pred = np.asarray(pred).astype(bool)
    tp = (y & pred).sum()
    denom = pred.sum() + y.sum()
    return float(2 * tp / denom) if denom else 0.0


def sign_metrics(y, s, calls, rule_columns):
    """
    How interpretable is the sign (or call) of a score?

    Parameters
    ----------
    y : np.ndarray or None
        Binary ground truth (1 = target region), or None for a null or absent signature.
    s : np.ndarray
        Score (ranking metrics).
    calls : dict
        Output of `sign_calls`.
    rule_columns : dict
        Rule -> (fraction, target, background, F1) column names.

    Returns
    -------
    dict
        Fraction called per rule, plus AUROC/AUPRC and target/background call rates and F1 when `y` is given.
    """
    out = {}
    for rule, c in calls.items():
        out[rule_columns[rule][0]] = float(np.mean(c))
    if y is not None:
        y = np.asarray(y).astype(bool)
        out["AUROC"] = roc_auc_score(y, s)
        out["AUPRC"] = average_precision_score(y, s)
        for rule, c in calls.items():
            _, tgt, bg, f1 = rule_columns[rule]
            out[tgt] = float(np.mean(c[y]))
            out[bg] = float(np.mean(c[~y]))
            out[f1] = f1_of_calls(y, c)
    return out


def calls_long(tab, keys, call_rules, rule_columns, rule_label):
    """
    One row per (keys, method, call rule) from a wide sign-metric table.

    Parameters
    ----------
    tab : pd.DataFrame
        Wide table with a "score" column holding the method label.
    keys : list
        Identifier columns kept (e.g. signature, domain, seed).
    call_rules : list
        (method, rule) pairs, in plotting order.
    rule_columns : dict
        Rule -> (fraction, target, background, F1) column names.
    rule_label : dict
        Rule -> label.

    Returns
    -------
    pd.DataFrame
        Columns keys + method, rule, call, fraction called, target called, background called, F1, AUROC.
    """
    rows = []
    for method, rule in call_rules:
        frac, tgt, bg, f1 = rule_columns[rule]
        d = tab[tab["score"] == method]
        for _, r in d.iterrows():
            rows.append(
                {
                    **{k: r[k] for k in keys},
                    "method": method,
                    "rule": rule,
                    "call": f"{method}: {rule_label[rule]}",
                    "fraction called": r.get(frac, np.nan),
                    "target called": r.get(tgt, np.nan),
                    "background called": r.get(bg, np.nan),
                    "F1": r.get(f1, np.nan),
                    "AUROC": r.get("AUROC", np.nan),
                }
            )
    return pd.DataFrame(rows)


def read_gsdensity_relevance(path, spots):
    """
    GSDensity per-spot output (spot_id, relevance, label) in spot order.

    Parameters
    ----------
    path : str
        CSV of revision/data/gsd_cache.
    spots : pd.Index
        Spot ids of the slide.

    Returns
    -------
    tuple
        (relevance x number of spots, antimode label as boolean or None where compute.cell.label found no antimode).
    """
    tab = pd.read_csv(path).set_index("spot_id").loc[list(spots)]
    label = None if tab["label"].isna().all() else (tab["label"] == "positive").values
    return tab["relevance"].values.astype(np.float64) * len(tab), label


def simulated_sign_metrics(seed, domain, variants, gsd_sim_dir, rule_columns):
    """
    Sign metrics of the five scores for one simulated slide (target and null signature).

    Parameters
    ----------
    seed : int
        Simulation seed.
    domain : str
        "voronoi" or "stripe".
    variants : dict
        EnrichMap variant ("null", "smooth", "gam") -> em.tl.score arguments.
    gsd_sim_dir : str
        Directory of the GSDensity results of the simulated slides.
    rule_columns : dict
        Rule -> (fraction, target, background, F1) column names.

    Returns
    -------
    list of dict
        One row per signature and method.
    """
    labels = {
        "null": "EnrichMap (null)",
        "smooth": "EnrichMap (smoothing-only)",
        "gam": "EnrichMap (GAM)",
    }
    files = {"target": "target", "null signature": "null_signature"}
    sim, sim_sigs, target_label = simulated_slide_with_null_signature(seed, domain)
    y = (sim.obs["cluster"] == target_label).values.astype(int)
    rows = []
    for signature, genes in sim_sigs.items():
        for variant, kwargs in variants.items():
            em.tl.score(
                sim, gene_set=genes, score_key=f"{signature}_{variant}", **kwargs
            )
        scores = {
            labels[v]: sim.obs[f"{signature}_{v}_score"].values.astype(np.float64)
            for v in variants
        }
        smooth = scores["EnrichMap (smoothing-only)"]
        scores["EnrichMap (slide z-score)"] = (smooth - smooth.mean()) / smooth.std()
        p_emp = sim.uns["null_model"][f"{signature}_null"]["p_empirical"].astype(
            np.float64
        )
        relevance, label = read_gsdensity_relevance(
            gsd_sim_dir + f"{domain}_seed{seed:03d}_{files[signature]}.csv",
            sim.obs_names,
        )
        scores["GSDensity"] = relevance
        truth = y if signature == "target" else None
        for method, s in scores.items():
            calls = sign_calls(method, s, p_emp=p_emp, label=label)
            row = {
                "domain": domain,
                "seed": seed,
                "signature": signature,
                "score": method,
                "n spots": sim.n_obs,
                "prevalence": y.mean(),
            }
            if method == "GSDensity":
                row["GSDensity label available"] = label is not None
            if method == "EnrichMap (null)":
                row["sign agreement with slide z-score"] = float(
                    np.mean((s > 0) == (scores["EnrichMap (slide z-score)"] > 0))
                )
            row.update(sign_metrics(truth, s, calls, rule_columns))
            rows.append(row)
    return rows


def mean_with_range(d, value):
    """
    Mean with min-max range over signatures.

    Parameters
    ----------
    d : pd.DataFrame
        Table holding `value`.
    value : str
        Column.

    Returns
    -------
    pd.Series
        value, value_lo, value_hi.
    """
    return pd.Series(
        {
            value: d[value].mean(),
            f"{value}_lo": d[value].min(),
            f"{value}_hi": d[value].max(),
        }
    )


def mean_with_t_interval(d, value):
    """
    Mean with a 95% t-interval over slides.

    Parameters
    ----------
    d : pd.DataFrame
        Table holding `value`.
    value : str
        Column.

    Returns
    -------
    pd.Series
        value, value_lo, value_hi.
    """
    v = d[value].dropna().values
    half = (
        student_t.ppf(0.975, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v))
        if len(v) > 1
        else 0.0
    )
    return pd.Series(
        {
            value: v.mean(),
            f"{value}_lo": v.mean() - half,
            f"{value}_hi": v.mean() + half,
        }
    )


def call_rule_pair_table(d, x, y, agg, call_rules, rule_markers):
    """
    One row per (method, rule) in plotting order, with mean and interval columns for `x` and `y`.

    Parameters
    ----------
    d : pd.DataFrame
        Long call table (calls_long).
    x, y : str
        Value columns.
    agg : callable
        mean_with_range or mean_with_t_interval.
    call_rules : list
        (method, rule) pairs.
    rule_markers : dict
        Rule -> marker.

    Returns
    -------
    pd.DataFrame
        Columns method, rule, marker and x, y with _lo/_hi.
    """
    rows = []
    for method, rule in call_rules:
        g = d[(d.method == method) & (d.rule == rule)]
        rows.append(
            {
                "method": method,
                "rule": rule,
                "marker": rule_markers[rule],
                **agg(g, x).to_dict(),
                **agg(g, y).to_dict(),
            }
        )
    return pd.DataFrame(rows)


def primary_call_mean(d, value, primary_rule):
    """
    Mean of `value` per method for its own call rule.

    Parameters
    ----------
    d : pd.DataFrame
        Long call table (calls_long).
    value : str
        Column.
    primary_rule : dict
        Method -> its call rule.

    Returns
    -------
    pd.Series
        Method -> mean.
    """
    return pd.Series(
        {
            m: d[(d.method == m) & (d.rule == r)][value].mean()
            for m, r in primary_rule.items()
        }
    )


def mean_ranking_metric(d, metric, methods):
    """
    Mean of a ranking metric (AUROC, AUPRC) per method from a wide sign-metric table.

    Parameters
    ----------
    d : pd.DataFrame
        Wide table with a "score" column holding the method label.
    metric : str
        Column.
    methods : list
        Method order.

    Returns
    -------
    pd.Series
        Method -> mean.
    """
    return d.groupby("score")[metric].mean().reindex(methods)


def mean_sd_over_slides(d, value):
    """
    Mean and SD over slides per (domain, method, rule).

    Parameters
    ----------
    d : pd.DataFrame
        Table holding `value` and some of the columns domain, method, rule.
    value : str
        Column.

    Returns
    -------
    pd.DataFrame
        Grouping columns, mean and std.
    """
    keys = [c for c in ["domain", "method", "rule"] if c in d.columns]
    return d.groupby(keys)[value].agg(["mean", "std"]).reset_index()


def plot_sign_metric_pair(table, x, y, ax, colours, ref_x=None, marker_col=None):
    """
    One metric against another, one point per method (and rule or task), with whiskers when columns
    f"{x}_lo"/f"{x}_hi" and f"{y}_lo"/f"{y}_hi" exist.

    Parameters
    ----------
    table : pd.DataFrame
        One row per point, with a "method" column and the metric columns.
    x, y : str
        Metric columns.
    ax : plt.Axes
        Target axes.
    colours : dict
        Method -> colour.
    ref_x : float or None
        Dotted vertical reference line (e.g. AUROC 0.5).
    marker_col : str or None
        Column holding a marker per row; default circles.

    Returns
    -------
    plt.Axes
    """
    for _, r in table.iterrows():
        xerr = yerr = None
        if f"{x}_lo" in table:
            xerr = [[r[x] - r[f"{x}_lo"]], [r[f"{x}_hi"] - r[x]]]
        if f"{y}_lo" in table:
            yerr = [[r[y] - r[f"{y}_lo"]], [r[f"{y}_hi"] - r[y]]]
        ax.errorbar(
            r[x],
            r[y],
            xerr=xerr,
            yerr=yerr,
            fmt=r[marker_col] if marker_col else "o",
            color=colours[r["method"]],
            ms=5,
            mec="white",
            mew=0.4,
            elinewidth=0.6,
            capsize=0,
            label=r["method"],
        )
    if ref_x is not None:
        ax.axvline(ref_x, ls=":", c="0.6", lw=0.6)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.grid(False)
    return ax


def plot_sign_summary(table, ax, colours, lower_is_better=()):
    """
    Summary panel (methods x metrics): circles sized and shaded by each metric scaled to 0-1 across methods
    (best = 1), raw values annotated, plus a mean-rank column.

    Parameters
    ----------
    table : pd.DataFrame
        Index = methods, columns = metrics.
    ax : plt.Axes
        Target axes.
    colours : dict
        Method -> colour (tick labels).
    lower_is_better : tuple of str
        Metrics where smaller is better.

    Returns
    -------
    plt.Axes
    """
    scaled = table.copy().astype(float)
    for c in table.columns:
        col = table[c].astype(float)
        rng_ = col.max() - col.min()
        s = (col - col.min()) / rng_ if rng_ > 0 else col * 0 + 1
        scaled[c] = 1 - s if c in lower_is_better else s
    ranks = scaled.rank(ascending=False).mean(axis=1)
    cmap = plt.get_cmap("Purples")
    for i, m in enumerate(table.index):
        for j, c in enumerate(table.columns):
            s = scaled.loc[m, c]
            ax.scatter(
                j, i, s=30 + 220 * s, color=cmap(0.25 + 0.75 * s), edgecolors="none"
            )
            ax.text(
                j,
                i,
                f"{table.loc[m, c]:.2g}",
                ha="center",
                va="center",
                fontsize=5,
                color="white" if s > 0.6 else "black",
            )
        ax.text(
            len(table.columns) + 0.2,
            i,
            f"{ranks[m]:.1f}",
            ha="left",
            va="center",
            fontsize=6,
        )
    ax.set_xticks(range(len(table.columns)))
    ax.set_xticklabels(table.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(table.index)))
    ax.set_yticklabels(table.index)
    for lab in ax.get_yticklabels():
        lab.set_color(colours[lab.get_text()])
    ax.text(len(table.columns) + 0.2, -0.9, "mean rank", fontsize=6, ha="left")
    ax.set_xlim(-0.6, len(table.columns) + 0.9)
    ax.set_ylim(len(table.index) - 0.4, -1.2)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    return ax


def grouped_call_bars(
    ax, d, group_col, groups, value, items, colours, rule_hatch, err=None, labels=None
):
    """
    Grouped bars: one bar per (method, rule) in each group; colour = method, hatch = rule.

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    d : pd.DataFrame
        Table with method, rule (optional), the group column and the value (and error) columns.
    group_col : str
        Column defining the x groups.
    groups : list
        Group order.
    value : str
        Bar height column.
    items : list
        (method, rule) pairs; (method, None) for tables without rules.
    colours : dict
        Method -> colour.
    rule_hatch : dict
        Rule -> hatch.
    err : str or None
        Error-bar column.
    labels : list or None
        Tick labels (default: the group names).

    Returns
    -------
    list
        (x position, method) of every bar whose value is 0 (invisible on log axes).
    """
    width = 0.84 / len(items)
    zeros = []
    for k, (method, rule) in enumerate(items):
        g = (
            d[d.method == method]
            if rule is None
            else d[(d.method == method) & (d.rule == rule)]
        )
        g = g.set_index(group_col).reindex(groups)
        x = np.arange(len(groups)) - 0.42 + width * (k + 0.5)
        ax.bar(
            x,
            g[value].values,
            width,
            color=colours[method],
            hatch=rule_hatch[rule] if rule else "",
            edgecolor="white",
            linewidth=0.3,
            yerr=g[err].values if err else None,
            error_kw=dict(lw=0.5, capsize=0, ecolor="0.25"),
        )
        zeros += [(xi, method) for xi, v in zip(x, g[value].values) if v == 0]
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(labels or groups, rotation=30, ha="right", fontsize=6)
    return zeros


# --------------------------------------------------------------------------------------------------
# R4_internal_audit
# --------------------------------------------------------------------------------------------------


def simulate_visium_fig1(
    seed=0,
    n_genes=3000,
    n_signature_genes=20,
    n_clusters=4,
    target_cluster=1,
    n_spots_range=(1500, 4000),
    spot_spacing_um=100,
    spot_diameter_um=55,
    mean_counts_per_spot=5000,
    dispersion_range=(0.1, 2.0),
    fold_change_range=(1.5, 5.0),
    dropout_rate=0.1,
    n_patches_per_cluster=100,
    smoothing_iterations=1,
):
    """
    Simulate a 10x Visium spatial transcriptomics dataset.

    The Fig. 1 generator of 05_00_simulated_dataset (`simulate_visium` there), copied verbatim and renamed
    so that it does not clash with the Sup. Fig. 1 generator of 05_01_simulated_benchmarking.

    The grid layout matches the Visium platform: a square capture area
    with hexagonally arranged spots at 100 um centre-to-centre spacing.
    The number of spots is drawn randomly from ``n_spots_range`` to
    mimic variable tissue coverage across real slides. This is achieved
    by varying the capture area dimensions per simulation.

    Clusters are formed by seeding ``n_patches_per_cluster`` random foci
    per cluster and then iteratively smoothing assignments by majority
    vote within a local neighbourhood, producing irregular, tissue-like
    patches instead of clean Voronoi regions.

    Returns (adata, sig_genes, fold_changes, target_label).
    """
    rng = np.random.default_rng(seed)

    # Variable spot count via variable capture area
    # Draw a target spot count, then compute the grid side that
    # produces approximately that many spots in a hexagonal grid.
    # For a hex grid: n_spots ≈ n_rows × n_cols ≈ side^2 / 0.866
    n_spots_target = rng.integers(
        n_spots_range[0],
        n_spots_range[1] + 1,
    )
    grid_side = int(np.ceil(np.sqrt(n_spots_target * 0.866)))
    capture_area_um = grid_side * spot_spacing_um

    # Hexagonal grid within square capture area
    n_cols = int(capture_area_um / spot_spacing_um)
    n_rows = int(capture_area_um / (spot_spacing_um * np.sqrt(3) / 2))

    x_coords, y_coords = [], []
    for row in range(n_rows):
        for col in range(n_cols):
            x = col * spot_spacing_um + (row % 2) * (spot_spacing_um / 2)
            y = row * spot_spacing_um * np.sqrt(3) / 2
            if x <= capture_area_um and y <= capture_area_um:
                x_coords.append(x)
                y_coords.append(y)

    x_coords = np.array(x_coords)
    y_coords = np.array(y_coords)
    n_spots = len(x_coords)
    coords = np.column_stack([x_coords, y_coords])

    # Patchy clusters: multiple seed foci per cluster + spatial smoothing
    # Each cluster gets n_patches_per_cluster random seed points so the
    # same label can appear as several disconnected blobs across the tissue,
    # mimicking infiltrating immune zones or tumour patches.
    margin = capture_area_um * 0.15
    total_seeds = n_clusters * n_patches_per_cluster
    seed_cluster_ids = np.repeat(np.arange(n_clusters), n_patches_per_cluster)

    seed_coords = np.column_stack(
        [
            rng.uniform(margin, capture_area_um - margin, total_seeds),
            rng.uniform(margin, capture_area_um - margin, total_seeds),
        ]
    )

    dists = cdist(coords, seed_coords)
    nearest_seed = np.argmin(dists, axis=1)
    cluster_ids = seed_cluster_ids[nearest_seed]

    # Iterative majority-vote smoothing: each spot adopts the most common
    # label among its spatial neighbours, creating organic, rounded boundaries.
    tree = KDTree(coords)
    smoothing_radius = spot_spacing_um * 2.5
    neighbors = tree.query_ball_point(coords, smoothing_radius)
    for _ in range(smoothing_iterations):
        cluster_ids = np.array(
            [
                np.bincount(cluster_ids[nbrs], minlength=n_clusters).argmax()
                for nbrs in neighbors
            ]
        )

    cluster_labels = np.array([f"cluster_{i + 1}" for i in cluster_ids])
    target_label = f"cluster_{target_cluster + 1}"

    # Gene names
    gene_names = [f"GENE_{i + 1:05d}" for i in range(n_genes)]

    # Signature genes with graded fold changes
    sig_indices = rng.choice(n_genes, n_signature_genes, replace=False)
    fold_changes = rng.uniform(
        fold_change_range[0],
        fold_change_range[1],
        n_signature_genes,
    )
    order = np.argsort(-fold_changes)
    sig_indices = sig_indices[order]
    fold_changes = fold_changes[order]
    sig_genes = [gene_names[i] for i in sig_indices]

    # Expression matrix (negative binomial)
    gene_means = np.exp(rng.normal(1.5, 1.0, n_genes))
    gene_means = np.clip(gene_means, 0.1, 100)

    dispersions = rng.uniform(
        dispersion_range[0],
        dispersion_range[1],
        n_genes,
    )

    lib_size_factors = np.exp(rng.normal(0, 0.3, n_spots))
    lib_size_factors *= mean_counts_per_spot / (
        gene_means.sum() * lib_size_factors.mean()
    )

    X = np.zeros((n_spots, n_genes), dtype=np.float32)
    target_mask = cluster_labels == target_label

    for g in range(n_genes):
        mu = gene_means[g] * lib_size_factors

        if g in sig_indices:
            fc_idx = np.where(sig_indices == g)[0][0]
            mu[target_mask] *= fold_changes[fc_idx]

        r = dispersions[g]
        p = r / (r + mu)
        counts = rng.negative_binomial(r, p)

        dropout_mask = rng.random(n_spots) < dropout_rate
        counts[dropout_mask] = 0

        X[:, g] = counts

    # Assemble AnnData
    adata = anndata.AnnData(X=csr_matrix(X))
    adata.var_names = gene_names
    adata.obs_names = [f"spot_{i}" for i in range(n_spots)]
    adata.obs["cluster"] = pd.Categorical(cluster_labels)
    adata.obsm["spatial"] = coords

    adata.raw = adata.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    return adata, sig_genes, fold_changes, target_label


def _make_enrichmap_scorer():
    def scorer(ad, gene_list):
        em.tl.score(
            ad,
            gene_set=gene_list,
            score_key="_bench",
            smoothing=True,
            null_model=False,
            correct_spatial_covariates=False,
        )
        return ad.obs["_bench_score"].values

    return scorer


def _make_enrichmap_gam_scorer():
    def scorer(ad, gene_list):
        em.tl.score(
            ad,
            gene_set=gene_list,
            score_key="_bench",
            smoothing=True,
            null_model=False,
            correct_spatial_covariates=True,
        )
        return ad.obs["_bench_score"].values

    return scorer


def _make_enrichmap_null_scorer():
    def scorer(ad, gene_list):
        em.tl.score(
            ad,
            gene_set=gene_list,
            score_key="_bench",
            smoothing=True,
            null_model=True,
            correct_spatial_covariates=False,
        )
        return ad.obs["_bench_score"].values

    return scorer


def _make_scanpy_scorer():
    def scorer(ad, gene_list):
        sc.tl.score_genes(ad, gene_list=gene_list, score_name="_bench", use_raw=False)
        return ad.obs["_bench"].values

    return scorer


def _make_aucell_scorer():
    def scorer(ad, gene_list):
        net = pd.DataFrame({"source": "_bench", "target": gene_list})
        dc.mt.aucell(ad, net=net)
        return (
            dc.pp.get_obsm(adata=ad, key="score_aucell")
            .obsm["score_aucell"]["_bench"]
            .values
        )

    return scorer


def _make_zscore_scorer():
    def scorer(ad, gene_list):
        X = ad.X.toarray() if issparse(ad.X) else ad.X
        expr_df = pd.DataFrame(X.T, index=ad.var_names, columns=ad.obs_names)
        common = [g for g in gene_list if g in expr_df.index]
        return (expr_df.loc[common].sum(axis=0) / np.sqrt(len(common))).values

    return scorer


def _make_ssgsea_scorer():
    def scorer(ad, gene_list):
        X = ad.X.toarray() if issparse(ad.X) else ad.X
        expr_df = pd.DataFrame(X.T, index=ad.var_names, columns=ad.obs_names)
        expr_df.index = expr_df.index.str.upper()
        gl_upper = [g.upper() for g in gene_list]
        result = gp.ssgsea(
            data=expr_df,
            gene_sets={"_bench": gl_upper},
            outdir=None,
            sample_norm_method="rank",
            no_plot=True,
            min_size=3,
        )
        nes = result.res2d.pivot(index="Term", columns="Name", values="NES")
        return nes.loc["_bench"].T.astype(float).reindex(ad.obs_names).values

    return scorer


def moran_geary_libpysal(score, coords, k=None):
    """
    Moran's I and Geary's C as em.tl.compute_spatial_metrics computes them.

    libpysal KNN weights, row-standardised, esda Moran / Geary; `k=None` is the package call
    `KNN.from_array(coords)` (libpysal default, k = 2).

    Parameters
    ----------
    score : np.ndarray
        Score per spot.
    coords : np.ndarray
        Spot coordinates.
    k : int or None
        Number of neighbours (None = libpysal default).

    Returns
    -------
    dict
        Moran's I, Geary's C and the number of neighbours used.
    """
    W = (
        KNN.from_array(np.asarray(coords))
        if k is None
        else KNN.from_array(np.asarray(coords), k=k)
    )
    W.transform = "r"
    s = np.asarray(score, dtype=np.float64)
    return {
        "Moran's I": float(Moran(s, W, permutations=0).I),
        "Geary's C": float(Geary(s, W, permutations=0).C),
        "k": int(W.max_neighbors),
    }


def coherence_f1_threshold(scores, labels, W):
    """
    'Spatial coherence' as run_classification (05_01) and _benchmark_pooled (17), clipped and unclipped.

    F1-optimal threshold -> quantile q (clipped to [0.01, 0.99]) -> predictions = scores > quantile(scores, q)
    -> Moran's I(prediction) / Moran's I(truth); the published value clips the ratio to [0, 1].

    Parameters
    ----------
    scores : np.ndarray
        Method scores.
    labels : np.ndarray
        Binary ground truth.
    W : scipy.sparse matrix
        Binary 6-NN spatial weights.

    Returns
    -------
    dict
        Clipped and unclipped ratio, Moran's I of prediction and truth, predicted prevalence.
    """
    y = np.asarray(labels).astype(int)
    s = np.asarray(scores, dtype=np.float64)
    gt = _morans_i_continuous(y.astype(float), W)
    prec, rec, thr = precision_recall_curve(y, s)
    prec, rec = prec[:-1], rec[:-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        f1 = np.nan_to_num(2 * prec * rec / (prec + rec), nan=0.0)
    q = float(np.clip(np.mean(s <= float(thr[int(np.argmax(f1))])), 0.01, 0.99))
    preds = (s > np.quantile(s, q)).astype(float)
    mi = _morans_i_continuous(preds, W) if preds.std() > 0 else 0.0
    ratio = mi / gt if gt > 0 else np.nan
    return {
        "coherence (clipped, as published)": float(np.clip(ratio, 0, 1))
        if preds.std() > 0
        else 0.0,
        "coherence ratio (unclipped)": float(ratio) if preds.std() > 0 else 0.0,
        "Moran's I prediction": float(mi),
        "Moran's I truth": float(gt),
        "predicted prevalence": float(preds.mean()),
    }


def coherence_of_calls(calls, labels, W):
    """
    Spatial coherence of a given binary call (e.g. GSDensity's antimode labels), without a threshold search.

    Parameters
    ----------
    calls : np.ndarray
        Binary predictions.
    labels : np.ndarray
        Binary ground truth.
    W : scipy.sparse matrix
        Binary 6-NN spatial weights.

    Returns
    -------
    dict
        Same keys as coherence_f1_threshold, plus the F1 of the calls.
    """
    pred = np.asarray(calls, dtype=float)
    gt = _morans_i_continuous(np.asarray(labels, dtype=float), W)
    mi = _morans_i_continuous(pred, W) if pred.std() > 0 else 0.0
    ratio = mi / gt if pred.std() > 0 else 0.0
    return {
        "coherence (clipped, as published)": float(np.clip(ratio, 0, 1)),
        "coherence ratio (unclipped)": float(ratio),
        "Moran's I prediction": float(mi),
        "Moran's I truth": float(gt),
        "predicted prevalence": float(pred.mean()),
        "F1 of the calls": float(
            f1_score(np.asarray(labels).astype(int), pred.astype(int), zero_division=0)
        ),
    }


def coherence_at_prevalence(scores, prevalence, W, gt_morans):
    """
    Unclipped coherence ratio when the top `prevalence` fraction of spots is called positive.

    Parameters
    ----------
    scores : np.ndarray
        Score (or signal-free field) per spot.
    prevalence : float
        Fraction of spots called positive (the target prevalence).
    W : scipy.sparse matrix
        Binary 6-NN spatial weights.
    gt_morans : float
        Moran's I of the ground truth.

    Returns
    -------
    float
        Moran's I(prediction) / Moran's I(truth).
    """
    s = np.asarray(scores, dtype=np.float64)
    preds = (s > np.quantile(s, 1 - prevalence)).astype(float)
    return _morans_i_continuous(preds, W) / gt_morans


def assign_weights(pos_genes, zero_genes, rng):
    """Random uniform weights in (0, 1] for pos_genes; weight=0 for zero_genes."""
    weights = {}
    pos_w = rng.uniform(0, 1, size=len(pos_genes))
    for gene, w in zip(pos_genes, pos_w):
        weights[gene] = round(float(w), 4)
    for gene in zero_genes:
        weights[gene] = 0
    return weights


def cohens_d(group1, group2):
    """Pooled-std Cohen's d (group1 is the positive class)."""
    n1, n2 = len(group1), len(group2)
    s_pooled = np.sqrt(
        ((n1 - 1) * np.var(group1, ddof=1) + (n2 - 1) * np.var(group2, ddof=1))
        / (n1 + n2 - 2)
    )
    return (np.mean(group1) - np.mean(group2)) / s_pooled if s_pooled > 0 else 0.0


def bootstrap_spearman(x, y, n_boot=1000, seed=0):
    rng_b = np.random.default_rng(seed)
    n = len(x)
    corrs = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng_b.integers(0, n, size=n)
        corrs[b] = spearmanr(x[idx], y[idx]).statistic
    return corrs


def bootstrap_spearman_ci_genes(x, y, n_boot=2000, seed=0):
    """
    Percentile bootstrap 95% CI of Spearman's rho over genes.

    Parameters
    ----------
    x, y : np.ndarray
        Paired gene-level values.
    n_boot : int
        Bootstrap resamples.
    seed : int
        Random seed.

    Returns
    -------
    tuple of float
        Lower and upper 95% limits.
    """
    rng = np.random.default_rng(seed)
    n = len(x)
    vals = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        if np.unique(x[i]).size > 1 and np.unique(y[i]).size > 1:
            vals.append(spearmanr(x[i], y[i]).correlation)
    return tuple(np.nanquantile(vals, [0.025, 0.975]))


def log2fc_target_vs_rest(adata, genes, mask_target):
    """
    log2 fold change of log1p means, target vs rest, as compute_gene_stats in 18_visium_weights.

    Parameters
    ----------
    adata : AnnData
        Log-normalised expression.
    genes : list
        Gene names.
    mask_target : np.ndarray
        Boolean mask of the target spots.

    Returns
    -------
    pd.Series
        log2((mean_target + 1e-6) / (mean_rest + 1e-6)) per gene.
    """
    X = adata[:, genes].X
    X = X.toarray() if issparse(X) else np.asarray(X)
    mt = X[mask_target].mean(axis=0)
    mr = X[~mask_target].mean(axis=0)
    return pd.Series(np.log2((mt + 1e-6) / (mr + 1e-6)), index=genes)


def partial_rank_correlation(x, y, z):
    """
    Rank-based partial correlation of x and y given z.

    Ranks of x and y are regressed linearly on the rank of z; Pearson correlation of the residuals.

    Parameters
    ----------
    x, y, z : np.ndarray
        Variables of the same length.

    Returns
    -------
    float
        Partial Spearman correlation.
    """
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    Z = np.column_stack([np.ones_like(rz), rz])
    ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
    ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])


def partial_rank_correlation_by_group(x, y, groups):
    """
    Rank correlation of x and y after removing the group (region) means of the ranks.

    Parameters
    ----------
    x, y : np.ndarray
        Variables of the same length.
    groups : np.ndarray
        Group label per observation.

    Returns
    -------
    float
        Region-adjusted Spearman correlation.
    """
    rx, ry = rankdata(x), rankdata(y)
    g = np.asarray(groups)
    ex = rx - pd.Series(rx).groupby(g).transform("mean").values
    ey = ry - pd.Series(ry).groupby(g).transform("mean").values
    return float(np.corrcoef(ex, ey)[0, 1])


def audit_simulation(
    seed, variants, variant_labels, scoring_functions, gsd_csv, n_null_fields=10
):
    """
    Items 1 and 4 of the internal audit on one Sup. Fig. 1 simulation (simulate_visium, 05_01).

    Item 1: the Sup. Fig. 1d metric (mean |gene_contributions| of the smoothing-only score, as
    run_gene_weight_validation in 05_01) against the CV weight, gene SD and baseline expression. Item 4:
    spatial coherence as run_classification (05_01), clipped and unclipped, for every method, GSDensity's
    antimode labels, the null model's fixed rules (z > 0, z > 1.96, p_empirical <= 0.05), and signal-free
    fields (iid noise, smoothed once or twice like EnrichMap) with the top-p spots called positive.

    Parameters
    ----------
    seed : int
        Simulation seed.
    variants : dict
        EnrichMap variants, key -> em.tl.score arguments; must contain "null" and "smooth".
    variant_labels : dict
        Variant key -> method label.
    scoring_functions : dict
        Other methods, label -> scorer(adata, genes) (17_benchmarking factories).
    gsd_csv : str
        GSDensity relevance and labels of this simulation (columns "relevance", "label").
    n_null_fields : int
        Signal-free draws per simulation.

    Returns
    -------
    dict
        "gw" (per-seed item 1 statistics), "genes" (per-gene table), "coherence" and "null" (item 4 rows),
        "meta" (seed, spots, prevalence).
    """
    adata, sig, fc, tgt = simulate_visium(seed=seed)
    y = (adata.obs["cluster"] == tgt).values.astype(int)
    for name, kwargs in variants.items():
        em.tl.score(adata, gene_set=sig, score_key=f"target_{name}", **kwargs)
    scores = {
        variant_labels[name]: adata.obs[f"target_{name}_score"].values.astype(
            np.float64
        )
        for name in variants
    }
    for method, fn in scoring_functions.items():
        scores[method] = np.asarray(fn(adata.copy(), sig), dtype=np.float64).ravel()
    gsd = pd.read_csv(gsd_csv, index_col=0).reindex(adata.obs_names)
    scores["GSDensity"] = gsd["relevance"].values.astype(np.float64)

    # item 1: Sup. Fig. 1d metric (smoothing-only contributions = expression x CV weight)
    contrib = adata.uns["gene_contributions"]["target_smooth"]
    common = [g for g in sig if g in contrib]
    metric = np.array([np.mean(np.abs(contrib[g])) for g in common], dtype=np.float64)
    fc_m = np.array([fc[sig.index(g)] for g in common])
    w = em.tl.infer_gene_weights(adata, sig)
    cv = np.array([w[g] for g in common], dtype=np.float64)
    Xs = adata[:, common].X.toarray()
    sd = Xs.std(axis=0).astype(np.float64)
    base = Xs[y == 0].mean(axis=0).astype(np.float64)
    genes = pd.DataFrame(
        {
            "seed": seed,
            "gene": common,
            "fold_change": fc_m,
            "log2_fold_change": np.log2(fc_m),
            "cv_weight": cv,
            "sd": sd,
            "mean": Xs.mean(axis=0).astype(np.float64),
            "baseline_mean": base,
            "original_metric": metric,
        }
    )
    m_all = np.asarray(adata.X.mean(axis=0)).ravel().astype(np.float64)
    sd_all = np.sqrt(
        np.maximum(
            np.asarray(adata.X.multiply(adata.X).mean(axis=0)).ravel() - m_all**2, 0
        )
    )
    cv_all = np.where(m_all > 0, sd_all / np.where(m_all > 0, m_all, 1), 0.0)
    is_sig = adata.var_names.isin(sig).astype(int)
    gw = {
        "seed": seed,
        "rho(FC, original metric)": spearmanr(metric, fc_m).correlation,
        "rho(FC, SD)": spearmanr(sd, fc_m).correlation,
        "rho(FC, CV weight)": spearmanr(cv, fc_m).correlation,
        "rho(CV weight, baseline mean)": spearmanr(cv, base).correlation,
        "rho(FC, baseline mean)": spearmanr(fc_m, base).correlation,
        "partial rho(FC, CV weight | baseline mean)": partial_rank_correlation(
            cv, fc_m, base
        ),
        "max |metric - SD|": float(np.max(np.abs(metric - sd))),
        "max relative |metric - SD|": float(np.max(np.abs(metric - sd) / sd)),
        "identical ranks (metric vs SD)": bool(
            np.array_equal(
                np.argsort(metric, kind="stable"), np.argsort(sd, kind="stable")
            )
        ),
        "AUROC CV (signature vs other genes)": roc_auc_score(is_sig, cv_all),
        "AUROC SD (signature vs other genes)": roc_auc_score(is_sig, sd_all),
        "AUROC mean (signature vs other genes)": roc_auc_score(is_sig, m_all),
    }

    # item 4: spatial coherence, as published and unclipped
    conn = spatial_knn_connectivity(np.asarray(adata.obsm["spatial"]))
    coherence = []
    for method, s in scores.items():
        r = coherence_f1_threshold(s, y, conn)
        r.update({"seed": seed, "method": method, "AUROC": roc_auc_score(y, s)})
        coherence.append(r)
    r = coherence_of_calls((gsd["label"] == "positive").values, y, conn)
    r.update({"seed": seed, "method": "GSDensity (antimode labels)", "AUROC": np.nan})
    coherence.append(r)
    z = scores[variant_labels["null"]]
    p_emp = adata.uns["null_model"]["target_null"]["p_empirical"]
    label = variant_labels["null"]
    for rule, calls in (
        (f"{label}, z > 0", z > 0),
        (f"{label}, z > 1.96", z > 1.96),
        (f"{label}, p_empirical ≤ 0.05", p_emp <= 0.05),
    ):
        r = coherence_of_calls(calls, y, conn)
        r.update({"seed": seed, "method": rule, "AUROC": np.nan})
        coherence.append(r)

    # signal-free fields at the target prevalence, and the same rule applied to the scores (draw = -1)
    rng = np.random.default_rng(50_000 + seed)
    prev = float(y.mean())
    gt = _morans_i_continuous(y.astype(float), conn)
    null = []
    for i in range(n_null_fields):
        z0 = rng.normal(size=len(y))
        zs = neighbour_average(z0, conn)
        for field, v in (
            ("noise (raw)", z0),
            ("noise, 6-NN smoothed", zs),
            ("noise, smoothed twice", neighbour_average(zs, conn)),
        ):
            null.append(
                {
                    "seed": seed,
                    "draw": i,
                    "field": field,
                    "coherence ratio (unclipped)": coherence_at_prevalence(
                        v, prev, conn, gt
                    ),
                }
            )
    for method, s in scores.items():
        null.append(
            {
                "seed": seed,
                "draw": -1,
                "field": method,
                "coherence ratio (unclipped)": coherence_at_prevalence(
                    s, prev, conn, gt
                ),
            }
        )
    meta = {"seed": seed, "n_spots": len(y), "prevalence": prev}
    return {
        "gw": gw,
        "genes": genes,
        "coherence": coherence,
        "null": null,
        "meta": meta,
    }


def audit_noise_robustness(
    adata,
    region,
    genes,
    scoring_functions,
    gsd_clean,
    gsd_task,
    noise_levels=(0.0, 0.5, 1.0, 1.5, 2.0, 3.0),
    n_repeats=3,
    n_stream_repeats=5,
    seed=0,
):
    """
    Noise arm of item 5: self-agreement (Spearman with the clean score) and accuracy (AUROC vs the region).

    Gaussian noise N(0, level x gene SD) is added to the signature genes with the random stream of
    benchmark_noise_robustness in 17_benchmarking (default_rng(seed); n_stream_repeats draws per level), and
    the first `n_repeats` repeats are scored, so repeat r uses exactly the noise of 17's repeat r. GSDensity
    relevance of the same noisy matrices is read from the REV_Fig2 runs.

    Parameters
    ----------
    adata : AnnData
        Mouse brain, log-normalised, with "cluster" in obs and coordinates in obsm["spatial"].
    region : str
        Annotated region (positive class) and the column of GSDensity's score files.
    genes : list
        Signature genes.
    scoring_functions : dict
        Method label -> scorer(adata, genes).
    gsd_clean : np.ndarray
        GSDensity relevance on the clean matrix.
    gsd_task : str
        GSDensity task directory of a noisy matrix, with fields {level} and {rep}.
    noise_levels : sequence of float
        Noise levels (x gene SD).
    n_repeats : int
        Scored repeats per level.
    n_stream_repeats : int
        Noise draws per level in 17's stream.
    seed : int
        Seed of the noise stream.

    Returns
    -------
    pd.DataFrame
        One row per method, level and repeat: "rho vs clean" and "AUROC vs truth".
    """
    X = adata.X.toarray().astype(np.float32)
    coords = np.asarray(adata.obsm["spatial"])
    y = (adata.obs["cluster"] == region).values.astype(int)
    common_genes = [g for g in genes if g in adata.var_names]
    gene_idx = [adata.var_names.get_loc(g) for g in common_genes]
    gene_stds = np.zeros(len(common_genes))
    for i, gene in enumerate(common_genes):
        x = adata[:, gene].X
        x = x.toarray().flatten() if issparse(x) else x.flatten()
        gene_stds[i] = x.std()

    def _score_all(M):
        out = {}
        for method, fn in scoring_functions.items():
            light = ad.AnnData(
                X=M,
                obs=pd.DataFrame(index=adata.obs_names),
                var=pd.DataFrame(index=adata.var_names),
            )
            light.obsm["spatial"] = coords
            out[method] = np.asarray(fn(light, common_genes), dtype=np.float64).ravel()
        return out

    clean = {**_score_all(X), "GSDensity": np.asarray(gsd_clean, dtype=np.float64)}
    rng = np.random.default_rng(seed)
    rows = []
    for level in noise_levels:
        for rep in range(n_stream_repeats):
            noise = [
                rng.normal(0, gene_stds[j] * level, X.shape[0])
                for j in range(len(gene_idx))
            ]
            if rep >= n_repeats:
                continue
            Xn = X.copy()
            for j, gi in enumerate(gene_idx):
                Xn[:, gi] += noise[j]
            noisy = _score_all(Xn)
            task = gsd_task.format(level=level, rep=0 if level == 0 else rep)
            noisy["GSDensity"] = (
                pd.read_csv(task + "gsd_scores.csv", index_col=0)
                .reindex(adata.obs_names)[region]
                .values.astype(np.float64)
            )
            for method, s in noisy.items():
                rows.append(
                    {
                        "region": region,
                        "method": method,
                        "noise_level": level,
                        "repeat": rep,
                        "rho vs clean": spearmanr(clean[method], s).correlation,
                        "AUROC vs truth": roc_auc_score(y, s),
                    }
                )
    return pd.DataFrame(rows)


def plot_metric_pair_by_method(
    table, x, y, ax, colours, ref_x=None, ref_y=None, invert_y=False, marker_col=None
):
    """
    One metric against another, one point per method and task, with 95% CI whiskers when
    f"{x}_lo"/f"{x}_hi" and f"{y}_lo"/f"{y}_hi" columns exist.

    Parameters
    ----------
    table : pd.DataFrame
        One row per method (and task), with a "method" column and the metric columns.
    x, y : str
        Metric columns.
    ax : matplotlib.axes.Axes
        Target axes.
    colours : dict
        Method label -> colour.
    ref_x, ref_y : float or None
        Dotted reference lines.
    invert_y : bool
        Invert the y axis for lower-is-better metrics (e.g. Geary's C).
    marker_col : str or None
        Column holding a marker per row (default circles).

    Returns
    -------
    matplotlib.axes.Axes
    """
    for _, r in table.iterrows():
        xerr = yerr = None
        if f"{x}_lo" in table:
            xerr = [[r[x] - r[f"{x}_lo"]], [r[f"{x}_hi"] - r[x]]]
        if f"{y}_lo" in table:
            yerr = [[r[y] - r[f"{y}_lo"]], [r[f"{y}_hi"] - r[y]]]
        ax.errorbar(
            r[x],
            r[y],
            xerr=xerr,
            yerr=yerr,
            fmt=r[marker_col] if marker_col else "o",
            color=colours[r["method"]],
            ms=5,
            mec="white",
            mew=0.4,
            elinewidth=0.6,
            capsize=0,
            label=r["method"],
        )
    if ref_x is not None:
        ax.axvline(ref_x, ls=":", c="0.6", lw=0.6)
    if ref_y is not None:
        ax.axhline(ref_y, ls=":", c="0.6", lw=0.6)
    if invert_y:
        ax.invert_yaxis()
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.grid(False)
    return ax


def plot_methods_by_metrics_summary(
    table, ax, colours, lower_is_better=(), annotate=True
):
    """
    Methods x metrics summary: circles sized and shaded by each metric scaled to 0-1 across methods
    (best = 1), raw values annotated, plus the mean rank of each method.

    Parameters
    ----------
    table : pd.DataFrame
        Index = method labels, columns = metrics.
    ax : matplotlib.axes.Axes
        Target axes.
    colours : dict
        Method label -> colour (tick labels).
    lower_is_better : tuple of str
        Metrics where smaller is better.
    annotate : bool
        Print the raw value in each cell.

    Returns
    -------
    matplotlib.axes.Axes
    """
    scaled = table.copy().astype(float)
    for c in table.columns:
        col = table[c].astype(float)
        span = col.max() - col.min()
        s = (col - col.min()) / span if span > 0 else col * 0 + 1
        scaled[c] = 1 - s if c in lower_is_better else s
    ranks = scaled.rank(ascending=False).mean(axis=1)
    cmap = plt.get_cmap("Purples")
    for i, m in enumerate(table.index):
        for j, c in enumerate(table.columns):
            s = scaled.loc[m, c]
            ax.scatter(
                j, i, s=30 + 220 * s, color=cmap(0.25 + 0.75 * s), edgecolors="none"
            )
            if annotate:
                ax.text(
                    j,
                    i,
                    f"{table.loc[m, c]:.2g}",
                    ha="center",
                    va="center",
                    fontsize=5,
                    color="white" if s > 0.6 else "black",
                )
        ax.text(
            len(table.columns) + 0.2,
            i,
            f"{ranks[m]:.1f}",
            ha="left",
            va="center",
            fontsize=6,
        )
    ax.set_xticks(range(len(table.columns)))
    ax.set_xticklabels(table.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(table.index)))
    ax.set_yticklabels(table.index)
    for lab in ax.get_yticklabels():
        lab.set_color(colours[lab.get_text()])
    ax.text(len(table.columns) + 0.2, -0.9, "mean rank", fontsize=6, ha="left")
    ax.set_xlim(-0.6, len(table.columns) + 0.9)
    ax.set_ylim(len(table.index) - 0.4, -1.2)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    return ax


# --------------------------------------------------------------------------------------------------
# REV_Fig1_simulated
# --------------------------------------------------------------------------------------------------


def plot_metric_pair_points(
    ax, table, x, y, colours, ref_x=None, ref_y=None, invert_y=False
):
    """
    One metric against another, one filled circle per method (whiskers when f"{x}_lo" and f"{x}_hi" exist).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    table : pd.DataFrame
        One row per method, with a "method" column and the metric columns.
    x, y : str
        Metric columns, e.g. AUROC (x) and AUPRC (y).
    colours : dict
        Method -> colour.
    ref_x, ref_y : float or None
        Dotted reference lines (e.g. AUROC 0.5, AUPRC = prevalence).
    invert_y : bool
        Invert the y axis for a lower-is-better metric (Geary's C), so that better is up and right.
    """
    for _, r in table.iterrows():
        xerr = (
            [[r[x] - r[f"{x}_lo"]], [r[f"{x}_hi"] - r[x]]]
            if f"{x}_lo" in table
            else None
        )
        yerr = (
            [[r[y] - r[f"{y}_lo"]], [r[f"{y}_hi"] - r[y]]]
            if f"{y}_lo" in table
            else None
        )
        ax.errorbar(
            r[x],
            r[y],
            xerr=xerr,
            yerr=yerr,
            fmt="o",
            color=colours[r["method"]],
            ms=5,
            mec="white",
            mew=0.4,
            elinewidth=0.6,
            capsize=0,
            label=r["method"],
        )
    if ref_x is not None:
        ax.axvline(ref_x, ls=":", c="0.6", lw=0.6)
    if ref_y is not None:
        ax.axhline(ref_y, ls=":", c="0.6", lw=0.6)
    if invert_y:
        ax.invert_yaxis()
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.grid(False)


def plot_scaled_metric_dot_summary(ax, table, colours, lower_is_better=()):
    """
    Methods x metrics: circles sized and shaded by each metric scaled to 0-1 across methods (best = 1), the raw
    value in each circle and the mean rank in the last column.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    table : pd.DataFrame
        Index = methods, columns = metrics.
    colours : dict
        Method -> colour (method tick labels).
    lower_is_better : tuple of str
        Metrics where smaller is better (e.g. "Geary's C").
    """
    scaled = ((table - table.min()) / (table.max() - table.min())).fillna(1.0)
    for c in lower_is_better:
        scaled[c] = 1 - scaled[c]
    ranks = scaled.rank(ascending=False).mean(axis=1)
    cmap = plt.get_cmap("Purples")
    for i, m in enumerate(table.index):
        for j, c in enumerate(table.columns):
            s = scaled.loc[m, c]
            ax.scatter(
                j, i, s=30 + 220 * s, color=cmap(0.25 + 0.75 * s), edgecolors="none"
            )
            ax.text(
                j,
                i,
                f"{table.loc[m, c]:.2g}",
                ha="center",
                va="center",
                fontsize=5,
                color="white" if s > 0.6 else "black",
            )
        ax.text(
            len(table.columns) + 0.2,
            i,
            f"{ranks[m]:.1f}",
            ha="left",
            va="center",
            fontsize=6,
        )
    ax.set_xticks(range(len(table.columns)))
    ax.set_xticklabels(table.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(table.index)))
    ax.set_yticklabels(table.index)
    for label in ax.get_yticklabels():
        label.set_color(colours[label.get_text()])
    ax.text(len(table.columns) + 0.2, -0.9, "mean rank", fontsize=6, ha="left")
    ax.set_xlim(-0.6, len(table.columns) + 0.9)
    ax.set_ylim(len(table.index) - 0.4, -1.2)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)


def plot_method_variograms(ax, variogram_df, methods, colours, max_lag, col, log=False):
    """
    Experimental variogram of every method on one axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    variogram_df : pd.DataFrame
        One row per method and lag bin ("method", "spatial_lag_upper_edge_um", "semivariance", "variance",
        "semivariance_over_variance", "inside_plotted_x_range").
    methods : list of str
        Methods in plotting order.
    colours : dict
        Method -> colour.
    max_lag : float
        x limit (um).
    col : str
        "semivariance" (with the variance of each score as a dashed line) or "semivariance_over_variance".
    log : bool
        Log y axis.
    """
    for m in methods:
        v = variogram_df[
            (variogram_df["method"] == m) & variogram_df["inside_plotted_x_range"]
        ]
        ax.plot(
            v["spatial_lag_upper_edge_um"],
            v[col],
            "o-",
            color=colours[m],
            lw=1.1,
            markersize=2.5,
            label=m,
        )
        if col == "semivariance":
            ax.axhline(
                v["variance"].iloc[0], color=colours[m], lw=0.6, ls="--", alpha=0.7
            )
    if col != "semivariance":
        ax.axhline(1.0, color="black", lw=0.6, ls="--")
    if log:
        ax.set_yscale("log")
    ax.set_xlim(0, max_lag)
    ax.set_xlabel("Spatial lag (µm)")
    ax.grid(False)
    sns.despine(ax=ax)


def contribution_heatmap_matrix(contrib, coords, top_n_genes, groupby, n_neighbors=6):
    """
    Matrix drawn by em.pl.gene_contributions_heatmap (order_spots="coherence"), spots grouped by cluster.

    Parameters
    ----------
    contrib : dict
        Gene -> contribution vector (adata.uns["gene_contributions"][signature]).
    coords : np.ndarray
        Spatial coordinates of the spots (adata.obsm["spatial"]).
    top_n_genes : int
        Genes with the largest mean |contribution|, in that order.
    groupby : np.ndarray
        Group label per spot (annotation bar); groups in order of first appearance.
    n_neighbors : int
        Neighbours for the coherence ordering.

    Returns
    -------
    tuple
        Row-scaled matrix (genes x ordered spots), genes, spot order, mean |contribution| per gene.
    """
    mean_abs = pd.Series(
        {gene: np.mean(np.abs(v)) for gene, v in contrib.items()}
    ).sort_values(ascending=False)
    genes = list(mean_abs.index[:top_n_genes])
    X = np.array([contrib[g] for g in genes])

    # row scaling (z-score per gene across spots)
    sd = X.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    X = (X - X.mean(axis=1, keepdims=True)) / sd

    def coherence_order(idx):
        # mean |scaled contribution| smoothed over the 6-NN graph, descending
        signal = np.mean(np.abs(X[:, idx]), axis=0)
        a = AnnData(X=np.zeros((len(idx), 1)))
        a.obsm["spatial"] = coords[idx]
        sq.gr.spatial_neighbors(a, coord_type="generic", n_neighs=n_neighbors)
        W = csr_matrix(a.obsp["spatial_connectivities"])
        row_sum = np.array(W.sum(axis=1)).flatten()
        row_sum[row_sum == 0] = 1.0
        smoothed = W.multiply(1 / row_sum[:, None]).dot(signal)
        return idx[np.argsort(smoothed)[::-1]]

    spot_order = np.concatenate(
        [coherence_order(np.where(groupby == grp)[0]) for grp in dict.fromkeys(groupby)]
    )
    return X[:, spot_order], genes, spot_order, mean_abs


# --------------------------------------------------------------------------------------------------
# REV_Fig2_brain_benchmarks
# --------------------------------------------------------------------------------------------------


def plot_benchmark_curve(
    df,
    x,
    y,
    palette,
    xlabel,
    ylabel,
    title,
    hue="method",
    ylim=(0, 1.05),
    errorbar="sd",
    figsize=(6, 4),
    save=None,
):
    """
    Benchmark curve per method (mean, with an SD band over repeats), one filled circle per level.

    The layout of the plots drawn by benchmark_contamination, benchmark_noise_robustness and
    benchmark_spatial_concordance (helpers.py), with marker "o" for every method instead of one marker per method.

    Parameters
    ----------
    df : pd.DataFrame
        Tidy results of a benchmark (one row per method, level and repeat).
    x, y : str
        Columns of the perturbation level (x) and of the metric (y).
    palette : dict
        Method -> colour.
    xlabel, ylabel, title : str
        Axis labels and title.
    hue : str
        Column of the method names.
    ylim : tuple
        y limits (None keeps that end automatic).
    errorbar : str or None
        Band around the mean ("sd" over repeats; None without repeats).
    figsize : tuple of float
        Figure size.
    save : str or None
        Full path of the figure file.

    Returns
    -------
    matplotlib.axes.Axes
    """
    fig, ax = plt.subplots(figsize=figsize)
    sns.lineplot(
        data=df,
        x=x,
        y=y,
        hue=hue,
        marker="o",
        dashes=False,
        err_style="band",
        errorbar=errorbar,
        palette=palette,
        ax=ax,
    )
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.set_ylim(*ylim)
    ax.tick_params(labelsize=8)
    ax.legend(title="", fontsize=8, frameon=False)
    ax.grid(False)
    sns.despine(ax=ax)
    plt.tight_layout()
    set_box_ratio(fig)
    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")
    plt.show()
    return ax


# --------------------------------------------------------------------------------------------------
# REV_Fig3_gene_weights
# --------------------------------------------------------------------------------------------------


def compute_cv_log2fc_gene_stats(
    adata,
    target_genes,
    noise_genes,
    target_cluster,
    known_weights,
    cluster_key="cluster",
):
    """CV, log2FC (target vs rest), detection rate and known weight per gene."""
    mask_target = adata.obs[cluster_key] == target_cluster
    mask_rest = ~mask_target
    records = []
    for gene_list, is_target in [(target_genes, True), (noise_genes, False)]:
        for gene in gene_list:
            if gene not in adata.var_names:
                continue
            expr = adata[:, gene].X
            expr = expr.toarray().flatten() if issparse(expr) else expr.flatten()
            mu, sigma = expr.mean(), expr.std()
            cv = sigma / mu if mu > 0 else 0.0
            log2fc = np.log2(
                (expr[mask_target].mean() + 1e-6) / (expr[mask_rest].mean() + 1e-6)
            )
            records.append(
                {
                    "gene": gene,
                    "CV": cv,
                    "log2FC": log2fc,
                    "detection": np.mean(expr > 0),
                    "known_weight": known_weights.get(gene, np.nan),
                    "is_target": is_target,
                }
            )
    return pd.DataFrame(records)


# --------------------------------------------------------------------------------------------------
# REV_Fig4_SupFig10_1_scores
# --------------------------------------------------------------------------------------------------


def present_in_data(genes, var_names):
    """
    Signature genes present in the data, de-duplicated, in input order.

    Parameters
    ----------
    genes : list
        Gene symbols.
    var_names : pd.Index
        Genes of the data.

    Returns
    -------
    list
        Genes found in `var_names`.
    """
    return [g for g in dict.fromkeys(genes) if g in var_names]


def gene_mean_sd(X):
    """
    Mean and population SD (ddof = 0) of every gene.

    Parameters
    ----------
    X : sparse matrix
        Spots x genes expression.

    Returns
    -------
    tuple of np.ndarray
        Gene means and SDs (float64).
    """
    X = csr_matrix(X)
    m = np.asarray(X.mean(axis=0)).ravel().astype(np.float64)
    m2 = np.asarray(X.multiply(X).mean(axis=0)).ravel().astype(np.float64)
    return m, np.sqrt(np.maximum(m2 - m**2, 0))


def mean_expression_bins_fig4(gene_means, n_bins=25):
    """
    Bins of mean-expression rank, as in Scanpy's score_genes.

    Parameters
    ----------
    gene_means : np.ndarray
        Mean expression per gene.
    n_bins : int
        Number of bins.

    Returns
    -------
    np.ndarray
        Bin per gene.
    """
    return pd.qcut(
        pd.Series(gene_means).rank(method="first"), n_bins, labels=False
    ).values


def expression_matched_gene_index(bins, sig_idx, n_sets, seed, exclude=()):
    """
    Random gene sets matched to the expression bins of a signature.

    Parameters
    ----------
    bins : np.ndarray
        Expression bin per gene.
    sig_idx : np.ndarray
        Indices of the signature genes.
    n_sets : int
        Number of random sets.
    seed : int
        Seed of the draws.
    exclude : array-like
        Gene indices that may not be drawn (the signature genes are always excluded).

    Returns
    -------
    np.ndarray
        (n_sets x signature size) gene indices; the gene for signature gene i is drawn from i's bin, without
        replacement within a set unless the pool is too small.
    """
    rng = np.random.default_rng(seed)
    sig_idx = np.asarray(sig_idx)
    excl = np.zeros(len(bins), bool)
    excl[sig_idx] = True
    excl[np.asarray(list(exclude), dtype=int)] = True
    out = np.empty((n_sets, len(sig_idx)), dtype=np.int64)
    sig_bins = bins[sig_idx]
    pools = {b: np.where((bins == b) & ~excl)[0] for b in np.unique(sig_bins)}
    for s in range(n_sets):
        for b, pool in pools.items():
            pos = np.where(sig_bins == b)[0]
            out[s, pos] = rng.choice(pool, size=len(pos), replace=len(pos) > len(pool))
    return out


def expression_matched_updown_sets(up, down, n_sets, seed, exclude, bins, var_names):
    """
    Random up/down gene sets matched to the expression bins of a signed signature (disjoint within a set).

    Parameters
    ----------
    up, down : list
        Signature genes (present in the data).
    n_sets : int
        Number of random sets.
    seed : int
        Seed of the draws.
    exclude : np.ndarray
        Gene indices that may not be drawn.
    bins : np.ndarray
        Expression bin per gene (`mean_expression_bins_fig4`).
    var_names : pd.Index
        Genes of the data, in the order of `bins`.

    Returns
    -------
    list of tuple
        (up genes, down genes) per random set.
    """
    up_sets = expression_matched_gene_index(
        bins, var_names.get_indexer(up), n_sets, seed=seed, exclude=exclude
    )
    out = []
    for k in range(n_sets):
        dn = expression_matched_gene_index(
            bins,
            var_names.get_indexer(down),
            1,
            seed=seed * 1000 + k,
            exclude=np.r_[exclude, up_sets[k]],
        )[0]
        out.append((list(var_names[up_sets[k]]), list(var_names[dn])))
    return out


def package_fitted_signed_weights(adata, up, down):
    """
    Signed weights as em.tl.score fits them for gene_set={name: {"up": up, "down": down}}.

    Parameters
    ----------
    adata : AnnData
        Log-normalised data.
    up, down : list
        Up and down genes.

    Returns
    -------
    dict
        Gene -> signed weight.
    """
    var_names = pd.Index(adata.var_names)
    up_p = present_in_data(up, var_names)
    dn_p = [g for g in present_in_data(down, var_names) if g not in up_p]
    genes = up_p + dn_p
    cv = em.tl.infer_gene_weights(adata, genes)
    cols = adata[:, genes].X
    cols = cols.toarray() if hasattr(cols, "toarray") else np.asarray(cols)
    is_up = np.array([True] * len(up_p) + [False] * len(dn_p))
    return dict(
        zip(
            genes,
            _signed_weights(
                cols.astype(np.float64), is_up, np.array([cv[g] for g in genes])
            ),
        )
    )


# --------------------------------------------------------------------------------------------------
# REV_Fig4_SupFig10_1b_calibration
# --------------------------------------------------------------------------------------------------


def expression_matched_updown_sets_fig4(
    up, down, n_sets, seed, exclude, bins, var_names
):
    """
    Random up/down gene sets matched to the expression bins of a signed signature (disjoint within a set).

    Parameters
    ----------
    up, down : list
        Signature genes (present in the data).
    n_sets : int
        Number of random sets.
    seed : int
        Seed of the draws.
    exclude : np.ndarray
        Gene indices that may not be drawn.
    bins : np.ndarray
        Expression bin per gene (`mean_expression_bins_fig4`).
    var_names : pd.Index
        Genes of the data, in the order of `bins`.

    Returns
    -------
    list of tuple
        (up genes, down genes) per random set.
    """
    up_sets = expression_matched_gene_index(
        bins, var_names.get_indexer(up), n_sets, seed=seed, exclude=exclude
    )
    out = []
    for k in range(n_sets):
        dn = expression_matched_gene_index(
            bins,
            var_names.get_indexer(down),
            1,
            seed=seed * 1000 + k,
            exclude=np.r_[exclude, up_sets[k]],
        )[0]
        out.append((list(var_names[up_sets[k]]), list(var_names[dn])))
    return out


def fitted_signed_weights(up, down, X_csc, var_names, cv_all):
    """
    Fitted signed weights (the build_signed_weights rule, as em.tl.score applies it to up/down lists).

    Parameters
    ----------
    up, down : list
        Up and down genes (present in the data, disjoint).
    X_csc : scipy.sparse.csc_matrix
        Spots x genes expression.
    var_names : pd.Index
        Genes of the data (columns of `X_csc`).
    cv_all : np.ndarray
        Coefficient of variation of every gene.

    Returns
    -------
    dict
        Gene -> signed weight.
    """
    idx = var_names.get_indexer(up + down)
    cols = X_csc[:, idx].toarray().astype(np.float64)
    is_up = np.arange(len(idx)) < len(up)
    return dict(zip(up + down, _signed_weights(cols, is_up, cv_all[idx])))


# --------------------------------------------------------------------------------------------------
# REV_Fig4_SupFig10_2_panels
# --------------------------------------------------------------------------------------------------


def steiger_z_fig4(r12, r13, r23, n):
    """
    Steiger's Z for two dependent correlations sharing one variable (09_01_HRD_in_breast_cancer).

    Parameters
    ----------
    r12, r13 : float
        Correlations of the reference with the two scores.
    r23 : float
        Correlation between the two scores.
    n : int
        Number of spots.

    Returns
    -------
    tuple of float
        Z and two-sided p.
    """
    r_bar = (r12 + r13) / 2
    f = (1 - r23) / (2 * (1 - r_bar**2))
    h = (1 - f * r_bar**2) / (1 - r_bar**2)
    z1, z2 = np.arctanh(r12), np.arctanh(r13)
    Z = (z1 - z2) * np.sqrt((n - 3) / (2 * (1 - r23) * h))
    p = 2 * (1 - stats.norm.cdf(abs(Z)))
    return Z, p


def bootstrap_spearman_fig4(x, y, n_boot=1000, seed=0):
    """
    Bootstrapped Spearman correlation of two per-spot scores (09_01_HRD_in_breast_cancer).

    Every call with the same seed draws the same resamples, so the replicates of different pairs are paired.

    Parameters
    ----------
    x, y : np.ndarray
        Scores per spot.
    n_boot : int
        Number of bootstrap replicates.
    seed : int
        Seed of the resampling.

    Returns
    -------
    np.ndarray
        Spearman rho per replicate.
    """
    rng_b = np.random.default_rng(seed)
    n = len(x)
    corrs = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng_b.integers(0, n, size=n)
        corrs[b] = spearmanr(x[idx], y[idx]).statistic
    return corrs


def _morans_i_continuous_fig4(values, W):
    """
    Moran's I of a continuous score with a binary spatial weight matrix (17_benchmarking).

    Parameters
    ----------
    values : np.ndarray
        Values per spot.
    W : scipy.sparse matrix
        Spatial weights.

    Returns
    -------
    float
        Moran's I.
    """
    values = np.asarray(values, dtype=np.float64)
    n = len(values)
    z = values - values.mean()
    z2 = z @ z
    if z2 == 0:
        return 0.0
    return float((n / W.sum()) * (z @ W.dot(z) / z2))


def gearys_c_continuous(values, W):
    """
    Geary's C of a continuous score with a binary spatial weight matrix (1 = no autocorrelation, < 1 = positive).

    Parameters
    ----------
    values : np.ndarray
        Values per spot.
    W : scipy.sparse matrix
        Spatial weights.

    Returns
    -------
    float
        Geary's C.
    """
    v = np.asarray(values, dtype=np.float64)
    z = v - v.mean()
    z2 = z @ z
    if z2 == 0:
        return 1.0
    Wc = W.tocoo()
    num = (Wc.data * (v[Wc.row] - v[Wc.col]) ** 2).sum()
    return float((len(v) - 1) * num / (2 * Wc.data.sum() * z2))


def zscore_per_slide(values, batch):
    """
    Z-score values within each slide (columns separately for a matrix; SD 0 left unscaled).

    Parameters
    ----------
    values : np.ndarray
        Values per spot (spots or spots x k).
    batch : np.ndarray
        Slide of every spot.

    Returns
    -------
    np.ndarray
        Per-slide standardised values.
    """
    values = np.asarray(values, dtype=np.float64)
    out = np.empty_like(values)
    for b in pd.unique(batch):
        m = batch == b
        sd = values[m].std(0)
        out[m] = (values[m] - values[m].mean(0)) / np.where(sd > 0, sd, 1.0)
    return out


def smooth_per_slide(values, batch, coords, n_neighbors=6):
    """
    Neighbourhood average over the kNN graph of each slide, the smoothing step of em.tl.score.

    Parameters
    ----------
    values : np.ndarray
        Values per spot (spots or spots x k).
    batch : np.ndarray
        Slide of every spot.
    coords : np.ndarray
        Spatial coordinates of every spot.
    n_neighbors : int
        Number of neighbours.

    Returns
    -------
    np.ndarray
        Smoothed values.
    """
    values = np.asarray(values, dtype=np.float64)
    out = np.empty_like(values)
    for b in pd.unique(batch):
        m = batch == b
        conn = spatial_knn_connectivity(coords[m], n_neighbors)
        d = np.maximum(np.asarray(conn.sum(axis=1)).ravel(), 1e-10)
        out[m] = conn.dot(values[m]) / (d[:, None] if values.ndim == 2 else d)
    return out


def knn_spatial_lag_matrix(coords, k=6):
    """
    Row-standardised libpysal KNN(k) weights as a sparse matrix, the spatial lag of em.pl.cross_moran_scatter.

    Parameters
    ----------
    coords : np.ndarray
        Spatial coordinates.
    k : int
        Number of neighbours.

    Returns
    -------
    scipy.sparse.csr_matrix
        Spatial lag operator.
    """
    W = KNN.from_array(np.asarray(coords), k=k)
    W.transform = "r"
    return csr_matrix(W.sparse)


def cross_moran_statistics(x, y, coords, n_neighbours=6):
    """
    Score x against the spatial lag of score y on one slide, as computed in em.pl.cross_moran_scatter.

    Parameters
    ----------
    x, y : np.ndarray
        Scores of the plotted spots (tumour spots of one slide).
    coords : np.ndarray
        Spatial coordinates of the same spots.
    n_neighbours : int
        Number of nearest neighbours of the row-standardised libpysal KNN weights.

    Returns
    -------
    tuple
        Spatial lag of y, quadrant per spot, category per spot and a dict of statistics.
    """
    W = KNN.from_array(coords, k=n_neighbours)
    W.transform = "r"
    y_lag = lag_spatial(W, y)
    r, p = pearsonr(x, y_lag)
    quadrant = np.full(len(x), "LL", dtype=object)
    quadrant[(x >= 0) & (y_lag >= 0)] = "HH"
    quadrant[(x >= 0) & (y_lag < 0)] = "HL"
    quadrant[(x < 0) & (y_lag >= 0)] = "LH"
    if r >= 0:
        label = {
            "HH": "Co-enriched",
            "LL": "Co-depleted",
            "LH": "Discordant",
            "HL": "Discordant",
        }
    else:
        label = {
            "LH": "Anti-correlated",
            "HL": "Anti-correlated",
            "HH": "Discordant",
            "LL": "Discordant",
        }
    fit = linregress(x, y_lag)
    statistics = {
        "n_spots": len(x),
        "pearson_r": r,
        "p_value": p,
        "p_value_as_printed": f"{p:.2g}",
        "ols_slope": fit.slope,
        "ols_intercept": fit.intercept,
        "n_anti_correlated": int(np.isin(quadrant, ["LH", "HL"]).sum()),
        "n_discordant": int(np.isin(quadrant, ["HH", "LL"]).sum()),
    }
    return y_lag, quadrant, np.array([label[q] for q in quadrant]), statistics


def pearson_with_columns_fig4(x, Y):
    """
    Pearson r between a vector and each column of a matrix.

    Parameters
    ----------
    x : np.ndarray
        Vector (n).
    Y : np.ndarray
        Matrix (n x k).

    Returns
    -------
    np.ndarray
        k correlations.
    """
    x = np.asarray(x, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64)
    if Y.ndim == 1:
        Y = Y[:, None]
    xc = x - x.mean()
    Yc = Y - Y.mean(0)
    return (xc @ Yc) / np.sqrt((xc @ xc) * (Yc**2).sum(0))


def dutilleul_distance_classes(coords, nclass=13):
    """
    Equal-width distance classes of one set of spots for repeated Dutilleul modified t-tests.

    Parameters
    ----------
    coords : np.ndarray
        Spatial coordinates of the spots.
    nclass : int
        Number of distance classes (13, as SpatialPack::modified.ttest).

    Returns
    -------
    dict
        Class of every spot pair, the pair indices and the number of pairs per class.
    """
    D = pdist(np.asarray(coords, dtype=np.float64))
    edges = np.linspace(0, D.max() * (1 + 1e-9), nclass + 1)
    i, j = np.triu_indices(len(coords), 1)
    cls = np.digitize(D, edges[1:-1]).astype(np.int32)
    return {
        "cls": cls,
        "i": i,
        "j": j,
        "Nk": np.bincount(cls, minlength=nclass).astype(np.float64),
        "n": len(coords),
        "nclass": nclass,
    }


def dutilleul_modified_ttest(x, y, classes):
    """
    Modified t-test for the correlation of two spatially autocorrelated variables (Clifford, Richardson & Hemon
    1989; Dutilleul 1993).

    Parameters
    ----------
    x, y : np.ndarray
        Values per spot.
    classes : dict
        Distance classes of the same spots (dutilleul_distance_classes).

    Returns
    -------
    dict
        r, n, effective sample size n_eff, t and two-sided p.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    cls, i, j, Nk, n, nclass = (
        classes[k] for k in ("cls", "i", "j", "Nk", "n", "nclass")
    )
    xc, yc = x - x.mean(), y - y.mean()
    Cx = np.bincount(cls, weights=xc[i] * xc[j], minlength=nclass) / np.maximum(Nk, 1)
    Cy = np.bincount(cls, weights=yc[i] * yc[j], minlength=nclass) / np.maximum(Nk, 1)
    c0x, c0y = (xc @ xc) / n, (yc @ yc) / n
    den = n * c0x * c0y + 2 * np.sum(Nk * Cx * Cy)
    n_eff = 1 + (n * c0x) * (n * c0y) / den if den > 0 else np.nan
    r = float((xc @ yc) / np.sqrt((xc @ xc) * (yc @ yc)))
    dfree = max(n_eff - 2, 1)
    t = r * np.sqrt(dfree / max(1 - r**2, 1e-12))
    return {
        "r": r,
        "n": n,
        "n_eff": float(n_eff),
        "t": float(t),
        "p": float(2 * stats.t.sf(abs(t), dfree)),
    }


def partial_spearman(x, y, Z):
    """
    Partial Spearman correlation of x and y given the columns of Z.

    Parameters
    ----------
    x, y : np.ndarray
        Variables to correlate.
    Z : np.ndarray
        Covariates (spots x k).

    Returns
    -------
    float
        Pearson r of the rank residuals.
    """
    rx, ry = rankdata(x), rankdata(y)
    Zr = np.column_stack([np.ones(len(x))] + [rankdata(z) for z in Z.T])
    ex = rx - Zr @ np.linalg.lstsq(Zr, rx, rcond=None)[0]
    ey = ry - Zr @ np.linalg.lstsq(Zr, ry, rcond=None)[0]
    return pearsonr(ex, ey)[0]


def correlation_cluster_order(corr):
    """
    Row and column order of a correlation heatmap, as in em.pl.signature_correlation_heatmap.

    Parameters
    ----------
    corr : pd.DataFrame
        Square correlation matrix.

    Returns
    -------
    np.ndarray
        Leaf order of average-linkage clustering on 1 - correlation.
    """
    dist = 1 - corr.values
    np.fill_diagonal(dist, 0)
    return leaves_list(linkage(squareform(dist, checks=False), method="average"))


def select_extreme_genes(contributions, n_top, n_bottom):
    """
    Top (and bottom) genes by mean absolute contribution, as in em.pl.gene_contributions_heatmap.

    Parameters
    ----------
    contributions : dict
        Gene -> per-spot weighted expression.
    n_top, n_bottom : int
        Number of genes from each end.

    Returns
    -------
    tuple
        Selected genes and the mean |contribution| of every gene.
    """
    mean_abs = {gene: np.mean(np.abs(scores)) for gene, scores in contributions.items()}
    sorted_genes = sorted(mean_abs, key=mean_abs.get, reverse=True)
    top = sorted_genes[:n_top]
    bottom = sorted_genes[-n_bottom:] if n_bottom > 0 else []
    bottom = [g for g in bottom if g not in top]
    return top + bottom, mean_abs


def contribution_cluster_order(matrix):
    """
    Gene order of the gene-contribution heatmap, as in em.pl.gene_contributions_heatmap.

    Parameters
    ----------
    matrix : np.ndarray
        Genes x spots contributions.

    Returns
    -------
    np.ndarray
        Leaf order of average-linkage clustering on correlation distance.
    """
    if matrix.shape[0] < 2:
        return np.arange(matrix.shape[0])
    dist = np.nan_to_num(pdist(matrix, metric="correlation"), nan=0.0)
    return leaves_list(linkage(dist, method="average"))


def scale_contribution_rows(X):
    """
    Z-score each gene's contribution across spots (SD 0 left unscaled).

    Parameters
    ----------
    X : np.ndarray
        Genes x spots contributions.

    Returns
    -------
    np.ndarray
        Row-scaled contributions.
    """
    mean = X.mean(axis=1, keepdims=True)
    std = X.std(axis=1, keepdims=True)
    std[std == 0] = 1.0
    return (X - mean) / std


def contribution_morans_i(contributions, genes, coords, n_neighbors=6):
    """
    Moran's I of each gene's contribution (tie-breaker of the gene sort of em.pl.gene_contributions_heatmap).

    Parameters
    ----------
    contributions : dict
        Gene -> per-spot weighted expression.
    genes : list
        Selected genes.
    coords : np.ndarray
        Spatial coordinates of all spots.
    n_neighbors : int
        Neighbours of the generic spatial graph.

    Returns
    -------
    dict
        Gene -> Moran's I.
    """
    a = AnnData(X=np.array([contributions[g] for g in genes]).T)
    a.obsm["spatial"] = coords
    a.var_names = genes
    sq.gr.spatial_neighbors(a, coord_type="generic", n_neighs=n_neighbors)
    sq.gr.spatial_autocorr(a, mode="moran")
    return a.uns["moranI"]["I"].to_dict()


def spot_coherence(data, coords, n_neighbors=6):
    """
    Spot order key of em.pl.gene_contributions_heatmap: kNN average of the mean |z| across the selected genes.

    Parameters
    ----------
    data : np.ndarray
        Genes x spots z-scores of one slide.
    coords : np.ndarray
        Spatial coordinates of the same spots.
    n_neighbors : int
        Neighbours of the generic spatial graph.

    Returns
    -------
    np.ndarray
        Coherence per spot.
    """
    signal = np.mean(np.abs(data), axis=0)
    a = AnnData(X=np.zeros((coords.shape[0], 1)))
    a.obsm["spatial"] = coords
    sq.gr.spatial_neighbors(a, coord_type="generic", n_neighs=n_neighbors)
    W = csr_matrix(a.obsp["spatial_connectivities"])
    row_sum = np.array(W.sum(axis=1)).flatten()
    row_sum[row_sum == 0] = 1.0
    W = W.multiply(1 / row_sum[:, None])
    return W.dot(signal)


def colour_shades(base_colour, n, lightest=0.86, darkest=0.30):
    """
    Shades of the hue of one colour, from dark to light.

    Parameters
    ----------
    base_colour : str
        Hex colour whose hue and saturation are kept.
    n : int
        Number of shades.
    lightest, darkest : float
        HLS lightness range.

    Returns
    -------
    list of str
        Hex colours.
    """
    h, _, s = colorsys.rgb_to_hls(
        *(int(base_colour[i : i + 2], 16) / 255 for i in (1, 3, 5))
    )
    return [
        "#%02x%02x%02x"
        % tuple(round(v * 255) for v in colorsys.hls_to_rgb(h, lightness, s))
        for lightness in np.linspace(darkest, lightest, n)
    ]


# --------------------------------------------------------------------------------------------------
# REV_Fig5_HD_SupFig11_visiumhd
# --------------------------------------------------------------------------------------------------


def visiumhd_morans_i(values, W):
    """
    Moran's I of a continuous value on a spatial weight matrix.

    Parameters
    ----------
    values : np.ndarray
        Value per bin.
    W : scipy.sparse matrix
        Spatial kNN connectivity.

    Returns
    -------
    float
    """
    values = np.asarray(values, dtype=np.float64)
    n = len(values)
    z = values - values.mean()
    z2 = z @ z
    if z2 == 0:
        return 0.0
    return float((n / W.sum()) * (z @ W.dot(z) / z2))


def visiumhd_gearys_c(values, W):
    """
    Geary's C with a binary spatial weight matrix (1 = no autocorrelation, < 1 = positive autocorrelation).

    Parameters
    ----------
    values : np.ndarray
        Value per bin.
    W : scipy.sparse matrix
        Binary spatial kNN connectivity.

    Returns
    -------
    float
    """
    z = np.asarray(values, dtype=np.float64)
    z = z - z.mean()
    z2 = z @ z
    if z2 == 0:
        return 1.0
    rs = np.asarray(W.sum(axis=1)).ravel()
    cs = np.asarray(W.sum(axis=0)).ravel()
    num = (rs * z**2).sum() + (cs * z**2).sum() - 2 * z @ W.dot(z)
    return float((len(z) - 1) * num / (2 * W.sum() * z2))


def visiumhd_f1(y, pred):
    """F1 of a binary call against a binary proxy."""
    y, pred = np.asarray(y).astype(bool), np.asarray(pred).astype(bool)
    denom = pred.sum() + y.sum()
    return float(2 * (y & pred).sum() / denom) if denom else 0.0


def visiumhd_f1_topk(y, s):
    """F1 when the k highest-scoring bins are called, k = number of proxy bins (precision = recall)."""
    pred = np.zeros(len(y), bool)
    pred[np.argsort(-np.asarray(s), kind="stable")[: int(np.sum(y))]] = True
    return visiumhd_f1(y, pred)


def visiumhd_otsu_threshold(s, nbins=256):
    """Label-free Otsu threshold of a score."""
    hist, edges = np.histogram(np.asarray(s, dtype=np.float64), bins=nbins)
    mids = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * mids) / np.maximum(w0, 1)
    m1 = (np.sum(hist * mids) - np.cumsum(hist * mids)) / np.maximum(w1, 1)
    return float(mids[int(np.argmax(w0 * w1 * (m0 - m1) ** 2))])


def visiumhd_ring_distance(mask, coords, bin_px, bin_um=16.0):
    """
    Distance (um) of every bin to the nearest bin of a region (0 inside).

    Parameters
    ----------
    mask : np.ndarray
        Boolean region mask.
    coords : np.ndarray
        Bin coordinates (n x 2, plotted pixel units).
    bin_px : float
        Bin spacing in pixel units.
    bin_um : float
        Bin spacing in um.

    Returns
    -------
    np.ndarray
    """
    d, _ = cKDTree(coords[mask]).query(coords, k=1)
    d = d / bin_px * bin_um
    d[mask] = 0.0
    return d


def visiumhd_marker_stretch(band, control, marker, coords, eps, k=25):
    """
    Stretch of a band of bins that expresses a marker: the marker averaged over the k nearest bins of the band, above its
    highest such average in a negative-control band, largest spatially connected part.

    Parameters
    ----------
    band, control : np.ndarray
        Boolean masks per bin: the band and the negative control.
    marker : np.ndarray
        Marker expression per bin.
    coords : np.ndarray
        Bin coordinates.
    eps : float
        DBSCAN radius that defines connected bins (coordinate units).
    k : int
        Number of nearest bins that are averaged.

    Returns
    -------
    np.ndarray
        Boolean mask per bin.
    """

    def local_mean(mask):
        idx = np.flatnonzero(mask)
        nn = cKDTree(coords[idx]).query(coords[idx], k=k)[1]
        return idx, marker[idx][nn].mean(axis=1)

    idx, local = local_mean(band)
    keep = idx[local > local_mean(control)[1].max()]
    labels = DBSCAN(eps=eps, min_samples=3).fit_predict(coords[keep])
    largest = pd.Series(labels[labels >= 0]).value_counts().index[0]
    out = np.zeros(len(marker), dtype=bool)
    out[keep[labels == largest]] = True
    return out


def visiumhd_readouts(
    s, in_region, dist, method, ct, adata, gsd_positive, W, ring_edges
):
    """
    Readouts of one map against a region proxy.

    Parameters
    ----------
    s : np.ndarray
        Score per bin.
    in_region : np.ndarray
        Boolean region proxy.
    dist : np.ndarray
        Distance (um) to the region.
    method : str
        Method label (null-model and GSDensity calls are added for those methods).
    ct : str
        Signature.
    adata : AnnData
        Object with the null-model output in `adata.uns["null_model"]`.
    gsd_positive : dict
        Signature -> boolean GSDensity antimode call per bin.
    W : scipy.sparse matrix
        Binary 6-NN spatial connectivity (Moran's I and Geary's C).
    ring_edges : list of int
        Distance bands outside each region (um).

    Returns
    -------
    dict
    """
    ring = (dist > 0) & (dist <= ring_edges[2])
    far = dist > ring_edges[-1]
    out = {
        "AUROC vs Leiden proxy": roc_auc_score(in_region, s),
        "AUPRC vs Leiden proxy": average_precision_score(in_region, s),
        "F1 top-k": visiumhd_f1_topk(in_region, s),
        "F1 > 0": visiumhd_f1(in_region, s > 0),
        "F1 Otsu": visiumhd_f1(in_region, s > visiumhd_otsu_threshold(s)),
        "F1 median": visiumhd_f1(in_region, s > np.median(s)),
        "Moran's I (6-NN)": visiumhd_morans_i(s, W),
        "Geary's C (6-NN)": visiumhd_gearys_c(s, W),
        "% bins > 0": 100 * (s > 0).mean(),
        "% region bins > 0": 100 * (s[in_region] > 0).mean(),
        "% outside bins > 0": 100 * (s[~in_region] > 0).mean(),
        "% ring bins (16-160 um) < 0": 100 * (s[ring] < 0).mean(),
        "% far bins (> 640 um) < 0": 100 * (s[far] < 0).mean(),
        "median region": float(np.median(s[in_region])),
        "median ring (16-160 um)": float(np.median(s[ring])),
        "median far (> 640 um)": float(np.median(s[far])),
        # < 0: the ring sits below the far background (a halo); > 0: scores fall off monotonically
        "ring dip (median ring - median far) / SD": float(
            (np.median(s[ring]) - np.median(s[far])) / s.std()
        ),
    }
    if method == "EnrichMap (null)":
        p_emp = adata.uns["null_model"][f"{ct}_null"]["p_empirical"]
        out.update(
            {
                "F1 z > 1.96": visiumhd_f1(in_region, s > 1.96),
                "F1 p_emp <= 0.05": visiumhd_f1(in_region, p_emp <= 0.05),
                "% bins z > 1.96": 100 * (s > 1.96).mean(),
                "% region bins z > 1.96": 100 * (s[in_region] > 1.96).mean(),
                "% outside bins z > 1.96": 100 * (s[~in_region] > 1.96).mean(),
                "% bins p_emp <= 0.05": 100 * (p_emp <= 0.05).mean(),
                "% region bins p_emp <= 0.05": 100 * (p_emp[in_region] <= 0.05).mean(),
                "% outside bins p_emp <= 0.05": 100
                * (p_emp[~in_region] <= 0.05).mean(),
                "% ring bins (16-160 um) z < -1.96": 100 * (s[ring] < -1.96).mean(),
                "% far bins (> 640 um) z < -1.96": 100 * (s[far] < -1.96).mean(),
            }
        )
    if method == "GSDensity":
        pos = gsd_positive[ct]
        out.update(
            {
                "F1 antimode label": visiumhd_f1(in_region, pos),
                "% bins antimode positive": 100 * pos.mean(),
                "% region bins antimode positive": 100 * pos[in_region].mean(),
                "% outside bins antimode positive": 100 * pos[~in_region].mean(),
            }
        )
    return out


def visiumhd_plot_metric_pair(ax, t, x, y, colours):
    """
    One filled circle per method, with 95% CI whiskers where `t` has f"{x}_lo" / f"{x}_hi".

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    t : pd.DataFrame
        Rows of `metrics` with a "method" column.
    x, y : str
        Metric columns.
    colours : dict
        Method -> colour.
    """
    for _, r in t.iterrows():
        xerr = (
            [[r[x] - r[f"{x}_lo"]], [r[f"{x}_hi"] - r[x]]] if f"{x}_lo" in t else None
        )
        yerr = (
            [[r[y] - r[f"{y}_lo"]], [r[f"{y}_hi"] - r[y]]] if f"{y}_lo" in t else None
        )
        ax.errorbar(
            r[x],
            r[y],
            xerr=xerr,
            yerr=yerr,
            fmt="o",
            color=colours[r["method"]],
            ms=5,
            mec="white",
            mew=0.4,
            elinewidth=0.6,
            capsize=0,
        )
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.grid(False)


def visiumhd_plot_summary(table, ax, colours, lower_is_better=()):
    """
    Methods x metrics summary: circles sized and shaded by each metric scaled to 0-1 across methods (best = 1),
    raw values annotated, and the mean rank over metrics.

    Parameters
    ----------
    table : pd.DataFrame
        Index = methods, columns = metrics.
    ax : plt.Axes
        Target axes.
    colours : dict
        Method -> colour (tick labels).
    lower_is_better : tuple of str
        Metrics where smaller is better.
    """
    scaled = table.astype(float).copy()
    for c in table.columns:
        col = table[c].astype(float)
        span = col.max() - col.min()
        s = (col - col.min()) / span if span > 0 else col * 0 + 1
        scaled[c] = 1 - s if c in lower_is_better else s
    ranks = scaled.rank(ascending=False).mean(axis=1)
    cmap = plt.get_cmap("Purples")
    for i, m in enumerate(table.index):
        for j, c in enumerate(table.columns):
            s = scaled.loc[m, c]
            ax.scatter(
                j, i, s=30 + 220 * s, color=cmap(0.25 + 0.75 * s), edgecolors="none"
            )
            ax.text(
                j,
                i,
                f"{table.loc[m, c]:.2g}",
                ha="center",
                va="center",
                fontsize=5,
                color="white" if s > 0.6 else "black",
            )
        ax.text(
            len(table.columns) + 0.2,
            i,
            f"{ranks[m]:.1f}",
            ha="left",
            va="center",
            fontsize=6,
        )
    ax.text(len(table.columns) + 0.2, -0.9, "mean rank", fontsize=6, ha="left")
    ax.set_xticks(range(len(table.columns)))
    ax.set_xticklabels(table.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(table.index)))
    ax.set_yticklabels(table.index)
    for lab in ax.get_yticklabels():
        lab.set_color(colours[lab.get_text()])
    ax.set_xlim(-0.6, len(table.columns) + 0.9)
    ax.set_ylim(len(table.index) - 0.4, -1.2)
    ax.grid(False)
    for side in ax.spines.values():
        side.set_visible(False)


# --------------------------------------------------------------------------------------------------
# REV_Fig5_platforms_xenium_merfish_imc
# --------------------------------------------------------------------------------------------------


def pseudobulk_deseq2_signature(
    adata, is_target, cell_type, padj, lfc, top_n, k=3, seed=0
):
    """
    Up-regulated genes of a cell type against all other cells from pseudo-bulk DESeq2.

    Parameters
    ----------
    adata : AnnData
        Integer counts in `adata.X`.
    is_target : np.ndarray
        Boolean mask of the target cells.
    cell_type : str
        Name of the target group.
    padj, lfc : float
        Thresholds on the adjusted p-value and the log2 fold change.
    top_n : int
        Number of genes kept (largest log2 fold change first).
    k : int
        Pseudo-replicates per group.
    seed : int
        Seed of the pseudo-replicate assignment.

    Returns
    -------
    list
        Signature genes.
    """
    group = pd.Series(np.where(is_target, cell_type, "Rest"), index=adata.obs_names)

    # Create pseudo-replicates
    # Randomly assign spots to k pseudo-replicates within each group.
    rng = np.random.default_rng(seed)
    rep_labels = []
    for grp in [cell_type, "Rest"]:
        mask = group == grp
        n = mask.sum()
        assignments = rng.integers(0, k, size=n)
        rep_labels.append(
            pd.Series(
                [f"{grp}_rep{i}" for i in assignments],
                index=adata.obs_names[mask],
            )
        )
    pseudo_rep = pd.concat(rep_labels).reindex(adata.obs_names)

    # Aggregate raw counts per pseudo-replicate
    counts_mat = adata.X

    # densify if sparse
    if issparse(counts_mat):
        counts_mat = counts_mat.toarray()

    counts_df = pd.DataFrame(counts_mat, index=adata.obs_names, columns=adata.var_names)
    counts_df["pseudo_rep"] = pseudo_rep.values

    pseudobulk = counts_df.groupby("pseudo_rep").sum()

    # Build sample metadata
    sample_meta = pd.DataFrame({"pseudo_rep": pseudobulk.index})
    sample_meta["group"] = sample_meta["pseudo_rep"].str.rsplit("_rep", n=1).str[0]
    sample_meta = sample_meta.set_index("pseudo_rep")

    # Pre-filter low-count genes
    gene_keep = pseudobulk.columns[pseudobulk.sum(axis=0) >= 10]
    pseudobulk = pseudobulk[gene_keep]

    dds = DeseqDataSet(
        counts=pseudobulk,
        metadata=sample_meta,
        design="~group",
        ref_level=["group", "Rest"],
        n_cpus=2,  # parallelism only; the results do not depend on it
    )

    dds.deseq2()

    stat = DeseqStats(dds, contrast=["group", cell_type, "Rest"], n_cpus=2)
    stat.summary()

    results = stat.results_df.copy()

    # Filter for upregulated genes
    up = results[
        (results["padj"] < padj) & (results["log2FoldChange"] > lfc)
    ].sort_values("log2FoldChange", ascending=False)

    print(f"Upregulated genes in {cell_type}: {len(up)}")
    return up.head(top_n).index.tolist()


def platform_run_score(runs, method_keys, run, method, sig, graph="delaunay"):
    """
    Score of one run, method and signature.

    Parameters
    ----------
    runs : dict
        Run -> {"adata", "graphs" (graph -> score-key suffix), "gsd" (GSDensity suffix), ...}.
    method_keys : dict
        Method -> score column template with {sig} and {suffix}.
    run, method, sig, graph : str
        Run, method label, signature and graph.

    Returns
    -------
    np.ndarray
    """
    r = runs[run]
    suffix = r["gsd"] if method == "GSDensity" else r["graphs"][graph]
    return (
        r["adata"]
        .obs[method_keys[method].format(sig=sig, suffix=suffix)]
        .values.astype(float)
    )


def platform_run_p_empirical(runs, run, sig, graph="delaunay"):
    """Null-model empirical p-value of one run and signature (`runs` as in platform_run_score)."""
    r = runs[run]
    return (
        r["adata"]
        .uns["null_model"][f"{sig}_null{r['graphs'][graph]}"]["p_empirical"]
        .astype(float)
    )


def platform_run_gsd_label(runs, run, sig):
    """GSDensity antimode label of one run and signature (`runs` as in platform_run_score)."""
    r = runs[run]
    return r["adata"].obs[f"{sig}_gsdensity{r['gsd']}_label"].values


def platform_binary_calls(method, s, p_empirical=None, gsd_label=None):
    """
    Binary calls of one method: the original rule score > 0, plus z > 1.96 and p_empirical <= 0.05 for the null
    model and the native antimode label for GSDensity.

    Parameters
    ----------
    method : str
        Method label.
    s : np.ndarray
        Scores.
    p_empirical : np.ndarray, optional
        Null-model empirical p-value per cell.
    gsd_label : np.ndarray, optional
        GSDensity antimode label per cell ("positive", "negative", NaN).

    Returns
    -------
    dict
        Rule -> boolean array.
    """
    out = {"score > 0": s > 0}
    if method == "EnrichMap (null)":
        out["z > 0"] = s > 0
        out["z > 1.96"] = s > 1.96
        out["p_empirical <= 0.05"] = p_empirical <= 0.05
    if method == "GSDensity":
        out["antimode positive"] = np.asarray(gsd_label).astype(str) == "positive"
    return out


def platform_kmeans_blocks(
    xy, sections=None, cells_per_block=1500, min_blocks=25, seed=0
):
    """
    K-means blocks of the coordinates (within each section), the resampling units of the spatial block bootstrap.

    Parameters
    ----------
    xy : np.ndarray
        Coordinates (n x 2).
    sections : np.ndarray, optional
        Section per cell.
    cells_per_block : int
        Target block size.
    min_blocks : int
        Minimum number of blocks (shared over sections).
    seed : int
        K-means seed.

    Returns
    -------
    np.ndarray
        Block per cell.
    """
    secs = np.zeros(len(xy), int) if sections is None else np.asarray(sections)
    usecs = pd.unique(secs)
    out = np.empty(len(xy), dtype=object)
    for b in usecs:
        idx = np.where(secs == b)[0]
        k = max(
            min_blocks // len(usecs) if len(usecs) > 1 else min_blocks,
            int(round(len(idx) / cells_per_block)),
            2,
        )
        km = KMeans(n_clusters=k, n_init=4, random_state=seed).fit(xy[idx])
        out[idx] = [f"{b}_{l}" for l in km.labels_]
    return pd.factorize(out)[0]


class WeightedRankMetrics:
    """
    AUROC and average precision of one score vector under cell weights (as sklearn, tied scores grouped), so that
    each bootstrap resample only changes the weights.
    """

    def __init__(self, y, s):
        o = np.argsort(s, kind="mergesort")
        self.o, self.y = o, np.asarray(y, float)[o]
        ss = np.asarray(s, float)[o]
        self.starts = np.r_[0, np.flatnonzero(np.diff(ss) != 0) + 1]

    def __call__(self, w=None):
        ws = np.ones(len(self.o)) if w is None else np.asarray(w, float)[self.o]
        gp = np.add.reduceat(ws * self.y, self.starts)
        gn = np.add.reduceat(ws * (1 - self.y), self.starts)
        auc = (gp * (np.cumsum(gn) - gn + 0.5 * gn)).sum() / (gp.sum() * gn.sum())
        tp, fp = np.cumsum(gp[::-1]), np.cumsum(gn[::-1])
        prec = np.divide(tp, tp + fp, out=np.zeros_like(tp), where=(tp + fp) > 0)
        ap = np.sum(np.diff(np.r_[0.0, tp / tp[-1]]) * prec)
        return auc, ap


def platform_block_bootstrap(y, scores, blocks, n_boot=200, seed=0):
    """
    Spatial block bootstrap of AUROC and AUPRC (blocks resampled with replacement), paired across score vectors.

    Parameters
    ----------
    y : np.ndarray
        Boolean target indicator.
    scores : dict
        Key -> score vector.
    blocks : np.ndarray
        Block per cell.
    n_boot : int
        Resamples.
    seed : int
        Seed.

    Returns
    -------
    dict
        Key -> {"AUROC", "AUPRC", "boot AUROC", "boot AUPRC"}.
    """
    y = np.asarray(y, bool)
    rankers = {k: WeightedRankMetrics(y, v) for k, v in scores.items()}
    blocks = pd.factorize(np.asarray(blocks))[0]
    nb = blocks.max() + 1
    rng = np.random.default_rng(seed)
    out = {}
    for k, r in rankers.items():
        auc, ap = r()
        out[k] = {"AUROC": auc, "AUPRC": ap, "boot AUROC": [], "boot AUPRC": []}
    for _ in range(n_boot):
        w = np.bincount(rng.integers(0, nb, nb), minlength=nb)[blocks].astype(float)
        if y[w > 0].all() or not y[w > 0].any():
            continue
        for k, r in rankers.items():
            auc, ap = r(w)
            out[k]["boot AUROC"].append(auc)
            out[k]["boot AUPRC"].append(ap)
    return out


def platform_ci95(a):
    """2.5 and 97.5 percentiles."""
    return float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))


def platform_tpr_at_fpr(y, s, sections=None, fpr=0.10):
    """
    % of target cells above the (1 - fpr) quantile of the non-target cells of their section.

    Parameters
    ----------
    y : np.ndarray
        Boolean target indicator.
    s : np.ndarray
        Score.
    sections : np.ndarray, optional
        Section per cell.
    fpr : float
        False-positive rate defining the threshold.

    Returns
    -------
    float
        TPR in %.
    """
    secs = np.zeros(len(s), int) if sections is None else np.asarray(sections)
    hit = np.zeros(len(s), bool)
    for b in pd.unique(secs):
        m = secs == b
        if (m & ~y).sum() == 0 or (m & y).sum() == 0:
            continue
        hit[m] = s[m] > np.quantile(s[m & ~y], 1 - fpr)
    return 100 * hit[y].mean()


def platform_morans_i(values, W):
    """Moran's I of a continuous value on a spatial weight matrix."""
    values = np.asarray(values, dtype=np.float64)
    n = len(values)
    z = values - values.mean()
    z2 = z @ z
    if z2 == 0:
        return 0.0
    return float((n / W.sum()) * (z @ W.dot(z) / z2))


def platform_gearys_c(values, W):
    """Geary's C with binary symmetric weights (1 = no autocorrelation, < 1 = positive)."""
    x = np.asarray(values, dtype=np.float64)
    z = x - x.mean()
    num = (x**2 * np.asarray(W.sum(axis=1)).ravel()).sum() + (
        x**2 * np.asarray(W.sum(axis=0)).ravel()
    ).sum()
    num -= 2 * x @ W.dot(x)
    return float((len(x) - 1) * num / (2 * W.sum() * (z @ z)))


def platform_delaunay_weights(xy, sections=None):
    """
    Binary Delaunay 1-ring weights (within each section), the common weights of Moran's I and Geary's C.

    Parameters
    ----------
    xy : np.ndarray
        Coordinates (n x 2).
    sections : np.ndarray, optional
        Section per cell.

    Returns
    -------
    scipy.sparse.csr_matrix
        n x n binary weights, block diagonal over sections.
    """
    secs = np.zeros(len(xy), int) if sections is None else np.asarray(sections)
    rows, cols = [], []
    for b in pd.unique(secs):
        idx = np.where(secs == b)[0]
        tmp = sc.AnnData(
            X=np.zeros((len(idx), 1), dtype=np.float32), obsm={"spatial": xy[idx]}
        )
        sq.gr.spatial_neighbors(tmp, coord_type="generic", delaunay=True)
        c = tmp.obsp["spatial_connectivities"].tocoo()
        rows.append(idx[c.row])
        cols.append(idx[c.col])
    rows, cols = np.concatenate(rows), np.concatenate(cols)
    return csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(xy), len(xy)))


def enrichmap_null_control_sets(adata, genes, n_sets=100, random_state=0):
    """
    Control gene sets of one unsigned signature, drawn as em.tl.score(null_model=True) draws them.

    Each signature gene is replaced by a random non-signature gene of similar mean expression and CV; the candidate
    pools are widened until the control sets share less than 10% of their genes (mean Jaccard).

    Parameters
    ----------
    adata : AnnData
        Log-normalised expression, as scored.
    genes : list
        Signature genes.
    n_sets : int
        Number of control sets (package default).
    random_state : int
        Seed (package default).

    Returns
    -------
    dict
        "idx" (signature genes, columns of adata), "ctrl" (control sets x signature genes), "overlap" (mean Jaccard
        between control sets) and "coarse" (the package's rule for a coarse null).
    """
    # gene mean and CV as the package computes them (float64, 1,024 columns at a time)
    X = adata.X.tocsc() if issparse(adata.X) else np.asarray(adata.X)
    n, n_genes = adata.n_obs, adata.n_vars
    gene_mean, gene_sd = np.empty(n_genes), np.empty(n_genes)
    for a0 in range(0, n_genes, 1024):
        sl = slice(a0, min(a0 + 1024, n_genes))
        Xc = X[:, sl].astype(np.float64)
        m = np.asarray(Xc.sum(axis=0)).ravel() / n
        e2 = (
            np.asarray(
                (Xc.multiply(Xc) if issparse(Xc) else np.square(Xc)).sum(axis=0)
            ).ravel()
            / n
        )
        gene_mean[sl], gene_sd[sl] = m, np.sqrt(np.maximum(e2 - m**2, 0))
    gene_cv = np.divide(gene_sd, gene_mean, out=np.zeros(n_genes), where=gene_mean > 0)
    mean_rank = rankdata(np.round(gene_mean, 12))
    cv_rank = rankdata(np.round(gene_cv, 12))
    name_rank = rankdata(np.asarray(adata.var_names, dtype=str), method="ordinal")

    # signature genes in name order; the seed comes from the sorted gene names
    var_index = {g: i for i, g in enumerate(adata.var_names)}
    idx = np.asarray([var_index[g] for g in dict.fromkeys(genes) if g in var_index])
    idx = idx[np.argsort(name_rank[idx], kind="stable")]
    rng = np.random.default_rng(
        [random_state, zlib.crc32(";".join(adata.var_names[idx].astype(str)).encode())]
    )
    free = np.flatnonzero(~adata.var_names.isin(adata.var_names[idx]))
    k = len(idx)
    n_mean, n_cv = _score._N_MEAN_WINDOW, _score._N_MATCHED_CANDIDATES
    coarse_size = len(free) < max(5 * k, 50)
    while True:
        cand = _score._matched_candidates(
            mean_rank, cv_rank, idx, free, n_mean, n_cv, name_rank
        )
        ctrl = cand[np.arange(k), rng.integers(0, cand.shape[1], size=(n_sets, k))]
        # redraw genes that repeat within a control set
        dups_left = True
        for _ in range(100):
            order = np.argsort(ctrl, axis=1, kind="stable")
            srt = np.take_along_axis(ctrl, order, axis=1)
            dup_sorted = np.zeros(srt.shape, dtype=bool)
            dup_sorted[:, 1:] = srt[:, 1:] == srt[:, :-1]
            if not dup_sorted.any():
                dups_left = False
                break
            dup = np.zeros_like(dup_sorted)
            np.put_along_axis(dup, order, dup_sorted, axis=1)
            r_, c_ = np.nonzero(dup)
            ctrl[r_, c_] = cand[c_, rng.integers(0, cand.shape[1], size=len(c_))]
        # mean Jaccard overlap between control sets; widen the pools while it exceeds 10%
        used, inv = np.unique(ctrl, return_inverse=True)
        inc = csr_matrix(
            (np.ones(ctrl.size), (np.repeat(np.arange(n_sets), k), inv.ravel())),
            shape=(n_sets, len(used)),
        )
        inc.data[:] = 1.0
        shared = (inc @ inc.T).toarray()
        size = np.diag(shared)
        jaccard = shared / (size[:, None] + size[None, :] - shared)
        overlap = float(jaccard[~np.eye(n_sets, dtype=bool)].mean())
        if overlap <= _score._MAX_CONTROL_OVERLAP or coarse_size or n_mean >= len(free):
            break
        n_mean, n_cv = 2 * n_mean, 2 * n_cv
    n_distinct = len({tuple(np.sort(r)) for r in ctrl})
    coarse = (
        coarse_size
        or n_distinct < n_sets
        or dups_left
        or overlap > _score._MAX_CONTROL_OVERLAP
    )
    return {"idx": idx, "ctrl": ctrl, "overlap": overlap, "coarse": coarse}


def platform_plot_violins(rows, colours, call_rule):
    """
    Violin and box plots of the scores of target and rest cells, one column per method.

    Parameters
    ----------
    rows : list of tuple
        One (title, target label, target mask, {method: scores}, null-model p_empirical, GSDensity label) per row.
    colours : dict
        Method -> colour (one column per method).
    call_rule : dict
        Method -> rule of platform_binary_calls whose % of cells is shown.

    Returns
    -------
    pd.DataFrame
        Plotted values per cell (GSDensity as log1p(relevance x n)).
    """
    fig, axes = plt.subplots(
        len(rows),
        len(colours),
        figsize=(1.6 * len(colours), 2.05 * len(rows)),
        squeeze=False,
    )
    plotted = {}
    for r, (title, target_label, is_target, scores, p_emp, gsd_label) in enumerate(
        rows
    ):
        plotted[f"{title} | group"] = np.where(is_target, target_label, "Rest")
        for k, (method, colour) in enumerate(colours.items()):
            ax = axes[r, k]
            s = scores[method]
            values = np.log1p(s) if method == "GSDensity" else s
            plotted[f"{title} | {method}"] = values
            groups = [values[is_target], values[~is_target]]

            # Violin
            vp = ax.violinplot(
                groups,
                positions=[0, 1],
                showmedians=False,
                showextrema=False,
                widths=0.75,
            )
            for pc, face, alpha in zip(vp["bodies"], [colour, "white"], [0.6, 1.0]):
                pc.set_facecolor(to_rgba(face, alpha))
                pc.set_edgecolor(colour)
                pc.set_alpha(None)
                pc.set_linewidth(0.6)

            # Box plot overlay (IQR box + whiskers, no fliers)
            ax.boxplot(
                groups,
                positions=[0, 1],
                widths=0.12,
                patch_artist=True,
                showfliers=False,
                medianprops=dict(color="black", linewidth=1),
                whiskerprops=dict(color="black", linewidth=0.6),
                capprops=dict(color="black", linewidth=0.6),
                boxprops=dict(linewidth=0.6, facecolor="white"),
            )
            if method != "GSDensity":
                ax.axhline(0, color="black", linestyle="--", linewidth=0.6)
            if method == "EnrichMap (null)":
                ax.axhline(1.96, color="black", linestyle=":", linewidth=0.6)

            called = platform_binary_calls(method, s, p_emp, gsd_label)[
                call_rule[method]
            ]
            ax.set_xticks([0, 1])
            ax.set_xticklabels(
                [
                    f"{target_label}\n({100 * called[is_target].mean():.1f}%)",
                    f"Rest\n({100 * called[~is_target].mean():.1f}%)",
                ],
                fontsize=6,
            )
            pval = mannwhitneyu(
                s[is_target], s[~is_target], alternative="greater"
            ).pvalue
            text = "Wilcoxon p < 0.0001" if pval < 1e-4 else f"Wilcoxon p = {pval:.2g}"
            text += (
                f"\nAUROC {roc_auc_score(is_target, s):.3f}\n% = {call_rule[method]}"
            )
            lo, hi = values.min(), values.max()
            ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.5 * (hi - lo))
            ax.text(
                0.5,
                0.98,
                text,
                transform=ax.transAxes,
                ha="center",
                va="top",
                fontsize=5,
            )
            ax.set_ylabel(
                {
                    "EnrichMap (null)": "z vs matched random genes",
                    "EnrichMap (smoothing-only)": "EnrichMap score (a.u.)",
                    "GSDensity": "log1p(GSDensity relevance x n)",
                }[method],
                fontsize=6,
            )
            ax.tick_params(axis="y", labelsize=6)
            ax.set_title(f"{title}\n{method}", fontsize=6.5, color=colour)
            ax.grid(False)
            sns.despine(ax=ax)
    plt.tight_layout()
    return pd.DataFrame(plotted)


def platform_plot_roc(rows, colours):
    """
    ROC curves and fold enrichment of ground-truth cells among the top X% of cells, one line per method.

    Parameters
    ----------
    rows : list of tuple
        One (title, target mask, {method: scores}) per row.
    colours : dict
        Method -> colour.

    Returns
    -------
    pd.DataFrame
        Plotted curves (target, method, panel, x, y).
    """
    fig, axes = plt.subplots(
        len(rows), 2, figsize=(5.0, 2.3 * len(rows)), squeeze=False
    )
    curves = []
    for r, (title, is_target, scores) in enumerate(rows):
        ax_roc, ax_fe = axes[r]
        baseline = is_target.mean()
        percentiles = np.arange(1, 101)
        for method, colour in colours.items():
            s = scores[method]
            fpr, tpr, _ = roc_curve(is_target.astype(int), s)
            auc = roc_auc_score(is_target.astype(int), s)
            ax_roc.plot(
                fpr, tpr, color=colour, linewidth=1, label=f"{method}: AUC = {auc:.3f}"
            )
            fold_enrichments = [
                is_target[s >= np.percentile(s, 100 - pct)].mean() / baseline
                for pct in percentiles
            ]
            ax_fe.plot(
                percentiles, fold_enrichments, color=colour, linewidth=1, label=method
            )
            curves.append(
                pd.DataFrame(
                    {
                        "target": title,
                        "method": method,
                        "panel": "ROC",
                        "x": fpr,
                        "y": tpr,
                    }
                )
            )
            curves.append(
                pd.DataFrame(
                    {
                        "target": title,
                        "method": method,
                        "panel": "fold enrichment",
                        "x": percentiles,
                        "y": fold_enrichments,
                    }
                )
            )
        ax_roc.plot([0, 1], [0, 1], color="black", linestyle="--", linewidth=0.6)
        ax_roc.set_xlim(0, 1)
        ax_roc.set_ylim(0, 1)
        ax_roc.set_aspect("equal")
        ax_roc.set_xlabel("False positive rate")
        ax_roc.set_ylabel("True positive rate")
        ax_roc.set_title(f"{title}: ROC", fontsize=7)
        ax_roc.legend(loc="lower right", fontsize=5, frameon=False)
        ax_fe.axhline(1, color="black", linestyle="--", linewidth=0.6)
        ax_fe.set_xlabel("Top X% of cells")
        ax_fe.set_ylabel("Fold enrichment")
        ax_fe.set_title(f"{title}: ground truth cells recovered", fontsize=7)
        for ax in (ax_roc, ax_fe):
            ax.grid(False)
            sns.despine(ax=ax)
    plt.tight_layout()
    return pd.concat(curves, ignore_index=True)


def platform_map_table(adata, keys, mask=None):
    """Coordinates and plotted values of a score map, one row per cell."""
    mask = np.ones(adata.n_obs, bool) if mask is None else mask
    df = adata.obs.loc[mask, keys].copy()
    df.insert(0, "y", adata.obsm["spatial"][mask, 1])
    df.insert(0, "x", adata.obsm["spatial"][mask, 0])
    return df


def platform_plot_metric_pair(ax, table, x, y, colours):
    """One filled circle per row, coloured by method, with 95% CI whiskers when the _lo / _hi columns exist."""
    for _, r in table.iterrows():
        xerr = (
            [[r[x] - r[f"{x}_lo"]], [r[f"{x}_hi"] - r[x]]]
            if f"{x}_lo" in table
            else None
        )
        yerr = (
            [[r[y] - r[f"{y}_lo"]], [r[f"{y}_hi"] - r[y]]]
            if f"{y}_lo" in table
            else None
        )
        ax.errorbar(
            r[x],
            r[y],
            xerr=xerr,
            yerr=yerr,
            fmt="o",
            color=colours[r["method"]],
            ms=5,
            mec="white",
            mew=0.4,
            elinewidth=0.6,
            capsize=0,
        )
    ax.set_xlabel(x, fontsize=6)
    ax.set_ylabel(y, fontsize=6)
    ax.tick_params(labelsize=5.5)
    ax.grid(False)
    sns.despine(ax=ax)


# --------------------------------------------------------------------------------------------------
# REV_GSD_base_brain
# --------------------------------------------------------------------------------------------------


def describe_brain_gsd_signature(r, size_sets, panglao, housekeeping_genes, hk_missing):
    """
    One-line description of a mouse brain signature for the GSDensity index.

    Parameters
    ----------
    r : pd.Series
        Row of the signature map (signature, gene_set, name_in_notebook).
    size_sets : dict
        Gene-set-size series of the signature pickles (key -> genes).
    panglao : dict
        PanglaoDB sets (cell type -> genes).
    housekeeping_genes : list
        Housekeeping symbols of 01_visium_platform.
    hk_missing : list
        Housekeeping symbols absent from the data.

    Returns
    -------
    str
        Description.
    """
    if r.gene_set in size_sets:
        ct = r.gene_set.split("_top")[0]
        size = int(r.gene_set.split("_top")[1])
        txt = f"top-{size} pseudobulk DESeq2 genes ({ct}_signatures.pkl, key {r.gene_set})"
        return txt + (
            "; same file content as " + ct if size == 10 and r.signature != ct else ""
        )
    if r.name_in_notebook in panglao:
        if r.name_in_notebook == "Excitatory neurons":
            return "Excitatory neurons (Allen mouse cortex reference, top 20 by detection in glutamatergic minus other cells; positive control)"
        return f"PanglaoDB {r.name_in_notebook} (first 20 mouse canonical markers in the data; absent cell type)"
    n_hk = len(housekeeping_genes) - len(hk_missing)
    if r.gene_set == "Housekeeping":
        return (
            f"{len(housekeeping_genes)} housekeeping symbols of 01 c57 (contamination level 0); "
            f"{', '.join(hk_missing)} not in the data, so {n_hk} are used"
        )
    ct = r.gene_set.split("_and_")[0]
    return f"{ct}_top10 + the {len(housekeeping_genes)} housekeeping symbols of 01 c57 ({n_hk} in the data)"


def export_gsdensity_task(adata, genesets, task_dir, params):
    """
    Write the inputs of run_gsdensity.R: log-normalised expression (genes x spots), names, coordinates, gene sets
    and parameters; outputs of an earlier run of the task are removed.

    Parameters
    ----------
    adata : AnnData
        Log-normalised expression with coordinates in obsm["spatial"].
    genesets : dict
        Gene-set name -> genes.
    task_dir : str
        Task directory.
    params : dict
        GSDensity parameters (nn_use, restart, dims, seed, spatial, spatial_n, spatial_ntimes).

    Returns
    -------
    str
        The task directory.
    """
    os.makedirs(task_dir, exist_ok=True)
    X = adata.X.T
    X = X.tocsr() if issparse(X) else csr_matrix(X)
    mmwrite(task_dir + "data.mtx", X.astype(np.float64))
    with open(task_dir + "genes.txt", "w") as f:
        f.write("\n".join(map(str, adata.var_names)) + "\n")
    with open(task_dir + "spots.txt", "w") as f:
        f.write("\n".join(map(str, adata.obs_names)) + "\n")
    xy = np.asarray(adata.obsm["spatial"])[:, :2]
    pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1]}, index=adata.obs_names).to_csv(
        task_dir + "coords.csv"
    )
    with open(task_dir + "genesets.json", "w") as f:
        json.dump({k: list(map(str, v)) for k, v in genesets.items()}, f)
    with open(task_dir + "params.json", "w") as f:
        json.dump(params, f)
    for fn in (
        "gsd_scores.csv",
        "gsd_labels.csv",
        "gsd_spatial_kld.csv",
        "gsd_timing.json",
        "gsd_error.txt",
    ):
        if os.path.exists(task_dir + fn):
            os.remove(task_dir + fn)
    return task_dir


def run_gsdensity(task_dirs, rscript, n_jobs=1, log=None):
    """
    Run GSDensity on task directories in `n_jobs` parallel R sessions of the gsdensity_env conda environment.

    Parameters
    ----------
    task_dirs : list of str
        Task directories written by `export_gsdensity_task`.
    rscript : str
        Path of run_gsdensity.R.
    n_jobs : int
        Parallel R sessions (the tasks are split between them).
    log : str or None
        Log file stem; session i writes to log with suffix _i.

    Returns
    -------
    list of int
        Exit code of each R session.
    """
    chunks = [task_dirs[i::n_jobs] for i in range(n_jobs) if task_dirs[i::n_jobs]]
    procs = []
    for i, chunk in enumerate(chunks):
        cmd = (
            GSDENSITY_ACTIVATE +
            f"Rscript {rscript} " + " ".join(chunk)
        )
        out = open(log.replace(".log", f"_{i}.log"), "w") if log else subprocess.DEVNULL
        procs.append(
            subprocess.Popen(["bash", "-c", cmd], stdout=out, stderr=subprocess.STDOUT)
        )
    for p in procs:
        p.wait()
    return [p.returncode for p in procs]


def read_gsdensity(task_dir, spots):
    """
    Read the GSDensity output of one task directory.

    Parameters
    ----------
    task_dir : str
        Task directory.
    spots : array-like
        Spot order of the slide.

    Returns
    -------
    dict
        scores (RWR probabilities, spots x gene sets), labels (antimode labels of compute.cell.label), spatial_kld
        (per gene set) and timing (seconds per step).
    """
    return {
        "scores": pd.read_csv(task_dir + "gsd_scores.csv", index_col=0).reindex(spots),
        "labels": pd.read_csv(task_dir + "gsd_labels.csv", index_col=0).reindex(spots),
        "spatial_kld": pd.read_csv(task_dir + "gsd_spatial_kld.csv", index_col=0),
        "timing": json.load(open(task_dir + "gsd_timing.json")),
    }


def gsdensity_signature_table(scores, labels, gene_set, spots):
    """
    Per-spot GSDensity output of one gene set, in spot order.

    Parameters
    ----------
    scores : pd.DataFrame
        RWR probabilities (spots x gene sets).
    labels : pd.DataFrame
        Antimode labels of compute.cell.label (spots x gene sets).
    gene_set : str
        Column to write.
    spots : pd.Index
        Spot ids of the slide.

    Returns
    -------
    pd.DataFrame
        spot_id, relevance (RWR probability, sums to 1), label (positive / negative) and relevance_scaled
        (relevance x number of spots; mean 1).
    """
    s = scores[gene_set].reindex(spots)
    return pd.DataFrame(
        {
            "spot_id": spots,
            "relevance": s.values,
            "label": labels[gene_set].reindex(spots).values,
            "relevance_scaled": s.values * len(s),
        }
    )


# --------------------------------------------------------------------------------------------------
# REV_GSD_base_breast
# --------------------------------------------------------------------------------------------------


def read_gsdensity_outputs(task_dir):
    """
    Read the GSDensity output of one task directory, in the spot order of its files.

    Parameters
    ----------
    task_dir : str
        Task directory.

    Returns
    -------
    dict
        scores (spots x gene sets, RWR probability), labels (antimode, compute.cell.label), spatial_kld (per gene set)
        and timing (seconds per step).
    """
    return {
        "scores": pd.read_csv(task_dir + "gsd_scores.csv", index_col=0),
        "labels": pd.read_csv(task_dir + "gsd_labels.csv", index_col=0),
        "spatial_kld": pd.read_csv(task_dir + "gsd_spatial_kld.csv", index_col=0),
        "timing": json.load(open(task_dir + "gsd_timing.json")),
    }


# --------------------------------------------------------------------------------------------------
# REV_GSD_base_merfish_imc
# --------------------------------------------------------------------------------------------------


def seurat_safe_feature_names(names):
    """
    Feature names that Seurat keeps unchanged (every character other than letters, digits and "-" becomes "-").

    Parameters
    ----------
    names : iterable of str
        Feature names (e.g. "Blank_1", "651779Pr141Di Cytoker").

    Returns
    -------
    list
        Sanitised names.
    """
    return [re.sub(r"[^A-Za-z0-9-]", "-", str(n)) for n in names]


def export_gsdensity_task_gsdbase(adata, gene_sets, task_dir, params, xy=None):
    """
    Write the inputs of run_gsdensity.R: log-normalised expression (genes x cells), coordinates and gene sets.

    Parameters
    ----------
    adata : AnnData
        Cells of one MCA embedding (log-normalised X).
    gene_sets : dict
        Gene set name -> genes.
    task_dir : str
        Task directory.
    params : dict
        GSDensity parameters.
    xy : np.ndarray, optional
        Coordinates for the spatial KLD (default obsm["spatial"]).
    """
    os.makedirs(task_dir, exist_ok=True)
    for f in (
        "gsd_scores.csv",
        "gsd_labels.csv",
        "gsd_spatial_kld.csv",
        "gsd_timing.json",
        "gsd_error.txt",
    ):
        if os.path.exists(task_dir + f):
            os.remove(task_dir + f)
    safe = dict(zip(adata.var_names, seurat_safe_feature_names(adata.var_names)))
    keys = {
        f"gs{i + 1}": [safe[g] for g in genes]
        for i, genes in enumerate(gene_sets.values())
    }
    mmwrite(task_dir + "data.mtx", csr_matrix(adata.X.T).astype(np.float64))
    with open(task_dir + "genes.txt", "w") as f:
        f.write("\n".join(safe[g] for g in adata.var_names) + "\n")
    with open(task_dir + "spots.txt", "w") as f:
        f.write("\n".join(map(str, adata.obs_names)) + "\n")
    xy = np.asarray(adata.obsm["spatial"] if xy is None else xy)[:, :2]
    pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1]}, index=adata.obs_names).to_csv(
        task_dir + "coords.csv"
    )
    with open(task_dir + "genesets.json", "w") as f:
        json.dump(keys, f)
    with open(task_dir + "params.json", "w") as f:
        json.dump(params, f)
    with open(task_dir + "geneset_names.json", "w") as f:
        json.dump(dict(zip(keys, gene_sets)), f)


def start_gsdensity_session(task_dirs, log, rscript):
    """
    Start one R session (gsdensity_env) that runs GSDensity on the task directories one after the other.

    Parameters
    ----------
    task_dirs : list of str
        Task directories.
    log : str
        Log file of the R session.
    rscript : str
        Path of scripts/run_gsdensity.R.

    Returns
    -------
    subprocess.Popen
        The running R session.
    """
    cmd = (
        GSDENSITY_ACTIVATE +
        f"Rscript {rscript} " + " ".join(task_dirs)
    )
    return subprocess.Popen(
        ["bash", "-c", cmd], stdout=open(log, "w"), stderr=subprocess.STDOUT
    )


def read_gsdensity_outputs_gsdbase(task_dir, cells):
    """
    GSDensity output of one task, with the gene set names and the cells in the given order.

    Parameters
    ----------
    task_dir : str
        Task directory.
    cells : array-like
        Cell names (order).

    Returns
    -------
    dict
        "scores" (relevance, cells x gene sets), "labels" (antimode label), "spatial_kld" and "timing".
    """
    names = json.load(open(task_dir + "geneset_names.json"))
    scores = pd.read_csv(task_dir + "gsd_scores.csv", index_col=0).rename(columns=names)
    labels = pd.read_csv(task_dir + "gsd_labels.csv", index_col=0).rename(columns=names)
    scores.index, labels.index = scores.index.astype(str), labels.index.astype(str)
    kld = None
    if os.path.exists(task_dir + "gsd_spatial_kld.csv"):
        kld = pd.read_csv(task_dir + "gsd_spatial_kld.csv", index_col=0).rename(
            index=names
        )
    return {
        "scores": scores.reindex(cells),
        "labels": labels.reindex(cells),
        "spatial_kld": kld,
        "timing": json.load(open(task_dir + "gsd_timing.json")),
    }


def nonconstant_features(adata):
    """
    Features with non-zero variance.

    Parameters
    ----------
    adata : AnnData
        Cells of one run.

    Returns
    -------
    set
        Names of the features that vary across the cells.
    """
    X = adata.X.toarray() if issparse(adata.X) else np.asarray(adata.X)
    return set(adata.var_names[X.var(axis=0) > 1e-12])


# --------------------------------------------------------------------------------------------------
# REV_SupFig1_simulation_benchmarks
# --------------------------------------------------------------------------------------------------


def supfig1_run_concordance(adata, score_dict, labels, n_neighbors=6, max_iter=10):
    """Neighbourhood enrichment concordance for pre-computed scores."""
    adata_tmp = adata.copy()
    sq.gr.spatial_neighbors(
        adata_tmp,
        n_neighs=n_neighbors,
        coord_type="generic",
        key_added="spatial",
    )
    W = adata_tmp.obsp["spatial_connectivities"]
    row_sums = np.array(W.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1.0
    W_norm = W.multiply(1.0 / row_sums[:, None])

    results = {}
    for method, scores in score_dict.items():
        results[method] = {}
        spatial_gt = labels.astype(np.float64).copy()
        for it in range(1, max_iter + 1):
            spatial_gt = W_norm.dot(spatial_gt)
            rho, _ = spearmanr(scores, spatial_gt)
            results[method][it] = rho if not np.isnan(rho) else 0.0
    return results


def supfig1_run_contamination(
    adata,
    sig_genes,
    scoring_fns,
    clean_scores,
    fractions,
    n_repeats=3,
    rng=None,
):
    """Signature contamination robustness."""
    if rng is None:
        rng = np.random.default_rng(0)

    common = [g for g in sig_genes if g in adata.var_names]
    non_sig = [g for g in adata.var_names if g not in common]
    n_sig = len(common)

    results = {m: {f: [] for f in fractions} for m in scoring_fns}

    for frac in fractions:
        n_contam = int(round(n_sig * frac))
        for _ in range(n_repeats):
            noise_genes = (
                list(
                    rng.choice(
                        non_sig,
                        n_contam,
                        replace=len(non_sig) < n_contam,
                    )
                )
                if n_contam > 0
                else []
            )
            contam_set = common + noise_genes

            for method, fn in scoring_fns.items():
                try:
                    cs = fn(adata.copy(), contam_set).flatten().astype(np.float64)
                except Exception:
                    cs = np.zeros(adata.n_obs)
                rho, _ = spearmanr(clean_scores[method], cs)
                results[method][frac].append(rho if not np.isnan(rho) else 0.0)

    return results


def supfig1_run_noise_robustness(
    adata,
    sig_genes,
    scoring_fns,
    clean_scores,
    noise_levels,
    n_repeats=3,
    rng=None,
):
    """Noise robustness benchmark."""
    if rng is None:
        rng = np.random.default_rng(0)

    common = [g for g in sig_genes if g in adata.var_names]
    gene_idx = [list(adata.var_names).index(g) for g in common]

    X_dense = adata.X.toarray() if issparse(adata.X) else adata.X
    gene_stds = np.array([X_dense[:, gi].std() for gi in gene_idx])

    results = {m: {lv: [] for lv in noise_levels} for m in scoring_fns}

    for level in noise_levels:
        for _ in range(n_repeats):
            ad_noisy = adata.copy()
            X_n = X_dense.copy()
            for j, gi in enumerate(gene_idx):
                X_n[:, gi] += rng.normal(
                    0,
                    gene_stds[j] * level,
                    X_n.shape[0],
                )
            ad_noisy.X = X_n
            for layer in ad_noisy.layers:
                ad_noisy.layers[layer] = X_n.copy()

            for method, fn in scoring_fns.items():
                try:
                    ns = fn(ad_noisy.copy(), common).flatten().astype(np.float64)
                except Exception:
                    ns = np.zeros(adata.n_obs)
                rho, _ = spearmanr(clean_scores[method], ns)
                results[method][level].append(
                    rho if not np.isnan(rho) else 0.0,
                )

    return results


def supfig1_run_gene_weight_validation(adata, sig_genes, fold_changes):
    """Correlate inferred weights with true fold changes."""
    if "gene_contributions" not in adata.uns:
        return np.nan

    score_keys = list(adata.uns["gene_contributions"].keys())
    if not score_keys:
        return np.nan

    contributions = adata.uns["gene_contributions"][score_keys[0]]
    common = [g for g in sig_genes if g in contributions]
    if len(common) < 3:
        return np.nan

    weights = [np.mean(np.abs(contributions[g])) for g in common]
    fc_matched = [fold_changes[sig_genes.index(g)] for g in common]

    rho, _ = spearmanr(weights, fc_matched)
    return rho if not np.isnan(rho) else 0.0


def supfig1_run_classification(adata, score_dict, labels, n_neighbors=6):
    """
    Classification benchmark: AUROC, AUPRC, F1 (optimal threshold),
    spatial coherence, boundary F1, KS statistic.
    """
    from sklearn.metrics import (
        roc_auc_score,
        average_precision_score,
        precision_recall_curve,
        f1_score,
    )
    from scipy.stats import ks_2samp

    # Spatial graph
    adata_tmp = adata.copy()
    sq.gr.spatial_neighbors(
        adata_tmp,
        n_neighs=n_neighbors,
        coord_type="generic",
        key_added="spatial",
    )
    W = adata_tmp.obsp["spatial_connectivities"]

    # Ground truth Moran's I for spatial coherence normalisation
    y = labels.astype(np.float64)
    n = len(y)
    z = y - y.mean()
    z2 = z @ z
    gt_morans = float((n / W.sum()) * (z @ W.dot(z) / z2)) if z2 > 0 else 0.0

    # Boundary mask
    boundary_mask = np.zeros(n, dtype=bool)
    for i in range(n):
        neighbours = W[i].nonzero()[1]
        if len(neighbours) > 0 and np.any(labels[neighbours] != labels[i]):
            boundary_mask[i] = True

    results = {}
    for method, scores in score_dict.items():
        s = scores.astype(np.float64)
        valid = ~np.isnan(s)
        sv, yv = s[valid], labels[valid]

        if yv.sum() == 0 or yv.sum() == len(yv):
            continue

        rec = {}

        # AUROC
        try:
            rec["AUROC"] = roc_auc_score(yv, sv)
        except ValueError:
            rec["AUROC"] = np.nan

        # AUPRC
        try:
            rec["AUPRC"] = average_precision_score(yv, sv)
        except ValueError:
            rec["AUPRC"] = np.nan

        # Optimal F1
        prec, recall, thresholds = precision_recall_curve(yv, sv)
        prec, recall = prec[:-1], recall[:-1]
        with np.errstate(divide="ignore", invalid="ignore"):
            f1_vals = 2 * (prec * recall) / (prec + recall)
        f1_vals = np.nan_to_num(f1_vals, nan=0.0)
        best_idx = np.argmax(f1_vals)
        rec["F1"] = float(f1_vals[best_idx])
        best_thr = float(thresholds[best_idx])

        # KS statistic
        in_domain = sv[yv == 1]
        out_domain = sv[yv == 0]
        rec["KS"] = ks_2samp(in_domain, out_domain).statistic

        # Spatial coherence (predictions at optimal threshold)
        q = float(np.clip(np.mean(sv <= best_thr), 0.01, 0.99))
        preds = (s > np.quantile(s[valid], q)).astype(float)
        zp = preds - preds.mean()
        zp2 = zp @ zp
        if zp2 > 0 and gt_morans > 0:
            pred_morans = float((n / W.sum()) * (zp @ W.dot(zp) / zp2))
            rec["Spatial coherence"] = float(np.clip(pred_morans / gt_morans, 0, 1))
        else:
            rec["Spatial coherence"] = 0.0

        # Boundary F1
        if boundary_mask.sum() > 0:
            preds_int = preds.astype(int)
            yb = labels[boundary_mask]
            pb = preds_int[boundary_mask]
            if yb.sum() > 0 and yb.sum() < len(yb):
                rec["Boundary F1"] = float(f1_score(yb, pb, zero_division=0))
            else:
                rec["Boundary F1"] = np.nan
        else:
            rec["Boundary F1"] = np.nan

        results[method] = rec

    return results


def supfig1_f1_of_calls(y, pred):
    """
    F1 of binary calls.

    Parameters
    ----------
    y : np.ndarray
        Binary labels.
    pred : np.ndarray
        Binary calls.

    Returns
    -------
    float
        2 TP / (calls + positives).
    """
    y, pred = np.asarray(y).astype(bool), np.asarray(pred).astype(bool)
    return 2 * (y & pred).sum() / (pred.sum() + y.sum())


def supfig1_otsu_threshold(s, nbins=256):
    """
    Otsu threshold of a score (label-free: maximises the between-class variance of the histogram).

    Parameters
    ----------
    s : np.ndarray
        Score per spot.
    nbins : int
        Histogram bins.

    Returns
    -------
    float
        Threshold.
    """
    hist, edges = np.histogram(s, bins=nbins)
    mids = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * mids) / np.maximum(w0, 1)
    m1 = (np.sum(hist * mids) - np.cumsum(hist * mids)) / np.maximum(w1, 1)
    return mids[np.argmax(w0 * w1 * (m0 - m1) ** 2)]


def supfig1_spatial_folds(coords, n_blocks=25, n_folds=5, seed=0):
    """
    Spatial-block cross-validation folds: k-means blocks of the coordinates, assigned to folds at random.

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates.
    n_blocks, n_folds : int
        Number of spatial blocks and of folds.
    seed : int
        Random seed.

    Returns
    -------
    np.ndarray
        Fold per spot.
    """
    km = KMeans(n_clusters=n_blocks, n_init=10, random_state=seed).fit(coords)
    block_fold = np.random.default_rng(seed).permutation(np.arange(n_blocks) % n_folds)
    return block_fold[km.labels_]


def supfig1_gearys_c(values, W):
    """
    Geary's C with a binary spatial weight matrix (1 = no autocorrelation, < 1 = positive autocorrelation).

    Parameters
    ----------
    values : np.ndarray
        Score per spot.
    W : scipy.sparse matrix
        6-NN spatial connectivities.

    Returns
    -------
    float
        Geary's C.
    """
    v = np.asarray(values, dtype=np.float64)
    z = v - v.mean()
    W = csr_matrix(W)
    rows = np.asarray(W.sum(axis=1)).ravel()
    cols = np.asarray(W.sum(axis=0)).ravel()
    sq_diff = rows @ v**2 + cols @ v**2 - 2 * v @ W.dot(v)
    return (len(v) - 1) * sq_diff / (2 * W.sum() * (z @ z))


def supfig1_threshold_rules(y, s, W, folds):
    """
    F1 of one score under the five threshold rules of R4_D1_F1_audit, and Moran's I and Geary's C of the score.

    Parameters
    ----------
    y : np.ndarray
        Binary target-domain labels.
    s : np.ndarray
        Score per spot.
    W : scipy.sparse matrix
        6-NN spatial connectivities.
    folds : np.ndarray
        Spatial-block CV fold per spot.

    Returns
    -------
    dict
        Metric -> value.
    """
    f1_oracle, _ = _optimal_f1(y, s)
    # spatial-block CV: F1-optimal threshold learnt on the other folds, applied to the held-out fold
    pred_cv = np.zeros(len(y), bool)
    for fold in np.unique(folds):
        test = folds == fold
        if y[~test].sum() > 0:
            _, thr = _optimal_f1(y[~test], s[~test])
            pred_cv[test] = s[test] >= thr
    # top-k: the k highest scores called positive, k = number of target spots
    pred_topk = np.zeros(len(y), bool)
    pred_topk[np.argsort(-s, kind="stable")[: y.sum()]] = True
    return {
        "F1 (oracle threshold)": f1_oracle,
        "F1 (spatial-block CV threshold)": supfig1_f1_of_calls(y, pred_cv),
        "F1 (top-k, prevalence-matched)": supfig1_f1_of_calls(y, pred_topk),
        "F1 (Otsu, label-free)": supfig1_f1_of_calls(y, s > supfig1_otsu_threshold(s)),
        "F1 (median threshold)": supfig1_f1_of_calls(y, s > np.median(s)),
        "Moran's I": _morans_i_continuous(s, W),
        "Geary's C": supfig1_gearys_c(s, W),
    }


def supfig1_call_metrics(y, pred, W, gt_morans):
    """
    Binary calls against the target domain: F1, precision, recall, fraction called and spatial coherence.

    Parameters
    ----------
    y : np.ndarray
        Binary target-domain labels.
    pred : np.ndarray
        Boolean calls.
    W : scipy.sparse matrix
        6-NN spatial connectivities.
    gt_morans : float
        Moran's I of the labels.

    Returns
    -------
    dict
        Metric -> value.
    """
    tp = (pred & (y == 1)).sum()
    return {
        "F1": supfig1_f1_of_calls(y, pred),
        "precision": tp / max(pred.sum(), 1),
        "recall": tp / y.sum(),
        "fraction positive": pred.mean(),
        "fraction of background spots positive": pred[y == 0].mean(),
        "Spatial coherence": _spatial_coherence(pred, W, gt_morans),
    }


def supfig1_run_simulation(
    seed,
    scoring_functions,
    gsd_dir,
    concordance_iters,
    contamination_fracs,
    noise_levels,
    contam_repeats=3,
    noise_repeats=3,
):
    """
    One simulated slide: scores of every method, the five 05_01 benchmarks, the threshold rules and the fixed,
    label-free calls of the null model and GSDensity.

    Parameters
    ----------
    seed : int
        Simulation seed (= iteration of 05_01).
    scoring_functions : dict
        Method -> scorer(adata, genes) of the Python methods (GSDensity is read from the cache).
    gsd_dir : str
        GSDensity cache of the simulated slides (SupFig1_seed{seed:03d}_cluster_2.csv).
    concordance_iters : list of int
        Smoothing iterations reported for the neighbourhood concordance.
    contamination_fracs : list of float
        Random genes added, as a multiple of the signature size.
    noise_levels : list of float
        Gaussian noise added to the signature genes, as a multiple of each gene's SD.
    contam_repeats, noise_repeats : int
        Random draws per contamination fraction and per noise level.

    Returns
    -------
    dict
        Long records of every benchmark, and the per-spot scores of the slide.
    """
    rng = np.random.default_rng(seed + 1000)

    # simulate
    adata, sig_genes, fold_changes, target_label = simulate_visium(seed=seed)
    labels = (adata.obs["cluster"] == target_label).astype(int).values

    # score with EnrichMap first (gene contributions and the null-model p-values)
    em.tl.score(
        adata,
        gene_set=sig_genes,
        score_key="_em",
        smoothing=True,
        null_model=True,
        correct_spatial_covariates=False,
        compute_q_vals=True,
    )
    null = adata.uns["null_model"]["_em"]

    # score with all methods; GSDensity relevance (x spots) from the cache
    clean_scores = {}
    for method, fn in scoring_functions.items():
        clean_scores[method] = fn(adata.copy(), sig_genes).flatten().astype(np.float64)
    gsd = pd.read_csv(
        gsd_dir + f"SupFig1_seed{seed:03d}_cluster_2.csv", index_col=0
    ).reindex(adata.obs_names)
    score_dict = {**clean_scores, "GSDensity": gsd["relevance_scaled"].values}
    gsd_positive = (gsd["label"] == "positive").values

    # benchmarks 1-5 (05_01)
    conc = supfig1_run_concordance(
        adata, score_dict, labels, max_iter=max(concordance_iters)
    )
    cont = supfig1_run_contamination(
        adata,
        sig_genes,
        scoring_functions,
        clean_scores,
        contamination_fracs,
        n_repeats=contam_repeats,
        rng=rng,
    )
    noise = supfig1_run_noise_robustness(
        adata,
        sig_genes,
        scoring_functions,
        clean_scores,
        noise_levels,
        n_repeats=noise_repeats,
        rng=rng,
    )
    gw_rho = supfig1_run_gene_weight_validation(adata, sig_genes, fold_changes)
    clf = supfig1_run_classification(adata, score_dict, labels)

    # gene weights (panel d): CV weight, gene SD and mean |contribution| of each signature gene
    weights = em.tl.infer_gene_weights(adata, sig_genes)
    contributions = adata.uns["gene_contributions"]["_em"]
    X = adata[:, sig_genes].X.toarray()
    gw = pd.DataFrame(
        {
            "seed": seed,
            "gene": sig_genes,
            "true_fold_change": fold_changes,
            "true_log2_fold_change": np.log2(fold_changes),
            "cv_weight": [weights[g] for g in sig_genes],
            "gene_sd": X.std(axis=0),
            "baseline_mean": X[labels == 0].mean(axis=0),
            "mean_abs_contribution": [
                np.mean(np.abs(contributions[g])) for g in sig_genes
            ],
        }
    )

    # threshold rules, Moran's I and Geary's C of each score (6-NN graph, 25 spatial blocks in 5 folds)
    coords = adata.obsm["spatial"]
    sq.gr.spatial_neighbors(
        adata, n_neighs=6, coord_type="generic", key_added="spatial"
    )
    W = adata.obsp["spatial_connectivities"]
    folds = supfig1_spatial_folds(coords, n_blocks=25, n_folds=5, seed=0)
    rules = {
        m: supfig1_threshold_rules(labels, s, W, folds) for m, s in score_dict.items()
    }
    gt_morans = _morans_i_continuous(labels, W)

    # fixed, label-free calls: null-model z and empirical p-value, GSDensity's antimode labels
    z = clean_scores["EnrichMap (null)"]
    calls = {
        ("EnrichMap (null)", "z > 0"): z > 0,
        ("EnrichMap (null)", "z > 1.96"): z > 1.96,
        ("EnrichMap (null)", "p_empirical <= 0.05"): null["p_empirical"] <= 0.05,
        ("GSDensity", "GSDensity antimode labels"): gsd_positive,
    }
    for (m, rule), pred in calls.items():
        rules[m].update(
            {
                f"{k} ({rule})": v
                for k, v in supfig1_call_metrics(labels, pred, W, gt_morans).items()
            }
        )

    return {
        "seed": seed,
        "n_spots": adata.n_obs,
        "prevalence": labels.mean(),
        "conc": [
            {"seed": seed, "method": m, "smoothing_iter": it, "value": conc[m][it]}
            for m in conc
            for it in concordance_iters
        ],
        "contam": [
            {
                "seed": seed,
                "method": m,
                "contamination": frac,
                "repeat": i,
                "value": rho,
            }
            for m in cont
            for frac, rhos in cont[m].items()
            for i, rho in enumerate(rhos)
        ],
        "noise": [
            {"seed": seed, "method": m, "noise_level": level, "repeat": i, "value": rho}
            for m in noise
            for level, rhos in noise[m].items()
            for i, rho in enumerate(rhos)
        ],
        "gw": gw,
        "gw_rho": gw_rho,
        "clf": [{"seed": seed, "method": m, **met} for m, met in clf.items()],
        "rules": [
            {"seed": seed, "method": m, "metric": k, "value": v}
            for m, met in rules.items()
            for k, v in met.items()
        ],
        "truth": {
            "seed": seed,
            "Moran's I": gt_morans,
            "Geary's C": supfig1_gearys_c(labels, W),
        },
        "null": {
            "seed": seed,
            "median null_sd": np.median(null["null_sd"]),
            "fraction p_empirical <= 0.05 (background)": np.mean(
                null["p_empirical"][labels == 0] <= 0.05
            ),
            "fraction z > 1.96 (background)": np.mean(z[labels == 0] > 1.96),
            "fraction z > 1.96 (target)": np.mean(z[labels == 1] > 1.96),
        },
        "maps": pd.DataFrame(
            {
                "x": coords[:, 0],
                "y": coords[:, 1],
                "target domain": labels,
                **score_dict,
                "EnrichMap (null) p_empirical": null["p_empirical"],
                "GSDensity antimode label": gsd_positive,
            },
            index=adata.obs_names,
        ),
    }


def supfig1_line_summary(df, x):
    """
    Mean ± SD of each method at each x, as drawn.

    Parameters
    ----------
    df : pd.DataFrame
        Long table with seed, method, x and value.
    x : str
        Column on the x axis.

    Returns
    -------
    pd.DataFrame
        One row per method and x value.
    """
    agg = (
        df.groupby(["method", x], observed=True)
        .agg(
            mean=("value", "mean"),
            SD=("value", "std"),
            n_values=("value", "size"),
            n_simulations=("seed", "nunique"),
        )
        .reset_index()
    )
    agg["band_lower (mean - SD)"] = agg["mean"] - agg["SD"]
    agg["band_upper (mean + SD)"] = agg["mean"] + agg["SD"]
    agg["method"] = agg["method"].astype(str)
    return agg


def supfig1_per_slide_area(df, x):
    """
    Normalised area under a robustness curve per slide: mean rho over the levels > 0 (draws averaged first).

    Parameters
    ----------
    df : pd.DataFrame
        Long table with seed, method, x and value.
    x : str
        Perturbation level column.

    Returns
    -------
    pd.Series
        Indexed by (method, seed).
    """
    d = df[df[x] > 0].groupby(["method", "seed", x], observed=True)["value"].mean()
    return d.groupby(["method", "seed"], observed=True).mean()


def supfig1_boot_ci(values, n_boot=2000, seed=0):
    """
    Percentile bootstrap 95% CI of the mean over slides.

    Parameters
    ----------
    values : np.ndarray
        One value per slide.
    n_boot : int
        Bootstrap resamples.
    seed : int
        Random seed.

    Returns
    -------
    tuple
        (lower, upper).
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    means = v[np.random.default_rng(seed).integers(0, len(v), (n_boot, len(v)))].mean(
        axis=1
    )
    return tuple(np.percentile(means, [2.5, 97.5]))


def supfig1_metric_table(per_slide, metrics, methods):
    """
    Mean and 95% bootstrap CI over slides of each metric, one row per method.

    Parameters
    ----------
    per_slide : pd.DataFrame
        Columns method, seed and the metrics.
    metrics : list of str
        Metric columns.
    methods : list of str
        Method order of the output.

    Returns
    -------
    pd.DataFrame
        method, metric, metric_lo, metric_hi and the number of slides for each metric.
    """
    summary_rows = []
    for m, g in per_slide.groupby("method", observed=True, sort=False):
        row = {"method": str(m)}
        for k in metrics:
            v = g[k].dropna().values
            row[k] = v.mean()
            row[f"{k}_lo"], row[f"{k}_hi"] = supfig1_boot_ci(v)
            row[f"{k} (n slides)"] = len(v)
        summary_rows.append(row)
    out = pd.DataFrame(summary_rows)
    order = [m for m in methods + ["Truth (target domain)"] if m in set(out["method"])]
    return out.set_index("method").loc[order].reset_index()


def supfig1_method_legend(ax, colours, n_note=None, kind="line", extra=None):
    """
    Legend of the scoring methods to the right of an axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes that carries the legend.
    colours : dict
        Method -> colour, in legend order.
    n_note : dict or None
        Method -> number of slides, appended to the label when it differs from 100.
    kind : str
        "line" (line and marker) or "dot" (marker only).
    extra : list or None
        Further legend handles.
    """
    n_note = n_note or {}
    methods = list(colours)
    handles = []
    for m in methods:
        label = m if n_note.get(m, 100) == 100 else f"{m} (n = {n_note[m]})"
        handles.append(
            Line2D(
                [0],
                [0],
                color=colours[m],
                lw=1 if kind == "line" else 0,
                ls="--" if m == "GSDensity" and kind == "line" else "-",
                marker="o",
                ms=3,
                label=label,
            )
        )
    handles += extra or []
    ax.legend(
        handles=handles,
        title="Scoring method",
        frameon=False,
        fontsize=6,
        title_fontsize=7,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
    )


def supfig1_line_panel(
    summary, x, xticks, xlabel, ylabel, title, colours, bands=True, figsize=(2.4, 2.2)
):
    """
    Mean (± SD) lines of every method (Sup. Fig. 1a-c layout).

    Parameters
    ----------
    summary : pd.DataFrame
        Output of line_summary.
    x : str
        x variable.
    xticks : list
        x ticks.
    xlabel, ylabel, title : str
        Labels.
    colours : dict
        Method -> colour, in plotting order.
    bands : bool
        Draw mean ± SD bands (full version).
    figsize : tuple
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure
        The figure.
    """
    methods = list(colours)
    fig, ax = plt.subplots(figsize=figsize)
    for m in methods:
        d = summary[summary["method"] == m].sort_values(x)
        # null model drawn wider underneath, so it stays visible where it overlaps the smoothing-only line
        null = m == "EnrichMap (null)"
        ax.plot(
            d[x],
            d["mean"],
            "--" if m == "GSDensity" else "-",
            color=colours[m],
            lw=2.4 if null else 1,
            marker="o",
            ms=3.5 if null else 2,
            zorder=2 if null else 3,
        )
        if bands:
            ax.fill_between(
                d[x],
                d["band_lower (mean - SD)"],
                d["band_upper (mean + SD)"],
                color=colours[m],
                alpha=0.12,
                lw=0,
            )
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{t:g}" for t in xticks])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=8)
    ax.grid(False)
    sns.despine(ax=ax)
    supfig1_method_legend(
        ax, colours, summary.groupby("method")["n_simulations"].max().to_dict()
    )
    return fig


def supfig1_metric_pair(
    ax, table, x, y, colours, ref_x=None, ref_y=None, invert_y=False
):
    """
    One metric against another, one point per method, with 95% CI whiskers (f"{x}_lo", f"{x}_hi").

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    table : pd.DataFrame
        One row per method, with a "method" column and the metric columns.
    x, y : str
        Metric columns.
    colours : dict
        Method -> colour.
    ref_x, ref_y : float or None
        Dotted reference lines.
    invert_y : bool
        Invert the y axis for a lower-is-better metric (Geary's C), so that better is up and right.
    """
    for _, r in table.iterrows():
        xerr = (
            [[r[x] - r[f"{x}_lo"]], [r[f"{x}_hi"] - r[x]]]
            if f"{x}_lo" in table
            else None
        )
        yerr = (
            [[r[y] - r[f"{y}_lo"]], [r[f"{y}_hi"] - r[y]]]
            if f"{y}_lo" in table
            else None
        )
        ax.errorbar(
            r[x],
            r[y],
            xerr=xerr,
            yerr=yerr,
            fmt="o",
            color=colours[r["method"]],
            ms=5,
            mec="white",
            mew=0.4,
            elinewidth=0.6,
            capsize=0,
            label=r["method"],
        )
    if ref_x is not None:
        ax.axvline(ref_x, ls=":", c="0.6", lw=0.6)
    if ref_y is not None:
        ax.axhline(ref_y, ls=":", c="0.6", lw=0.6)
    if invert_y:
        ax.invert_yaxis()
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.grid(False)


def supfig1_benchmark_summary(ax, table, colours, lower_is_better=()):
    """
    Methods x metrics: circles sized and shaded by each metric scaled to 0-1 across methods (best = 1), the raw
    value in each circle and the mean rank in the last column.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    table : pd.DataFrame
        Index = methods, columns = metrics.
    colours : dict
        Method -> colour of the method labels.
    lower_is_better : tuple of str
        Metrics where smaller is better (e.g. "Geary's C").
    """
    scaled = ((table - table.min()) / (table.max() - table.min())).fillna(1.0)
    for c in lower_is_better:
        scaled[c] = 1 - scaled[c]
    ranks = scaled.rank(ascending=False).mean(axis=1)
    cmap = plt.get_cmap("Purples")
    for i, m in enumerate(table.index):
        for j, c in enumerate(table.columns):
            s = scaled.loc[m, c]
            ax.scatter(
                j, i, s=30 + 220 * s, color=cmap(0.25 + 0.75 * s), edgecolors="none"
            )
            ax.text(
                j,
                i,
                f"{table.loc[m, c]:.2g}",
                ha="center",
                va="center",
                fontsize=5,
                color="white" if s > 0.6 else "black",
            )
        ax.text(
            len(table.columns) + 0.2,
            i,
            f"{ranks[m]:.1f}",
            ha="left",
            va="center",
            fontsize=6,
        )
    ax.set_xticks(range(len(table.columns)))
    ax.set_xticklabels(table.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(table.index)))
    ax.set_yticklabels(table.index)
    for label in ax.get_yticklabels():
        label.set_color(colours[label.get_text()])
    ax.text(len(table.columns) + 0.2, -0.9, "mean rank", fontsize=6, ha="left")
    ax.set_xlim(-0.6, len(table.columns) + 0.9)
    ax.set_ylim(len(table.index) - 0.4, -1.2)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)


def supfig1_mean_at(summary, x, value):
    """
    Mean of each method at one x value.

    Parameters
    ----------
    summary : pd.DataFrame
        Output of line_summary.
    x : str
        x variable.
    value : float
        x value.

    Returns
    -------
    pd.Series
        Mean per method.
    """
    return summary[summary[x] == value].set_index("method")["mean"]


# --------------------------------------------------------------------------------------------------
# REV_SupFig2_brain_spatial_fidelity
# --------------------------------------------------------------------------------------------------


def ring_distance_to_region(conn, mask, max_hops=12):
    """
    Graph distance (number of neighbour rings) of every spot to the nearest spot of a region.

    Parameters
    ----------
    conn : scipy.sparse matrix
        Spatial connectivity (made symmetric here).
    mask : np.ndarray
        Boolean region membership per spot.
    max_hops : int
        Largest ring computed; spots further away stay at inf.

    Returns
    -------
    np.ndarray
        Ring index per spot (0 inside the region).
    """
    mask = np.asarray(mask, bool)
    A = csr_matrix((conn + conn.T) > 0).astype(np.int8)
    dist = np.full(len(mask), np.inf)
    dist[mask] = 0
    frontier = mask.copy()
    seen = mask.copy()
    for h in range(1, max_hops + 1):
        nxt = (A.dot(frontier.astype(np.int8)) > 0) & ~seen
        dist[nxt] = h
        seen |= nxt
        frontier = nxt
        if not nxt.any():
            break
    return dist


def ring_leakage_profile(score, mask, dist, max_ring=6, far_from=6):
    """
    Standardised score by ring outside a region:
    L(d) = (mean score at ring d - mean far background) / (mean in region - mean far background).

    Parameters
    ----------
    score : np.ndarray
        Score per spot.
    mask : np.ndarray
        Boolean region membership per spot.
    dist : np.ndarray
        Ring index per spot (ring_distance_to_region).
    max_ring : int
        Last ring reported.
    far_from : int
        Spots at least this many rings away form the far background.

    Returns
    -------
    dict
        Ring -> L(d); 0 = as low as the far background, 1 = as high as inside, < 0 = negative halo.
    """
    score = np.asarray(score, dtype=np.float64)
    inside = score[np.asarray(mask, bool)].mean()
    far = score[dist >= far_from].mean()
    span = inside - far
    prof = {}
    for d in range(1, max_ring + 1):
        sel = dist == d
        prof[d] = (
            float((score[sel].mean() - far) / span)
            if sel.any() and span != 0
            else np.nan
        )
    return prof


def f1_median_threshold(y, s):
    """
    F1 when the score is thresholded at its median (the rule of the manuscript Methods).

    Parameters
    ----------
    y : np.ndarray
        Binary labels.
    s : np.ndarray
        Score per spot.

    Returns
    -------
    float
    """
    return f1_of_calls(y, s > np.median(s))


def f1_prevalence_top_k(y, s):
    """
    F1 when the k highest-scoring spots are called positive, k = number of positives (precision = recall = F1).

    Parameters
    ----------
    y : np.ndarray
        Binary labels.
    s : np.ndarray
        Score per spot.

    Returns
    -------
    float
    """
    y = np.asarray(y).astype(bool)
    k = int(y.sum())
    order = np.argsort(-np.asarray(s), kind="stable")
    pred = np.zeros(len(y), bool)
    pred[order[:k]] = True
    return f1_of_calls(y, pred)


def otsu_threshold_1d(s, nbins=256):
    """
    Otsu threshold of a score (maximal between-class variance of a 256-bin histogram).

    Parameters
    ----------
    s : np.ndarray
        Score per spot.
    nbins : int
        Histogram bins.

    Returns
    -------
    float
        Threshold (bin centre).
    """
    s = np.asarray(s, dtype=np.float64)
    hist, edges = np.histogram(s, bins=nbins)
    mids = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * mids) / np.maximum(w0, 1)
    m1 = (np.sum(hist * mids) - np.cumsum(hist * mids)) / np.maximum(w1, 1)
    between = w0 * w1 * (m0 - m1) ** 2
    return float(mids[int(np.argmax(between))])


def f1_otsu_threshold(y, s):
    """
    F1 at the label-free Otsu threshold.

    Parameters
    ----------
    y : np.ndarray
        Binary labels.
    s : np.ndarray
        Score per spot.

    Returns
    -------
    float
    """
    return f1_of_calls(y, np.asarray(s) > otsu_threshold_1d(s))


def spatial_block_partition(coords, n_blocks=25, n_folds=5, seed=0):
    """
    Spatial blocks (k-means of the coordinates) assigned at random to folds.

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates.
    n_blocks : int
        Number of k-means blocks.
    n_folds : int
        Number of folds (n_folds = n_blocks gives one label per block).
    seed : int
        Seed of k-means and of the block-to-fold assignment.

    Returns
    -------
    np.ndarray
        Fold (or block) label per spot.
    """
    km = KMeans(n_clusters=n_blocks, n_init=10, random_state=seed).fit(
        np.asarray(coords)
    )
    rng = np.random.default_rng(seed)
    block_fold = rng.permutation(np.arange(n_blocks) % n_folds)
    return block_fold[km.labels_]


def f1_spatial_block_cv(y, s, folds):
    """
    F1 with the F1-optimal threshold learnt on the training blocks and applied to the held-out blocks
    (calls pooled over folds).

    Parameters
    ----------
    y : np.ndarray
        Binary labels.
    s : np.ndarray
        Score per spot.
    folds : np.ndarray
        Fold label per spot (spatial_block_partition).

    Returns
    -------
    float
    """
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=np.float64)
    pred = np.zeros(len(y), bool)
    for f in np.unique(folds):
        test = folds == f
        train = ~test
        if y[train].sum() == 0:
            continue
        _, thr = _optimal_f1(y[train], s[train])
        pred[test] = s[test] >= thr
    return f1_of_calls(y, pred)


def threshold_rule_metrics(y, s, W, folds):
    """
    AUROC, AUPRC and F1 under the threshold rules of R4_D1_F1_audit, plus the spatial coherence of the oracle calls
    and Moran's I of the score.

    Parameters
    ----------
    y : np.ndarray
        Binary labels.
    s : np.ndarray
        Score per spot.
    W : scipy.sparse matrix
        Binary spatial connectivity (spatial_knn_connectivity).
    folds : np.ndarray
        Fold label per spot for the spatial-block CV rule.

    Returns
    -------
    dict
        Metric name -> value.
    """
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=np.float64)
    f1o, thr = _optimal_f1(y, s)
    gt = _morans_i_continuous(y.astype(float), W)
    return {
        "AUROC": roc_auc_score(y, s),
        "AUPRC": average_precision_score(y, s),
        "F1 (oracle threshold)": f1o,
        "F1 (median threshold)": f1_median_threshold(y, s),
        "F1 (top-k, prevalence-matched)": f1_prevalence_top_k(y, s),
        "F1 (Otsu, label-free)": f1_otsu_threshold(y, s),
        "F1 (spatial-block CV threshold)": f1_spatial_block_cv(y, s, folds),
        "Spatial coherence (oracle thr)": _spatial_coherence(s >= thr, W, gt),
        "Moran's I": _morans_i_continuous(s, W),
    }


def null_model_rule_calls(z, p):
    """
    Binary calls of the null model under its fixed rules.

    Parameters
    ----------
    z : np.ndarray
        Null-model score (z vs matched random gene sets).
    p : np.ndarray
        Empirical p-value per spot.

    Returns
    -------
    dict
        Rule name -> boolean calls.
    """
    return {"z > 0": z > 0, "z > 1.96": z > 1.96, "p_empirical <= 0.05": p <= 0.05}


def paired_block_bootstrap_rules(
    y,
    score_dict,
    p_null,
    blocks,
    partitions,
    metric_fns,
    null_method="EnrichMap (null)",
    n_boot=1000,
    seed=7,
    min_pos=3,
    W=None,
):
    """
    Paired spatial block bootstrap of R4_D1_F1_audit: every method is evaluated on the same resampled blocks, with
    every threshold re-derived per replicate. Adds the spatial-block CV rule (partition r mod len(partitions) in
    replicate r) and the fixed rules of the null model.

    Parameters
    ----------
    y : np.ndarray
        Binary region labels.
    score_dict : dict
        Method -> score per spot.
    p_null : np.ndarray
        Empirical p-value of the null model per spot.
    blocks : np.ndarray
        Bootstrap block per spot.
    partitions : list of np.ndarray
        CV fold labels per spot, one array per partition.
    metric_fns : dict
        Metric name -> fn(y, s).
    null_method : str
        Key of the null model in `score_dict`.
    n_boot : int
        Number of replicates.
    seed : int
        Seed of the block resampling.
    min_pos : int
        Replicates with fewer region spots (or no background spot) are skipped.
    W : scipy.sparse matrix or None
        Spatial connectivity. If given, adds the spatial coherence of the calls at the F1-optimal threshold, on the
        graph restricted to the resampled spots.

    Returns
    -------
    pd.DataFrame
        replicate, method, metric, value.
    """
    rng = np.random.default_rng(seed)
    y = np.asarray(y).astype(int)
    ub = np.unique(blocks)
    idx_by_block = {b: np.flatnonzero(blocks == b) for b in ub}
    records = []
    for r in range(n_boot):
        pick = rng.choice(ub, size=len(ub), replace=True)
        idx = np.concatenate([idx_by_block[b] for b in pick])
        yy = y[idx]
        if yy.sum() < min_pos or yy.sum() == len(yy):
            continue
        fold = partitions[r % len(partitions)][idx]
        if W is not None:
            W_boot = W[idx][:, idx]
            gt_morans = _morans_i_continuous(yy.astype(float), W_boot)
        for m, s in score_dict.items():
            ss = s[idx]
            for name, fn in metric_fns.items():
                records.append((r, m, name, float(fn(yy, ss))))
            records.append(
                (
                    r,
                    m,
                    "F1 (spatial-block CV threshold)",
                    f1_spatial_block_cv(yy, ss, fold),
                )
            )
            if W is not None:
                calls = ss >= _optimal_f1(yy, ss)[1]
                records.append(
                    (
                        r,
                        m,
                        "Spatial coherence (oracle thr)",
                        _spatial_coherence(calls, W_boot, gt_morans),
                    )
                )
        for rule, calls in null_model_rule_calls(
            score_dict[null_method][idx], p_null[idx]
        ).items():
            records.append((r, null_method, f"F1 ({rule})", f1_of_calls(yy, calls)))
    return pd.DataFrame(records, columns=["replicate", "method", "metric", "value"])


def percentile_interval(x, level=0.95):
    """
    Percentile interval of bootstrap replicates.

    Parameters
    ----------
    x : array-like
        Replicate values.
    level : float
        Coverage.

    Returns
    -------
    tuple of float
        (lower, upper).
    """
    a = (1 - level) / 2
    return float(np.quantile(x, a)), float(np.quantile(x, 1 - a))


def median_threshold_f1_ceiling(prevalence):
    """
    Largest F1 when exactly half of the spots are called positive: 2 min(p, 1/2) / (1/2 + p).

    Parameters
    ----------
    prevalence : float or np.ndarray
        Share of positive spots.

    Returns
    -------
    float or np.ndarray
    """
    p = np.asarray(prevalence, dtype=np.float64)
    return 2 * np.minimum(p, 0.5) / (0.5 + p)


def format_estimate_ci(r):
    """
    Format an estimate with its bootstrap interval.

    Parameters
    ----------
    r : pd.Series
        Row with `value`, `lo` and `hi`.

    Returns
    -------
    str
        "value [lo, hi]", or "value" when no interval exists.
    """
    return (
        f"{r.value:.3f} [{r.lo:.3f}, {r.hi:.3f}]"
        if pd.notna(r.lo)
        else f"{r.value:.3f}"
    )


def moran_scatter_quadrants(score, coords, n_neighbours=6, seed=0):
    """
    Score, spatial lag, quadrant and Moran's I as computed by em.pl.morans_correlogram; quadrants are taken relative
    to the score mean (the centring of Moran's I).

    Parameters
    ----------
    score : np.ndarray
        Score per spot.
    coords : np.ndarray
        Spot coordinates.
    n_neighbours : int
        k of the row-standardised KNN weights.
    seed : int
        Seed of the 999 permutations of the pseudo p-value.

    Returns
    -------
    tuple
        (DataFrame with score, spatial lag and quadrant; esda Moran object).
    """
    np.random.seed(seed)
    W = KNN.from_array(coords, k=n_neighbours)
    W.transform = "r"
    moran = Moran(score, W)
    score_lag = lag_spatial(W, score)
    centre = score.mean()
    high, high_lag = score >= centre, score_lag >= centre
    quadrant = np.where(
        high & high_lag,
        "Co-enriched",
        np.where(~high & ~high_lag, "Co-depleted", "Discordant"),
    )
    df = pd.DataFrame({"score": score, "spatial lag": score_lag, "quadrant": quadrant})
    return df, moran


def plot_spatial_metric_bars(spatial_stats, metric, ylabel, methods, regions, colours):
    """
    Bar chart of one spatial statistic per method, one axes per region.

    Parameters
    ----------
    spatial_stats : pd.DataFrame
        Columns region, method and the metric.
    metric : str
        Metric column.
    ylabel : str
        y-axis label.
    methods : list of str
        Methods in plotting order.
    regions : list of str
        Region labels (one axes each).
    colours : dict
        Method -> colour.

    Returns
    -------
    pd.DataFrame
        The plotted values (region, method, metric).
    """
    fig, axes = plt.subplots(1, len(regions), figsize=(5.2, 2.3), sharey=True)
    plotted = []
    for ax, region in zip(axes, regions):
        d = (
            spatial_stats[spatial_stats.region == region]
            .set_index("method")
            .loc[methods, metric]
        )
        ax.bar(
            range(len(methods)),
            d.values,
            color=[colours[m] for m in methods],
            width=0.75,
        )
        ax.set_xticks(range(len(methods)))
        ax.set_xticklabels(methods, rotation=60, ha="right", fontsize=6)
        ax.tick_params(axis="y", labelsize=6)
        ax.set_title(region, fontsize=8)
        ax.grid(False)
        sns.despine(ax=ax)
        plotted.append(
            pd.DataFrame({"region": region, "method": methods, metric: d.values})
        )
    axes[0].set_ylabel(ylabel, fontsize=7)
    plt.tight_layout()
    return pd.concat(plotted, ignore_index=True)


def plot_benchmark_dot_summary(table, ax, colours, lower_is_better=(), annotate=True):
    """
    Methods x metrics summary: circles sized and shaded by each metric scaled to 0-1 across methods (best = 1), raw
    values annotated, plus the mean rank across metrics.

    Parameters
    ----------
    table : pd.DataFrame
        Index = methods, columns = metrics.
    ax : plt.Axes
        Target axes.
    colours : dict
        Method -> colour of the tick label.
    lower_is_better : tuple of str
        Metrics where smaller is better (e.g. "Geary's C").
    annotate : bool
        Print the raw value in each cell.

    Returns
    -------
    plt.Axes
    """
    scaled = table.copy().astype(float)
    for c in table.columns:
        col = table[c].astype(float)
        rng_ = col.max() - col.min()
        s = (col - col.min()) / rng_ if rng_ > 0 else col * 0 + 1
        scaled[c] = 1 - s if c in lower_is_better else s
    ranks = scaled.rank(ascending=False).mean(axis=1)
    cmap = plt.get_cmap("Purples")
    for i, m in enumerate(table.index):
        for j, c in enumerate(table.columns):
            s = scaled.loc[m, c]
            ax.scatter(
                j, i, s=30 + 220 * s, color=cmap(0.25 + 0.75 * s), edgecolors="none"
            )
            if annotate:
                ax.text(
                    j,
                    i,
                    f"{table.loc[m, c]:.2g}",
                    ha="center",
                    va="center",
                    fontsize=5,
                    color="white" if s > 0.6 else "black",
                )
        ax.text(
            len(table.columns) + 0.2,
            i,
            f"{ranks[m]:.1f}",
            ha="left",
            va="center",
            fontsize=6,
        )
    ax.set_xticks(range(len(table.columns)))
    ax.set_xticklabels(table.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(table.index)))
    ax.set_yticklabels(table.index)
    for lab in ax.get_yticklabels():
        lab.set_color(colours[lab.get_text()])
    ax.text(len(table.columns) + 0.2, -0.9, "mean rank", fontsize=6, ha="left")
    ax.set_xlim(-0.6, len(table.columns) + 0.9)
    ax.set_ylim(len(table.index) - 0.4, -1.2)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)
    return ax


def plot_estimate_ci(ax, k, v, lo, hi, colour, marker="o"):
    """
    One estimate with its interval, drawn horizontally at row k.

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    k : float
        Row position.
    v, lo, hi : float
        Estimate and interval limits.
    colour : str
        Method colour.
    marker : str
        Marker of the estimate.
    """
    ax.plot([lo, hi], [k, k], color=colour, lw=0.9, solid_capstyle="butt")
    ax.plot(v, k, marker, ms=3.2, color=colour, mec="white", mew=0.3, zorder=3)


def format_benchmark_metric_axis(ax, met, n_rows):
    """
    x range, ticks and grid of one metric axes of the classification panel.

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    met : str
        Metric name ("AUROC" gets the 0.5-1 range, every other metric 0-1).
    n_rows : int
        Number of rows (methods or rules).
    """
    if met == "AUROC":
        ax.set_xlim(0.5, 1.02)
        ax.set_xticks([0.5, 0.75, 1.0])
        ax.set_xticklabels(["0.5", "", "1"])
    else:
        ax.set_xlim(0, 1.04)
        ax.set_xticks([0, 0.5, 1.0])
        ax.set_xticklabels(["0", "0.5", "1"])
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.tick_params(labelsize=6)
    ax.grid(False)
    ax.grid(axis="x", lw=0.3, color="0.88")
    ax.set_axisbelow(True)
    sns.despine(ax=ax)


def plot_fixed_rule_f1(
    ax, summary, region, fixed_rules, fixed_labels, colours, plotted=None
):
    """
    F1 of the fixed rules of the null model (z > 0, z > 1.96, p_empirical <= 0.05) with bootstrap intervals; a
    vertical tick marks the oracle F1 of the same method.

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    summary : pd.DataFrame
        Estimates with intervals (region, method, metric, value, lo, hi, n_rep).
    region : str
        Region label.
    fixed_rules : list of tuple
        (method, rule) pairs in plotting order.
    fixed_labels : dict
        Rule -> tick label.
    colours : dict
        Method -> colour.
    plotted : list or None
        Rows of plotted values are appended here.
    """
    for k, (m, rule) in enumerate(fixed_rules):
        d = summary[(summary.region == region) & (summary.method == m)].set_index(
            "metric"
        )
        r = d.loc[f"F1 ({rule})"]
        plot_estimate_ci(ax, k, r["value"], r["lo"], r["hi"], colours[m])
        oracle = d.loc["F1 (oracle threshold)", "value"]
        ax.plot([oracle, oracle], [k - 0.22, k + 0.22], color=colours[m], lw=0.8)
        if plotted is not None:
            plotted.append(
                {
                    "region": region,
                    "metric": f"F1 ({rule})",
                    "method / rule": f"{m}, {rule}",
                    "estimate": r["value"],
                    "95% CI lo": r["lo"],
                    "95% CI hi": r["hi"],
                    "n replicates": int(r["n_rep"]),
                    "oracle F1 of the same method (tick)": oracle,
                }
            )
    format_benchmark_metric_axis(ax, "F1", len(fixed_rules))
    ax.set_yticks(range(len(fixed_rules)))
    ax.set_yticklabels([fixed_labels[r] for _, r in fixed_rules], fontsize=6)
    for lab, (m, _) in zip(ax.get_yticklabels(), fixed_rules):
        lab.set_color(colours[m])


# --------------------------------------------------------------------------------------------------
# REV_SupFig3_4_controls
# --------------------------------------------------------------------------------------------------


def _get_dense_gene_expr(adata):
    """Return a (genes x cells) DataFrame from adata.X."""
    X_dense = adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X
    gene_expr = pd.DataFrame(X_dense.T, index=adata.var_names, columns=adata.obs_names)
    gene_expr.index = gene_expr.index.str.upper()
    return gene_expr


def score_signatures(
    adata,
    signatures_dict,
    methods=None,
    enrichmap_kwargs=None,
    ssgsea_kwargs=None,
    gsva_kwargs=None,
    min_size=5,
):
    """
    Score cells/spots for multiple gene signatures using several methods.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix (with spatial coordinates if needed downstream).
    signatures_dict : dict[str, list[str]]
        Mapping of signature names to gene lists,
        e.g. {"Pyramidal_layer": ["GENE1", "GENE2", ...], "Astrocyte": [...]}.
    methods : list[str] or None
        Subset of methods to run. Choose from:
        "enrichmap", "aucell", "scanpy", "ssgsea", "gsva", "zscore".
        If None, all methods are run.
    enrichmap_kwargs : dict or None
        Extra keyword arguments forwarded to em.tl.score().
    ssgsea_kwargs : dict or None
        Extra keyword arguments forwarded to gp.ssgsea().
    gsva_kwargs : dict or None
        Extra keyword arguments forwarded to gp.gsva().
    min_size : int
        Minimum gene set size for ssGSEA and GSVA (default 5).

    Returns
    -------
    adata : AnnData
        The input object with scores added to adata.obs as
        "{signature_name}_{method}" columns.
    """
    all_methods = ["enrichmap", "aucell", "scanpy", "ssgsea", "gsva", "zscore"]
    if methods is None:
        methods = all_methods
    else:
        invalid = set(methods) - set(all_methods)
        if invalid:
            raise ValueError(f"Unknown methods: {invalid}. Choose from {all_methods}")

    if enrichmap_kwargs is None:
        enrichmap_kwargs = {}
    if ssgsea_kwargs is None:
        ssgsea_kwargs = {}
    if gsva_kwargs is None:
        gsva_kwargs = {}

    # Upper-case version of signatures for gseapy methods
    sigs_upper = {k: [g.upper() for g in v] for k, v in signatures_dict.items()}

    # Build decoupler-format net DataFrame
    net_df = pd.DataFrame(
        [(key, gene) for key, genes in signatures_dict.items() for gene in genes],
        columns=["source", "target"],
    )

    # Pre-compute dense gene expression matrix (genes x cells) for gseapy/zscore
    need_dense = bool(set(methods) & {"ssgsea", "gsva", "zscore"})
    gene_expr = _get_dense_gene_expr(adata) if need_dense else None

    if "enrichmap" in methods:
        for sig_name, gene_set in signatures_dict.items():
            em.tl.score(adata, gene_set=gene_set, **enrichmap_kwargs)
            col = "enrichmap_score"
            if col in adata.obs.columns:
                adata.obs[f"{sig_name}_enrichmap"] = adata.obs[col].copy()

    if "aucell" in methods:
        dc.mt.aucell(adata, net=net_df)
        score_aucell = dc.pp.get_obsm(adata=adata, key="score_aucell")
        for sig_name in signatures_dict:
            if sig_name in score_aucell.obsm["score_aucell"].columns:
                adata.obs[f"{sig_name}_aucell"] = score_aucell.obsm["score_aucell"][
                    sig_name
                ].values
            elif "score_aucell" in adata.obsm:
                adata.obs[f"{sig_name}_aucell"] = adata.obsm[
                    "score_aucell"
                ].values.ravel()

    if "scanpy" in methods:
        for sig_name, gene_set in signatures_dict.items():
            sc.tl.score_genes(
                adata, gene_set, score_name=f"{sig_name}_scanpy", use_raw=False
            )

    if "ssgsea" in methods:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ssgsea_res = gp.ssgsea(
                data=gene_expr,
                gene_sets=sigs_upper,
                outdir=None,
                sample_norm_method="rank",
                no_plot=True,
                min_size=min_size,
                **ssgsea_kwargs,
            )
        nes = ssgsea_res.res2d.pivot(index="Term", columns="Name", values="NES")
        for sig_name in sigs_upper:
            if sig_name in nes.index:
                adata.obs[f"{sig_name}_ssgsea"] = (
                    nes.loc[sig_name].T.astype(float).reindex(adata.obs_names).values
                )

    if "gsva" in methods:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gsva_res = gp.gsva(
                data=gene_expr,
                gene_sets=sigs_upper,
                outdir=None,
                min_size=min_size,
                **gsva_kwargs,
            )
        es = gsva_res.res2d.pivot(index="Term", columns="Name", values="ES")
        for sig_name in sigs_upper:
            if sig_name in es.index:
                adata.obs[f"{sig_name}_gsva"] = (
                    es.loc[sig_name].T.astype(float).reindex(adata.obs_names).values
                )

    if "zscore" in methods:
        for sig_name, gene_set in sigs_upper.items():
            present = [g for g in gene_set if g in gene_expr.index]
            if present:
                z = gene_expr.loc[present].sum(axis=0) / np.sqrt(len(present))
                adata.obs[f"{sig_name}_zscore"] = z.values
            else:
                adata.obs[f"{sig_name}_zscore"] = np.nan

    return adata


def benchmark_rank_consistency(
    adata: AnnData,
    score_keys: Sequence[str],
    plot: bool = True,
    figsize: tuple[float, float] = (6, 5),
    cmap: str = "coolwarm",
    save: str | None = None,
) -> pd.DataFrame:
    """
    Pairwise Spearman rank correlation between all scoring methods.

    Tests whether EnrichMap produces genuinely different information or
    is simply a smoothed version of another method. If EnrichMap's
    correlation with, say, scanpy is ~1.0, the spatial processing adds
    no new information beyond post-hoc smoothing. A moderate correlation
    (e.g. 0.6-0.8) demonstrates that gene weighting and spatial
    smoothing produce a qualitatively different ranking of spots.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with scores in ``adata.obs``.
    score_keys : sequence of str
        Column names for each scoring method.
    plot : bool, default True
    figsize, cmap, save
        Plotting options.

    Returns
    -------
    pd.DataFrame
        Square correlation matrix indexed by method name.
    """
    # Collect valid scores
    score_data = {}
    for method in score_keys:
        if method in adata.obs.columns:
            score_data[method] = adata.obs[method].values.astype(np.float64)

    methods = list(score_data.keys())
    n = len(methods)
    corr_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            s1, s2 = score_data[methods[i]], score_data[methods[j]]
            valid = ~(np.isnan(s1) | np.isnan(s2))
            rho, _ = spearmanr(s1[valid], s2[valid])
            corr_matrix[i, j] = rho if not np.isnan(rho) else 0.0

    # Clean up labels
    clean_labels = [m.replace("_score", "") for m in methods]
    result = pd.DataFrame(corr_matrix, index=clean_labels, columns=clean_labels)

    if plot and not result.empty:
        _plot_rank_consistency(result, figsize=figsize, cmap=cmap, save=save)

    return result


def _plot_rank_consistency(df, figsize=(6, 5), cmap="coolwarm", save=None):
    """Heatmap of pairwise Spearman correlations between methods."""
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        df,
        annot=True,
        fmt=".2f",
        cmap=cmap,
        vmin=-1,
        vmax=1,
        square=True,
        linewidths=0,
        ax=ax,
        annot_kws={"size": 8},
    )
    ax.set_title("Rank consistency across methods", fontsize=10)
    ax.tick_params(axis="x", rotation=45, labelsize=8)
    ax.tick_params(axis="y", rotation=0, labelsize=8)
    ax.grid(False)
    plt.tight_layout()

    if save:
        os.makedirs(os.path.dirname(save) or "figures", exist_ok=True)
        path = os.path.join("figures", save) if not os.path.dirname(save) else save
        plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.show()
    return ax


def collect_method_scores(adata, gene_sets, method_keys):
    """
    Per-method score tables of one panel, taken from adata.obs.

    Parameters
    ----------
    adata : AnnData
        Data with the score columns in adata.obs.
    gene_sets : dict
        Signature name -> genes (only the names are used).
    method_keys : dict
        Method label -> obs column pattern with a "{gs}" placeholder, e.g. "{gs}_null_score".

    Returns
    -------
    dict
        Method label -> DataFrame (spots x signatures).
    """
    return {
        m: pd.DataFrame(
            {gs: adata.obs[key.format(gs=gs)].values for gs in gene_sets},
            index=adata.obs_names,
        )
        for m, key in method_keys.items()
    }


def gearys_c_binary(values, W):
    """
    Geary's C with a binary spatial weight matrix (lower = stronger positive autocorrelation; 1 = none).

    Parameters
    ----------
    values : np.ndarray
        Value per spot.
    W : scipy.sparse matrix
        Spatial weights (the 6-NN graph of em.tl.score).

    Returns
    -------
    float
        Geary's C.
    """
    v = np.asarray(values, dtype=np.float64)
    z = v - v.mean()
    z2 = z @ z
    if z2 == 0:
        return np.nan
    W = W.tocoo()
    num = (W.data * (v[W.row] - v[W.col]) ** 2).sum()
    return float((len(v) - 1) * num / (2 * W.data.sum() * z2))


def stratified_bootstrap_indices(y, n_boot=1000, seed=0):
    """
    Stratified spot bootstrap (positives and negatives resampled separately), shared by every method (paired).

    Parameters
    ----------
    y : np.ndarray
        Binary labels.
    n_boot : int
        Number of bootstrap samples.
    seed : int
        Random seed.

    Returns
    -------
    np.ndarray
        n_boot x n spot indices.
    """
    rng = np.random.default_rng(seed)
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    return np.array(
        [
            np.concatenate([rng.choice(pos, len(pos)), rng.choice(neg, len(neg))])
            for _ in range(n_boot)
        ]
    )


def auroc_auprc_with_ci(y, s, boot):
    """
    AUROC and AUPRC with percentile 95% intervals over bootstrap samples.

    Parameters
    ----------
    y : np.ndarray
        Binary labels.
    s : np.ndarray
        Score.
    boot : np.ndarray
        Bootstrap indices (stratified_bootstrap_indices).

    Returns
    -------
    dict
        AUROC, AUROC_lo, AUROC_hi, AUPRC, AUPRC_lo, AUPRC_hi.
    """
    roc = np.array([roc_auc_score(y[b], s[b]) for b in boot])
    prc = np.array([average_precision_score(y[b], s[b]) for b in boot])
    return {
        "AUROC": roc_auc_score(y, s),
        "AUROC_lo": np.percentile(roc, 2.5),
        "AUROC_hi": np.percentile(roc, 97.5),
        "AUPRC": average_precision_score(y, s),
        "AUPRC_lo": np.percentile(prc, 2.5),
        "AUPRC_hi": np.percentile(prc, 97.5),
    }


def spot_call_fractions(z_null, p_null, s_smooth, gsd_positive):
    """
    Fraction of spots called positive by the null-model rules, smoothing-only score > 0 and GSDensity antimode.

    Parameters
    ----------
    z_null, p_null : np.ndarray
        Null-model z and empirical p-value per spot.
    s_smooth : np.ndarray
        EnrichMap (smoothing-only) score.
    gsd_positive : np.ndarray
        GSDensity antimode label == positive.

    Returns
    -------
    dict
        "<method>: <rule>" -> fraction of spots.
    """
    return {
        "EnrichMap (null): z > 0": float(np.mean(z_null > 0)),
        "EnrichMap (null): z > 1.96": float(np.mean(z_null > 1.96)),
        "EnrichMap (null): p_empirical <= 0.05": float(np.mean(p_null <= 0.05)),
        "EnrichMap (smoothing-only): score > 0": float(np.mean(s_smooth > 0)),
        "GSDensity: antimode positive": float(np.mean(gsd_positive)),
    }


def within_region_spearman(s, covariate, clusters, min_spots=30):
    """
    Library-size dependence that is not explained by regional (biological) differences.

    Parameters
    ----------
    s : np.ndarray
        Score per spot.
    covariate : np.ndarray
        Library-size measure per spot.
    clusters : np.ndarray
        Region annotation per spot.
    min_spots : int
        Regions with fewer spots are left out.

    Returns
    -------
    float
        n-weighted mean of the Spearman ρ computed within each annotated region with >= min_spots spots.
    """
    rhos, weights = [], []
    for k in np.unique(clusters):
        m = clusters == k
        if m.sum() < min_spots:
            continue
        rho = spearmanr(s[m], covariate[m]).correlation
        if not np.isnan(rho):
            rhos.append(rho)
            weights.append(m.sum())
    return float(np.average(rhos, weights=weights))


def gene_contributions_pca_table(contrib_dict, housekeeping_genes, top_n_genes=10):
    """
    PCA coordinates of per-gene contribution vectors, as plotted by em.pl.gene_contributions_pca.

    Parameters
    ----------
    contrib_dict : dict
        Gene -> contribution vector of one signature (adata.uns["gene_contributions"][signature]).
    housekeeping_genes : list of str
        Genes labelled "housekeeping" (all others "signature").
    top_n_genes : int
        Number of genes highlighted (largest contribution-vector norm).

    Returns
    -------
    pd.DataFrame
        One row per gene: gene, PC1, PC2, contribution_norm, top_n (highlighted and labelled), gene_group.
    np.ndarray
        Explained variance ratio of PC1 and PC2.
    """
    gene_names = list(contrib_dict.keys())
    X = np.vstack([contrib_dict[gene] for gene in gene_names])
    pca = PCA(n_components=2, random_state=0)
    pcs = pca.fit_transform(X)
    norms = np.linalg.norm(X, axis=1)
    top_idx = np.argsort(norms)[-top_n_genes:]
    out = pd.DataFrame(
        {
            "gene": gene_names,
            "PC1": pcs[:, 0],
            "PC2": pcs[:, 1],
            "contribution_norm": norms,
            "top_n": np.isin(np.arange(len(gene_names)), top_idx),
        }
    )
    out["gene_group"] = np.where(
        out["gene"].isin(housekeeping_genes), "housekeeping", "signature"
    )
    return out.sort_values("contribution_norm", ascending=False).reset_index(
        drop=True
    ), pca.explained_variance_ratio_


def plot_gene_contributions_pca(ax, df, evr, title, group_colours):
    """
    Gene-contribution PCA: signature vs housekeeping genes, the 10 largest contributions labelled.

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    df : pd.DataFrame
        Output of gene_contributions_pca_table.
    evr : np.ndarray
        Explained variance ratio of PC1 and PC2.
    title : str
        Panel title.
    group_colours : dict
        Gene group ("signature", "housekeeping") -> colour.
    """
    ax.axhline(0, ls="--", lw=0.4, c="k", alpha=0.4)
    ax.axvline(0, ls="--", lw=0.4, c="k", alpha=0.4)
    for group in ["housekeeping", "signature"]:
        d = df[df["gene_group"] == group]
        ax.scatter(
            d["PC1"],
            d["PC2"],
            s=10,
            c=group_colours[group],
            edgecolors="k",
            linewidths=0.3,
            label=f"{group.capitalize()} genes",
            zorder=3,
        )
    texts = [
        ax.text(r.PC1, r.PC2, r.gene, fontsize=5.5, style="italic")
        for r in df[df["top_n"]].itertuples()
    ]
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", lw=0.3, color="k"))
    ax.set_xlabel(f"PC1 ({100 * evr[0]:.0f}%)")
    ax.set_ylabel(f"PC2 ({100 * evr[1]:.0f}%)")
    ax.set_title(title, fontsize=8)
    ax.grid(False)
    sns.despine(ax=ax)


def plot_metric_pair_ci(
    table, x, y, ax, colours, ref_x=None, ref_y=None, invert_y=False, marker_col=None
):
    """
    One metric against another, one point per method (and gene set), with 95% CI whiskers when columns
    f"{x}_lo"/f"{x}_hi" and f"{y}_lo"/f"{y}_hi" exist.

    Parameters
    ----------
    table : pd.DataFrame
        One row per method (and gene set), with a "method" column and the metric columns.
    x, y : str
        Metric columns, e.g. "AUROC" and "AUPRC".
    ax : plt.Axes
        Target axes.
    colours : dict
        Method label -> colour.
    ref_x, ref_y : float or None
        Dotted reference lines (e.g. AUPRC = prevalence).
    invert_y : bool
        Invert the y axis for lower-is-better metrics (e.g. Geary's C), so better is always up and right.
    marker_col : str or None
        Column holding a marker per row; default circles.
    """
    for _, r in table.iterrows():
        xerr = yerr = None
        if f"{x}_lo" in table:
            xerr = [[r[x] - r[f"{x}_lo"]], [r[f"{x}_hi"] - r[x]]]
        if f"{y}_lo" in table:
            yerr = [[r[y] - r[f"{y}_lo"]], [r[f"{y}_hi"] - r[y]]]
        ax.errorbar(
            r[x],
            r[y],
            xerr=xerr,
            yerr=yerr,
            fmt=r[marker_col] if marker_col else "o",
            color=colours[r["method"]],
            ms=5,
            mec="white",
            mew=0.4,
            elinewidth=0.6,
            capsize=0,
            label=r["method"],
        )
    if ref_x is not None:
        ax.axvline(ref_x, ls=":", c="0.6", lw=0.6)
    if ref_y is not None:
        ax.axhline(ref_y, ls=":", c="0.6", lw=0.6)
    if invert_y:
        ax.invert_yaxis()
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.grid(False)


def plot_control_score_histograms(scores, cell_types, titles, colours, calls_table, reference=None):
    """
    Histograms of every Sup. Fig. 3c map (rows = methods, columns = cell types; 50 bins, bins and x axis shared per method),
    with the fractions of spots called under the null-model and GSDensity rows.

    Parameters
    ----------
    scores : dict
        Method label -> DataFrame (spots x cell types).
    cell_types : list of str
        Columns of the score tables, in plotting order.
    titles : list of str
        Column titles.
    colours : dict
        Method label -> colour (also sets the row order).
    calls_table : pd.DataFrame
        Fraction of spots called, rows "<method>: <rule>" (spot_call_fractions), columns = cell types.
    reference : str or None
        Cell type whose range sets the x axis of each method (default: the range of all its cell types).

    Returns
    -------
    plt.Figure
    """
    null = "EnrichMap (null)"
    methods = list(colours)
    fig, axes = plt.subplots(
        len(methods), len(cell_types), figsize=(5.4, 0.62 * len(methods) + 0.3)
    )
    for i, m in enumerate(methods):
        ref = scores[m].values if reference is None else scores[m][reference].values
        lo, hi = np.nanmin(ref), np.nanmax(ref)
        edges = np.linspace(lo - 0.02 * (hi - lo), hi + 0.02 * (hi - lo), 51)
        for j, ct in enumerate(cell_types):
            ax = axes[i, j]
            ax.hist(scores[m][ct].values, bins=edges, color=colours[m], linewidth=0)
            ax.set_xlim(edges[0], edges[-1])
            ax.set_yticks([])
            ax.grid(False)
            ax.tick_params(axis="x", labelsize=5, length=1.5, width=0.4, pad=1)
            if i == 0:
                ax.set_title(titles[j], fontsize=7, pad=14)
            if j == 0:
                ax.set_ylabel(
                    m.replace(" (", "\n("),
                    fontsize=6,
                    color=colours[m],
                    rotation=0,
                    ha="right",
                    va="center",
                )
            sns.despine(ax=ax, left=True)
            if m == null:
                ax.axvline(1.96, ls="--", lw=0.4, c="k")
                ax.text(
                    1.96,
                    1.0,
                    "z > 1.96",
                    transform=ax.get_xaxis_transform(),
                    ha="center",
                    va="bottom",
                    fontsize=4.2,
                )
    plt.tight_layout()
    return fig


def plot_library_size_dependence(
    axes,
    lib,
    lib_random,
    regions,
    region_labels,
    colours,
    value,
    within,
    random_value,
    xlabel,
    truth="Annotation (truth)",
):
    """
    Library-size panels: one bar chart per region (bars = ρ over all spots, diamonds = within-region ρ) and the
    random-set strip plot.

    Parameters
    ----------
    axes : np.ndarray
        Three axes (one per region, then random sets).
    lib : pd.DataFrame
        One row per region and method with the ρ columns.
    lib_random : pd.DataFrame
        One row per random gene set and method.
    regions : list of str
        Region keys, in plotting order.
    region_labels : dict
        Region key -> label used in lib["region"].
    colours : dict
        Method label -> colour (also sets the order).
    value, within : str
        Columns of `lib` for the bars (all spots) and diamonds (within-region ρ).
    random_value : str
        Column of `lib_random`.
    xlabel : str
        x-axis label.
    truth : str
        Method label of the region indicator (drawn in black).
    """
    methods = list(colours)
    bar_methods = methods + [truth]
    y = np.arange(len(bar_methods))
    for ax, ct in zip(axes[:2], regions):
        d = lib[lib["region"] == region_labels[ct]].set_index("method").loc[bar_methods]
        ax.barh(
            y,
            d[value],
            color=[colours.get(m, "k") for m in bar_methods],
            height=0.7,
            linewidth=0,
        )
        ax.scatter(
            d[within],
            y,
            marker="D",
            s=9,
            facecolors="white",
            edgecolors="k",
            linewidths=0.5,
            zorder=3,
        )
        ax.axvline(0, ls="--", lw=0.5, c="k")
        ax.set_yticks(y)
        ax.set_yticklabels(bar_methods if ct == regions[0] else [])
        ax.invert_yaxis()
        ax.set_xlabel(xlabel)
        ax.set_title(region_labels[ct], fontsize=8)
    ax = axes[2]
    for k, m in enumerate(methods):
        v = lib_random.loc[lib_random["method"] == m, random_value].values
        jitter = np.random.default_rng(k).uniform(-0.25, 0.25, len(v))
        ax.scatter(v, k + jitter, s=4, color=colours[m], linewidths=0)
        ax.vlines(np.median(v), k - 0.35, k + 0.35, colors="k", lw=0.8)
    ax.axvline(0, ls="--", lw=0.5, c="k")
    ax.set_yticks(np.arange(len(methods)))
    ax.set_yticklabels(methods)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(f"{lib_random['gene set'].nunique()} random 10-gene sets", fontsize=8)
    for ax in axes:
        ax.grid(False)
        ax.tick_params(labelsize=6)
        for lab in ax.get_yticklabels():
            lab.set_color(colours.get(lab.get_text(), "k"))
        sns.despine(ax=ax)
    handles = [
        Line2D(
            [],
            [],
            marker="D",
            ls="",
            mfc="white",
            mec="k",
            mew=0.5,
            ms=3,
            label="within-region ρ",
        )
    ]
    axes[1].legend(handles=handles, frameon=False, fontsize=5.5, loc="upper left")


def plot_method_summary_dots(table, ax, colours, lower_is_better=()):
    """
    Compact benchmark summary (methods x metrics): circles sized and shaded by each metric scaled to 0-1 across
    methods (best = 1), raw values annotated, plus a mean-rank column.

    Parameters
    ----------
    table : pd.DataFrame
        Index = method labels, columns = metrics.
    ax : plt.Axes
        Target axes.
    colours : dict
        Method label -> colour of the row label.
    lower_is_better : tuple of str
        Metrics where smaller is better.
    """
    scaled = table.copy().astype(float)
    for c in table.columns:
        col = table[c].astype(float)
        rng_ = col.max() - col.min()
        s = (col - col.min()) / rng_ if rng_ > 0 else col * 0 + 1
        scaled[c] = 1 - s if c in lower_is_better else s
    ranks = scaled.rank(ascending=False).mean(axis=1)
    cmap = plt.get_cmap("Purples")
    for i, m in enumerate(table.index):
        for j, c in enumerate(table.columns):
            s = scaled.loc[m, c]
            ax.scatter(
                j, i, s=30 + 220 * s, color=cmap(0.25 + 0.75 * s), edgecolors="none"
            )
            ax.text(
                j,
                i,
                f"{table.loc[m, c]:.2g}",
                ha="center",
                va="center",
                fontsize=5,
                color="white" if s > 0.6 else "black",
            )
        ax.text(
            len(table.columns) + 0.2,
            i,
            f"{ranks[m]:.1f}",
            ha="left",
            va="center",
            fontsize=6,
        )
    ax.set_xticks(range(len(table.columns)))
    ax.set_xticklabels(table.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(table.index)))
    ax.set_yticklabels(table.index)
    for lab in ax.get_yticklabels():
        lab.set_color(colours[lab.get_text()])
    ax.text(len(table.columns) + 0.2, -0.9, "mean rank", fontsize=6, ha="left")
    ax.set_xlim(-0.6, len(table.columns) + 0.9)
    ax.set_ylim(len(table.index) - 0.4, -1.2)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)


# --------------------------------------------------------------------------------------------------
# REV_SupFig5_6_7_ed4_layout
# --------------------------------------------------------------------------------------------------


def ed4_wide_auc(auc_table, ct):
    """
    AUROC/AUPRC with 95% intervals of one region, one row per method and size.

    Parameters
    ----------
    auc_table : pd.DataFrame
        Long table with 'region', 'method', 'gene set size', 'metric', 'value', 'lo', 'hi'.
    ct : str
        Region.

    Returns
    -------
    pd.DataFrame
        method, gene set size, AUROC, AUROC_lo, AUROC_hi, AUPRC, AUPRC_lo, AUPRC_hi.
    """
    d = auc_table[auc_table.region == ct]
    out = d.pivot_table(
        index=["method", "gene set size"],
        columns="metric",
        values=["value", "lo", "hi"],
    )
    t = pd.DataFrame(index=out.index)
    for met in ["AUROC", "AUPRC"]:
        t[met], t[f"{met}_lo"], t[f"{met}_hi"] = (
            out[("value", met)],
            out[("lo", met)],
            out[("hi", met)],
        )
    return t.reset_index()


def ed4_style(ax, title=""):
    """
    Final-size styling: 5 pt ticks, 6 pt labels and titles, no grid, despined.

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    title : str
        Axes title (none if empty).
    """
    ax.tick_params(labelsize=5, length=1.5, width=0.4, pad=1)
    ax.xaxis.label.set_fontsize(6)
    ax.yaxis.label.set_fontsize(6)
    if title:
        ax.set_title(title, fontsize=6, pad=2)
    ax.grid(False)
    sns.despine(ax=ax)
    for s in ax.spines.values():
        s.set_linewidth(0.4)


def ed4_metric_pair(
    ax, table, x, y, methods, colours, ref_y=None, invert_y=False, label_fontsize=4
):
    """
    Metric `x` against metric `y`, one filled circle per method (colour) and gene set size, 95% interval whiskers
    where f"{x}_lo" / f"{x}_hi" columns exist; the sizes of one method are joined and each point is labelled with
    its gene set size.

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    table : pd.DataFrame
        Rows of one region with 'method', 'gene set size', `x`, `y` (optional *_lo / *_hi columns).
    x, y : str
        Metric columns.
    methods : list of str
        Methods.
    colours : dict
        Method -> colour.
    ref_y : float or None
        Dotted reference line (AUPRC = prevalence).
    invert_y : bool
        Invert the y axis for a lower-is-better metric (Geary's C).
    label_fontsize : float
        Font size of the gene set size labels.
    """
    for m in methods:
        d = table[table.method == m].sort_values("gene set size")
        ax.plot(d[x], d[y], "-", color=colours[m], lw=0.5, alpha=0.6, zorder=0)
        for _, r in d.iterrows():
            xerr = (
                [[r[x] - r[f"{x}_lo"]], [r[f"{x}_hi"] - r[x]]]
                if f"{x}_lo" in d
                else None
            )
            yerr = (
                [[r[y] - r[f"{y}_lo"]], [r[f"{y}_hi"] - r[y]]]
                if f"{y}_lo" in d
                else None
            )
            ax.errorbar(
                r[x],
                r[y],
                xerr=xerr,
                yerr=yerr,
                fmt="o",
                color=colours[m],
                ms=2.5,
                mec="white",
                mew=0.3,
                elinewidth=0.4,
                capsize=0,
            )
            ax.annotate(
                str(int(r["gene set size"])),
                (r[x], r[y]),
                xytext=(1.5, 1.5),
                textcoords="offset points",
                fontsize=label_fontsize,
                color=colours[m],
                ha="left",
                va="bottom",
            )
    if ref_y is not None:
        ax.axhline(ref_y, ls=":", c="0.6", lw=0.6)
    if invert_y:
        ax.invert_yaxis()
    ax.set_xlabel(x)
    ax.set_ylabel(y)


def ed4_statistic_vs_size(ax, corr, col, ct, all_sizes_ct, methods, colours):
    """
    Statistic `col` of Sup. Fig. 6a-d against gene set size (categorical axis), one line per method.

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    corr : pd.DataFrame
        'region', 'method', 'gene set size' and `col`.
    col : str
        Statistic ('pearson r' or 'spearman rho').
    ct : str
        Region.
    all_sizes_ct : list of int
        Gene set sizes of the region (x order).
    methods : list of str
        Methods.
    colours : dict
        Method -> colour.
    """
    for m in methods:
        d = (
            corr[(corr.region == ct) & (corr.method == m)]
            .dropna(subset=[col])
            .sort_values("gene set size")
        )
        x = [all_sizes_ct.index(n) for n in d["gene set size"]]
        ax.plot(x, d[col], "o-", ms=2, lw=0.6, color=colours[m], mec="white", mew=0.2)
    ax.set_xticks(range(len(all_sizes_ct)))
    ax.set_xticklabels([str(n) for n in all_sizes_ct], rotation=90)
    ax.set_xlabel("Gene set size (genes)")
    ax.axhline(0, color="0.75", lw=0.4, ls=":")


def ed4_runtime_axis(
    ax, runtime, ct, all_sizes_ct, runtime_methods, published_methods, colours
):
    """
    Runtime against gene set size (log-log): this server (mean, min-max whiskers, solid) and the published single
    run (dashed); every point is a filled circle in the method colour.

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    runtime : pd.DataFrame
        'region', 'method', 'gene set size', 'mean', 'min', 'max', 'published time (s)'.
    ct : str
        Region.
    all_sizes_ct : list of int
        Gene set sizes of the region.
    runtime_methods : list of str
        Methods.
    published_methods : list of str
        Methods with a published run.
    colours : dict
        Method -> colour.
    """
    for m in runtime_methods:
        d = runtime[(runtime.region == ct) & (runtime.method == m)].sort_values(
            "gene set size"
        )
        ax.errorbar(
            d["gene set size"],
            d["mean"],
            yerr=[d["mean"] - d["min"], d["max"] - d["mean"]],
            fmt="o-",
            ms=1.8,
            lw=0.6,
            elinewidth=0.4,
            capsize=0,
            color=colours[m],
            mec="white",
            mew=0.2,
        )
        if m in published_methods:
            p = d.dropna(subset=["published time (s)"])
            ax.plot(
                p["gene set size"],
                p["published time (s)"],
                "o--",
                ms=1.8,
                lw=0.4,
                mec="white",
                mew=0.2,
                color=colours[m],
                alpha=0.8,
            )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(all_sizes_ct)
    ax.set_xticklabels([str(n) for n in all_sizes_ct], rotation=90)
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_xlabel("Gene set size (genes)")


def ed4_panel_letter(fig, ax, letter, dx=-0.045, dy=0.012):
    """
    Bold panel letter (8 pt) at the top left of an axes; the offset is in figure fractions from the axes corner and
    follows the axes box when its aspect changes at draw time.

    Parameters
    ----------
    fig : plt.Figure
        Figure.
    ax : plt.Axes
        Axes the letter belongs to.
    letter : str
        Panel letter.
    dx, dy : float
        Offset from the top-left corner of the axes (fractions of the figure width and height).
    """
    w, h = fig.get_size_inches()
    trans = mtransforms.offset_copy(
        ax.transAxes, fig=fig, x=dx * w, y=dy * h, units="inches"
    )
    ax.text(
        0,
        1,
        letter,
        transform=trans,
        fontsize=8,
        fontweight="bold",
        va="bottom",
        ha="left",
        clip_on=False,
    )


def ed4_draw_metrics(
    autocorr,
    auc_table,
    corr,
    runtime,
    msigdb,
    maps,
    regions,
    all_sizes,
    methods,
    runtime_methods,
    published_methods,
    colours,
):
    """
    Extended Data Fig. 4b-f: b, Moran's I vs Geary's C; c, AUROC vs AUPRC (95% intervals); d, cross-Moran r and
    Spearman rho against gene set size; e, runtime; f, MSigDB gene set sizes.

    Parameters
    ----------
    autocorr : pd.DataFrame
        Moran's I and Geary's C per region, method and gene set size.
    auc_table : pd.DataFrame
        AUROC/AUPRC with 95% intervals (long).
    corr : pd.DataFrame
        Cross-Moran r and Spearman rho per region, method and gene set size.
    runtime : pd.DataFrame
        Runtime statistics per region, method and gene set size.
    msigdb : pd.DataFrame
        MSigDB gene set sizes ('collection', 'collection_label', 'x_order', 'gene_set_size').
    maps : pd.DataFrame
        Spot table with 'cluster' (region prevalence).
    regions : list of str
        Regions.
    all_sizes : dict
        Region -> the eight gene set sizes.
    methods : list of str
        Methods of panels b-d.
    runtime_methods : list of str
        Methods of panel e.
    published_methods : list of str
        Methods with a published runtime.
    colours : dict
        Method -> colour.

    Returns
    -------
    plt.Figure
    """
    fig = plt.figure(figsize=(180 / 25.4, 136 / 25.4))
    outer = fig.add_gridspec(
        3,
        1,
        height_ratios=[1, 1, 1.1],
        hspace=0.6,
        left=0.07,
        right=0.84,
        top=0.95,
        bottom=0.1,
    )
    row1 = outer[0].subgridspec(1, 4, wspace=0.62)
    row2 = outer[1].subgridspec(1, 4, wspace=0.62)
    row3 = outer[2].subgridspec(1, 4, wspace=0.62, width_ratios=[1, 1, 0.55, 1.7])
    axes = {}
    # b: Moran's I vs Geary's C (inverted: better is up and right)
    for j, ct in enumerate(regions):
        ax = fig.add_subplot(row1[0, j])
        ed4_metric_pair(
            ax,
            autocorr[autocorr.region == ct],
            "Moran's I",
            "Geary's C",
            methods,
            colours,
            invert_y=True,
        )
        ed4_style(ax, ct.replace("_", " "))
        axes[f"b{j}"] = ax
    # c: AUROC vs AUPRC with 95% paired spatial block bootstrap intervals
    for j, ct in enumerate(regions):
        ax = fig.add_subplot(row1[0, 2 + j])
        prevalence = float((maps["cluster"] == ct).mean())
        ed4_metric_pair(
            ax,
            ed4_wide_auc(auc_table, ct),
            "AUROC",
            "AUPRC",
            methods,
            colours,
            ref_y=prevalence,
        )
        ax.set_ylim(0, 1.02)
        ed4_style(ax, f"{ct.replace('_', ' ')} ({100 * prevalence:.1f}%)")
        axes[f"c{j}"] = ax
    # d: cross-Moran r and Spearman rho against the n = 50 score of the same method
    for j, ct in enumerate(regions):
        ax = fig.add_subplot(row2[0, j])
        ed4_statistic_vs_size(
            ax, corr, "pearson r", ct, all_sizes[ct], methods, colours
        )
        ax.set_ylabel("Cross-Moran r (vs n = 50)")
        ax.set_ylim(0.35, 1.02)
        ed4_style(ax, ct.replace("_", " "))
        axes[f"d{j}"] = ax
        ax = fig.add_subplot(row2[0, 2 + j])
        ed4_statistic_vs_size(
            ax, corr, "spearman rho", ct, all_sizes[ct], methods, colours
        )
        ax.set_ylabel("Spearman ρ (vs n = 50)")
        ax.set_ylim(-0.1, 1.02)
        ed4_style(ax, ct.replace("_", " "))
        axes[f"d{2 + j}"] = ax
    # e: runtime on this server and published
    for j, ct in enumerate(regions):
        ax = fig.add_subplot(row3[0, j], sharey=axes.get("e0"))
        ed4_runtime_axis(
            ax, runtime, ct, all_sizes[ct], runtime_methods, published_methods, colours
        )
        ax.set_ylabel("Time (s)" if j == 0 else "")
        ed4_style(ax, ct.replace("_", " "))
        axes[f"e{j}"] = ax
    # f: MSigDB gene set sizes (horizontal violins, log scale; grey band 10-100 genes)
    ax = fig.add_subplot(row3[0, 3])
    order = msigdb.drop_duplicates("collection").sort_values("x_order")
    hue_order = sorted(msigdb["collection"].unique())
    ax.axvspan(10, 100, color="0.92", lw=0, zorder=0)
    sns.violinplot(
        data=msigdb,
        y="collection",
        x="gene_set_size",
        orient="h",
        order=order["collection"],
        log_scale=True,
        hue="collection",
        hue_order=hue_order,
        palette=sns.color_palette("husl", len(hue_order)),
        inner="box",
        linewidth=0.3,
        width=0.85,
        legend=False,
        ax=ax,
    )
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order["collection_label"])
    ax.set_ylabel("")
    ax.set_xlabel("Gene set size (genes)")
    ed4_style(ax, "MSigDB collections")
    axes["f"] = ax
    # legends (x in figure fractions, y following the top of panels c and e at draw time): methods of b-d (point
    # labels = number of genes), runtime methods and line styles of e
    w, h = fig.get_size_inches()
    handles = [
        Line2D([], [], ls="-", marker="o", ms=3, lw=0.6, color=colours[m], label=m)
        for m in methods
    ]
    trans = mtransforms.offset_copy(
        mtransforms.blended_transform_factory(fig.transFigure, axes["c1"].transAxes),
        fig=fig,
        y=0.01 * h,
        units="inches",
    )
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.855, 1.0),
        bbox_transform=trans,
        frameon=False,
        fontsize=5,
        handletextpad=0.3,
        labelspacing=0.35,
        title="Point labels: number of genes",
        title_fontsize=5,
    )
    handles = [
        Line2D([], [], ls="-", marker="o", ms=2.5, lw=0.6, color=colours[m], label=m)
        for m in runtime_methods
    ]
    handles += [
        Line2D(
            [],
            [],
            ls="-",
            marker="o",
            ms=2.5,
            lw=0.6,
            color="0.3",
            label="this server, 3 repeats",
        ),
        Line2D(
            [],
            [],
            ls="--",
            marker="o",
            ms=2.5,
            lw=0.4,
            color="0.3",
            label="published, 1 run",
        ),
    ]
    trans = mtransforms.offset_copy(
        mtransforms.blended_transform_factory(fig.transFigure, axes["e1"].transAxes),
        fig=fig,
        y=0.03 * h,
        units="inches",
    )
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.855, 1.0),
        bbox_transform=trans,
        frameon=False,
        fontsize=5,
        handletextpad=0.3,
        labelspacing=0.35,
    )
    for key, letter in [("b0", "b"), ("c0", "c"), ("d0", "d"), ("e0", "e"), ("f", "f")]:
        ed4_panel_letter(fig, axes[key], letter)
    return fig


# --------------------------------------------------------------------------------------------------
# REV_SupFig5_6_7_gene_set_size
# --------------------------------------------------------------------------------------------------


def gss_size_ticks(sizes_ct):
    """
    Tick labels 'n = ...' of one region's gene set sizes.

    Parameters
    ----------
    sizes_ct : list of int
        Gene set sizes of one region.

    Returns
    -------
    list of str
    """
    return [f"n = {n}" for n in sizes_ct]


def gss_grouped_bars(ax, df, value, sizes_ct, methods, colours):
    """
    Bars of one metric per gene set size, one bar per method.

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    df : pd.DataFrame
        Rows of one region with columns 'method', 'gene set size' and `value`.
    value : str
        Column to plot.
    sizes_ct : list of int
        Gene set sizes of the region (x order).
    methods : list of str
        Methods (bar order).
    colours : dict
        Method -> colour.
    """
    width = 0.8 / len(methods)
    for j, m in enumerate(methods):
        d = df[df.method == m].set_index("gene set size").reindex(sizes_ct)
        x = np.arange(len(sizes_ct)) + (j - (len(methods) - 1) / 2) * width
        ax.bar(x, d[value].values, width=width, color=colours[m], linewidth=0, label=m)
    ax.set_xticks(np.arange(len(sizes_ct)))
    ax.set_xticklabels(gss_size_ticks(sizes_ct), rotation=90, fontsize=6)
    ax.tick_params(axis="y", labelsize=6)
    ax.grid(False)
    sns.despine(ax=ax)


def gss_bar_panel(df, value, ylabel, ylim, regions, sizes, methods, colours):
    """
    Sup. Fig. 5b/c as published: one axes per region (stacked), grouped bars per method.

    Parameters
    ----------
    df : pd.DataFrame
        Long table with 'region', 'method', 'gene set size' and `value`.
    value : str
        Column to plot.
    ylabel : str
        y-axis label.
    ylim : tuple
        y-axis limits.
    regions : list of str
        Regions (one axes each).
    sizes : dict
        Region -> gene set sizes.
    methods : list of str
        Methods.
    colours : dict
        Method -> colour.

    Returns
    -------
    plt.Figure
    """
    fig, axes = plt.subplots(2, 1, figsize=(3.6, 4.2))
    for ax, ct in zip(axes, regions):
        gss_grouped_bars(ax, df[df.region == ct], value, sizes[ct], methods, colours)
        ax.set_ylim(*ylim)
        ax.set_ylabel(ylabel, fontsize=7)
        ax.set_title(ct.replace("_", " "), fontsize=8)
    axes[0].legend(
        frameon=False, fontsize=6, loc="upper left", bbox_to_anchor=(1.01, 1.0)
    )
    plt.tight_layout()
    return fig


def gss_pair_panel(
    ax,
    table,
    x,
    y,
    methods,
    colours,
    title="",
    ref_y=None,
    invert_y=False,
    label_fontsize=4.5,
):
    """
    Compact panel: metric `x` against metric `y`, one filled circle per method (colour) and gene set size, with
    95% interval whiskers where f"{x}_lo" / f"{x}_hi" columns exist; the points of one method are joined by size and
    each point is labelled with its gene set size.

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    table : pd.DataFrame
        Rows of one region with 'method', 'gene set size', `x`, `y` (and optional *_lo / *_hi columns).
    x, y : str
        Metric columns.
    methods : list of str
        Methods.
    colours : dict
        Method -> colour.
    title : str
        Axes title.
    ref_y : float or None
        Dotted reference line (AUPRC = prevalence).
    invert_y : bool
        Invert the y axis for a lower-is-better metric (Geary's C), so better is up and right.
    label_fontsize : float
        Font size of the gene set size labels.
    """
    for m in methods:
        d = table[table.method == m].sort_values("gene set size")
        ax.plot(d[x], d[y], "-", color=colours[m], lw=0.6, alpha=0.6, zorder=0)
        for _, r in d.iterrows():
            xerr = (
                [[r[x] - r[f"{x}_lo"]], [r[f"{x}_hi"] - r[x]]]
                if f"{x}_lo" in d
                else None
            )
            yerr = (
                [[r[y] - r[f"{y}_lo"]], [r[f"{y}_hi"] - r[y]]]
                if f"{y}_lo" in d
                else None
            )
            ax.errorbar(
                r[x],
                r[y],
                xerr=xerr,
                yerr=yerr,
                fmt="o",
                color=colours[m],
                ms=4,
                mec="white",
                mew=0.4,
                elinewidth=0.6,
                capsize=0,
            )
            ax.annotate(
                str(int(r["gene set size"])),
                (r[x], r[y]),
                xytext=(2.5, 2.5),
                textcoords="offset points",
                fontsize=label_fontsize,
                color=colours[m],
                ha="left",
                va="bottom",
            )
    if ref_y is not None:
        ax.axhline(ref_y, ls=":", c="0.6", lw=0.6)
    if invert_y:
        ax.invert_yaxis()
    ax.set_xlabel(x, fontsize=7)
    ax.set_ylabel(y, fontsize=7)
    ax.set_title(title, fontsize=8)
    ax.tick_params(labelsize=6)
    ax.grid(False)
    sns.despine(ax=ax)


def gss_method_legend(ax, methods, colours, fontsize=5.5):
    """
    Legend of the compact panels: one filled circle per method; the point labels are the gene set sizes.

    Parameters
    ----------
    ax : plt.Axes
        Axes the legend is attached to (placed to its right).
    methods : list of str
        Methods.
    colours : dict
        Method -> colour.
    fontsize : float
        Legend font size.
    """
    handles = [
        Line2D([], [], ls="", marker="o", ms=4, color=colours[m], label=m)
        for m in methods
    ]
    ax.legend(
        handles=handles,
        frameon=False,
        fontsize=fontsize,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        handletextpad=0.2,
        title="Point labels: number of genes",
        title_fontsize=fontsize,
    )


def gss_plot_summary(table, ax, colours, lower_is_better=()):
    """
    Methods x metrics summary: circles sized and shaded by each metric scaled to 0-1 across methods (best = 1),
    raw values annotated and the mean rank on the right.

    Parameters
    ----------
    table : pd.DataFrame
        Index = methods, columns = metrics.
    ax : plt.Axes
        Target axes.
    colours : dict
        Method -> colour of the tick label.
    lower_is_better : tuple of str
        Metrics where smaller is better (e.g. Geary's C).
    """
    scaled = table.astype(float).copy()
    for c in table.columns:
        col = table[c].astype(float)
        s = (
            (col - col.min()) / (col.max() - col.min())
            if col.max() > col.min()
            else col * 0 + 1
        )
        scaled[c] = 1 - s if c in lower_is_better else s
    ranks = scaled.rank(ascending=False).mean(axis=1)
    cmap = plt.get_cmap("Purples")
    for i, m in enumerate(table.index):
        for j, c in enumerate(table.columns):
            s = scaled.loc[m, c]
            ax.scatter(
                j, i, s=30 + 220 * s, color=cmap(0.25 + 0.75 * s), edgecolors="none"
            )
            ax.text(
                j,
                i,
                f"{table.loc[m, c]:.2g}",
                ha="center",
                va="center",
                fontsize=5,
                color="white" if s > 0.6 else "black",
            )
        ax.text(
            len(table.columns) + 0.2,
            i,
            f"{ranks[m]:.1f}",
            ha="left",
            va="center",
            fontsize=6,
        )
    ax.set_xticks(range(len(table.columns)))
    ax.set_xticklabels(table.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(table.index)))
    ax.set_yticklabels(table.index)
    for lab in ax.get_yticklabels():
        lab.set_color(colours.get(lab.get_text(), "black"))
    ax.text(len(table.columns) + 0.2, -0.9, "mean rank", fontsize=6, ha="left")
    ax.set_xlim(-0.6, len(table.columns) + 0.9)
    ax.set_ylim(len(table.index) - 0.4, -1.2)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)


def gss_draw_variograms(
    variograms,
    normalised,
    methods,
    regions,
    sizes,
    colours,
    size_shades,
    max_lag,
    row_label,
    semivariance_unit,
):
    """
    Sup. Fig. 5d: one axes per method (rows) and region (columns), one line per gene set size.

    Parameters
    ----------
    variograms : pd.DataFrame
        Experimental variogram bins ('region', 'method', 'gene set size', 'spatial lag (px, bin upper edge)',
        'semivariance', 'score variance (dashed line)', 'semivariance / variance', 'within plotted x range').
    normalised : bool
        True: semivariance divided by the score variance (unit-free, comparable across methods; compact panel);
        False: semivariance with the score variance as dashed lines (as published; full panel).
    methods : list of str
        Methods (rows).
    regions : list of str
        Regions (columns).
    sizes : dict
        Region -> gene set sizes.
    colours : dict
        Method -> colour (row labels).
    size_shades : array-like
        One colour per gene set size (light = small, dark = large).
    max_lag : float
        x-axis limit (95th percentile of the pairwise spot distances).
    row_label : dict
        Method -> row label.
    semivariance_unit : dict
        Method -> semivariance unit of the full panel.

    Returns
    -------
    plt.Figure
    """
    fig, axes = plt.subplots(
        len(methods),
        2,
        figsize=(4.6, 1.35 * len(methods)) if normalised else (5.6, 2.07 * len(methods)),
        sharex=True,
        sharey=normalised,
        squeeze=False,
    )
    y_top = (
        1.05
        * variograms.loc[
            variograms["within plotted x range"], "semivariance / variance"
        ].max()
    )
    for i, m in enumerate(methods):
        for j, ct in enumerate(regions):
            ax = axes[i, j]
            d = variograms[
                (variograms.region == ct)
                & (variograms.method == m)
                & variograms["within plotted x range"]
            ]
            for k, n in enumerate(sizes[ct]):
                dn = d[d["gene set size"] == n]
                yv = dn["semivariance / variance"] if normalised else dn["semivariance"]
                ax.plot(
                    dn["spatial lag (px, bin upper edge)"],
                    yv,
                    "o-",
                    ms=1.8,
                    lw=0.7,
                    color=size_shades[k],
                    label=f"n = {n}" if k < 7 else "n = 500 / 434",
                )
                if not normalised:
                    ax.axhline(
                        dn["score variance (dashed line)"].iloc[0],
                        color=size_shades[k],
                        ls="--",
                        lw=0.4,
                        alpha=0.7,
                    )
            if normalised:
                ax.axhline(1, color="0.6", ls="--", lw=0.4)
            else:
                ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 3))
                ax.yaxis.get_offset_text().set_fontsize(5.5)
            ax.set_xlim(0, max_lag)
            ax.set_ylim(0, y_top if normalised else None)
            ax.tick_params(labelsize=5.5)
            if i == 0:
                ax.set_title(ct.replace("_", " "), fontsize=8)
            if i == len(methods) - 1:
                ax.set_xlabel("Spatial lag (px)", fontsize=6.5)
            if j == 0:
                unit = (
                    "Semivariance / variance"
                    if normalised
                    else f"Semivariance ({semivariance_unit[m]})"
                )
                ax.set_ylabel(f"{row_label[m]}\n{unit}", fontsize=5.5, color=colours[m])
            ax.grid(False)
            sns.despine(ax=ax)
    axes[0, 1].legend(
        frameon=False,
        fontsize=5.5,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        title="Number of genes",
        title_fontsize=6,
    )
    plt.tight_layout()
    return fig


def gss_spatial_folds(coords, n_blocks=25, n_folds=5, seed=0):
    """
    Spatial blocks: k-means clusters of the spot coordinates, assigned to folds in a random order.

    Parameters
    ----------
    coords : np.ndarray
        Spot coordinates (n x 2).
    n_blocks : int
        Number of k-means blocks.
    n_folds : int
        Number of folds (n_folds = n_blocks gives one label per block).
    seed : int
        Seed of k-means and of the block order.

    Returns
    -------
    np.ndarray
        Block (fold) label per spot.
    """
    km = KMeans(n_clusters=n_blocks, n_init=10, random_state=seed).fit(
        np.asarray(coords)
    )
    rng = np.random.default_rng(seed)
    block_fold = rng.permutation(np.arange(n_blocks) % n_folds)
    return block_fold[km.labels_]


def gss_threshold_rule_bootstrap(y, s, p_empirical, blocks, n_boot=1000, seed=7, min_pos=3):
    """F1 of each threshold rule on a spatial block bootstrap. The same blocks and seed give the same replicates
    for every gene set size; rules that depend on the data (oracle, top-k, Otsu, median) are refitted per replicate."""
    rng = np.random.default_rng(seed)
    y = np.asarray(y).astype(int)
    s, p_empirical = np.asarray(s, dtype=np.float64), np.asarray(p_empirical, dtype=np.float64)
    ub = np.unique(blocks)
    idx_by_block = {b: np.flatnonzero(blocks == b) for b in ub}
    records = []
    for r in range(n_boot):
        pick = rng.choice(ub, size=len(ub), replace=True)
        idx = np.concatenate([idx_by_block[b] for b in pick])
        yy = y[idx]
        if yy.sum() < min_pos or yy.sum() == len(yy):
            continue
        ss = s[idx]
        calls = {
            "Top-k": gss_topk_calls(yy, ss),
            "Otsu": ss > gss_otsu_threshold(ss),
            "Median": ss > np.median(ss),
            "z > 0": ss > 0,
            "z > 1.96": ss > 1.96,
            "p_emp ≤ 0.05": p_empirical[idx] <= 0.05,
        }
        records.append((r, "Oracle", _optimal_f1(yy, ss)[0]))
        records.extend((r, rule, f1_of_calls(yy, pred)) for rule, pred in calls.items())
    return pd.DataFrame(records, columns=["replicate", "rule", "F1"])


def gss_paired_block_bootstrap(y, score_dict, blocks, n_boot=1000, seed=7, min_pos=3):
    """
    Resample spatial blocks with replacement and evaluate every score on the same resample (paired).

    Parameters
    ----------
    y : np.ndarray
        Binary region labels.
    score_dict : dict
        (method, gene set size) -> score per spot.
    blocks : np.ndarray
        Block label per spot.
    n_boot : int
        Number of replicates.
    seed : int
        Seed of the block resampling.
    min_pos : int
        Replicates with fewer region spots are skipped.

    Returns
    -------
    pd.DataFrame
        replicate, method, gene set size, metric, value.
    """
    rng = np.random.default_rng(seed)
    y = np.asarray(y).astype(int)
    ub = np.unique(blocks)
    idx_by_block = {b: np.flatnonzero(blocks == b) for b in ub}
    records = []
    for r in range(n_boot):
        pick = rng.choice(ub, size=len(ub), replace=True)
        idx = np.concatenate([idx_by_block[b] for b in pick])
        yy = y[idx]
        if yy.sum() < min_pos or yy.sum() == len(yy):
            continue
        for (m, n), s in score_dict.items():
            ss = s[idx]
            records.append((r, m, n, "AUROC", roc_auc_score(yy, ss)))
            records.append((r, m, n, "AUPRC", average_precision_score(yy, ss)))
    return pd.DataFrame(
        records, columns=["replicate", "method", "gene set size", "metric", "value"]
    )


def gss_wide_auc(auc_table, ct):
    """
    AUROC/AUPRC with 95% intervals of one region, one row per method and size (gss_pair_panel columns).

    Parameters
    ----------
    auc_table : pd.DataFrame
        Long table with 'region', 'method', 'gene set size', 'metric', 'value', 'lo', 'hi'.
    ct : str
        Region.

    Returns
    -------
    pd.DataFrame
        method, gene set size, AUROC, AUROC_lo, AUROC_hi, AUPRC, AUPRC_lo, AUPRC_hi.
    """
    d = auc_table[auc_table.region == ct]
    out = d.pivot_table(
        index=["method", "gene set size"],
        columns="metric",
        values=["value", "lo", "hi"],
    )
    t = pd.DataFrame(index=out.index)
    for met in ["AUROC", "AUPRC"]:
        t[met], t[f"{met}_lo"], t[f"{met}_hi"] = (
            out[("value", met)],
            out[("lo", met)],
            out[("hi", met)],
        )
    return t.reset_index()


def gss_otsu_threshold(s, nbins=256):
    """
    Label-free histogram threshold (Otsu): maximises the between-class variance of the scores.

    Parameters
    ----------
    s : np.ndarray
        Score per spot.
    nbins : int
        Number of histogram bins.

    Returns
    -------
    float
        Threshold (bin midpoint).
    """
    s = np.asarray(s, dtype=np.float64)
    hist, edges = np.histogram(s, bins=nbins)
    mids = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * mids) / np.maximum(w0, 1)
    m1 = (np.sum(hist * mids) - np.cumsum(hist * mids)) / np.maximum(w1, 1)
    between = w0 * w1 * (m0 - m1) ** 2
    return float(mids[int(np.argmax(between))])


def gss_topk_calls(y, s):
    """
    Prevalence-matched rule: the k highest-scoring spots, k = number of region spots.

    Parameters
    ----------
    y : np.ndarray
        Binary region labels (only their sum is used).
    s : np.ndarray
        Score per spot.

    Returns
    -------
    np.ndarray
        Boolean calls (ties broken by spot order).
    """
    order = np.argsort(-np.asarray(s), kind="stable")
    pred = np.zeros(len(s), bool)
    pred[order[: int(np.sum(y))]] = True
    return pred


def gss_fixed_rules_axis(ax, rule_means, ct, fixed_rules, fixed_labels, colours):
    """
    Mean F1 over the gene set sizes of the fixed native rules (whiskers = range over sizes); a vertical tick marks
    the same method's mean oracle F1.

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    rule_means : pd.DataFrame
        'region', 'method', 'rule', 'mean', 'min', 'max' of F1 over the gene set sizes.
    ct : str
        Region.
    fixed_rules : list of tuple
        (method, rule) pairs, one row each.
    fixed_labels : dict
        Rule -> tick label.
    colours : dict
        Method -> colour.
    """
    d = rule_means[rule_means.region == ct].set_index(["method", "rule"])
    for k, (m, rule) in enumerate(fixed_rules):
        r = d.loc[(m, rule)]
        ax.errorbar(
            r["mean"],
            k,
            xerr=[[r["mean"] - r["min"]], [r["max"] - r["mean"]]],
            fmt="o",
            ms=4,
            color=colours[m],
            mec="white",
            mew=0.4,
            elinewidth=0.6,
            capsize=0,
        )
        oracle = d.loc[(m, "Oracle"), "mean"]
        ax.plot([oracle, oracle], [k - 0.25, k + 0.25], color=colours[m], lw=0.9)
    ax.set_yticks(range(len(fixed_rules)))
    ax.set_yticklabels([fixed_labels[r] for _, r in fixed_rules], fontsize=6)
    for lab, (m, _) in zip(ax.get_yticklabels(), fixed_rules):
        lab.set_color(colours[m])
    ax.set_ylim(len(fixed_rules) - 0.5, -0.5)
    ax.set_xlim(0, 1)
    ax.set_xlabel("F1, mean over sizes (range; tick: oracle F1)", fontsize=6.5)
    ax.tick_params(axis="x", labelsize=6)
    ax.grid(False)
    sns.despine(ax=ax)


def gss_label_line_ends(
    ax, x_ends, y_ends, labels, colour, x_text, y_max=1.0, min_gap=0.05, fontsize=5
):
    """
    Text labels to the right of a set of lines, spread vertically so that they do not overlap; a thin connector
    joins each label to the last point of its line.

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    x_ends, y_ends : list of float
        Last point of each line.
    labels : list of str
        One label per line.
    colour : str
        Text and connector colour.
    x_text : float
        x of the labels (data units, right of every line end).
    y_max : float
        Highest label position (data units).
    min_gap : float
        Minimum vertical distance between labels (data units).
    fontsize : float
        Label font size.
    """
    if not len(labels):
        return
    order = np.argsort(y_ends)
    placed = np.asarray(y_ends, dtype=float)[order].copy()
    for i in range(1, len(placed)):
        placed[i] = max(placed[i], placed[i - 1] + min_gap)
    placed -= max(0.0, placed[-1] - y_max)
    for i in range(len(placed) - 2, -1, -1):
        placed[i] = min(placed[i], placed[i + 1] - min_gap)
    for k, yl in zip(order, placed):
        ax.plot(
            [x_ends[k] + 0.1, x_text - 0.05],
            [y_ends[k], yl],
            color=colour,
            lw=0.3,
            alpha=0.6,
            clip_on=False,
        )
        ax.text(
            x_text,
            yl,
            labels[k],
            color=colour,
            fontsize=fontsize,
            ha="left",
            va="center",
        )


def gss_log10_p(r, n):
    """
    log10 of the two-sided p-value of a correlation coefficient (exact Student t tail, arbitrary precision).

    Parameters
    ----------
    r : float
        Pearson or Spearman correlation.
    n : int
        Number of observations.

    Returns
    -------
    float
        log10(p); finite where the double-precision p-value underflows to 0.
    """
    df = n - 2
    t2 = r**2 * df / (1 - r**2)
    return float(
        mpmath.log10(mpmath.betainc(df / 2, 0.5, 0, df / (df + t2), regularized=True))
    )


def gss_p_label(p):
    """
    p-value as printed in the published panels ('p < 0.0001' below 1e-4).

    Parameters
    ----------
    p : float
        p-value.

    Returns
    -------
    str
    """
    return "p < 0.0001" if p < 1e-4 else f"p = {p:.2g}"


def gss_draw_cross_moran(
    cm_points_ct,
    cm_stats,
    ct,
    methods,
    sizes_ct,
    colours,
    quadrant_colours,
    axis_unit,
    row_label,
):
    """
    Sup. Fig. 6a/b: one row per method, one column per gene set size (score against the spatial lag of the n = 50
    score of the same method; spots coloured by quadrant).

    Parameters
    ----------
    cm_points_ct : pd.DataFrame
        Spots of one region: 'method', 'gene set size', 'score', 'spatial lag of n = 50 score', 'quadrant'.
    cm_stats : pd.DataFrame
        'region', 'method', 'gene set size', 'pearson r', 'p value', 'ols slope', 'ols intercept'.
    ct : str
        Region.
    methods : list of str
        Methods (rows).
    sizes_ct : list of int
        Gene set sizes of the region (columns).
    colours : dict
        Method -> colour (row labels).
    quadrant_colours : dict
        Quadrant -> colour.
    axis_unit : dict
        Method -> unit of the score axis.
    row_label : dict
        Method -> row label.

    Returns
    -------
    plt.Figure
    """
    fig, axes = plt.subplots(
        len(methods), len(sizes_ct), figsize=(13.0, 1.75 * len(methods)),
        squeeze=False,
    )
    for i, m in enumerate(methods):
        d = cm_points_ct[cm_points_ct.method == m]
        for k, n in enumerate(sizes_ct):
            ax = axes[i, k]
            dn = d[d["gene set size"] == n]
            for q, colour in quadrant_colours.items():
                dq = dn[dn.quadrant == q]
                ax.scatter(
                    dq["score"],
                    dq["spatial lag of n = 50 score"],
                    s=1.5,
                    color=colour,
                    alpha=0.5,
                    linewidths=0,
                )
            st = cm_stats[
                (cm_stats.region == ct)
                & (cm_stats.method == m)
                & (cm_stats["gene set size"] == n)
            ].iloc[0]
            xs = np.array([dn["score"].min(), dn["score"].max()])
            ax.plot(
                xs, st["ols intercept"] + st["ols slope"] * xs, color="black", lw=0.7
            )
            ax.axhline(0, color="0.6", lw=0.4, ls="--")
            ax.axvline(0, color="0.6", lw=0.4, ls="--")
            ax.text(
                0.97,
                0.04,
                f"r = {st['pearson r']:.2f}\n{gss_p_label(st['p value'])}",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=5.5,
            )
            ax.tick_params(labelsize=5, length=2, pad=1)
            ax.ticklabel_format(style="sci", scilimits=(-2, 3))
            ax.xaxis.get_offset_text().set_fontsize(5)
            ax.yaxis.get_offset_text().set_fontsize(5)
            if i == 0:
                ax.set_title(f"n = {n}", fontsize=7)
            ax.set_xlabel(f"{axis_unit[m]}, n = {n}", fontsize=5.5, labelpad=1)
            if k == 0:
                ax.set_ylabel(
                    f"{row_label[m]}\nlag of n = 50", fontsize=5.5, color=colours[m]
                )
            ax.grid(False)
            sns.despine(ax=ax)
    fig.suptitle(
        f"{ct.replace('_', ' ')}: cross-Moran scatter (score at each size vs spatial lag of the n = 50 score)",
        fontsize=8,
    )
    handles = [
        Line2D([], [], ls="", marker="o", ms=4, color=c, label=q)
        for q, c in reversed(quadrant_colours.items())
    ]
    axes[0, -1].legend(
        handles=handles,
        frameon=False,
        fontsize=6,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
    )
    plt.tight_layout()
    return fig


def gss_draw_spearman(
    sp_points_ct, sp_stats, ct, methods, sizes_ct, colours, axis_unit, row_label
):
    """
    Sup. Fig. 6c/d: one row per method, one column per gene set size (score of n = 50 on x).

    Parameters
    ----------
    sp_points_ct : pd.DataFrame
        Spots of one region: 'method', 'gene set size', 'score n = 50', 'score n'.
    sp_stats : pd.DataFrame
        'region', 'method', 'gene set size', 'spearman rho', 'p value', 'log10 p', 'ols slope', 'ols intercept'.
    ct : str
        Region.
    methods : list of str
        Methods (rows).
    sizes_ct : list of int
        Gene set sizes of the region (n = 50 is the reference and is left out).
    colours : dict
        Method -> colour.
    axis_unit : dict
        Method -> unit of the score axes.
    row_label : dict
        Method -> row label.

    Returns
    -------
    plt.Figure
    """
    n_other = [k for k in sizes_ct if k != 50]
    fig, axes = plt.subplots(
        len(methods), len(n_other), figsize=(11.6, 1.75 * len(methods)),
        squeeze=False,
    )
    for i, m in enumerate(methods):
        d = sp_points_ct[sp_points_ct.method == m]
        for k, n in enumerate(n_other):
            ax = axes[i, k]
            dn = d[d["gene set size"] == n]
            ax.scatter(
                dn["score n = 50"],
                dn["score n"],
                s=1.5,
                color=colours[m],
                alpha=0.4,
                linewidths=0,
            )
            st = sp_stats[
                (sp_stats.region == ct)
                & (sp_stats.method == m)
                & (sp_stats["gene set size"] == n)
            ].iloc[0]
            xs = np.array([dn["score n = 50"].min(), dn["score n = 50"].max()])
            ax.plot(
                xs, st["ols intercept"] + st["ols slope"] * xs, color="black", lw=0.7
            )
            p_text = (
                f"p = {st['p value']:.1e}"
                if st["p value"] > 0
                else f"log10 p = {st['log10 p']:.0f}"
            )
            ax.text(
                0.04,
                0.97,
                f"ρ = {st['spearman rho']:.4f}\n{p_text}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=5.5,
            )
            ax.tick_params(labelsize=5, length=2, pad=1)
            ax.ticklabel_format(style="sci", scilimits=(-2, 3))
            ax.xaxis.get_offset_text().set_fontsize(5)
            ax.yaxis.get_offset_text().set_fontsize(5)
            if i == 0:
                ax.set_title(f"n = {n}", fontsize=7)
            ax.set_xlabel(f"{axis_unit[m]}, n = 50", fontsize=5.5, labelpad=1)
            if k == 0:
                ax.set_ylabel(
                    f"{row_label[m]}\n{axis_unit[m]} at n",
                    fontsize=5.5,
                    color=colours[m],
                )
            ax.grid(False)
            sns.despine(ax=ax)
    fig.suptitle(
        f"Spearman correlation with the n = 50 score ({ct.replace('_', ' ')})",
        fontsize=8,
    )
    plt.tight_layout()
    return fig


def gss_fetch_msigdb(annotations_url="https://omnipathdb.org/annotations"):
    """
    Download the MSigDB annotations from OmniPath (academic licence).

    Parameters
    ----------
    annotations_url : str
        OmniPath annotations endpoint.

    Returns
    -------
    pd.DataFrame
        Long annotation table (record_id, label, value, genesymbol, ...).
    """
    url = (
        annotations_url
        + "?"
        + urllib.parse.urlencode({"databases": "MSigDB", "license": "academic"})
    )
    return pd.read_csv(urllib.request.urlopen(url), sep="\t", low_memory=False)


def gss_parse_msigdb(df):
    """
    One row per (collection, gene set, gene) from the OmniPath MSigDB annotations.

    Parameters
    ----------
    df : pd.DataFrame
        Output of gss_fetch_msigdb.

    Returns
    -------
    pd.DataFrame
        collection, geneset, genesymbol.
    """
    df = df[df["label"].isin(["collection", "geneset", "genesymbol"])]
    pivot = df.pivot_table(
        index=["record_id", "genesymbol"],
        columns="label",
        values="value",
        aggfunc="first",
    ).reset_index()
    pivot = pivot.dropna(subset=["collection", "geneset", "genesymbol"])
    return pivot[["collection", "geneset", "genesymbol"]]


def gss_msigdb_gene_set_sizes(pivot):
    """
    Number of distinct genes of every MSigDB gene set.

    Parameters
    ----------
    pivot : pd.DataFrame
        Output of gss_parse_msigdb.

    Returns
    -------
    pd.DataFrame
        collection, geneset, gene_set_size.
    """
    return (
        pivot.groupby(["collection", "geneset"])["genesymbol"]
        .nunique()
        .reset_index(name="gene_set_size")
    )


# --------------------------------------------------------------------------------------------------
# REV_SupFig5_6_7_runtime
# --------------------------------------------------------------------------------------------------


def runtime_time_method(
    adata,
    gene_expr,
    method,
    cell_type,
    gene_set,
    variants,
    enrichmap_variants,
    threads=2,
):
    """
    Wall-clock time of one scoring call, measured as in 06_performance.

    Parameters
    ----------
    adata : AnnData
        Log-normalised data (copied inside the timed block for the AnnData methods).
    gene_expr : pd.DataFrame
        Genes (upper case) x spots expression of `adata`, built outside the timing.
    method : str
        "EnrichMap (null)", "EnrichMap (smoothing-only)", "AUCell", "Scanpy", "ssGSEA", "GSVA" or "Z score".
    cell_type : str
        Region (signature name).
    gene_set : list
        Signature genes.
    variants : dict
        EnrichMap variant -> em.tl.score keyword arguments.
    enrichmap_variants : dict
        EnrichMap method name -> variant.
    threads : int
        gseapy threads (ssGSEA, GSVA).

    Returns
    -------
    float
        Seconds.
    """
    # prepare signatures for AUCell, ssGSEA, GSVA
    signatures_dict = {cell_type: gene_set}
    signatures = pd.DataFrame(
        [(key, gene) for key, genes in signatures_dict.items() for gene in genes],
        columns=["source", "target"],
    )
    signatures_dict_upper = {
        k: [g.upper() for g in v] for k, v in signatures_dict.items()
    }

    if method in enrichmap_variants:
        variant = enrichmap_variants[method]
        t1 = time.time()
        adata_temp = adata.copy()
        em.tl.score(
            adata_temp,
            gene_set=gene_set,
            score_key=f"{cell_type}_{variant}",
            smoothing=True,
            batch_key=None,
            n_neighbors=6,
            **variants[variant],
        )
        t2 = time.time()
    elif method == "AUCell":
        t1 = time.time()
        adata_temp = adata.copy()
        dc.mt.aucell(adata_temp, net=signatures, tmin=0)
        score_aucell = dc.pp.get_obsm(adata_temp, key="score_aucell")
        adata_temp.obs["aucell_score"] = score_aucell.obsm["score_aucell"]
        t2 = time.time()
    elif method == "Scanpy":
        t1 = time.time()
        adata_temp = adata.copy()
        sc.tl.score_genes(
            adata_temp, gene_set, score_name="scanpy_score", use_raw=False
        )
        t2 = time.time()
    elif method == "ssGSEA":
        t1 = time.time()
        ssgsea = gp.ssgsea(
            data=gene_expr,
            gene_sets=signatures_dict_upper,
            outdir=None,
            sample_norm_method="rank",
            no_plot=True,
            min_size=2,
            threads=threads,
        )
        nes = ssgsea.res2d.pivot(index="Term", columns="Name", values="NES")
        ssgsea_score = nes.loc[cell_type].T.astype(float).reindex(adata.obs_names)
        t2 = time.time()
    elif method == "GSVA":
        t1 = time.time()
        gsva = gp.gsva(
            data=gene_expr,
            gene_sets=signatures_dict_upper,
            outdir=None,
            min_size=2,
            threads=threads,
        )
        es = gsva.res2d.pivot(index="Term", columns="Name", values="ES")
        gsva_score = es.loc[cell_type].T.astype(float).reindex(adata.obs_names)
        t2 = time.time()
    elif method == "Z score":
        gene_set_upper = signatures_dict_upper[cell_type]
        t1 = time.time()
        z_score = gene_expr.loc[gene_set_upper].sum(axis=0) / np.sqrt(
            len(gene_set_upper)
        )
        t2 = time.time()
    return t2 - t1


def runtime_vs_size_axis(
    ax,
    stats,
    published,
    ct,
    sizes_ct,
    methods,
    published_methods,
    colours,
    max_size=None,
    log_x=True,
):
    """
    Runtime against gene set size for one region: mean of the three repeats on this server (solid lines, min-max
    whiskers) and the published single run (dashed lines); every point is a filled circle in the method colour.

    Parameters
    ----------
    ax : plt.Axes
        Target axes.
    stats : pd.DataFrame
        'region', 'method', 'gene set size', 'mean', 'min', 'max' of the repeats on this server.
    published : pd.DataFrame
        'region', 'method', 'gene set size', 'published time (s)'.
    ct : str
        Region.
    sizes_ct : list of int
        Gene set sizes of the region.
    methods : list of str
        Methods.
    published_methods : list of str
        Methods with a published run.
    colours : dict
        Method -> colour.
    max_size : int or None
        Show only gene set sizes up to this value (the published zoom panels: 50).
    log_x : bool
        Logarithmic x axis.
    """
    keep = [n for n in sizes_ct if max_size is None or n <= max_size]
    for m in methods:
        d = stats[
            (stats.region == ct)
            & (stats.method == m)
            & stats["gene set size"].isin(keep)
        ]
        if not d.empty:
            ax.errorbar(
                d["gene set size"],
                d["mean"],
                yerr=[d["mean"] - d["min"], d["max"] - d["mean"]],
                fmt="o-",
                ms=2.6,
                lw=0.9,
                elinewidth=0.6,
                capsize=0,
                color=colours[m],
                label=m,
                mec="white",
                mew=0.3,
            )
        if m in published_methods:
            p = published[
                (published.region == ct)
                & (published.method == m)
                & published["gene set size"].isin(keep)
            ]
            ax.plot(
                p["gene set size"],
                p["published time (s)"],
                "o--",
                ms=2.6,
                lw=0.6,
                color=colours[m],
                mec="white",
                mew=0.3,
                alpha=0.8,
            )
    ax.set_yscale("log")
    if log_x:
        ax.set_xscale("log")
    ax.set_xticks(keep)
    ax.set_xticklabels([str(n) for n in keep], fontsize=6)
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_xlabel("Gene set size (genes)", fontsize=7)
    ax.set_ylabel("Time (seconds, log scale)", fontsize=7)
    ax.tick_params(axis="y", labelsize=6)
    ax.grid(False)
    sns.despine(ax=ax)


def runtime_method_legend(ax, methods, colours):
    """
    Legend: one entry per method plus the line styles (this server: solid, published: dashed).

    Parameters
    ----------
    ax : plt.Axes
        Axes the legend is attached to (placed to its right).
    methods : list of str
        Methods.
    colours : dict
        Method -> colour.
    """
    handles = [
        Line2D([], [], marker="o", ms=3, lw=0.9, color=colours[m], label=m)
        for m in methods
    ]
    handles += [
        Line2D(
            [],
            [],
            marker="o",
            ms=3,
            lw=0.9,
            color="0.3",
            label="this server: mean of 3 repeats (range)",
        ),
        Line2D(
            [],
            [],
            marker="o",
            ms=3,
            lw=0.6,
            ls="--",
            color="0.3",
            label="published: one run",
        ),
    ]
    ax.legend(
        handles=handles,
        frameon=False,
        fontsize=5.5,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
    )


def runtime_plot_summary(table, ax, colours, lower_is_better=()):
    """
    Methods x metrics summary: circles sized and shaded by each metric scaled to 0-1 across methods (best = 1),
    raw values annotated and the mean rank on the right.

    Parameters
    ----------
    table : pd.DataFrame
        Index = methods, columns = metrics.
    ax : plt.Axes
        Target axes.
    colours : dict
        Method -> colour of the tick label.
    lower_is_better : tuple of str
        Metrics where smaller is better (here every runtime column).
    """
    scaled = table.astype(float).copy()
    for c in table.columns:
        col = table[c].astype(float)
        s = (
            (col - col.min()) / (col.max() - col.min())
            if col.max() > col.min()
            else col * 0 + 1
        )
        scaled[c] = 1 - s if c in lower_is_better else s
    ranks = scaled.rank(ascending=False).mean(axis=1)
    cmap = plt.get_cmap("Purples")
    for i, m in enumerate(table.index):
        for j, c in enumerate(table.columns):
            s = scaled.loc[m, c]
            ax.scatter(
                j, i, s=30 + 220 * s, color=cmap(0.25 + 0.75 * s), edgecolors="none"
            )
            ax.text(
                j,
                i,
                f"{table.loc[m, c]:.2g}",
                ha="center",
                va="center",
                fontsize=5,
                color="white" if s > 0.6 else "black",
            )
        ax.text(
            len(table.columns) + 0.2,
            i,
            f"{ranks[m]:.1f}",
            ha="left",
            va="center",
            fontsize=6,
        )
    ax.set_xticks(range(len(table.columns)))
    ax.set_xticklabels(table.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(table.index)))
    ax.set_yticklabels(table.index)
    for lab in ax.get_yticklabels():
        lab.set_color(colours.get(lab.get_text(), "black"))
    ax.text(len(table.columns) + 0.2, -0.9, "mean rank", fontsize=6, ha="left")
    ax.set_xlim(-0.6, len(table.columns) + 0.9)
    ax.set_ylim(len(table.index) - 0.4, -1.2)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)


# --------------------------------------------------------------------------------------------------
# REV_GSD_base_visiumhd, REV_GSD_base_xenium
# --------------------------------------------------------------------------------------------------


def spatial_grid_subsample(xy, n, n_bins=50, seed=0):
    """
    Spatially stratified subsample: the section is cut into an n_bins x n_bins grid and every occupied bin keeps
    cells in proportion to its size (largest-remainder rounding), drawn without replacement.

    Parameters
    ----------
    xy : np.ndarray
        Cell coordinates (n_cells x 2).
    n : int
        Number of cells to keep.
    n_bins : int
        Grid bins per axis.
    seed : int
        Random seed.

    Returns
    -------
    np.ndarray
        Sorted integer indices of the kept cells.
    """
    rng = np.random.default_rng(seed)
    span = np.ptp(xy, axis=0)
    ij = np.clip(((xy - xy.min(axis=0)) / span * n_bins).astype(int), 0, n_bins - 1)
    b = ij[:, 0] * n_bins + ij[:, 1]
    bins, inv, counts = np.unique(b, return_inverse=True, return_counts=True)
    quota = counts * n / len(xy)
    k = np.floor(quota).astype(int)
    order = np.argsort(-(quota - k), kind="stable")
    k[order[: n - k.sum()]] += 1
    members = np.split(np.argsort(inv, kind="stable"), np.cumsum(counts)[:-1])
    keep = [rng.choice(m, kj, replace=False) for m, kj in zip(members, k) if kj > 0]
    return np.sort(np.concatenate(keep))


def write_gsdensity_cache(
    adata,
    gsd,
    gene_sets,
    out_dir,
    index_path,
    dataset,
    params,
    note,
    n_total,
    working_dir,
    extra_columns=None,
):
    """
    Write one GSDensity CSV per gene set and the dataset index.

    Per gene set: spot_id (adata.obs_names), any `extra_columns`, relevance (RWR probability; sums to 1 over the
    spots scored), label (compute.cell.label antimode call) and relevance_scaled (relevance x n spots; mean 1).

    Parameters
    ----------
    adata : AnnData
        The spots or cells that were scored.
    gsd : dict
        Output of read_gsdensity_outputs (scores, labels, spatial_kld, timing).
    gene_sets : dict
        Gene-set name -> genes given.
    out_dir : str
        Folder of the per-gene-set CSVs.
    index_path : str
        Path of the index CSV (one row per gene set).
    dataset : str
        Dataset name written to the index.
    params : dict
        GSDensity parameters written to the index.
    note : str
        Subsample note written to the index.
    n_total : int
        Spots or cells in the whole dataset.
    working_dir : str
        Paths in the index are relative to it.
    extra_columns : dict or None
        Column name -> values (in adata order) written after spot_id.

    Returns
    -------
    pd.DataFrame
        The index.
    """
    os.makedirs(out_dir, exist_ok=True)
    spots = adata.obs_names.astype(str)
    scores, labels = gsd["scores"].copy(), gsd["labels"].copy()
    scores.index, labels.index = scores.index.astype(str), labels.index.astype(str)
    scores, labels = scores.reindex(spots), labels.reindex(spots)
    sk = gsd["spatial_kld"]
    var = set(adata.var_names)
    rows = []
    for s, genes in gene_sets.items():
        rel = scores[s].values.astype(float)
        df = pd.DataFrame({"spot_id": spots})
        for k, v in (extra_columns or {}).items():
            df[k] = v
        df["relevance"] = rel
        df["label"] = labels[s].values
        df["relevance_scaled"] = rel * len(rel)
        path = out_dir + f"{s}.csv"
        df.to_csv(path, index=False)
        used = [g for g in genes if g in var]
        rows.append(
            {
                "dataset": dataset,
                "signature": s,
                "genes_used": len(used),
                "genes_used_of_given": f"{len(used)}/{len(genes)}",
                "genes_given": len(genes),
                "genes_used_list": ";".join(used),
                "genes_missing": ";".join(g for g in genes if g not in var),
                "n_spots": adata.n_obs,
                "n_spots_in_dataset": n_total,
                "n_positive_label": int((df["label"] == "positive").sum()),
                "n_label_NA": int(df["label"].isna().sum()),
                "n_zero_relevance": int((rel == 0).sum()),
                "subsample_note": note,
                "runtime_s": round(float(gsd["timing"]["total"]), 1),
                "params": json.dumps(params),
                "spatial_kld": float(sk.loc[s, "spatial.kld"])
                if s in sk.index
                else np.nan,
                "delta_kld": float(sk.loc[s, "delta.kld"]) if s in sk.index else np.nan,
                "spatial_kld_p": float(sk.loc[s, "p.value"])
                if s in sk.index
                else np.nan,
                "file": os.path.relpath(path, working_dir),
            }
        )
    index = pd.DataFrame(rows)
    index.to_csv(index_path, index=False)
    return index


# --------------------------------------------------------------------------------------------------
# REV_SupFig3_4_controls
# --------------------------------------------------------------------------------------------------


def plot_maps_shared_limits(adata, keys_by_method, reference, vcenter=0, size=20, panel_size=(2.4, 2.1), plot_fdr=()):
    """
    Spatial maps, one row per method, whose panels share the colour limits of the method's reference map.

    Parameters
    ----------
    adata : AnnData
        Object with the plotted columns in obs.
    keys_by_method : dict
        Method -> obs columns of its row.
    reference : dict
        Method -> obs column that sets the limits of the row: symmetric around `vcenter` (vcenter +/- the largest
        distance from it) or, with vcenter=None, its minimum and maximum. With None, each panel sets its own limits.
    vcenter : float or None
        Centre of the colour scale.
    size : float
        Spot size.
    panel_size : tuple
        Width and height of one panel (inches).
    plot_fdr : collection of str
        Methods whose maps get a black stroke around spots with q < 0.05 (obs column "<key>_q").

    Returns
    -------
    plt.Figure
    """
    from matplotlib.colors import CenteredNorm, Normalize

    n_cols = max(len(keys) for keys in keys_by_method.values())
    fig, axes = plt.subplots(
        len(keys_by_method),
        n_cols,
        figsize=(panel_size[0] * n_cols, panel_size[1] * len(keys_by_method)),
        squeeze=False,
    )
    for row, (method, keys) in zip(axes, keys_by_method.items()):
        centre = vcenter[method] if isinstance(vcenter, dict) else vcenter
        panels = [(keys, reference[method], list(row[: len(keys)]) if len(keys) > 1 else row[0])]
        if reference[method] is None:
            panels = [([key], key, ax) for key, ax in zip(keys, row)]
        for panel_keys, ref_key, panel_ax in panels:
            ref = np.asarray(adata.obs[ref_key], dtype=float)
            if centre is None:
                norm = Normalize(vmin=np.nanmin(ref), vmax=np.nanmax(ref))
            else:
                norm = CenteredNorm(vcenter=centre, halfrange=np.nanmax(np.abs(ref - centre)))
            em.pl.spatial_enrichmap(
                adata,
                score_key=panel_keys,
                vcenter=None,
                norm=norm,
                ncols=n_cols,
                shape=None,
                size=size,
                plot_fdr=method in plot_fdr,
                fig=fig,
                ax=panel_ax,
            )
    return fig
