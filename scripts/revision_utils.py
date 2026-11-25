"""
Shared utilities for the Reviewer #4 revision analyses (EnrichMap, NG-A69290R1).

Every revision notebook imports from here so that scorers, spatial operators and
metrics are identical across analyses. Nothing in the original project notebooks
is modified; the original scoring/benchmark code is re-implemented verbatim where
noted so that numbers can be reproduced and then extended.

Design
------
EnrichMap's pipeline is a composition of linear operators acting on a raw score
vector r (manuscript Eq. 1-5):

    smoothed   = S r                 S = D^-1 A (row-normalised 6-NN spatial graph)
    corrected  = (I - H) S r         H = B (B'B + P)^-1 B'  (pyGAM te(0,1) hat matrix)

Because S and H depend only on spot coordinates, the same operators can be
applied to the raw score of *any* method (AUCell, GSVA, ...). This lets us run a
full factorial ablation {scorer} x {smoothing on/off} x {GAM on/off} with a
single raw scoring per method.
"""
from __future__ import annotations

import os
import time
import warnings

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import squidpy as sq
from scipy.sparse import csr_matrix, issparse
from scipy.sparse.csgraph import shortest_path
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score
from pygam import LinearGAM, s, te

warnings.filterwarnings("ignore")

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/"  # repository root
REV = PROJECT
RES = REV + "results/"
FIG = REV + "figures/"
DATA = REV + "data/"
for _d in (RES, FIG, DATA):
    os.makedirs(_d, exist_ok=True)

METHOD_ORDER = ["EnrichMap", "z-score", "scanpy", "AUCell", "GSVA", "ssGSEA", "GSDensity"]


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def read_h5ad_no_uns(path, keep_uns=("spatial",)):
    """Read an h5ad written by a newer anndata (whose 'null' uns encoding anndata
    0.10 cannot parse) by reading every slot except uns."""
    import h5py
    from anndata.experimental import read_elem

    try:
        return ad.read_h5ad(path)
    except Exception:
        pass
    with h5py.File(path, "r") as f:
        kw = {k: read_elem(f[k]) for k in ("X", "obs", "var", "obsm", "obsp", "layers", "varm") if k in f}
        uns = {}
        for k in keep_uns:
            if "uns" in f and k in f["uns"]:
                try:
                    uns[k] = read_elem(f["uns"][k])
                except Exception:
                    pass
    return ad.AnnData(uns=uns, **kw)


def load_brain():
    """Mouse brain Visium (squidpy visium_hne) as preprocessed for the manuscript, plus the
    top-10 pyramidal layer / striatum signatures used in Fig. 2 and Sup. Fig. 2."""
    import pickle

    adata = read_h5ad_no_uns(PROJECT + "processed_data/visium_brain_adata.h5ad")
    keep = [c for c in adata.obs.columns if not c.endswith("_score")]
    adata.obs = adata.obs[keep].copy()
    sigs = {}
    for ct in ("Pyramidal_layer", "Striatum"):
        with open(PROJECT + f"{ct}_signatures.pkl", "rb") as f:
            sigs[ct] = pickle.load(f)
    return adata, sigs


def dense(X):
    return X.toarray() if issparse(X) else np.asarray(X)


# ---------------------------------------------------------------------------
# Spatial operators
# ---------------------------------------------------------------------------

def spatial_conn(coords, n_neighbors=6):
    """Binary spatial kNN connectivity exactly as built inside em.tl.score()."""
    tmp = ad.AnnData(X=np.zeros((len(coords), 1), dtype=np.float32))
    tmp.obsm["spatial"] = np.asarray(coords)
    sq.gr.spatial_neighbors(tmp, n_neighs=n_neighbors, coord_type="generic", key_added="spatial")
    return csr_matrix(tmp.obsp["spatial_connectivities"])


def smooth(Y, conn, n_iter=1):
    """EnrichMap Eq. 3: neighbourhood average (row-normalised adjacency), optionally iterated."""
    Y = np.asarray(Y, dtype=np.float64)
    d = np.maximum(np.asarray(conn.sum(axis=1)).ravel(), 1e-10)
    out = Y
    for _ in range(n_iter):
        out = conn.dot(out) / (d[:, None] if out.ndim == 2 else d)
    return out


class SpatialGAM:
    """The spatial-covariate GAM of EnrichMap written as a linear smoother.

    fitted = H y, residual = (I - H) y with H = B (B'B + P)^-1 B'. With
    basis='te', n_splines=None and lam=None this is identical to
    LinearGAM(te(0, 1)) as called by em.tl.score (pyGAM defaults: 10 x 10 cubic
    B-splines, lam = 0.6 per margin, intercept). basis='additive' gives
    s(0) + s(1) (pyGAM default 20 splines per margin). `edof` is trace(H).
    """

    def __init__(self, coords, basis="te", n_splines=None, lam=None):
        coords = np.asarray(coords, dtype=np.float64)
        kw = {}
        if n_splines is not None:
            kw["n_splines"] = n_splines
        if lam is not None:
            kw["lam"] = lam
        if basis == "te":
            term = te(0, 1, **kw)
        elif basis == "additive":
            term = s(0, **kw) + s(1, **kw)
        else:
            raise ValueError(basis)
        gam = LinearGAM(term).fit(coords, coords[:, 0])
        B = gam._modelmat(coords)
        B = B.toarray() if hasattr(B, "toarray") else np.asarray(B)
        P = gam._P()
        P = P.toarray() if hasattr(P, "toarray") else np.asarray(P)
        self.B, self.P = B, P
        # te(0, 1) is a partition of unity, so with the intercept B is rank-deficient and B'B + P is
        # singular along one direction; add the same sqrt(machine-eps) ridge pyGAM uses before solving.
        # Fitted values were already correct (~1e-9); this makes `edof` match pyGAM's statistics_['edof'].
        self._A = B.T @ B + P + np.sqrt(np.finfo(float).eps) * np.eye(B.shape[1])
        self.edof = float(np.trace(np.linalg.solve(self._A, B.T @ B)))
        self.basis, self.n_splines, self.lam = basis, n_splines, lam
        self.terms = str(gam.terms)

    def fitted(self, Y):
        Y = np.asarray(Y, dtype=np.float64)
        return self.B @ np.linalg.solve(self._A, self.B.T @ Y)

    def residual(self, Y):
        return np.asarray(Y, dtype=np.float64) - self.fitted(Y)

    def r2(self, y):
        """Fraction of the (centred) variance of y that the GAM absorbs."""
        y = np.asarray(y, dtype=np.float64)
        yc = y - y.mean()
        return float(1 - ((y - self.fitted(y)) ** 2).sum() / (yc ** 2).sum())

    @property
    def label(self):
        if self.basis == "additive":
            return f"additive s(x)+s(y) (edf={self.edof:.0f})"
        ns = "10x10" if self.n_splines is None else f"{self.n_splines}x{self.n_splines}"
        lam = 0.6 if self.lam is None else self.lam
        return f"te {ns}, lam={lam:g} (edf={self.edof:.0f})"


def statsmodels_additive_gam_residual(coords, Y, df=10, degree=3):
    """GAM used in the first submission (statsmodels GLMGam, additive cubic B-splines,
    df=[10, 10], unpenalised) -- kept to test the 'reduces leakage' statement."""
    from statsmodels.gam.api import BSplines, GLMGam

    coords = np.asarray(coords, dtype=np.float64)
    bs = BSplines(coords, df=[df, df], degree=[degree, degree])
    Y = np.atleast_2d(np.asarray(Y, dtype=np.float64).T).T
    out = np.empty_like(Y)
    for j in range(Y.shape[1]):
        res = GLMGam.from_formula("y ~ 1", data={"y": Y[:, j]}, smoother=bs).fit()
        out[:, j] = Y[:, j] - res.fittedvalues
    return out.ravel() if out.shape[1] == 1 else out


def configurations(raw, conn, gam):
    """Factorial ablation of one raw score vector: {smoothing} x {GAM}."""
    raw = np.asarray(raw, dtype=np.float64)
    sm = smooth(raw, conn)
    return {
        "raw": raw,
        "smooth": sm,
        "GAM": gam.residual(raw),
        "smooth+GAM": gam.residual(sm),
    }


# ---------------------------------------------------------------------------
# Raw scorers (as in 17_benchmarking.ipynb / 05_01_simulated_benchmarking.py)
# ---------------------------------------------------------------------------

def raw_enrichmap(adata, genes, weights=None):
    """EnrichMap Eq. 1 (CV-weighted mean) before batch z-scoring, smoothing and GAM."""
    import enrichmap as em

    genes = [g for g in dict.fromkeys(genes) if g in adata.var_names]
    w = weights if weights is not None else em.tl.infer_gene_weights(adata, genes)
    wv = np.array([w.get(g, 1.0) for g in genes], dtype=np.float64)
    X = dense(adata[:, genes].X)
    return X @ wv / np.abs(wv).sum()


def raw_zscore(adata, genes):
    genes = [g for g in dict.fromkeys(genes) if g in adata.var_names]
    return dense(adata[:, genes].X).sum(axis=1) / np.sqrt(len(genes))


def raw_scanpy(adata, genes):
    a = ad.AnnData(X=adata.X, obs=pd.DataFrame(index=adata.obs_names), var=pd.DataFrame(index=adata.var_names))
    sc.tl.score_genes(a, gene_list=list(genes), score_name="_bench", use_raw=False)
    return a.obs["_bench"].values.astype(np.float64)


def raw_aucell(adata, genes):
    import decoupler as dc

    a = ad.AnnData(X=adata.X, obs=pd.DataFrame(index=adata.obs_names), var=pd.DataFrame(index=adata.var_names))
    net = pd.DataFrame({"source": "_bench", "target": list(genes)})
    dc.mt.aucell(a, net=net)
    return dc.pp.get_obsm(adata=a, key="score_aucell").obsm["score_aucell"]["_bench"].values.astype(np.float64)


def _expr_df(adata, upper=False):
    df = pd.DataFrame(dense(adata.X).T, index=adata.var_names, columns=adata.obs_names)
    if upper:
        df.index = df.index.str.upper()
    return df


def raw_gsva(adata, genes, threads=1):
    import gseapy as gp

    res = gp.gsva(data=_expr_df(adata), gene_sets={"_bench": list(genes)}, outdir=None, min_size=3, threads=threads)
    es = res.res2d.pivot(index="Term", columns="Name", values="ES")
    return es.loc["_bench"].T.astype(float).reindex(adata.obs_names).values


def raw_ssgsea(adata, genes, threads=1):
    import gseapy as gp

    res = gp.ssgsea(
        data=_expr_df(adata, upper=True), gene_sets={"_bench": [g.upper() for g in genes]}, outdir=None,
        sample_norm_method="rank", no_plot=True, min_size=3, threads=threads,
    )
    nes = res.res2d.pivot(index="Term", columns="Name", values="NES")
    return nes.loc["_bench"].T.astype(float).reindex(adata.obs_names).values


RAW_SCORERS = {
    "EnrichMap": raw_enrichmap,
    "z-score": raw_zscore,
    "scanpy": raw_scanpy,
    "AUCell": raw_aucell,
    "GSVA": raw_gsva,
    "ssGSEA": raw_ssgsea,
}
FAST_SCORERS = {k: RAW_SCORERS[k] for k in ("EnrichMap", "z-score", "scanpy", "AUCell")}


def score_all(adata, genes, scorers=RAW_SCORERS, timings=None):
    out = {}
    for name, fn in scorers.items():
        t0 = time.time()
        try:
            out[name] = np.asarray(fn(adata, genes), dtype=np.float64).ravel()
        except Exception as e:  # keep the benchmark running; record failure as NaN
            print(f"[score_all] {name} failed: {e}")
            out[name] = np.full(adata.n_obs, np.nan)
        if timings is not None:
            timings[name] = time.time() - t0
    return out


# ---------------------------------------------------------------------------
# Classification metrics (one binary task = one region vs rest)
# ---------------------------------------------------------------------------

def _f1(y, pred):
    y = np.asarray(y).astype(bool)
    pred = np.asarray(pred).astype(bool)
    tp = (y & pred).sum()
    denom = pred.sum() + y.sum()
    return float(2 * tp / denom) if denom else 0.0


def f1_oracle(y, s):
    """F1 at the F1-maximising threshold chosen on the evaluated data (what Sup. Figs 1e/2e
    report as 'F1'); an optimistic upper bound."""
    p, r, t = precision_recall_curve(y, s)
    p, r = p[:-1], r[:-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.nan_to_num(2 * p * r / (p + r))
    i = int(np.argmax(f))
    return float(f[i]), float(t[i])


def f1_median(y, s):
    """F1 when the score is thresholded at its median (as stated in the manuscript Methods)."""
    return _f1(y, s > np.median(s))


def f1_topk(y, s):
    """Prevalence-matched threshold: call the k highest-scoring spots positive, k = #positives.
    Then precision = recall = F1 (R-precision); no oracle threshold search."""
    y = np.asarray(y).astype(bool)
    k = int(y.sum())
    order = np.argsort(-np.asarray(s), kind="stable")
    pred = np.zeros(len(y), bool)
    pred[order[:k]] = True
    return _f1(y, pred)


def otsu_threshold(s, nbins=256):
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
    """Label-free threshold (Otsu) -- what a user without annotations could do."""
    return _f1(y, np.asarray(s) > otsu_threshold(s))


def spatial_folds(coords, n_blocks=25, n_folds=5, seed=0):
    km = KMeans(n_clusters=n_blocks, n_init=10, random_state=seed).fit(np.asarray(coords))
    rng = np.random.default_rng(seed)
    block_fold = rng.permutation(np.arange(n_blocks) % n_folds)
    return block_fold[km.labels_]


def f1_spatial_cv(y, s, folds):
    """Threshold learnt (F1-optimal) on training spatial blocks, applied to held-out blocks;
    predictions pooled over folds. Avoids selecting the threshold on the evaluated spots."""
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=np.float64)
    pred = np.zeros(len(y), bool)
    for f in np.unique(folds):
        te_ = folds == f
        tr = ~te_
        if y[tr].sum() == 0:
            continue
        _, thr = f1_oracle(y[tr], s[tr])
        pred[te_] = s[te_] >= thr
    return _f1(y, pred)


def morans_i(values, W):
    """Moran's I with a binary spatial weight matrix (as _morans_i_continuous in 17_benchmarking)."""
    v = np.asarray(values, dtype=np.float64)
    z = v - v.mean()
    z2 = z @ z
    if z2 == 0:
        return 0.0
    return float((len(v) / W.sum()) * (z @ W.dot(z)) / z2)


def spatial_coherence(pred, W, gt_morans):
    """As defined in 17_benchmarking: Moran's I of the binary prediction / Moran's I of the truth."""
    mi = morans_i(np.asarray(pred, dtype=float), W)
    return float(np.clip(mi / gt_morans, 0, 1)) if gt_morans > 0 else max(mi, 0.0)


def classification_metrics(y, s, W=None, folds=None):
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
    """Max F1 achievable when exactly half the spots are called positive: all positives
    recovered (recall 1) and precision 2p -> F1 = 4p / (1 + 2p) (for p <= 0.5)."""
    p = np.asarray(prevalence, dtype=np.float64)
    return 4 * p / (1 + 2 * p)


# ---------------------------------------------------------------------------
# Spatial benchmarks used in the manuscript (re-implemented verbatim)
# ---------------------------------------------------------------------------

def neighbourhood_concordance(scores, labels, conn, iters=range(1, 11)):
    """Fig. 2c / Sup. Fig. 1a as implemented in 05_01_simulated_benchmarking.run_concordance:
    the *ground-truth labels* are smoothed k times with the row-normalised graph and
    Spearman-correlated with the (unchanged) method scores."""
    gt = np.asarray(labels, dtype=np.float64)
    out = {}
    cur = gt.copy()
    d = np.maximum(np.asarray(conn.sum(axis=1)).ravel(), 1e-10)
    for k in range(1, max(iters) + 1):
        cur = conn.dot(cur) / d
        if k in iters:
            rho = spearmanr(scores, cur).correlation
            out[k] = 0.0 if np.isnan(rho) else float(rho)
    return out


# ---------------------------------------------------------------------------
# Leakage
# ---------------------------------------------------------------------------

def hop_distance(conn, mask, max_hops=12):
    """Graph distance (in neighbour rings) of every spot to the nearest spot in `mask`."""
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
    """Standardised score by ring outside the region:
    L(d) = (mean score at ring d - mean far background) / (mean in region - mean far background),
    far background = spots >= `far_from` rings away. L(d)=0: no leakage; L(d)=1: as high as inside."""
    score = np.asarray(score, dtype=np.float64)
    inside = score[np.asarray(mask, bool)].mean()
    far = score[dist >= far_from].mean()
    span = inside - far
    prof = {}
    for d in range(1, max_ring + 1):
        sel = dist == d
        prof[d] = float((score[sel].mean() - far) / span) if sel.any() and span != 0 else np.nan
    return prof


def leakage_summary(score, mask, dist, halo=(1, 2, 3)):
    """Scalar leakage metrics:
    - halo index: mean L(d) over the halo rings (lower = less leakage);
    - halo AUROC: separation of in-region spots from halo spots (higher = sharper boundary);
    - halo false-positive share: among the top-k spots (k = region size), the share that lies in the halo."""
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


# ---------------------------------------------------------------------------
# Simulations
# ---------------------------------------------------------------------------

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
    """Copy of simulate_visium() from 05_01_simulated_benchmarking.py (Sup. Fig. 1), extended.

    With the default arguments the random stream and the output are identical to the
    original function. Extensions (all drawn from an independent RNG so they do not
    perturb the original stream):

    domain : 'voronoi' (original, n_clusters Voronoi cells, target = cluster_2),
             'stripe' (a thin sinusoidal band ~2*domain_width_spots spots wide, layer-like),
             'patches' (n_patches small discs of radius patch_radius_spots, niche-like).
    tech_field : None | 'gradient' | 'edge' | 'grf'  -- smooth multiplicative technical
             field F(x, y) acting on a random fraction `tech_gene_frac` of genes before
             count sampling (log F has range `tech_strength` in natural-log units).
             'gradient': linear trend in a random direction; 'edge': capture loss towards
             the capture-area border; 'grf': smooth Gaussian random field with length scale
             `tech_scale` x side.  tech_on_signature: None (signature genes affected at the
             same rate as other genes), True (all affected) or False (none affected).
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
    seed_coords = np.column_stack([
        rng.uniform(margin, capture_area_um - margin, n_clusters),
        rng.uniform(margin, capture_area_um - margin, n_clusters),
    ])
    dists = cdist(coords, seed_coords)
    dists += rng.normal(0, spot_spacing_um * 0.8, dists.shape)
    cluster_ids = np.argmin(dists, axis=1)
    cluster_labels = np.array([f"cluster_{i+1}" for i in cluster_ids])
    target_label = f"cluster_{target_cluster + 1}"

    if domain == "stripe":
        L = capture_area_um
        phase = rng_ext.uniform(0, 2 * np.pi)
        amp = rng_ext.uniform(0.08, 0.15) * L
        y0 = rng_ext.uniform(0.35, 0.65) * L
        curve = y0 + amp * np.sin(2 * np.pi * coords[:, 0] / (0.8 * L) + phase)
        in_target = np.abs(coords[:, 1] - curve) <= domain_width_spots * spot_spacing_um
        cluster_labels = np.where(in_target, target_label, np.where(cluster_labels == target_label, "cluster_bg", cluster_labels))
    elif domain == "patches":
        L = capture_area_um
        centres = rng_ext.uniform(0.08 * L, 0.92 * L, size=(n_patches, 2))
        dd = cdist(coords, centres).min(axis=1)
        in_target = dd <= patch_radius_spots * spot_spacing_um
        cluster_labels = np.where(in_target, target_label, np.where(cluster_labels == target_label, "cluster_bg", cluster_labels))
    elif domain != "voronoi":
        raise ValueError(domain)

    gene_names = [f"GENE_{i+1:05d}" for i in range(n_genes)]
    sig_indices = rng.choice(n_genes, n_signature_genes, replace=False)
    fold_changes = rng.uniform(fold_change_range[0], fold_change_range[1], n_signature_genes)
    order = np.argsort(-fold_changes)
    sig_indices = sig_indices[order]
    fold_changes = fold_changes[order]
    sig_genes = [gene_names[i] for i in sig_indices]

    gene_means = np.exp(rng.normal(1.5, 1.0, n_genes))
    gene_means = np.clip(gene_means, 0.1, 100)
    dispersions = rng.uniform(dispersion_range[0], dispersion_range[1], n_genes)
    lib_size_factors = np.exp(rng.normal(0, 0.3, n_spots))
    lib_size_factors *= mean_counts_per_spot / (gene_means.sum() * lib_size_factors.mean())

    # --- technical field (independent RNG) -------------------------------------
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
        return adata, sig_genes, fold_changes, target_label, {"log_field": log_field, "affected": affected}
    return adata, sig_genes, fold_changes, target_label


# ---------------------------------------------------------------------------
# Reference-based (control gene set) spatial correction -- a proposed alternative
# ---------------------------------------------------------------------------

def matched_control_sets(adata, genes, n_sets=20, n_bins=25, seed=0):
    """Expression-matched random control gene sets (same size as `genes`, drawn from the
    same mean-expression bins, signature genes excluded), as in Scanpy/Seurat module scores."""
    rng = np.random.default_rng(seed)
    genes = [g for g in dict.fromkeys(genes) if g in adata.var_names]
    means = pd.Series(np.asarray(dense(adata.X).mean(axis=0)).ravel(), index=adata.var_names)
    bins = pd.qcut(means.rank(method="first"), n_bins, labels=False)
    sig = set(genes)
    pool = {b: [g for g in bins.index[bins == b] if g not in sig] for b in range(n_bins)}
    need = pd.Series([bins[g] for g in genes]).value_counts()
    sets = []
    for _ in range(n_sets):
        cs = []
        for b, k in need.items():
            cs += list(rng.choice(pool[b], size=min(k, len(pool[b])), replace=False))
        sets.append(cs)
    return sets


def reference_corrected(adata, genes, conn, gam, n_sets=20, seed=0, smooth_first=True, scorer=raw_enrichmap):
    """Remove only the spatial trend that is *shared with expression-matched control gene sets*:
    corrected = S r_sig - H (mean_c S r_c). Biology specific to the signature is not absorbed by the
    GAM because the GAM is fitted to the control score, not to the signature score."""
    r = scorer(adata, genes)
    ctrl = np.mean([scorer(adata, cs) for cs in matched_control_sets(adata, genes, n_sets=n_sets, seed=seed)], axis=0)
    if smooth_first:
        r, ctrl = smooth(r, conn), smooth(ctrl, conn)
    return r - gam.fitted(ctrl)


# ---------------------------------------------------------------------------
# GSDensity bridge (R package; run in the isolated conda env `gsdensity_env`)
# ---------------------------------------------------------------------------

GSD_RSCRIPT = REV + "scripts/run_gsdensity.R"


def export_gsd_task(adata, genesets, task_dir, params=None):
    """Write the inputs run_gsdensity.R expects (log-normalised adata.X, genes x spots)."""
    import json
    from scipy.io import mmwrite

    os.makedirs(task_dir, exist_ok=True)
    X = adata.X.T
    X = csr_matrix(X) if not issparse(X) else X.tocsr()
    mmwrite(os.path.join(task_dir, "data.mtx"), X.astype(np.float64))
    with open(os.path.join(task_dir, "genes.txt"), "w") as f:
        f.write("\n".join(map(str, adata.var_names)) + "\n")
    with open(os.path.join(task_dir, "spots.txt"), "w") as f:
        f.write("\n".join(map(str, adata.obs_names)) + "\n")
    xy = np.asarray(adata.obsm["spatial"])[:, :2]
    pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1]}, index=adata.obs_names).to_csv(os.path.join(task_dir, "coords.csv"))
    with open(os.path.join(task_dir, "genesets.json"), "w") as f:
        json.dump({k: list(map(str, v)) for k, v in genesets.items()}, f)
    if params:
        with open(os.path.join(task_dir, "params.json"), "w") as f:
            json.dump(params, f)
    for fn in ("gsd_scores.csv", "gsd_labels.csv", "gsd_spatial_kld.csv", "gsd_timing.json", "gsd_error.txt"):
        if os.path.exists(os.path.join(task_dir, fn)):
            os.remove(os.path.join(task_dir, fn))
    return task_dir


def run_gsd(task_dirs, n_jobs=1, log=None):
    """Run GSDensity on task dirs in `n_jobs` parallel R sessions (activating gsdensity_env)."""
    import subprocess

    task_dirs = list(task_dirs)
    chunks = [task_dirs[i::n_jobs] for i in range(n_jobs) if task_dirs[i::n_jobs]]
    procs = []
    for i, ch in enumerate(chunks):
        cmd = (os.environ.get("GSDENSITY_ACTIVATE", 'eval "$(conda shell.bash hook)" && conda activate gsdensity_env && ') +
               f"Rscript {GSD_RSCRIPT} " + " ".join(ch))
        out = open(log.replace(".log", f"_{i}.log"), "w") if log else subprocess.DEVNULL
        procs.append(subprocess.Popen(["bash", "-c", cmd], stdout=out, stderr=subprocess.STDOUT))
    for p in procs:
        p.wait()
    return [p.returncode for p in procs]


def read_gsd(task_dir, spots=None):
    import json

    err = os.path.join(task_dir, "gsd_error.txt")
    if os.path.exists(err):
        raise RuntimeError(open(err).read())
    sc_ = pd.read_csv(os.path.join(task_dir, "gsd_scores.csv"), index_col=0)
    lab = pd.read_csv(os.path.join(task_dir, "gsd_labels.csv"), index_col=0)
    skp = os.path.join(task_dir, "gsd_spatial_kld.csv")
    sk = pd.read_csv(skp, index_col=0) if os.path.exists(skp) else None
    tim = json.load(open(os.path.join(task_dir, "gsd_timing.json")))
    if spots is not None:
        sc_, lab = sc_.reindex(spots), lab.reindex(spots)
    return {"scores": sc_, "labels": lab, "spatial_kld": sk, "timing": tim}


# ---------------------------------------------------------------------------
# Plot style
# ---------------------------------------------------------------------------

def set_style():
    import matplotlib as mpl

    mpl.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 7, "axes.titlesize": 8,
        "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
        "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150, "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })


CONFIG_COLOURS = {"raw": "#bdbdbd", "smooth": "#3182bd", "GAM": "#fdae6b", "smooth+GAM": "#e6550d"}
