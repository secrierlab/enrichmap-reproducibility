#!/usr/bin/env python
"""
F2b -- factorial ablation on the ORIGINAL Sup. Fig. 1 simulations (Reviewer 4, points F1/F2).

For every simulated slide (ru.simulate_visium(seed), default arguments = original generator) and every scorer
{EnrichMap, z-score, scanpy, AUCell, GSVA, ssGSEA, GSDensity}, one raw score vector r is computed and the four
pipeline configurations are obtained with the same linear operators used inside em.tl.score:

    raw          r
    smooth       S r                     (6-NN neighbourhood mean, Eq. 3)
    GAM          (I - H) r               (pyGAM te(0,1) residual, Eq. 4-5, default lam / splines)
    smooth+GAM   (I - H) S r             (= EnrichMap default when r is the CV-weighted mean)

Benchmarks (as in 05_01_simulated_benchmarking.py, extended):
  * classification vs ground truth (AUROC, AUPRC, F1 oracle/median/top-k/Otsu/spatial-CV, Moran's I, boundary F1, KS)
  * neighbourhood concordance (labels smoothed k times, as implemented; plus the Methods-text variant)
  * contamination robustness (fractions 0, .5, 1, 2, 4 x signature size; 3 repeats) and noise robustness
    (levels 0, .5, 1, 2, 3 x gene SD on signature genes; 3 repeats), using the ORIGINAL random stream
    (np.random.default_rng(seed + 1000), contamination draws first, then noise draws). For each perturbation we
    record self-consistency (Spearman with the clean score of the same configuration; the original metric) and
    accuracy vs truth (AUROC / AUPRC of the perturbed score).
  * Sup. Fig. 1d gene-weight analysis (true CV weights vs the original mean|contribution| metric).

Parts (run as separate processes so slow methods never block fast ones):
  fast    EnrichMap, z-score, scanpy, AUCell (+ oracle 'Truth' and alternative EnrichMap weightings), all seeds
  ssgsea  ssGSEA (clean + contamination batched in one call, all seeds; noise for seed < N_SSGSEA_NOISE)
  gsva    GSVA   (clean + contamination batched in one call, all seeds)
  gsva_noise  GSVA noise, one repeat per level, seed < N_GSVA_NOISE (needs part gsva for the seed)
  gsd     GSDensity (R; clean + spatial KLD all seeds; contamination seed < N_GSD_CONTAM; noise seed < N_GSD_NOISE;
          parameter-sensitivity arm seed < N_GSD_PARAM)

Usage: python F2b_compute.py --part fast --seeds 0-99 --workers 2
Outputs: data/F2b/<part>/seed<NNN>.pkl (metrics, long format) and .npz (raw score vectors). Resumable.
"""
import os
import sys

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_v] = "2"
os.environ["PYTHONWARNINGS"] = "ignore"

import argparse
import json
import pickle
import shutil
import time
import traceback
import warnings

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import anndata as ad
from scipy.io import mmwrite
from scipy.sparse import csr_matrix, issparse
from scipy.stats import spearmanr, ks_2samp
from sklearn.metrics import roc_auc_score, average_precision_score, f1_score

import revision_utils as ru

OUT = ru.DATA + "F2b/"
LOGS = ru.REV + "logs/"
FRACS = [0.0, 0.5, 1.0, 2.0, 4.0]
LEVELS = [0.0, 0.5, 1.0, 2.0, 3.0]
REPS = 3
CONC_K = [1, 3, 5, 7, 10]
CONFIGS = ["raw", "smooth", "GAM", "smooth+GAM"]

# seed limits for the expensive arms (set after timing seeds 0-1 under the current server load; see notebook Section 0)
# GSVA/ssGSEA: contamination sets are scored in the same call as the clean set (the kernel/rank step dominates the
# cost, so this is ~free) -> contamination on all seeds. Noise needs one call per perturbed matrix.
N_SSGSEA_CONTAM = 100
N_SSGSEA_NOISE = 20     # full noise grid (4 levels x 3 repeats)
N_GSVA_CONTAM = 100
N_GSVA_NOISE = 10       # GSVA noise arm reduced to ONE repeat per level (part gsva_noise), seeds 0-9
GSVA_NOISE_LEVELS = [0.5, 1.0, 2.0, 3.0]
N_GSD_CONTAM = 20
N_GSD_NOISE = 10
N_GSD_PARAM = 20
GSD_PARAMS = {"dims30": {"dims": 30}, "dims50": {"dims": 50}, "restart0.5": {"restart": 0.5}, "nn100": {"nn_use": 100}}


# ---------------------------------------------------------------------------
# per-seed context: simulation, operators, perturbations (original RNG stream)
# ---------------------------------------------------------------------------

def boundary_mask(conn, y):
    """Spots with >= 1 neighbour of the other label (as run_classification in the original script)."""
    A = csr_matrix(conn)
    pos_nb = np.asarray(A.dot(y.astype(float))).ravel()
    deg = np.asarray(A.sum(axis=1)).ravel()
    return ((y == 1) & (pos_nb < deg)) | ((y == 0) & (pos_nb > 0))


def make_context(seed):
    adata, sig, fc, tgt = ru.simulate_visium(seed=seed)
    y = (adata.obs["cluster"] == tgt).values.astype(int)
    coords = np.asarray(adata.obsm["spatial"])
    conn = ru.spatial_conn(coords)
    gam = ru.SpatialGAM(coords)
    folds = ru.spatial_folds(coords, seed=seed)

    # --- original perturbation stream (05_01_simulated_benchmarking.py) ---
    rng = np.random.default_rng(seed + 1000)
    var = list(adata.var_names)
    common = [g for g in sig if g in adata.var_names]
    cs = set(common)
    non_sig = [g for g in var if g not in cs]
    contam = {}
    for frac in FRACS:
        n_contam = int(round(len(common) * frac))
        for rep in range(REPS):
            noise_genes = list(rng.choice(non_sig, n_contam, replace=len(non_sig) < n_contam)) if n_contam > 0 else []
            contam[(frac, rep)] = common + noise_genes
    Xd = ru.dense(adata.X)
    gidx = [var.index(g) for g in common]
    sds = np.array([Xd[:, gi].std() for gi in gidx])
    noise = {}
    for level in LEVELS:
        for rep in range(REPS):
            noise[(level, rep)] = np.column_stack([rng.normal(0, sds[j] * level, Xd.shape[0]) for j in range(len(gidx))])
    rng_ctrl = np.random.default_rng(seed + 777_777)  # independent stream for a random control set (GSDensity KLD)
    rand_set = list(rng_ctrl.choice(non_sig, len(common), replace=False))
    return dict(seed=seed, adata=adata, sig=sig, common=common, fc=fc, tgt=tgt, y=y, coords=coords, conn=conn,
                gam=gam, folds=folds, bmask=boundary_mask(conn, y), contam=contam, noise=noise, Xd=Xd, gidx=gidx,
                rand_set=rand_set)


def noisy_adata(ctx, level, rep):
    """X_n = X + N(0, (sd_g * level)^2) on signature genes, applied column-wise exactly as the original."""
    Xn = ctx["Xd"].copy()
    E = ctx["noise"][(level, rep)]
    for j, gi in enumerate(ctx["gidx"]):
        Xn[:, gi] += E[:, j]
    a = ctx["adata"]
    return ad.AnnData(X=Xn, obs=pd.DataFrame(index=a.obs_names), var=pd.DataFrame(index=a.var_names))


# ---------------------------------------------------------------------------
# batched scorers (one call for several gene sets; identical per-set results, verified in the notebook)
# ---------------------------------------------------------------------------

def aucell_multi(adata, sets):
    import decoupler as dc
    a = ad.AnnData(X=adata.X, obs=pd.DataFrame(index=adata.obs_names), var=pd.DataFrame(index=adata.var_names))
    net = pd.concat([pd.DataFrame({"source": k, "target": list(v)}) for k, v in sets.items()], ignore_index=True)
    dc.mt.aucell(a, net=net)
    df = dc.pp.get_obsm(adata=a, key="score_aucell").obsm["score_aucell"]
    return {k: df[k].values.astype(np.float64) for k in sets}


def gsva_multi(adata, sets, threads=2):
    import gseapy as gp
    res = gp.gsva(data=ru._expr_df(adata), gene_sets={k: list(v) for k, v in sets.items()}, outdir=None,
                  min_size=3, threads=threads)
    es = res.res2d.pivot(index="Term", columns="Name", values="ES")
    return {k: es.loc[k].T.astype(float).reindex(adata.obs_names).values for k in sets}


def ssgsea_multi(adata, sets, threads=2):
    """NOTE: gseapy's NES may be rescaled jointly over all sets of one call; every metric used here is invariant
    to a positive rescaling of the score (ranks, AUROC, Spearman, Moran's I; S and H are linear)."""
    import gseapy as gp
    res = gp.ssgsea(data=ru._expr_df(adata, upper=True), gene_sets={k: [g.upper() for g in v] for k, v in sets.items()},
                    outdir=None, sample_norm_method="rank", no_plot=True, min_size=3, threads=threads)
    nes = res.res2d.pivot(index="Term", columns="Name", values="NES")
    return {k: nes.loc[k].T.astype(float).reindex(adata.obs_names).values for k in sets}


def contam_sets(ctx, fracs=None):
    return {f"c{frac}_{rep}": g for (frac, rep), g in ctx["contam"].items() if frac > 0 and (fracs is None or frac in fracs)}


# ---------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------

def clf_metrics(y, s, ctx):
    s = np.asarray(s, dtype=np.float64)
    m = ru.classification_metrics(y, s, W=ctx["conn"], folds=ctx["folds"])
    if not m:
        return m
    # original definitions (05_01 run_classification) for the reproduction check
    _, thr = ru.f1_oracle(y, s)
    q = float(np.clip(np.mean(s <= thr), 0.01, 0.99))
    preds = (s > np.quantile(s, q)).astype(int)
    bm = ctx["bmask"]
    m["Boundary F1 (oracle thr; original)"] = float(f1_score(y[bm], preds[bm], zero_division=0))
    gt = ru.morans_i(y.astype(float), ctx["conn"])
    m["Spatial coherence (original)"] = ru.spatial_coherence(preds, ctx["conn"], gt)
    # prevalence-matched threshold (k = #positives): no threshold is tuned on the labels
    k = int(y.sum())
    top = np.zeros(len(y), int)
    top[np.argsort(-s, kind="stable")[:k]] = 1
    m["Boundary F1 (top-k)"] = float(f1_score(y[bm], top[bm], zero_division=0))
    m["Spatial coherence (top-k)"] = ru.spatial_coherence(top, ctx["conn"], gt)
    m["KS"] = float(ks_2samp(s[y == 1], s[y == 0]).statistic)
    return m


def concordance_scores_smoothed(scores, labels, conn, iters):
    """Methods-text variant: the SCORE is smoothed k times and correlated with the binary truth."""
    out, cur = {}, np.asarray(scores, dtype=np.float64)
    d = np.maximum(np.asarray(conn.sum(axis=1)).ravel(), 1e-10)
    for k in range(1, max(iters) + 1):
        cur = conn.dot(cur) / d
        if k in iters:
            rho = spearmanr(cur, labels).correlation
            out[k] = 0.0 if np.isnan(rho) else float(rho)
    return out


def configs(raw, ctx):
    return ru.configurations(raw, ctx["conn"], ctx["gam"])


def eval_clean(method, raw, ctx, rec, extra_cfgs=None):
    """Classification + concordance rows for every configuration of one raw score."""
    y, seed = ctx["y"], ctx["seed"]
    cf = configs(raw, ctx)
    if extra_cfgs:
        cf.update(extra_cfgs)
    for cfg, s in cf.items():
        if np.isnan(s).any():
            continue
        for metric, v in clf_metrics(y, s, ctx).items():
            rec["clf"].append(dict(seed=seed, method=method, config=cfg, metric=metric, value=float(v)))
        for k, rho in ru.neighbourhood_concordance(s, y, ctx["conn"], iters=CONC_K).items():
            rec["conc"].append(dict(seed=seed, method=method, config=cfg, k=k, variant="labels smoothed (as implemented)", rho=rho))
        for k, rho in concordance_scores_smoothed(s, y, ctx["conn"], CONC_K).items():
            rec["conc"].append(dict(seed=seed, method=method, config=cfg, k=k, variant="scores smoothed (Methods text)", rho=rho))
    return cf


def eval_pert(kind, method, level, rep, raw_pert, clean_cfgs, ctx, rec):
    y, seed = ctx["y"], ctx["seed"]
    raw_pert = np.asarray(raw_pert, dtype=np.float64)
    pc = configs(raw_pert, ctx) if not np.isnan(raw_pert).any() else {c: None for c in CONFIGS}
    for cfg in CONFIGS:
        s = pc[cfg]
        if s is None or np.nanstd(s) == 0:
            rho, auc, ap = np.nan, np.nan, np.nan
        else:
            rho = spearmanr(clean_cfgs[cfg], s).correlation
            auc, ap = roc_auc_score(y, s), average_precision_score(y, s)
        rec[kind].append(dict(seed=seed, method=method, config=cfg, level=level, rep=rep,
                              rho_self=float(rho), auroc=float(auc), auprc=float(ap)))


def new_rec():
    return {"clf": [], "conc": [], "contam": [], "noise": [], "gw": [], "diag": [], "gsd": [], "timing": {}}


def save(part, seed, rec, scores):
    os.makedirs(OUT + part, exist_ok=True)
    np.savez_compressed(OUT + f"{part}/seed{seed:03d}_scores.npz", **{k: np.asarray(v, np.float64) for k, v in scores.items()})
    with open(OUT + f"{part}/seed{seed:03d}.pkl.tmp", "wb") as f:
        pickle.dump(rec, f)
    os.replace(OUT + f"{part}/seed{seed:03d}.pkl.tmp", OUT + f"{part}/seed{seed:03d}.pkl")


# ---------------------------------------------------------------------------
# part: fast scorers
# ---------------------------------------------------------------------------

def run_fast(seed):
    import enrichmap as em
    t0 = time.time()
    ctx = make_context(seed)
    a, y, sig, conn, gam = ctx["adata"], ctx["y"], ctx["common"], ctx["conn"], ctx["gam"]
    rec, scores = new_rec(), {}
    rec["timing"]["setup"] = time.time() - t0

    # --- diagnostics of the GAM operator on this slide ---
    Sy = ru.smooth(y.astype(float), conn)
    rec["diag"].append(dict(seed=seed, n_spots=a.n_obs, prevalence=float(y.mean()), gam_edof=gam.edof,
                            gam_r2_labels=gam.r2(y.astype(float)), gam_r2_smoothed_labels=gam.r2(Sy),
                            morans_truth=ru.morans_i(y.astype(float), conn), n_boundary=int(ctx["bmask"].sum())))

    # --- clean scores ---
    fast = {"EnrichMap": ru.raw_enrichmap, "z-score": ru.raw_zscore, "scanpy": ru.raw_scanpy}
    clean = {}
    for name, fn in fast.items():
        t1 = time.time()
        clean[name] = np.asarray(fn(a, sig), dtype=np.float64)
        rec["timing"][name] = time.time() - t1
    t1 = time.time()
    clean["AUCell"] = aucell_multi(a, {"target": sig})["target"]
    rec["timing"]["AUCell"] = time.time() - t1
    clean["Truth (oracle)"] = y.astype(np.float64)

    # alternative EnrichMap weightings (Sup. Fig. 1d follow-up): SD weights and oracle log-FC weights
    w_cv = em.tl.infer_gene_weights(a, sig)
    Xs = ru.dense(a[:, sig].X)
    sd, mu = Xs.std(axis=0), Xs.mean(axis=0)
    clean["EnrichMap (SD weights)"] = ru.raw_enrichmap(a, sig, weights=dict(zip(sig, sd)))
    clean["EnrichMap (oracle log-FC weights)"] = ru.raw_enrichmap(a, sig, weights=dict(zip(ctx["sig"], np.log(ctx["fc"]))))

    clean_cfgs = {}
    for name, r in clean.items():
        extra = None
        if name in ("EnrichMap", "z-score"):
            scorer = ru.raw_enrichmap if name == "EnrichMap" else ru.raw_zscore
            t1 = time.time()
            extra = {"smooth+refGAM": ru.reference_corrected(a, sig, conn, gam, n_sets=20, seed=seed, smooth_first=True, scorer=scorer)}
            rec["timing"][name + " refGAM"] = time.time() - t1
        clean_cfgs[name] = eval_clean(name, r, ctx, rec, extra_cfgs=extra)
        scores[f"clean|{name}"] = r

    # --- reproduction check against em.tl.score (smoothing on, GAM off / on) ---
    a2 = a.copy()
    em.tl.score(a2, gene_set=list(ctx["sig"]), score_key="_em", smoothing=True, correct_spatial_covariates=False)
    em_nogam = a2.obs["_em_score"].values.copy()
    contributions = a2.uns["gene_contributions"]["_em"]
    em.tl.score(a2, gene_set=list(ctx["sig"]), score_key="_emg", smoothing=True, correct_spatial_covariates=True)
    em_gam = a2.obs["_emg_score"].values.copy()
    rec["diag"][0]["max_abs_diff_smooth_vs_em"] = float(np.max(np.abs(em_nogam - clean_cfgs["EnrichMap"]["smooth"])))
    rec["diag"][0]["max_abs_diff_smoothGAM_vs_em"] = float(np.max(np.abs(em_gam - clean_cfgs["EnrichMap"]["smooth+GAM"])))
    rec["diag"][0]["sd_em_gam"] = float(np.std(em_gam))

    # --- gene weights (Sup. Fig. 1d) ---
    common_c = [g for g in ctx["sig"] if g in contributions]
    mean_abs_contrib = np.array([np.mean(np.abs(contributions[g])) for g in common_c])
    fc = np.array([ctx["fc"][list(ctx["sig"]).index(g)] for g in common_c])
    cv = np.array([w_cv[g] for g in common_c])
    ii = [sig.index(g) for g in common_c]
    in_t = y == 1
    emp_lfc = Xs[in_t][:, ii].mean(axis=0) - Xs[~in_t][:, ii].mean(axis=0)
    base_mean = Xs[~in_t][:, ii].mean(axis=0)
    gene_auc = np.array([roc_auc_score(y, Xs[:, i]) for i in ii])
    for g, c, s_, m_, f_, e_, b_, au in zip(common_c, cv, sd[ii], mean_abs_contrib, fc, emp_lfc, base_mean, gene_auc):
        rec["gw"].append(dict(seed=seed, gene=g, cv_weight=float(c), sd=float(s_), mean_abs_contribution=float(m_),
                              true_fc=float(f_), empirical_dlog=float(e_), background_mean=float(b_), gene_auroc=float(au)))

    # --- contamination (original stream) ---
    t1 = time.time()
    cs = contam_sets(ctx)
    au = aucell_multi(a, cs)
    for (frac, rep), genes in ctx["contam"].items():
        key = f"c{frac}_{rep}"
        for name in ("EnrichMap", "z-score", "scanpy", "AUCell"):
            if frac == 0:
                rp = clean[name]  # identical gene set -> identical score
            elif name == "AUCell":
                rp = au[key]
            else:
                rp = fast[name](a, genes)
            eval_pert("contam", name, frac, rep, rp, clean_cfgs[name], ctx, rec)
            if frac > 0 and rep == 0:
                scores[f"contam|{name}|{key}"] = rp
    rec["timing"]["contam"] = time.time() - t1

    # --- noise (original stream) ---
    t1 = time.time()
    for (level, rep) in ctx["noise"]:
        an = noisy_adata(ctx, level, rep)
        for name in ("EnrichMap", "z-score", "scanpy", "AUCell"):
            if level == 0:
                rp = clean[name]
            elif name == "AUCell":
                rp = aucell_multi(an, {"target": sig})["target"]
            else:
                rp = fast[name](an, sig)
            eval_pert("noise", name, level, rep, rp, clean_cfgs[name], ctx, rec)
    rec["timing"]["noise"] = time.time() - t1
    rec["timing"]["total"] = time.time() - t0
    save("fast", seed, rec, scores)
    return seed, rec["timing"]["total"]


# ---------------------------------------------------------------------------
# parts: ssGSEA, GSVA
# ---------------------------------------------------------------------------

def run_rank_method(seed, method):
    t0 = time.time()
    ctx = make_context(seed)
    a, sig = ctx["adata"], ctx["common"]
    rec, scores = new_rec(), {}
    multi = gsva_multi if method == "GSVA" else ssgsea_multi
    do_contam = seed < (N_GSVA_CONTAM if method == "GSVA" else N_SSGSEA_CONTAM)
    if method == "GSVA":
        noise_list = []  # GSVA noise is computed by part gsva_noise (one repeat per level)
    else:
        noise_list = [(lv, rp) for (lv, rp) in ctx["noise"] if lv > 0] if seed < N_SSGSEA_NOISE else []

    sets = {"target": sig}
    if do_contam:
        sets.update(contam_sets(ctx))
    t1 = time.time()
    res = multi(a, sets)
    rec["timing"]["clean+contam call"] = time.time() - t1
    rec["timing"]["n_sets"] = len(sets)
    clean = res["target"]
    clean_cfgs = eval_clean(method, clean, ctx, rec)
    scores[f"clean|{method}"] = clean
    if do_contam:
        for (frac, rep) in ctx["contam"]:
            rp = clean if frac == 0 else res[f"c{frac}_{rep}"]
            eval_pert("contam", method, frac, rep, rp, clean_cfgs, ctx, rec)
    if noise_list:
        for rep in range(REPS):  # level 0 == clean input
            eval_pert("noise", method, 0.0, rep, clean, clean_cfgs, ctx, rec)
        for (level, rep) in noise_list:
            t1 = time.time()
            rp = multi(noisy_adata(ctx, level, rep), {"target": sig})["target"]
            rec["timing"][f"noise {level} {rep}"] = time.time() - t1
            eval_pert("noise", method, level, rep, rp, clean_cfgs, ctx, rec)
    rec["timing"]["total"] = time.time() - t0
    save(method.lower(), seed, rec, scores)
    return seed, rec["timing"]["total"]


def run_gsva_noise(seed):
    """GSVA noise arm, one repeat (rep 0) per level in GSVA_NOISE_LEVELS. The clean GSVA score is taken from part
    'gsva' of the same seed (computed there in the batched clean+contamination call). Levels already present in the
    'gsva' pickle (seeds 0-1 were timed with an earlier setting that included 1x and 3x) are not recomputed."""
    t0 = time.time()
    gp_ = OUT + f"gsva/seed{seed:03d}.pkl"
    if not os.path.exists(gp_):
        raise RuntimeError("run part gsva for this seed first")
    have = {(r["level"], r["rep"]) for r in pickle.load(open(gp_, "rb"))["noise"] if r["method"] == "GSVA"}
    ctx = make_context(seed)
    sig = ctx["common"]
    clean = np.load(OUT + f"gsva/seed{seed:03d}_scores.npz")["clean|GSVA"].astype(np.float64)
    clean_cfgs = configs(clean, ctx)
    rec, scores = new_rec(), {}
    if not any(lv == 0 for lv, _ in have):
        for rep in range(REPS):
            eval_pert("noise", "GSVA", 0.0, rep, clean, clean_cfgs, ctx, rec)
    for level in GSVA_NOISE_LEVELS:
        if (level, 0) in have:
            continue
        t1 = time.time()
        rp = gsva_multi(noisy_adata(ctx, level, 0), {"target": sig})["target"]
        rec["timing"][f"noise {level} 0"] = time.time() - t1
        eval_pert("noise", "GSVA", level, 0, rp, clean_cfgs, ctx, rec)
        scores[f"noise|GSVA|{level}_0"] = rp
    rec["timing"]["total"] = time.time() - t0
    save("gsva_noise", seed, rec, scores)
    return seed, rec["timing"]["total"]


# ---------------------------------------------------------------------------
# part: GSDensity (R, isolated env) via run_gsdensity.R
# ---------------------------------------------------------------------------

def _write_task(td, X_gxs, adata, genesets, params, data_from=None):
    """Same file layout as ru.export_gsd_task; data.mtx can be hard-linked from another task (parameter arms)."""
    os.makedirs(td, exist_ok=True)
    for fn in ("gsd_scores.csv", "gsd_labels.csv", "gsd_spatial_kld.csv", "gsd_timing.json", "gsd_error.txt", "data.mtx"):
        if os.path.exists(os.path.join(td, fn)):
            os.remove(os.path.join(td, fn))
    if data_from is not None:
        os.link(os.path.join(data_from, "data.mtx"), os.path.join(td, "data.mtx"))
    else:
        mmwrite(os.path.join(td, "data.mtx"), csr_matrix(X_gxs).astype(np.float64))
    with open(os.path.join(td, "genes.txt"), "w") as f:
        f.write("\n".join(map(str, adata.var_names)) + "\n")
    with open(os.path.join(td, "spots.txt"), "w") as f:
        f.write("\n".join(map(str, adata.obs_names)) + "\n")
    xy = np.asarray(adata.obsm["spatial"])[:, :2] if "spatial" in adata.obsm else None
    pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1]}, index=adata.obs_names).to_csv(os.path.join(td, "coords.csv"))
    with open(os.path.join(td, "genesets.json"), "w") as f:
        json.dump({k: list(map(str, v)) for k, v in genesets.items()}, f)
    with open(os.path.join(td, "params.json"), "w") as f:
        json.dump(params, f)
    return td


def gsd_transforms(p, n):
    """GSDensity's per-spot output is an RWR probability vector (>= 0, sums to 1). Primary: n * p (density relative
    to uniform; a positive rescaling, so every metric is unchanged and S, H act on the native quantity, as GSDensity's
    own weighted-KDE spatial step does). Sensitivity: log1p(n * p) (variance-stabilised)."""
    p = np.asarray(p, dtype=np.float64)
    return {"GSDensity": n * p, "GSDensity (log1p)": np.log1p(n * p)}


def run_gsd(seed):
    t0 = time.time()
    ctx = make_context(seed)
    a, sig, y = ctx["adata"], ctx["common"], ctx["y"]
    n = a.n_obs
    rec, scores = new_rec(), {}
    root = OUT + f"gsd_tasks/seed{seed:03d}/"
    tasks = {}
    main_sets = {"target": sig, "random": ctx["rand_set"]}
    if seed < N_GSD_CONTAM:
        main_sets.update(contam_sets(ctx))
    XT = a.X.T
    tasks["main"] = _write_task(root + "main", XT, a, main_sets, {"spatial": True})
    if seed < N_GSD_PARAM:
        for nm, pr in GSD_PARAMS.items():
            tasks["param_" + nm] = _write_task(root + "param_" + nm, None, a, {"target": sig}, dict(pr, spatial=False),
                                               data_from=root + "main")
    if seed < N_GSD_NOISE:
        for (level, rep) in ctx["noise"]:
            if level == 0:
                continue
            an = noisy_adata(ctx, level, rep)
            an.obsm["spatial"] = ctx["coords"]
            tasks[f"noise_{level}_{rep}"] = _write_task(root + f"noise_{level}_{rep}", an.X.T, an, {"target": sig}, {"spatial": False})
    rec["timing"]["export"] = time.time() - t0
    t1 = time.time()
    rc = ru.run_gsd(list(tasks.values()), n_jobs=1, log=LOGS + f"F2b_gsd_seed{seed:03d}.log")
    rec["timing"]["R"] = time.time() - t1
    rec["timing"]["n_tasks"] = len(tasks)

    out = {k: ru.read_gsd(td, spots=a.obs_names) for k, td in tasks.items()}
    for k, td in tasks.items():
        rec["timing"]["R " + k] = out[k]["timing"]
        os.remove(os.path.join(td, "data.mtx"))  # large; outputs (csv) are kept

    m = out["main"]
    clean_cfgs = {}
    for nm, sv in gsd_transforms(m["scores"]["target"].values, n).items():
        clean_cfgs[nm] = eval_clean(nm, sv, ctx, rec)
        scores[f"clean|{nm}"] = sv
    lab = (m["labels"]["target"] == "positive").values
    p_t = m["scores"]["target"].values
    rec["gsd"].append(dict(seed=seed, setting="default", native_label_F1=ru._f1(y, lab), native_frac_positive=float(lab.mean()),
                           frac_zero=float((p_t == 0).mean()), **{f"kld_{c}": float(m["spatial_kld"].loc["target", c])
                                                                  for c in ("spatial.kld", "rspatial.kld", "p.value")}))
    sk = m["spatial_kld"].copy()
    sk["seed"] = seed
    rec["gsd_kld"] = sk.reset_index().rename(columns={"index": "set"}).to_dict("records")
    # parameter-sensitivity arm
    for k in tasks:
        if k.startswith("param_"):
            nm = k[len("param_"):]
            pv = out[k]["scores"]["target"].values
            for tn, sv in gsd_transforms(pv, n).items():
                eval_clean(f"{tn} [{nm}]", sv, ctx, rec)
            labk = (out[k]["labels"]["target"] == "positive").values
            rec["gsd"].append(dict(seed=seed, setting=nm, native_label_F1=ru._f1(y, labk), native_frac_positive=float(labk.mean()),
                                   frac_zero=float((pv == 0).mean())))
            scores[f"param|{nm}"] = pv * n
    # contamination
    if seed < N_GSD_CONTAM:
        for (frac, rep) in ctx["contam"]:
            pv = p_t if frac == 0 else m["scores"][f"c{frac}_{rep}"].values
            for tn, sv in gsd_transforms(pv, n).items():
                eval_pert("contam", tn, frac, rep, sv, clean_cfgs[tn], ctx, rec)
    if seed < N_GSD_NOISE:
        for rep in range(REPS):
            for tn, sv in gsd_transforms(p_t, n).items():
                eval_pert("noise", tn, 0.0, rep, sv, clean_cfgs[tn], ctx, rec)
        for (level, rep) in ctx["noise"]:
            if level == 0:
                continue
            pv = out[f"noise_{level}_{rep}"]["scores"]["target"].values
            for tn, sv in gsd_transforms(pv, n).items():
                eval_pert("noise", tn, level, rep, sv, clean_cfgs[tn], ctx, rec)
    rec["timing"]["total"] = time.time() - t0
    save("gsd", seed, rec, scores)
    return seed, rec["timing"]["total"]


# ---------------------------------------------------------------------------
# part: noiseclip -- sensitivity of the noise arm to the (unphysical) negative values created by additive Gaussian
# noise on log-normalised data. Same perturbation as the original (rep 0), then clipped at 0, for the fast scorers and
# GSDensity (whose MCA fuzzy-codes each gene by min-max scaling). Seeds < N_CLIP; levels CLIP_LEVELS.
# ---------------------------------------------------------------------------
N_CLIP = 5
CLIP_LEVELS = [0.5, 1.0, 3.0]


def run_noiseclip(seed):
    t0 = time.time()
    ctx = make_context(seed)
    a, sig, n = ctx["adata"], ctx["common"], ctx["adata"].n_obs
    rec, scores = new_rec(), {}
    fast = {"EnrichMap": ru.raw_enrichmap, "z-score": ru.raw_zscore, "scanpy": ru.raw_scanpy}
    clean = {k: np.asarray(f(a, sig), dtype=np.float64) for k, f in fast.items()}
    clean["AUCell"] = aucell_multi(a, {"target": sig})["target"]
    clean["GSDensity"] = np.load(OUT + f"gsd/seed{seed:03d}_scores.npz")["clean|GSDensity"].astype(np.float64)
    clean_cfgs = {k: configs(v, ctx) for k, v in clean.items()}
    root = OUT + f"gsd_tasks/seed{seed:03d}/"
    tasks, noisy = {}, {}
    for level in CLIP_LEVELS:
        an = noisy_adata(ctx, level, 0)
        an.X = np.maximum(an.X, 0)
        an.obsm["spatial"] = ctx["coords"]
        noisy[level] = an
        tasks[level] = _write_task(root + f"noiseclip_{level}_0", an.X.T, an, {"target": sig}, {"spatial": False})
    ru.run_gsd(list(tasks.values()), n_jobs=1, log=LOGS + f"F2b_noiseclip_seed{seed:03d}.log")
    for level, an in noisy.items():
        for k, f in fast.items():
            eval_pert("noise", k, level, 0, f(an, sig), clean_cfgs[k], ctx, rec)
        eval_pert("noise", "AUCell", level, 0, aucell_multi(an, {"target": sig})["target"], clean_cfgs["AUCell"], ctx, rec)
        pv = ru.read_gsd(tasks[level], spots=a.obs_names)["scores"]["target"].values
        eval_pert("noise", "GSDensity", level, 0, n * pv, clean_cfgs["GSDensity"], ctx, rec)
        os.remove(os.path.join(tasks[level], "data.mtx"))
    rec["timing"]["total"] = time.time() - t0
    save("noiseclip", seed, rec, scores)
    return seed, rec["timing"]["total"]


# ---------------------------------------------------------------------------

PARTS = {"fast": run_fast, "ssgsea": lambda s: run_rank_method(s, "ssGSEA"), "gsva": lambda s: run_rank_method(s, "GSVA"),
         "gsva_noise": run_gsva_noise, "gsd": run_gsd, "noiseclip": run_noiseclip}


def _safe(part, seed):
    if os.path.exists(OUT + f"{part}/seed{seed:03d}.pkl"):
        return seed, 0.0
    try:
        return PARTS[part](seed)
    except Exception:
        msg = traceback.format_exc()
        os.makedirs(OUT + part, exist_ok=True)
        with open(OUT + f"{part}/seed{seed:03d}.error", "w") as f:
            f.write(msg)
        print(f"[{part}] seed {seed} FAILED\n{msg}", flush=True)
        return seed, np.nan


def parse_seeds(s):
    out = []
    for tok in s.split(","):
        if "-" in tok:
            lo, hi = tok.split("-")
            out += list(range(int(lo), int(hi) + 1))
        else:
            out.append(int(tok))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", required=True, choices=list(PARTS))
    ap.add_argument("--seeds", default="0-99")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--order", default="asc", choices=["asc", "desc"])
    args = ap.parse_args()
    seeds = parse_seeds(args.seeds)
    if args.order == "desc":
        seeds = seeds[::-1]
    t0 = time.time()
    print(f"[F2b {args.part}] {len(seeds)} seeds, {args.workers} workers", flush=True)
    if args.workers == 1:
        for sd in seeds:
            s_, tt = _safe(args.part, sd)
            print(f"[F2b {args.part}] seed {s_} {tt:.1f}s (elapsed {time.time() - t0:.0f}s)", flush=True)
    else:
        from joblib import Parallel, delayed
        for s_, tt in Parallel(n_jobs=args.workers, backend="loky", return_as="generator_unordered")(
                delayed(_safe)(args.part, sd) for sd in seeds):
            print(f"[F2b {args.part}] seed {s_} {tt:.1f}s (elapsed {time.time() - t0:.0f}s)", flush=True)
    print(f"[F2b {args.part}] finished in {time.time() - t0:.0f}s", flush=True)
