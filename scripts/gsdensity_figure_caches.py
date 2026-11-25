#!/usr/bin/env python
"""
GSDensity results that the figure notebooks read from `data/` but that the REV_GSD_base_* notebooks do not write.

Each part is the code that wrote the cache for the published figures, with paths changed to this repository:

  f1    data/F1/F1_genesets.json and data/F1/gsd_main/
        50 random 10-gene sets and one GSDensity run on the clean section, which holds them
        (read by ExtendedDataFig3_ExtendedDataFig4a-b; the clean matrix is also reused by part fig2).
  fig2  data/gsd_cache/brain/REV_Fig2/ and gsdensity_robustness.csv
        GSDensity on the Fig. 2d contaminated gene sets and the Fig. 2e noise matrices (read by Fig2).
  fig3  data/gsd_cache/brain/fig3/
        GSDensity on the noisy (10 signal + 10 noise genes) and signal-only gene sets of Fig. 3 (read by Fig3).
  ed5d  data/gsd_cache/brain/SupFig5_6_7/runtime_20261029/
        timed GSDensity runs, one per gene set size and repeat (read by ExtendedDataFig5d).

Run after 00_gene_set, 00_visium_brain_preprocessing and REV_GSD_base_brain, and run part f1 before part fig2.
Every part is resumable: GSDensity tasks that already finished are skipped.

Usage:  python scripts/gsdensity_figure_caches.py f1 fig2 fig3 ed5d
        python scripts/gsdensity_figure_caches.py all

GSDensity runs in the conda environment gsdensity_env (see README); set GSDENSITY_ACTIVATE to activate it another way.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Sequence

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import revision_utils as ru  # noqa: E402
from benchmark_functions import benchmark_contamination, benchmark_noise_robustness  # noqa: E402

DATA = ru.REV + "data/"
LOGS = ru.REV + "logs/"
os.makedirs(LOGS, exist_ok=True)
GSD_ACTIVATE = os.environ.get(
    "GSDENSITY_ACTIVATE", 'eval "$(conda shell.bash hook)" && conda activate gsdensity_env && '
)
REGIONS = ("Pyramidal_layer", "Striatum")


def gsd_task_done(task_dir: str) -> bool:
    """True when a GSDensity task directory holds scores and timing and no error."""
    return (
        os.path.exists(os.path.join(task_dir, "gsd_scores.csv"))
        and os.path.exists(os.path.join(task_dir, "gsd_timing.json"))
        and not os.path.exists(os.path.join(task_dir, "gsd_error.txt"))
    )


# ---------------------------------------------------------------------------
# f1: random gene sets and the GSDensity run on the clean section (scripts/F1_compute.py, F1_helpers.py)
# ---------------------------------------------------------------------------

F1 = DATA + "F1/"
PANGLAO_TYPES = ("Neurons", "Cardiomyocytes", "Embryonic stem cells", "Hepatocytes")
CONTAM_FRACS = (0.25, 0.5, 1.0, 2.0, 4.0)


def safe(name):
    return name.replace(" ", "_")


def panglao_sets(cache=F1 + "F1_panglao_sets.json"):
    """PanglaoDB control sets as in 14_unit_tests.ipynb (cells 10-11): mouse, canonical
    markers, mouse sensitivity > 0.5, duplicates removed, first 10 per cell type. Cached to JSON."""
    if os.path.exists(cache):
        return json.load(open(cache))
    import decoupler as dc

    m = dc.op.resource("PanglaoDB", organism="mouse")
    m = m[m["mouse"].astype(bool) & m["canonical_marker"].astype(bool) & (m["mouse_sensitivity"].astype(float) > 0.5)]
    m = m[~m.duplicated(["cell_type", "genesymbol"])]
    m = m.rename(columns={"cell_type": "source", "genesymbol": "target"})[["source", "target"]]
    m = m[m["source"].isin(list(PANGLAO_TYPES))]
    d = m.groupby("source")["target"].apply(list).to_dict()
    out = {k: v[:10] for k, v in d.items()}
    json.dump(out, open(cache, "w"), indent=1)
    return out


def contamination_set(sig, var_names, frac, rep, region_idx):
    """Signature + round(n*frac) random non-signature genes (independent seed per region, level and repeat)."""
    rng = np.random.default_rng(10_000 + 1000 * region_idx + int(round(100 * frac)) + 7 * rep)
    non_sig = [g for g in var_names if g not in set(sig)]
    n_contam = int(round(len(sig) * frac))
    return list(sig) + list(rng.choice(non_sig, n_contam, replace=False))


def build_gene_sets(adata, sigs, n_random=50, n_matched=50, n_contam_rep=5, seed=0):
    """All gene sets scored on the clean brain section, with a metadata table.
    category: region | panglao_positive | panglao_absent | random | matched | contaminated."""
    var = list(adata.var_names)
    vset = set(var)
    sets, meta = {}, []

    def add(name, genes, category, base=None, level=np.nan, rep=np.nan):
        genes = [g for g in dict.fromkeys(genes) if g in vset]
        sets[name] = genes
        meta.append(dict(set=name, category=category, base=base, level=level, rep=rep, n_genes=len(genes)))

    region_sigs = {r: list(sigs[r][f"{r}_top10"]) for r in REGIONS}
    for r in REGIONS:
        add(r, region_sigs[r], "region", base=r, level=0.0)
    pg = panglao_sets()
    for ct in PANGLAO_TYPES:
        add(safe(ct), pg[ct], "panglao_positive" if ct == "Neurons" else "panglao_absent", base=ct)
    # uniformly random 10-gene sets (region-signature genes excluded)
    rng = np.random.default_rng(seed)
    excl = set(region_sigs["Pyramidal_layer"]) | set(region_sigs["Striatum"])
    pool = [g for g in var if g not in excl]
    for i in range(n_random):
        add(f"random_{i:02d}", list(rng.choice(pool, 10, replace=False)), "random")
    # expression-matched control sets (same mean-expression bins as each region signature)
    for ri, r in enumerate(REGIONS):
        for i, cs in enumerate(ru.matched_control_sets(adata, region_sigs[r], n_sets=n_matched, seed=seed + 1 + ri)):
            add(f"matched_{r}_{i:02d}", cs, "matched", base=r)
    # contaminated signatures
    for ri, r in enumerate(REGIONS):
        for frac in CONTAM_FRACS:
            for rep in range(n_contam_rep):
                add(f"contam_{r}_{frac:g}_{rep}", contamination_set(region_sigs[r], var, frac, rep, ri),
                    "contaminated", base=r, level=frac, rep=rep)
    return sets, pd.DataFrame(meta).set_index("set")


def load_sets(adata, sigs):
    """Build (once) and cache all gene sets; later runs reload the cached JSON."""
    fj, fm = F1 + "F1_genesets.json", F1 + "F1_geneset_meta.csv"
    if os.path.exists(fj) and os.path.exists(fm):
        return json.load(open(fj)), pd.read_csv(fm, index_col=0)
    sets, meta = build_gene_sets(adata, sigs, n_random=50, n_matched=50, n_contam_rep=5, seed=0)
    # expression-matched nulls for the four PanglaoDB sets as well
    extra = []
    for i, ct in enumerate(PANGLAO_TYPES):
        base = sets[safe(ct)]
        for j, cs in enumerate(ru.matched_control_sets(adata, base, n_sets=50, seed=10 + i)):
            nm = f"matched_{safe(ct)}_{j:02d}"
            sets[nm] = [g for g in dict.fromkeys(cs)]
            extra.append(dict(set=nm, category="matched", base=safe(ct), level=np.nan, rep=np.nan, n_genes=len(sets[nm])))
    meta = pd.concat([meta, pd.DataFrame(extra).set_index("set")])
    json.dump(sets, open(fj, "w"))
    meta.to_csv(fm)
    return sets, meta


def part_f1():
    os.makedirs(F1, exist_ok=True)
    adata, sigs = ru.load_brain()
    sets, meta = load_sets(adata, sigs)
    cat = meta["category"]
    named = [k for k in meta.index if cat[k] in ("region", "panglao_positive", "panglao_absent")]
    main_sets = named + [k for k in meta.index if cat[k] == "random" or (cat[k] == "matched" and meta.base[k] in REGIONS)]
    main_dir = F1 + "gsd_main/"
    if gsd_task_done(main_dir):
        print(f"f1: GSDensity already done ({main_dir})")
        return
    ru.export_gsd_task(adata, {k: sets[k] for k in main_sets}, main_dir)
    rc = ru.run_gsd([main_dir], n_jobs=1, log=LOGS + "F1_gsd_main.log")
    assert rc == [0] and gsd_task_done(main_dir), f"GSDensity failed (see {LOGS}F1_gsd_main_0.log)"
    print(f"f1: {len(main_sets)} gene sets in {main_dir}")


# ---------------------------------------------------------------------------
# fig2: GSDensity on the Fig. 2d/2e perturbations (scripts/rev_fig2_helpers.py, scratch build_fig2.py)
# ---------------------------------------------------------------------------

CONTAMINATION_FRACTIONS = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0]  # 17_benchmarking cell 19
NOISE_LEVELS = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0]  # 17_benchmarking cell 16
GSD_BASE_DIR = DATA + "gsd_cache/brain/"  # base relevance per signature (REV_GSD_base_brain)
GSD_FIG2_DIR = GSD_BASE_DIR + "REV_Fig2/"
GSD_MAIN_DIR = F1 + "gsd_main/"  # clean matrix of the base run (part f1), read via symlinks
GSD_FIG2_PARAMS = {"spatial": False}  # relevance only; nn_use 300, restart 0.75, dims 10, seed 0


def record_contamination_sets(adata, sigs) -> dict:
    """
    Gene lists of the Fig. 2d contamination benchmark in call order (clean first): benchmark_contamination
    (seed 0, 10 repeats) replayed with a scorer that records each list instead of scoring it.
    """
    out = {}
    for ct in REGIONS:
        seen = []

        def record(a, gene_list, seen=seen):
            seen.append([str(g) for g in gene_list])
            return np.zeros(a.n_obs)

        benchmark_contamination(
            adata,
            gene_set=sigs[ct][f"{ct}_top10"],
            scoring_functions={"record": record},
            contamination_fractions=CONTAMINATION_FRACTIONS,
            n_repeats=10,
        )
        out[ct] = seen
    return out


def read_gsd_base(ct: str, spots) -> pd.DataFrame:
    """Base (unperturbed) GSDensity relevance and antimode label of one brain signature."""
    df = pd.read_csv(GSD_BASE_DIR + f"{ct}.csv", index_col=0)
    return df.reindex(spots)[["relevance", "label"]]


def _link_main_data(task_dir: str, genesets: dict) -> None:
    """Task on the clean matrix of the base run (symlinked inputs; only gene sets and parameters are written)."""
    os.makedirs(task_dir, exist_ok=True)
    for fn in ("data.mtx", "genes.txt", "spots.txt", "coords.csv"):
        os.symlink(GSD_MAIN_DIR + fn, os.path.join(task_dir, fn))
    with open(os.path.join(task_dir, "genesets.json"), "w") as f:
        json.dump({k: list(map(str, v)) for k, v in genesets.items()}, f)
    with open(os.path.join(task_dir, "params.json"), "w") as f:
        json.dump(GSD_FIG2_PARAMS, f)


def _run_gsd_task(name: str) -> str:
    """Run run_gsdensity.R on the prepared temporary task `_tmp_<name>` and rename it to `<name>` when complete."""
    final = GSD_FIG2_DIR + name + "/"
    tmp = GSD_FIG2_DIR + "_tmp_" + name + "/"
    os.makedirs(GSD_FIG2_DIR + "logs/", exist_ok=True)
    cmd = (
        GSD_ACTIVATE
        + "export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 && "
        + f"Rscript {ru.GSD_RSCRIPT} {tmp}"
    )
    with open(GSD_FIG2_DIR + f"logs/{name}.log", "w") as log:
        rc = subprocess.call(["bash", "-c", cmd], stdout=log, stderr=subprocess.STDOUT)
    if rc != 0 or not gsd_task_done(tmp):
        raise RuntimeError(f"GSDensity task {name} failed (see {GSD_FIG2_DIR}logs/{name}.log)")
    mtx = os.path.join(tmp, "data.mtx")
    if os.path.exists(mtx) and not os.path.islink(mtx):
        os.remove(mtx)  # perturbed matrix (~0.5 GB); reproducible from the seed
    os.rename(tmp, final)
    return final


def contamination_task_sets(contamination_sets: dict) -> tuple:
    """Unique gene lists of the contamination benchmark, named for one multi-set GSDensity task."""
    genesets, key = {}, {}
    for ct, sets in contamination_sets.items():
        for g in sets:
            k = (ct, tuple(g))
            if k not in key:
                key[k] = f"{ct}__{len([1 for c, _ in key if c == ct]):02d}"
                genesets[key[k]] = list(g)
    return genesets, key


def run_gsd_robustness(contamination_sets: dict, n_jobs: int = 2, regions: Sequence[str] = REGIONS) -> dict:
    """
    GSDensity relevance for the contamination gene lists and noise matrices (resumable).

    Contamination: one multi-set task on the clean matrix (the random walk of a set does not depend on the other
    sets). Noise: benchmark_noise_robustness (seed 0, 5 repeats) is replayed with a scorer that exports each noisy
    matrix, so GSDensity sees exactly the matrices the other methods saw; one task per level and repeat. Level 0
    leaves the matrix bit-identical to the clean one (asserted): repeat 0 is run on the clean matrix and stands for
    the 5 identical repeats.
    """
    os.makedirs(GSD_FIG2_DIR, exist_ok=True)
    pool = ThreadPoolExecutor(max_workers=n_jobs)
    futures = {}
    pending = threading.Semaphore(n_jobs + 1)  # at most one exported task waiting for an R session

    def _submit(name, build):
        if gsd_task_done(GSD_FIG2_DIR + name + "/"):
            return
        tmp = GSD_FIG2_DIR + "_tmp_" + name + "/"
        if os.path.exists(tmp):
            shutil.rmtree(tmp)  # leftover of an interrupted run
        pending.acquire()
        build(tmp)

        def _job():
            try:
                return _run_gsd_task(name)
            finally:
                pending.release()

        futures[name] = pool.submit(_job)

    # contamination: one task, every unique list of both regions
    genesets, _ = contamination_task_sets({ct: contamination_sets[ct] for ct in regions})
    _submit("contamination", lambda tmp: _link_main_data(tmp, genesets))

    # noise: replay the benchmark's random stream and export each noisy matrix
    ad, sigs = ru.load_brain()
    clean_X = ru.dense(ad.X)
    for ct in regions:
        genes = [g for g in sigs[ct][f"{ct}_top10"] if g in ad.var_names]
        calls = [0]

        def exporter(a, gene_list, ct=ct, genes=genes, calls=calls):
            k = calls[0]
            calls[0] += 1
            if k == 0:  # clean reference call
                return np.zeros(a.n_obs)
            level, rep = NOISE_LEVELS[(k - 1) // 5], (k - 1) % 5
            name = f"noise_{ct}_{level:g}_{rep}"
            if level == 0:
                assert np.array_equal(ru.dense(a.X), clean_X)
                if rep == 0:
                    _submit(name, lambda tmp: _link_main_data(tmp, {ct: genes}))
            else:
                _submit(name, lambda tmp: ru.export_gsd_task(a, {ct: genes}, tmp, params=GSD_FIG2_PARAMS))
            return np.zeros(a.n_obs)

        benchmark_noise_robustness(
            ad,
            gene_set=sigs[ct][f"{ct}_top10"],
            scoring_functions={"export": exporter},
            noise_levels=NOISE_LEVELS,
            n_repeats=5,
        )
    for name, fut in futures.items():
        fut.result()
    pool.shutdown()
    names = ["contamination"] + [
        f"noise_{ct}_{lvl:g}_{rep}" for ct in regions for lvl in NOISE_LEVELS for rep in range(5)
        if lvl > 0 or rep == 0
    ]
    return {n: GSD_FIG2_DIR + n + "/" for n in names}


def gsd_contamination_long(contamination_sets: dict, spots) -> tuple:
    """Spearman rho of GSDensity relevance between contaminated and clean gene lists (benchmark_contamination replayed)."""
    ad, sigs = ru.load_brain()
    genesets, key = contamination_task_sets(contamination_sets)
    scores = ru.read_gsd(GSD_FIG2_DIR + "contamination/", spots=spots)["scores"]
    out, check = [], {}
    for ct in contamination_sets:
        base = read_gsd_base(ct, spots)["relevance"].values.astype(np.float64)
        clean = tuple(contamination_sets[ct][0])
        check[ct] = float(np.abs(scores[key[(ct, clean)]].values - base).max())
        table = {g: scores[name].values.astype(np.float64) for (c, g), name in key.items() if c == ct}
        table[clean] = base  # reference = base relevance (data/gsd_cache/brain)

        def lookup(a, gene_list, table=table):
            return table[tuple(gene_list)]

        res = benchmark_contamination(
            ad,
            gene_set=sigs[ct][f"{ct}_top10"],
            scoring_functions={"GSDensity": lookup},
            contamination_fractions=CONTAMINATION_FRACTIONS,
            n_repeats=10,
        )
        res["region"] = ct
        out.append(res)
    return pd.concat(out, ignore_index=True), check


def gsd_noise_long(regions: Sequence[str], spots) -> tuple:
    """Spearman rho of GSDensity relevance between noisy and clean data (same noise matrices as the other methods)."""
    rows, check = [], {}
    for ct in regions:
        base = read_gsd_base(ct, spots)["relevance"].values.astype(np.float64)
        for lvl in NOISE_LEVELS:
            for rep in range(5):
                name = f"noise_{ct}_{lvl:g}_{0 if lvl == 0 else rep}"
                v = ru.read_gsd(GSD_FIG2_DIR + name + "/", spots=spots)["scores"][ct].values.astype(np.float64)
                if lvl == 0:
                    check[ct] = float(np.abs(v - base).max())
                rho, _ = spearmanr(base, v)
                rows.append({"method": "GSDensity", "noise_level": lvl, "repeat": rep, "correlation": rho,
                             "region": ct})
    return pd.DataFrame(rows), check


def part_fig2(n_jobs=2):
    assert gsd_task_done(GSD_MAIN_DIR), "run part f1 first (the contamination task uses its clean matrix)"
    adata, sigs = ru.load_brain()
    contamination_sets = record_contamination_sets(adata, sigs)
    tasks = run_gsd_robustness(contamination_sets, n_jobs=n_jobs)
    contamination, contamination_check = gsd_contamination_long(contamination_sets, adata.obs_names)
    noise, noise_check = gsd_noise_long(REGIONS, adata.obs_names)
    print(f"fig2: {len(tasks)} GSDensity tasks")
    print(f"max |clean relevance in the contamination task - base relevance|: {contamination_check}")
    print(f"max |noise level 0 relevance - base relevance|: {noise_check}")
    # the table Fig2 reads: experiment, region, x (contamination fraction or noise level), repeat, correlation
    rows = []
    for experiment, df, x in (("contamination", contamination, "contamination"), ("noise", noise, "noise_level")):
        d = df.rename(columns={x: "x"})[["region", "x", "repeat", "correlation"]].copy()
        d.insert(0, "experiment", experiment)
        rows.append(d)
    pd.concat(rows).to_csv(GSD_FIG2_DIR + "gsdensity_robustness.csv", index=False)
    print(f"fig2: wrote {GSD_FIG2_DIR}gsdensity_robustness.csv")


# ---------------------------------------------------------------------------
# fig3: GSDensity on the Fig. 3 gene sets (REV_Fig3_gene_weights, steps 2, 7 and 8, 28 Oct 2026)
# ---------------------------------------------------------------------------

def fig3_noise_genes(adata, signature_1, signature_2):
    """The 10 noise genes of Fig. 3 (18_visium_weights cell 7)."""
    from scipy.sparse import issparse

    if issparse(adata.X):
        gene_means_arr = np.asarray(adata.X.mean(axis=0)).flatten()
        gene_sq_means = np.asarray(adata.X.multiply(adata.X).mean(axis=0)).flatten()
        gene_vars_arr = gene_sq_means - gene_means_arr ** 2
    else:
        gene_means_arr = adata.X.mean(axis=0)
        gene_vars_arr = adata.X.var(axis=0)

    gene_mean_series = pd.Series(gene_means_arr, index=adata.var_names)
    gene_cv_series = pd.Series(
        np.where(gene_means_arr > 0, np.sqrt(np.maximum(gene_vars_arr, 0)) / gene_means_arr, 0.0),
        index=adata.var_names,
    )
    all_sig_genes = set(list(signature_1) + list(signature_2))
    sig_cvs = gene_cv_series[[g for g in all_sig_genes if g in adata.var_names]]
    min_sig_cv = sig_cvs.min()
    bg_expressed = [g for g in adata.var_names if g not in all_sig_genes and gene_mean_series[g] > 0]
    expr_floor = gene_mean_series[bg_expressed].quantile(0.25)
    noise_candidates = [
        g for g in bg_expressed if gene_cv_series[g] < min_sig_cv and gene_mean_series[g] >= expr_floor
    ]
    rng_noise = np.random.default_rng(seed=42)
    noise_genes = sorted(rng_noise.choice(noise_candidates, size=10, replace=False).tolist())
    assert noise_genes == [
        "Anp32e", "Apba2", "Clptm1l", "Kcnh3", "Mdk", "Prkra", "Rbm6", "Sh3gl3", "Ssb", "Trappc1",
    ], noise_genes
    return noise_genes


def part_fig3():
    gsd_base_dir = GSD_BASE_DIR
    gsd_dir = gsd_base_dir + "fig3/"
    os.makedirs(gsd_dir, exist_ok=True)
    adata, sigs = ru.load_brain()
    cell_type_1, cell_type_2 = REGIONS
    signature_1 = sigs[cell_type_1][f"{cell_type_1}_top10"]
    signature_2 = sigs[cell_type_2][f"{cell_type_2}_top10"]
    signatures = {cell_type_1: list(signature_1), cell_type_2: list(signature_2)}
    noise_genes = fig3_noise_genes(adata, signature_1, signature_2)

    # parameters of the base cache (run_gsdensity.R defaults); one MCA embedding for the section
    GSD_PARAMS = {"nn_use": 300, "restart": 0.75, "dims": 10, "seed": 0, "spatial": True, "spatial_n": 10,
                  "spatial_ntimes": 100}
    gsd_sets = {}
    for ct in [cell_type_1, cell_type_2]:
        gsd_sets[f"{ct}_noisy"] = signatures[ct] + noise_genes
        gsd_sets[f"{ct}_top10"] = signatures[ct]
    assert all(g in adata.var_names for v in gsd_sets.values() for g in v)

    def _task_done(task_dir: str, genesets: dict, params: dict) -> bool:
        """True if scores, labels and timing exist and the stored inputs match."""
        need = ["gsd_scores.csv", "gsd_labels.csv", "gsd_timing.json", "genesets.json", "params.json", "spots.txt"]
        if not all(os.path.exists(task_dir + f) for f in need) or os.path.exists(task_dir + "gsd_error.txt"):
            return False
        same_sets = json.load(open(task_dir + "genesets.json")) == genesets
        same_params = json.load(open(task_dir + "params.json")) == params
        same_spots = open(task_dir + "spots.txt").read().split("\n")[:-1] == list(adata.obs_names)
        return same_sets and same_params and same_spots

    task_dir = gsd_dir + "_task/"
    if _task_done(task_dir, gsd_sets, GSD_PARAMS):
        print(f"fig3: GSDensity cached ({task_dir})")
    else:
        t0 = time.time()
        ru.export_gsd_task(adata, gsd_sets, task_dir, params=GSD_PARAMS)
        rc = ru.run_gsd([task_dir], n_jobs=1, log=LOGS + "REV_Fig3_gsd_R.log")
        print(f"fig3: GSDensity run: exit codes {rc}, {time.time() - t0:.1f} s")

    gsd = ru.read_gsd(task_dir, spots=adata.obs_names)
    assert not gsd["scores"].isna().any().any()
    assert np.allclose(gsd["scores"][list(gsd_sets)].sum(axis=0), 1.0)

    # the signal-only runs must equal the base cache (same section and parameters)
    for ct in [cell_type_1, cell_type_2]:
        base = pd.read_csv(gsd_base_dir + f"{ct}_top10.csv").set_index("spot_id").reindex(adata.obs_names)
        diff = float(np.abs(gsd["scores"][f"{ct}_top10"].values - base["relevance"].values).max())
        print(f"fig3: {ct}_top10 max |new - base cache| relevance = {diff:.2e}")

    # one file per gene set (temporary name, then rename) and an index of this key
    index_rows = []
    for name, genes in gsd_sets.items():
        tab = pd.DataFrame({
            "spot_id": adata.obs_names,
            "relevance": gsd["scores"][name].values,
            "label": gsd["labels"][name].reindex(adata.obs_names).values,
        })
        tab["relevance_scaled"] = tab["relevance"] * len(tab)
        tmp = gsd_dir + f"{name}.csv.tmp"
        tab.to_csv(tmp, index=False)
        os.replace(tmp, gsd_dir + f"{name}.csv")
        index_rows.append({
            "dataset": "brain", "signature": name, "used_in": "Fig. 3",
            "genes_used": len(genes), "genes": ";".join(genes), "n_spots": adata.n_obs,
            "n_positive_label": int((tab["label"] == "positive").sum()), "n_label_NA": int(tab["label"].isna().sum()),
            "subsample_note": "none: all spots, one MCA embedding for the section",
            "runtime_s": gsd["timing"]["total"], "params": json.dumps(GSD_PARAMS),
            "file": f"data/gsd_cache/brain/fig3/{name}.csv",
        })
    pd.DataFrame(index_rows).to_csv(gsd_dir + "INDEX_fig3.csv", index=False)
    print(f"fig3: wrote {len(gsd_sets)} gene sets to {gsd_dir}")


# ---------------------------------------------------------------------------
# ed5d: GSDensity runtime (REV_SupFig5_6_7_runtime, step 3, first run on 29 Oct 2026)
# ---------------------------------------------------------------------------

def part_ed5d(run_tag="20261029", n_repeats=3):
    import pickle

    gsd_dir = GSD_BASE_DIR + f"SupFig5_6_7/runtime_{run_tag}/"
    os.makedirs(gsd_dir, exist_ok=True)
    adata, _ = ru.load_brain()  # log-normalised as in 06_performance (no gene filter)
    regions = list(REGIONS)
    signature_dicts = {}
    for ct in regions:
        with open(ru.PROJECT + f"{ct}_signatures.pkl", "rb") as f:
            signature_dicts[ct] = pickle.load(f)
    sizes = {ct: [len(v) for v in signature_dicts[ct].values()] for ct in regions}
    gene_sets = {ct: {len(v): list(v) for v in signature_dicts[ct].values()} for ct in regions}

    gsd_params = {"nn_use": 300, "restart": 0.75, "dims": 10, "seed": 0, "spatial": True, "spatial_n": 10,
                  "spatial_ntimes": 100}
    base_task = gsd_dir + "_data/"
    if not os.path.exists(base_task + "data.mtx"):
        tmp = base_task.rstrip("/") + "_tmp/"
        ru.export_gsd_task(adata, {}, tmp, gsd_params)
        os.rename(tmp, base_task)

    gsd_tasks = []
    for rep in range(1, n_repeats + 1):
        for ct in regions:
            for n in sizes[ct]:
                td = gsd_dir + f"rep{rep}_{ct}_top{n}/"
                os.makedirs(td, exist_ok=True)
                for fn in ("data.mtx", "genes.txt", "spots.txt", "coords.csv"):
                    if not os.path.exists(td + fn):
                        os.symlink(base_task + fn, td + fn)
                with open(td + "genesets.json", "w") as f:
                    json.dump({f"{ct}_top{n}": gene_sets[ct][n]}, f)
                with open(td + "params.json", "w") as f:
                    json.dump(gsd_params, f)
                gsd_tasks.append(td)
    todo = [t for t in gsd_tasks if not os.path.exists(t + "gsd_timing.json")]
    print(f"ed5d: GSDensity tasks: {len(gsd_tasks)}, to run: {len(todo)}")
    if todo:
        # one R process runs the tasks one after another; each task times its own steps (gsd_timing.json)
        rc = ru.run_gsd(todo, n_jobs=1, log=LOGS + f"gsd_runtime_{run_tag}.log")
        assert rc == [0], f"GSDensity failed (see {LOGS}gsd_runtime_{run_tag}_0.log)"
    for t in gsd_tasks:
        assert not os.path.exists(t + "gsd_error.txt"), open(t + "gsd_error.txt").read()
        assert os.path.exists(t + "gsd_timing.json"), t
    print(f"ed5d: timings in {gsd_dir}")


PARTS = {"f1": part_f1, "fig2": part_fig2, "fig3": part_fig3, "ed5d": part_ed5d}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("parts", nargs="+", choices=list(PARTS) + ["all"])
    args = ap.parse_args()
    for p in (list(PARTS) if "all" in args.parts else args.parts):
        PARTS[p]()
