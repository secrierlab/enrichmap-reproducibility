import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
import gseapy as gp
import enrichmap as em
import warnings


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
