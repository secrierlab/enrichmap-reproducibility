#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
os.environ["PYTHONWARNINGS"] = "ignore"


# In[2]:


import warnings
warnings.filterwarnings("ignore")


# In[ ]:


import scanpy as sc
import squidpy as sq
import anndata
import pandas as pd
import numpy as np
import gseapy as gp
import decoupler as dc
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.colors import ListedColormap
from scipy.stats import pearsonr
from skimage.draw import disk
from scipy.spatial import distance


# In[4]:


sc.set_figure_params(dpi=100, facecolor=None, color_map="seismic", frameon=False, vector_friendly=True)
sc.settings._vector_friendly = False


# In[5]:


# set working directory
project_dir = "/Users/cenkcelik/Cenk_scoring/"
working_dir = project_dir + ""
os.chdir(working_dir)

# set figure directory
figure_dir = working_dir + "figures/"

# processed data directory
processed_data = working_dir + "processed_data/"


# In[6]:


# import a local package
import sys
sys.path.append("/Users/cenkcelik/Documents/GitHub/EnrichMap/")
import enrichmap as em


# In[7]:


np.random.seed(0)


# In[ ]:


# Load human gene names from a reference file or database
def load_human_genes(n_genes=15000):
    gene_list = pd.read_csv("https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym&format=text", sep="\t")
    return np.random.choice(gene_list["Approved symbol"].values, n_genes, replace=False)

# Generate hexagonal grid for 10X Visium-like spatial coordinates
def generate_hexagonal_grid(rows=70, cols=70, spacing=50):
    x_coords, y_coords = [], []
    for row in range(rows):
        for col in range(cols):
            x = col * spacing + (row % 2) * (spacing / 2)
            y = row * spacing * np.sqrt(3) / 2
            x_coords.append(x)
            y_coords.append(y)
    return np.array(x_coords), np.array(y_coords)

# Assign biological labels in spatially coherent clusters
def assign_spatial_clusters(x_coords, y_coords, labels):
    num_spots = len(x_coords)
    label_assignment = np.random.choice(labels, num_spots, replace=True)
    
    for i in range(num_spots):
        nearest_neighbours = distance.cdist([(x_coords[i], y_coords[i])], np.column_stack((x_coords, y_coords))).flatten()
        close_spots = np.argsort(nearest_neighbours)[1:4]
        label_assignment[close_spots] = label_assignment[i] # Increase local clustering
    
    return label_assignment

# Simulate expression matrix
def generate_expression_matrix(n_spots, n_genes, labels, gene_signature):
    expression_matrix = np.random.poisson(5, (n_spots, n_genes))
    
    # Boost state gene expression in made up-labelled spots
    made_up_spots = np.where(labels == "cluster_2")[0]
    expression_matrix[made_up_spots[:, None], gene_signature] *= 5
    
    return expression_matrix

# generate a number between 0-255
def random_color():
    return np.random.randint(0, 256)

# Generate a dummy tissue image
def generate_tissue_image(image_size=(2000, 2000)):
    image = np.ones((*image_size, 3), dtype=np.uint8) * 252
    rr, cc = disk((500, 500), 300, shape=image.shape[:2])
    image[rr, cc] = [random_color(), random_color(), random_color()]
    return image

# Main function to generate simulated spatial AnnData
def generate_spatial_cancer_data(seed=0):
    np.random.seed(seed)
    n_genes = 15000
    x_coords, y_coords = generate_hexagonal_grid()
    n_spots = len(x_coords)
    gene_names = load_human_genes(n_genes)
    
    # Define labels
    labels = np.array(["cluster_1", "cluster_2", "cluster_3", "cluster_4"])
    assigned_labels = assign_spatial_clusters(x_coords, y_coords, labels)
    
    # Pick Made-up spots-related genes
    spot_gene_indices = np.random.choice(n_genes, size=100, replace=False)
    
    # Generate expression data
    expression_matrix = generate_expression_matrix(n_spots, n_genes, assigned_labels, spot_gene_indices)
    
    # Create AnnData object
    adata = anndata.AnnData(X=expression_matrix)
    adata.var_names = gene_names
    adata.obs["cluster"] = assigned_labels
    adata.obsm["spatial"] = np.column_stack((x_coords, y_coords))
    
    # Add Visium-like metadata
    adata.uns["spatial"] = {
        "image_1": {
            "images": {
                "hires": generate_tissue_image(),
                "lowres": generate_tissue_image((500, 500))
            },
            "metadata": {
                "chemistry_description": "Spatial 3' v1",
                "software_version": "spaceranger-1.2.2"
            },
            "scalefactors": {
                "fiducial_diameter_fullres": 28.89725002221966,
                "spot_diameter_fullres": 17.888773823278836,
                "tissue_hires_scalef": 0.856531,
                "tissue_lowres_scalef": 0.25695932
            }
        }
    }
    
    return adata


# In[ ]:


adata = generate_spatial_cancer_data(seed=0)


# In[10]:


cluster_palettes = {
    "cluster_1": "#F6C141",
    "cluster_2": "#882E72",
    "cluster_3": "#7BAFDE",
    "cluster_4": "#F4A582"
}

# Convert 'cluster' column to categorical dtype
adata.obs["cluster"] = pd.Categorical(adata.obs["cluster"])

# Now you can access the categories
clusters = adata.obs["cluster"].cat.categories

# Reorder the palette list to match the category order
palette_list = [cluster_palettes[c] for c in clusters]

# Create a ListedColormap
cmap = ListedColormap(palette_list)


# In[11]:


sq.pl.spatial_scatter(adata, color="cluster", img_alpha=0, size=3, palette=cmap, save=figure_dir + "spatial_scatter_simulated_clusters.pdf")


# In[12]:


sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)


# In[13]:


sc.tl.rank_genes_groups(adata, groupby="cluster", method="t-test_overestim_var", use_raw=False)
sc.tl.filter_rank_genes_groups(adata, min_in_group_fraction=0.5, min_fold_change=2, max_out_group_fraction=0.5)


# In[14]:


gene_set_10 = adata.uns["rank_genes_groups"]["names"]["cluster_2"][:10].tolist()


# ### EnrichMap

# In[ ]:


em.tl.score(
    adata,
    gene_set=gene_set_10,
    smoothing=True,
    correct_spatial_covariates=True,
    batch_key=None
)


# In[16]:


sq.pl.spatial_scatter(adata, color=["enrichmap_score"], img_alpha=0, size=3, vcenter=0, save=figure_dir + "spatial_scatter_enrichmap_score.pdf")


# ### AUCell

# In[ ]:


# convert the signature dict to a long format
signatures_dict_10 = {
    "cluster_2": gene_set_10
}

signatures_10 = pd.DataFrame(
    [(key, gene) for key, genes in signatures_dict_10.items() for gene in genes],
    columns=["geneset", "genesymbol"]
)


dc.run_aucell(
    adata,
    net=signatures_10,
    source="geneset",
    target="genesymbol",
    seed=0,
    verbose=True,
    use_raw=False
)



acts = dc.get_acts(adata, obsm_key="aucell_estimate")

adata.obs["aucell_score"] = acts.obsm["aucell_estimate"]["cluster_2"]


# ## `score_genes()`/`AddModuleScore()`

# In[ ]:


sc.tl.score_genes(adata, gene_set_10, score_name="scanpy_score", use_raw=False)


# ## ssGSEA

# In[20]:


# If adata.X is sparse, convert to dense first
X_dense = adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X

# Create DataFrame using var names as columns and obs names as index
gene_expr = pd.DataFrame(X_dense.T, index=adata.var_names, columns=adata.obs_names)


# In[21]:


gene_expr.index = gene_expr.index.str.upper()
signatures_dict_10 = {k: [g.upper() for g in v] for k, v in signatures_dict_10.items()}


# In[ ]:


ssgsea = gp.ssgsea(
    data=gene_expr,
    gene_sets=signatures_dict_10,
    outdir=None,
    sample_norm_method="rank",
    no_plot=True,
    min_size=5
)


# In[23]:


nes = ssgsea.res2d.pivot(index="Term", columns="Name", values="NES")
adata.obs["ssgsea_score"] = nes.loc["cluster_2"].T.astype(float).reindex(adata.obs_names)


# ### GSVA

# In[ ]:


gsva = gp.gsva(
    data=gene_expr,
    gene_sets=signatures_dict_10,
    outdir=None,
    min_size=5
)


# In[25]:


es = gsva.res2d.pivot(index="Term", columns="Name", values="ES")
adata.obs["gsva_score"] = es.loc["cluster_2"].T.astype(float).reindex(adata.obs_names)


# ### Z score

# In[ ]:


gene_set = signatures_dict_10["cluster_2"]
z_score = gene_expr.loc[gene_set].sum(axis=0) / np.sqrt(len(gene_set))
adata.obs["z_score"] = z_score


# In[27]:


methods = ["enrichmap_score", "aucell_score", "z_score", "scanpy_score", "ssgsea_score", "gsva_score"]


# In[29]:


em.pl.spatial_metrics(adata, methods, metric="Moran's I", figsize=(3, 2), save=figure_dir + "spatial_metrics_simulated_moransI.pdf")


# In[30]:


em.pl.spatial_metrics(adata, methods, metric="Geary's C", figsize=(3, 2), save=figure_dir + "spatial_metrics_simulated_gearysC.pdf")


# In[31]:


em.pl.variogram(
    adata,
    score_keys=methods,
    save="variogram_simulated_method_comparison.pdf",
)


# In[32]:


em.pl.morans_correlogram(
    adata,
    score_key="enrichmap_score",
    save="morans_correlogram_enrichmap_score.pdf",
)


# ## Gene set sensitivity

# In[33]:


gene_set_2 = adata.uns["rank_genes_groups"]["names"]["cluster_2"][:2].tolist()
gene_set_5 = adata.uns["rank_genes_groups"]["names"]["cluster_2"][:5].tolist()
gene_set_10 = adata.uns["rank_genes_groups"]["names"]["cluster_2"][:10].tolist()
gene_set_20 = adata.uns["rank_genes_groups"]["names"]["cluster_2"][:20].tolist()
gene_set_50 = adata.uns["rank_genes_groups"]["names"]["cluster_2"][:50].tolist()
gene_set_100 = adata.uns["rank_genes_groups"]["names"]["cluster_2"][:100].tolist()
gene_set_200 = adata.uns["rank_genes_groups"]["names"]["cluster_2"][:200].tolist()
gene_set_500 = adata.uns["rank_genes_groups"]["names"]["cluster_2"][:500].tolist()


# In[34]:


# generate a dict
gene_set_sensitivity_dict = {
    "N_genes=2": gene_set_2,
    "N_genes=5": gene_set_5,
    "N_genes=10": gene_set_10,
    "N_genes=20": gene_set_20,
    "N_genes=50": gene_set_50,
    "N_genes=100": gene_set_100,
    "N_genes=200": gene_set_200,
    "N_genes=500": gene_set_500
}


# In[35]:


score_keys = [
    "N_genes=2_score",
    "N_genes=5_score",
    "N_genes=10_score",
    "N_genes=20_score",
    "N_genes=50_score",
    "N_genes=100_score",
    "N_genes=200_score",
    "N_genes=500_score"
]


# In[36]:


em.tl.score(
    adata,
    gene_set=gene_set_sensitivity_dict,
    smoothing=True,
    correct_spatial_covariates=True,
    batch_key=None
)


# In[37]:


em.pl.spatial_enrichmap(
    adata,
    score_key=score_keys,
    img_alpha=0,
    size=3,
    ncols=4,
    save=figure_dir + "spatial_enrichmap_gene_set_sensitivity_score_comparison.pdf",
)


# In[38]:


em.pl.spatial_metrics(adata, methods, metric="Moran's I", figsize=(3,2), save=figure_dir + "spatial_metrics_gene_set_sensitivity_moransI.pdf")


# In[39]:


em.pl.spatial_metrics(adata, methods, metric="Geary's C", figsize=(3,2), save=figure_dir + "spatial_metrics_gene_set_sensitivity_gearysC.pdf")


# In[40]:


# List of gene score labels to compare with the reference
scores = ["2", "5", "10", "20", "100", "200", "500"]
x_label = "N_genes=50_score"

# Set up the figure with subplots
n_plots = len(scores)
n_cols = 3
n_rows = (n_plots + n_cols - 1) // n_cols
fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 4))
axes = axes.flatten()

for i, score in enumerate(scores):
    y_label = f"N_genes={score}_score"
    df = adata.obs[[x_label, y_label]].dropna()

    corr, pval = pearsonr(df[x_label], df[y_label])

    ax = axes[i]
    sns.regplot(x=x_label, y=y_label, data=df,
                scatter_kws={"s": 10, "color": "lightblue"},
                line_kws={"color": "black"},
                ax=ax)

    ax.set_title(f"{y_label} vs {x_label}")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.text(0.05, 0.95,
            f"r = {corr:.4f}\np = {pval:.1e}",
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment="top",
            bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.3"))
    ax.grid(False)

# Hide any unused subplots
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.savefig(figure_dir + "spatial_enrichmap_gene_set_sensitivity_score_correlation.pdf", bbox_inches="tight")
plt.show()


# In[41]:


em.pl.gene_contributions_heatmap(
    adata,
    score_key="N_genes=20",
    top_n_genes=5,
    figsize=(2, 2),
    save=figure_dir + "spatial_enrichmap_gene_set_sensitivity_top_5_contributions.pdf"
)


# In[42]:


em.pl.gene_contributions_pca(
    adata,
    score_key="N_genes=20",
    top_n_genes=3,
    figsize=(3, 3),
    save=figure_dir + "spatial_enrichmap_gene_contributions_pca.pdf"
)


# In[43]:


sq.pl.spatial_scatter(
    adata,
    color=["RNU6-935P", "CFAP47", "AOC3"],
    size=3,
    use_raw=False,
    img_alpha=0,
    ncols=3,
    cmap="PuOr",
    save=figure_dir + "spatial_scatter_simulated_genes.pdf",
)


# In[ ]:


em.pl.variogram(
    adata,
    score_keys=score_keys,
    save="variogram_simulated_geneset_comparison.pdf",
)


# In[ ]:


em.pl.morans_correlogram(
    adata,
    score_key="N_genes=20_score",
    save="morans_correlogram_simulated_geneset_enrichmap_score.pdf",
)


# In[ ]:


signature_name = "enrichmap"
binary_labels = em.tl.generate_binary_labels(adata, signature_name)
methods = [
    "enrichmap_score",
    "aucell_score",
    "z_score",
    "scanpy_score",
    "ssgsea_score",
    "gsva_score"
]
f1_scores = em.tl.compute_f1_scores(adata, methods, binary_labels)


# In[ ]:


f1_scores = em.tl.compute_f1_scores(adata, score_keys, binary_labels)


# In[ ]:


em.pl.variogram_all(
    adata,
    score_keys=methods,
    save="variogram_all_simulated_method_comparison.pdf",
)

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt

from anndata import AnnData
from typing import List

from scipy.spatial.distance import pdist
from skgstat import Variogram

plt.rcParams["pdf.fonttype"] = "truetype"

def variogram(
    adata: AnnData,
    score_keys: List[str],
    save: None | str = None,
    max_lag: float | None = None,
    lag_percentile: float = 95
) -> None:
    """
    Compute empirical variograms to assess spatial dependence for any score key.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with spatial coordinates in `adata.obsm["spatial"]`.
    score_keys : list of str
        List of keys in `adata.obs` to compute variograms for.
    save : str or None, optional
        If provided, path to save the figure as a PDF file.
    max_lag : float or None, optional
        If set, limits the x-axis of the variogram plot to this value.
        If None, computed from pairwise distances using `lag_percentile`.
    lag_percentile : float, optional
        Percentile of pairwise distances to set as max_lag when max_lag is None (default: 95).
    
    Returns
    -------
    variograms : list of Variogram
        List of computed Variogram objects for each score key.
    """
    coords = adata.obsm["spatial"]

    if max_lag is None:
        dists = pdist(coords)
        max_lag = np.percentile(dists, lag_percentile)

    n = len(score_keys)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 4), constrained_layout=True)
    if n == 1:
        axes = [axes]

    variograms = []
    for ax, key in zip(axes, score_keys):
        values = adata.obs[key].values
        V = Variogram(coords, values, method="cressie", model="gaussian", verbose=True)
        variograms.append(V)

        ax.plot(V.bins, V.experimental, "o-", label="Experimental variogram")
        ax.axhline(np.var(values), color="r", linestyle="--", label="Variance")
        ax.set_xlabel("Spatial lag")
        ax.set_ylabel("Semivariance")
        ax.set_title(f"{key}")
        ax.legend()
        ax.grid(False)
        ax.set_xlim(0, max_lag)


    if save:
        # Ensure 'figures/' directory exists
        os.makedirs("figures", exist_ok=True)

        # If 'save' has no directory path, prepend 'figures/'
        if not os.path.dirname(save):
            save = os.path.join("figures", save)
        plt.savefig(save, dpi=300, bbox_inches="tight")
    plt.show()


from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt

from anndata import AnnData
from typing import List
from scipy.spatial.distance import pdist
from skgstat import Variogram

plt.rcParams["pdf.fonttype"] = "truetype"

def variogram_all(
    adata: AnnData,
    score_keys: List[str],
    save: None | str = None,
    max_lag: float | None = None,
    lag_percentile: float = 95
) -> None:
    """
    Compute and plot empirical variograms for multiple score keys on the same plot.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with spatial coordinates in `adata.obsm["spatial"]`.
    score_keys : list of str
        List of keys in `adata.obs` to compute variograms for.
    save : str or None, optional
        If provided, path to save the figure as a PDF file.
    max_lag : float or None, optional
        If set, limits the x-axis of the variogram plot to this value.
        If None, computed from pairwise distances using `lag_percentile`.
    lag_percentile : float, optional
        Percentile of pairwise distances to set as max_lag when max_lag is None (default: 95).
    """
    coords = adata.obsm["spatial"]

    if max_lag is None:
        dists = pdist(coords)
        max_lag = np.percentile(dists, lag_percentile)

    fig, ax = plt.subplots(figsize=(4, 4), constrained_layout=True)

    colours = plt.cm.tab10(np.linspace(0, 1, len(score_keys)))
    variograms = []

    for key, colour in zip(score_keys, colours):
        values = adata.obs[key].values
        V = Variogram(coords, values, method="cressie", model="gaussian", verbose=False)
        variograms.append(V)

        ax.plot(V.bins, V.experimental, "o-", label=key, color=colour)
        ax.axhline(np.var(values), color=colour, linestyle="--", alpha=0.5)

    ax.set_xlabel("Spatial lag")
    ax.set_ylabel("Semivariance")
    ax.set_title("Empirical variograms")
    ax.legend(title="Signatures")
    ax.grid(False)
    ax.set_xlim(0, max_lag)

    if save:
        # Ensure 'figures/' directory exists
        os.makedirs("figures", exist_ok=True)

        # If 'save' has no directory path, prepend 'figures/'
        if not os.path.dirname(save):
            save = os.path.join("figures", save)
        plt.savefig(save, dpi=300, bbox_inches="tight")
    plt.show()


    from __future__ import annotations

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from ..tools._compute_spatial_metrics import compute_spatial_metrics
from anndata import AnnData
from collections.abc import Sequence

plt.rcParams["pdf.fonttype"] = "truetype"


def spatial_metrics(
    adata: AnnData,
    score_keys: Sequence[str],
    metric: str = "Moran's I" or "Geary's C",
    n_neighbors: int = 6,
    figsize: tuple[int, int] = (4, 4),
    save=None
) -> None:
    """
    Compute and visualise spatial metrics for different scoring methods in a given dataset.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing spatial and gene expression information.
    score_keys : sequence of str
        A list of method names for which to compute the spatial metric. These methods should correspond to
        columns in `adata.obs`.
    metric : str, optional
        The spatial metric to compute, e.g., "Moran's I" or "Geary's C". Defaults to "Moran's I".
    n_neighs : int, optional
        Number of neighbours to use when computing spatial weights. Defaults to 6.
    figsize : tuple of int, optional
        The size of the figure to be generated for the bar plot. Defaults to (4, 4).
    save : str or None, optional
        If specified, the plot will be saved to the file with the given filename. If None, the plot will not be saved.

    Returns
    -------
    None
    """
    results = []

    for method in score_keys:
        if method in adata.obs.columns:
            print(f"Computing {metric} for {method} with {n_neighbors} neighbours...")
            metrics, _ = compute_spatial_metrics(adata, score_key=method, n_neighbors=n_neighbors)
            results.append((method, metrics[metric]))

    df = pd.DataFrame(results, columns=["Method", metric])

    plt.figure(figsize=figsize)
    sns.barplot(data=df, x="Method", y=metric, palette="muted", edgecolor=None)

    plt.axhline(0, linestyle="--", color="gray", linewidth=1)
    plt.title(f"{metric}", fontsize=8)
    plt.ylabel("Continuity score", fontsize=6)
    plt.xlabel("Scoring methods", fontsize=6)
    plt.yticks(fontsize=6)
    plt.xticks(rotation=90, fontsize=6)
    plt.grid(False)

    if save:
        os.makedirs("figures", exist_ok=True)
        if not os.path.dirname(save):
            save = os.path.join("figures", save)
        plt.savefig(save, dpi=300, bbox_inches="tight")
    plt.show()

from __future__ import annotations

import squidpy as sq
import numpy as np

from anndata import AnnData
from libpysal.weights import KNN
from esda.moran import Moran, Moran_Local
from esda.geary import Geary
from esda.getisord import G_Local


def compute_spatial_metrics(
    adata: AnnData,
    score_key: str = "enrichment_score",
    n_neighbors: int = 6,
):
    """
    Compute global and local spatial autocorrelation metrics for a spatial score.

    This function calculates several spatial statistics to assess spatial autocorrelation in gene set
    enrichment or similar scores stored in `adata.obs[score_key]`. It returns Moran's I (global),
    Geary's C (local variance), Local Moran's I (spatial clusters), and Getis-Ord G (hotspots).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with spatial coordinates in `adata.obsm['spatial']`.
    score_key : str, default "enrichment_score"
        Column in `adata.obs` containing the score vector for which spatial metrics should be computed.
    n_neighbors : int, default 6
        Number of spatial neighbours to use for spatial weights calculation.

    Returns
    -------
    metrics : dict
        Dictionary containing:
            - "Moran's I": float
            - "Geary's C": float
            - "Local Moran I": np.ndarray of local Moran’s I values

    W : libpysal.weights.W
        Spatial weights matrix used for computation.

    local_moran : esda.Moran_Local
        Full Local Moran’s I result object.
    """
    # Compute spatial neighbors and weights
    sq.gr.spatial_neighbors(adata, n_neighs=n_neighbors, coord_type="generic", key_added="spatial")

    # Extract spatial weights from coordinates
    W = KNN.from_array(adata.obsm["spatial"])
    W.transform = "r"

    # Extract score vector
    score = adata.obs[score_key].values.astype(np.float64)

    # Compute global Moran's I
    moran = Moran(score, W)

    # Compute Geary's C
    geary = Geary(score, W)

    # Compute Local Moran's I
    local_moran = Moran_Local(score, W)

    metrics = {
        "Moran's I": moran.I,
        "Geary's C": geary.C,
        "Local Moran I": local_moran.Is,
    }

    return metrics, W


from __future__ import annotations

from sklearn.metrics import f1_score
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from typing import Sequence
from anndata import AnnData

def compute_f1_scores(
    adata: AnnData,
    score_keys: Sequence[str],
    binary_labels: np.ndarray,
    figsize: tuple[int, int] = (3, 3),
    save: str | None = None
) -> dict[str, float]:
    """
    Compute and plot F1 scores for multiple gene set scoring methods.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with scoring results in `adata.obs`.
    score_keys : list of str
        List of method score keys to benchmark.
    binary_labels : np.ndarray
        Ground truth binary labels.
    figsize : tuple, optional
        Size of the figure.
    save : str or None, optional
        If given, path to save the figure.

    Returns
    -------
    f1_scores : dict
        Dictionary of F1 scores per `score_keys` provided.
    """
    f1_scores = {}

    for method in score_keys:
        if method in adata.obs.columns:
            preds = adata.obs[method] > np.median(adata.obs[method])
            f1_scores[method] = f1_score(binary_labels, preds)

    # Plot
    plt.figure(figsize=figsize)
    sns.barplot(x=list(f1_scores.keys()), y=list(f1_scores.values()), palette="muted", edgecolor=None)
    plt.title("F1 scores", fontsize=10)
    plt.ylabel("F1 score", fontsize=8)
    plt.xlabel("Method", fontsize=8)
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.yticks(fontsize=8)
    plt.grid(False)

    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")
    plt.show()

    return f1_scores
