# %%
import os, sys
import logging
import warnings
warnings.filterwarnings("ignore")

DATASETS = {'chicken':{'annotation': "newannotation_high", 'cluster': "seurat_clusters",
                       'scientific_name': 'Ggal', 'root_cell_type': 'PROG', 'time': 'Stage',
                       'start_time': 6.0}, 
            'salamander':{'annotation': "newannotation_final_hires", 'cluster': "seurat_clusters",
                          'scientific_name': 'Pwal', 'root_cell_type': 'LGE', 'time': 'Stage',
                          'start_time': 30.0}}

# ========== 1. Specify dataset and Configure paths ==========
dataset = 'salamander'
data_dir = "./Data"
output_dir = f"./{dataset}_output"
os.makedirs(output_dir, exist_ok=True)
FILE_NAME = DATASETS[dataset]['scientific_name']
ANNOTATION_KEY = DATASETS[dataset]["annotation"]
CLUSTER_KEY = DATASETS[dataset]["cluster"]
TIME_KEY = DATASETS[dataset]['time']
PROG = DATASETS[dataset]['root_cell_type']

# ========== 2. Set up logging ==========
log_file = os.path.join(output_dir, f"{dataset}_traj.log")
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logging.info(f"❗Starting {dataset} trajectory inference script...")

# ========== 3. Read in dataset ==========
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
np.random.seed(42)

# %%
# ========== PAGA ==========
adata_hvg = adata[:, adata.var['highly_variable']==True].copy()
sc.pp.normalize_total(adata_hvg, target_sum=1e4, layer='counts')
sc.pp.log1p(adata_hvg)
sc.pp.scale(adata_hvg)

logging.info("Running PAGA...")
sc.pp.neighbors(adata_hvg, n_neighbors=5, n_pcs=30, use_rep='X_pca')
sc.tl.draw_graph(adata_hvg, init_pos='X_umap')
sc.pl.draw_graph(adata_hvg, color=CLUSTER_KEY, legend_loc='on data', legend_fontsize = 'xx-small')
plt.savefig(os.path.join(output_dir, "PAGA_annot_sc_FA.png"), dpi=300, bbox_inches="tight")
sc.tl.paga(adata_hvg, groups=CLUSTER_KEY)
sc.pl.paga(adata_hvg, color=[ANNOTATION_KEY, CLUSTER_KEY], edge_width_scale=0.2)
plt.savefig(os.path.join(output_dir, "PAGA_annot_nodes_FA.png"), dpi=300, bbox_inches="tight")

sc.pl.umap(adata_hvg, edges=True, color = ANNOTATION_KEY, legend_loc='on data', legend_fontsize= 'xx-small')
plt.savefig(os.path.join(output_dir, "PAGA_annot_sc_UMAP.png"), dpi=300, bbox_inches="tight")

adata_hvg.write(os.path.join(data_dir, f"{FILE_NAME}_hvg_PAGA.h5ad"))

logging.info("Algorithm PAGA complete!")

# %%
# ========== CellRank: Combining MOSCOT with DPT ==========
import pickle
import palantir
import scvelo as scv
import cellrank as cr
from sklearn.metrics import silhouette_score
from cellrank.estimators import GPCCA
from cellrank.kernels import RealTimeKernel
from cellrank.kernels import PseudotimeKernel
from cellrank.kernels import CytoTRACEKernel

# %%
# ========== MOSCOT ==========
from moscot.problems.time import TemporalProblem

adata = sc.read(os.path.join(data_dir, f"{FILE_NAME}.h5ad"))
sc.pp.normalize_total(adata, layer='counts')
sc.pp.log1p(adata)
sc.pp.neighbors(adata, use_rep='X_pca')

logging.info("Running CellRank...")
if dataset == 'chicken':
    adata.obs[TIME_KEY] = adata.obs[TIME_KEY].str.replace('E', '', regex=False)
adata.obs[TIME_KEY] = adata.obs[TIME_KEY].astype(float).astype("category")
adata.obs["Stage_numerical"] = adata.obs[TIME_KEY].astype(float)

logging.info("Computing RealTimeKernel...")
tp = TemporalProblem(adata)
tp = tp.score_genes_for_marginals(gene_set_proliferation="human", gene_set_apoptosis="human")
sc.pl.umap(adata, color=["proliferation", "apoptosis"])
plt.savefig(os.path.join(output_dir, "CellRank_birth_death_UMAP.png"), dpi=300, bbox_inches="tight")
tp = tp.prepare(time_key=TIME_KEY)
tp = tp.solve(epsilon=1e-3, tau_a=0.95, scale_cost="mean", device="cpu")
tmk = RealTimeKernel.from_moscot(tp)
tmk.compute_transition_matrix(self_transitions="all", conn_weight=0.5, threshold="auto")
tmk.plot_random_walks(
    max_iter=500,
    start_ixs={TIME_KEY: 30.0},
    basis="umap",
    seed=0,
    dpi=150,
    size=30,
    save=os.path.join(output_dir, "CellRank_random_walks_stage_UMAP.png")
)
tmk.write_to_adata(key="moscot_T", copy=False)
adata.write(os.path.join(data_dir, f"{FILE_NAME}_moscot.h5ad"))

# %%
# ========== DPT ==========
adata = sc.read(os.path.join(data_dir, f"{FILE_NAME}_moscot.h5ad"))
logging.info("Computing PseudotimeKernel...")
root_cell = np.random.choice(adata.obs_names[(adata.obs[ANNOTATION_KEY] == PROG)])
adata.uns['iroot'] = np.where(adata.obs_names == root_cell)[0][0]
sc.tl.diffmap(adata)
sc.tl.dpt(adata)
sc.pl.umap(adata, color = 'dpt_pseudotime')
plt.savefig(os.path.join(output_dir, "dpt_pseudotime_UMAP.png"), dpi=300, bbox_inches="tight")
pk = PseudotimeKernel(adata, time_key="dpt_pseudotime")
pk.compute_transition_matrix(threshold_scheme="soft")
pk.write_to_adata(key="pdt_T", copy=False)

# %%
# # ========== CytoTRACE ==========
# adata = sc.read(os.path.join(data_dir, f"{FILE_NAME}_CellRank.h5ad"))
# adata.layers["spliced"] = adata.X
# adata.layers["unspliced"] = adata.X
# scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
# logging.info("Computing CytoTRACE Kernel...")
# ctk = CytoTRACEKernel(adata).compute_cytotrace()
# sc.pl.umap(adata, color="ct_score", color_map="gnuplot2")
# plt.savefig(os.path.join(output_dir, "CytoTRACE_score_UMAP.png"), dpi=300, bbox_inches="tight")
# adata.write(os.path.join(data_dir, f"{FILE_NAME}_CellRank.h5ad"))

# %%
# ========== Combine kernels ==========
adata = sc.read(os.path.join(data_dir, f"{FILE_NAME}_CellRank.h5ad"))
pk = PseudotimeKernel.from_adata(adata, key="pdt_T")
tmk = RealTimeKernel.from_adata(adata, key="moscot_T")
combined_kernel = 0.7 * tmk + 0.3 * pk
combined_kernel.write_to_adata(key="combined_kernel_T", copy=False)

# %%
# ========== Computing initial and terminal states ==========
logging.info("Computing macrostates...")
g = cr.estimators.GPCCA(combined_kernel) # Use the Generalized Perron Cluster Cluster Analysis estimator
g.compute_schur(n_components=50, method='krylov', verbose=True)
g.plot_spectrum(real_only=True, 
                save=os.path.join(output_dir, "CellRank_real_eigenvals.png"), dpi=300)
g.fit(cluster_key=ANNOTATION_KEY, n_states=40)
g.plot_macrostates(which="all", legend_loc="on data", size=100, 
                   save=os.path.join(output_dir, "CellRank_all_states.png"), dpi=300)
g.plot_macrostate_composition(key=ANNOTATION_KEY, figsize=(10, 4), 
                              save=os.path.join(output_dir, "CellRank_macrostate_composition.png"), 
                              dpi=300)
g.plot_coarse_T(annotate=False, 
                save=os.path.join(output_dir, "CellRank_coarse_transitions.png"), dpi=300)

logging.info("Computing fate probabilities and driver genes...")
g.set_terminal_states(states=["OB_1", "OB_2", "OB_3", "OB_4", "OB_5", "OB_6", "OB_7", "OB_8", "OB_9", 
                              "OB_12", "PALL_2", "PALL_3", "PALL_4", "IN-MGE_1", "IN-MGE_2", "IN-MGE_3", 
                              "IN-MGE_4", "STR-MSN2_1", "STR-MSN2_2", "SEP_1", "SEP_2", "SEP_3", "SEP_4", 
                              "SEP_5", "IN-CGE", "STR-MSN1_1", "STR-MSN1_2", "STR-MSN1_3", "STR-MSN1_4", 
                              "STR-MSN1_5"])
g.set_initial_states(states=["CGE", "LGE_1", "LGE_2", "PREC_1"])
g.plot_macrostates(which="terminal", legend_loc="on data", s=100, 
                   save=os.path.join(output_dir, "CellRank_terminal_states.png"), dpi=300)
g.plot_macrostates(which="initial", legend_loc="on data", s=100, 
                   save=os.path.join(output_dir, "CellRank_initial_states.png"), dpi=300)
g.compute_fate_probabilities()
g.plot_fate_probabilities(save=os.path.join(output_dir, "CellRank_fate_probabilities.png"), dpi=300)
with open(os.path.join(data_dir, "gpcca_model.pkl"), "wb") as f:
    pickle.dump(g, f)

# %%
# ========== Visualizing and clustering gene expression trends ==========
dm_res = palantir.utils.run_diffusion_maps(adata, n_components=10)
ms_data = palantir.utils.determine_multiscale_space(adata)
imputed_X = palantir.utils.run_magic_imputation(adata)
palantir.plot.plot_diffusion_components(adata)
plt.savefig(os.path.join(output_dir, "Palantir_diffusion_components.png"), dpi=300, bbox_inches="tight")
model = cr.models.GAM(adata)
gene_sets = {}
for terminal_state in g.terminal_states.unique().dropna():
    lineage_drivers = g.compute_lineage_drivers(lineages=terminal_state)
    driver_genes = lineage_drivers.head(40).index.to_list()
    gene_sets[terminal_state] = driver_genes
    cr.pl.heatmap(
        adata,
        model=model,
        lineages=terminal_state,
        cluster_key=ANNOTATION_KEY,
        show_fate_probabilities=True,
        data_key= "MAGIC_imputed_data",
        genes=driver_genes,
        time_key="ct_score",
        figsize=(12, 10),
        show_all_genes=True,
        weight_threshold=(1e-3, 1e-3),
    )
    plt.savefig(os.path.join(output_dir, f"CellRank_{terminal_state}_gene_trend.png"),
                dpi=300, bbox_inches="tight")
pd.DataFrame(gene_sets).to_excel(os.path.join(output_dir, "CellRank_fates_top10_dexp.xlsx"), index=False)

adata.write(os.path.join(data_dir, f"{FILE_NAME}_CellRank.h5ad"))
logging.info("CellRank inference complete!")

# %%
# ========== Mapping progenitor and precursor states to terminal states ==========
adata = sc.read(os.path.join(data_dir, f"{FILE_NAME}_CellRank.h5ad"))
terminal_states = adata.obs['term_states_fwd'].unique().dropna().to_list()

prec_clusters = adata[adata.obs[ANNOTATION_KEY]=='PREC'].obs['SCT_snn_res.2'].value_counts() > 30
prec_clusters = prec_clusters[prec_clusters].index.to_list()
prec_clusters.remove('NaN')
sc.pl.umap(adata, color='SCT_snn_res.2', groups=prec_clusters, legend_loc="on data", 
           legend_fontsize = 'x-small', frameon=True, size=10)
plt.savefig(os.path.join(output_dir, "prec_clusters_UMAP.png"), dpi=300, bbox_inches="tight")
cr.pl.aggregate_fate_probabilities(adata, mode="bar", lineages=terminal_states, 
                                   cluster_key='SCT_snn_res.2', clusters=prec_clusters, ncols = 4)
plt.savefig(os.path.join(output_dir, "CellRank_prec_fate_probs.png"), dpi=300, bbox_inches="tight")

prog_clusters = adata[adata.obs['origin']=='progenitors'].obs['SCT_snn_res.2'].value_counts() > 30
prog_clusters = prog_clusters[prog_clusters].index.to_list()
prog_clusters.remove('NaN')
sc.pl.umap(adata, color='SCT_snn_res.2', groups=prog_clusters, legend_fontsize = 'xx-small', frameon=True, size=10)
plt.savefig(os.path.join(output_dir, "prog_clusters_UMAP.png"), dpi=300, bbox_inches="tight")
cr.pl.aggregate_fate_probabilities(adata, mode="bar", lineages=terminal_states, 
                                   cluster_key='SCT_snn_res.2', clusters=prog_clusters, ncols = 4)
plt.savefig(os.path.join(output_dir, "CellRank_prog_fate_probs.png"), dpi=300, bbox_inches="tight")

# %%
# ========== scFates ==========
os.environ['R_HOME'] = sys.exec_prefix+"/lib/R/"
import seaborn
import palantir
import scFates as scf
from anndata import AnnData
seaborn.reset_orig()
sc.set_figure_params()
scf.set_figure_pubready()
sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout

adata = sc.read(os.path.join(data_dir, f"{FILE_NAME}_CellRank.h5ad"))
scf.tl.cellrank_to_tree(adata, method="ppt", Nodes=400, time="dpt_pseudotime",
                        device="gpu", seed=1, ppt_lambda=1, ppt_sigma=0.025, ppt_nsteps=400)
scf.pl.graph(adata)
plt.savefig(os.path.join(output_dir, "scFates_Cellrank_tree_rep_UMAP.png"), dpi=300, bbox_inches="tight")

scf.tl.root(adata, 130)
scf.tl.pseudotime(adata, n_jobs=16, n_map=200, seed=42)
scf.pl.trajectory(adata, color_cells=ANNOTATION_KEY)
plt.savefig(os.path.join(output_dir, "scFates_Cellrank_trajectory_UMAP.png"), dpi=300, bbox_inches="tight")
sc.pl.umap(adata, color=[ANNOTATION_KEY, "seg"])
plt.savefig(os.path.join(output_dir, "scFates_Cellrank_annot_seg_UMAP.png"), dpi=300, bbox_inches="tight")

scf.tl.test_association(adata, n_jobs=32)
sc.set_figure_params()
scf.pl.test_association(adata)
plt.savefig(os.path.join(output_dir, "genes_scFates_tree_associations.png"), dpi=300, bbox_inches="tight")
scf.tl.fit(adata, n_jobs=32)

scf.tl.dendrogram(adata, n_jobs=8, crowdedness=0.5)
scf.pl.dendrogram(adata, color="macrostates_fwd")
plt.savefig(os.path.join(output_dir, "scFates_Cellrank_annot_dendrogram.png"), dpi=300, bbox_inches="tight")
scf.pl.dendrogram(adata, color="milestones", legend_loc="on data", color_milestones=True, legend_fontoutline=True)
plt.savefig(os.path.join(output_dir, "scFates_Cellrank_milestones_dendrogram.png"), dpi=300, bbox_inches="tight")

adata.write(os.path.join(data_dir, f"{FILE_NAME}_scFates.h5ad"))
logging.info("scFates inference complete!")