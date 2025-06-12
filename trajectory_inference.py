# %%
import os, sys
import logging
import warnings
warnings.filterwarnings("ignore")

DATASETS = {'chicken':{'annotation': "newannotation_high", 'cluster': "seurat_clusters",
                       'scientific_name': 'Ggal', 'root_cell_type': 'PROG', 'start_time': 6.0}, 
            'salamander':{'annotation': "newannotation_final_hires", 'cluster': "seurat_clusters",
                          'scientific_name': 'Pwal', 'root_cell_type': 'LGE', 'start_time': 30.0}}
GENES = ['ZIC2','ZIC1','OTX2','PENK','NR2F2','MAF','TSHZ1','ETV1',
              'ZEB2','ERBB4','PAX6','FOXP2','FOXP1','TAC1','DRD2','DRD1',
              'PPP1R1B','RASD2','BRS3','ZNF503','ISL1','MEIS2','EGR3',
              'GBX2','GBX1','LHX8','LHX6','NKX2-1','DLX5','DLX1','FOXG1']

# ========== 1. Specify dataset and Configure paths ==========
dataset = 'chicken'  # Choose between 'salamander' and 'chicken'
data_dir = "./Data"
output_dir = f"./{dataset}_output"
os.makedirs(output_dir, exist_ok=True)
TIME_KEY = 'Stage'
PROG = DATASETS[dataset]['root_cell_type']
CLUSTER_KEY = DATASETS[dataset]["cluster"]
START_TIME = DATASETS[dataset]['start_time']
FILE_NAME = DATASETS[dataset]['scientific_name']
ANNOTATION_KEY = DATASETS[dataset]["annotation"]

if dataset == 'salamander':
    PATHS = {'LGE trajectory': ['LGE','prec_LGE','imm_STR','STR_MSN1','STR_MSN1_2',
                                'STR_MSN1_FOXP2','STR_MSN1_SOX8','STR_MSN2',
                                'STR_MSN2_PAX6','STR_MSN2_PENK','STR_NODRD'],
            'dLGE trajectory': ['dLGE','prec_dLGE','imm_OB','OB_GC13','OB_GC2-10-11-12',
                                'OB_GC3-7','OB_GC1-6','OB_GC4-9-5-6','OB_GC8',
                                'OB_PGC_FOXP2','OB_PGC_TH'],
            'MGE+CGE trajectory': ['MGE','CGE','prec_CGE','prec_MGE','prec_MGE/POA',
                                'imm_PALL','imm_INs','PALL1','PALL2','IN_CGE',
                                'IN_MGE_P_LAMP5','IN_MGE_P_SST-','IN_MGE_P_SST+',
                                'IN_MGE_SP_SST-','IN_MGE_SP_SST+'],
            'SEP+POA trajectory': ['eSEP','POA','eSEP-POA','prec_MGE/POA','prec_SEP',
                                'imm_SEP','imm_PALL/SEP','SEP_LAT_PAX6','SEP_LAT1-3',
                                'SEP_LAT2','SEP_LAT4','SEP_MED1','SEP1','SEP2'],
            'Other': ['AMY_CEA_LA', 'AMY_MEA']}
elif dataset == 'chicken':
    PATHS = {'LGE trajectory': ['PROG','PREC','STR','STR_MSN1','STR_MSN2','FOXP2_INs']}

# ========== 2. Import necessary packages and some setups ==========
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

np.random.seed(42)
def saveas(filename, output_dir=output_dir, format='png'):
    plt.savefig(os.path.join(output_dir, filename), dpi=300, bbox_inches='tight',
                format=format)
    plt.close()

# %%
# ========== 3. Set up logging ==========
log_file = os.path.join(output_dir, f"{dataset}_traj.log")
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logging.info(f"❗Starting {dataset} trajectory inference script...")

# %%
# # ========== CONCORD ==========
# import concord as ccd
# import scanpy as sc
# import torch
# torch.manual_seed(42)

# adata = sc.read(os.path.join(data_dir, f"{FILE_NAME}_final.h5ad"))
# logging.info("Running Concord...")
# device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
# feature_list = adata.var_names[adata.var["highly_variable"] == True].tolist()
# logging.info("Initialize Concord with an AnnData object:")
# cur_ccd = ccd.Concord(adata=adata, input_feature=feature_list, device=device) 
# logging.info("Encode data, saving the latent embedding in adata.obsm['Concord']:")
# cur_ccd.encode_adata(output_key='Concord')
# logging.info("Running Concord UMAP...")
# ccd.ul.run_umap(adata, source_key='Concord', result_key='Concord_UMAP', n_components=2, n_neighbors=15, min_dist=0.1, metric='euclidean')
# color_by = [ANNOTATION_KEY, CLUSTER_KEY] # Choose which variables you want to visualize
# ccd.pl.plot_embedding(
#     adata, basis='Concord_UMAP', color_by=color_by, figsize=(10, 5), dpi=600, ncols=2, font_size=6, point_size=10, legend_loc='on data',
#     save_path=os.path.join(output_dir,'Concord_UMAP.png')
# )

# ccd.ul.run_umap(adata, source_key='Concord', result_key='Concord_UMAP_3D', n_components=3, n_neighbors=15, min_dist=0.1, metric='euclidean')

# # Plot the 3D UMAP embeddings
# fig = ccd.pl.plot_embedding_3d(
#     adata, basis='Concord_UMAP_3D', color_by=ANNOTATION_KEY, 
#     save_path=os.path.join(output_dir,'Concord_UMAP_3D.html'),
#     point_size=3, opacity=0.8, width=1500, height=1000
# )

# logging.info("Saving final object...")
# adata.write(os.path.join(data_dir, f"{FILE_NAME}_concord.h5ad"), compression='gzip')
# logging.info("Concord embedding complete!")

# %%
# # ========== GeneTrajectory ==========
# logging.info("3D visualization...")
# from gene_trajectory.add_gene_bin_score import add_gene_bin_score
# from gene_trajectory.coarse_grain import select_top_genes, coarse_grain_adata
# from gene_trajectory.extract_gene_trajectory import get_gene_embedding, extract_gene_trajectory
# from gene_trajectory.get_graph_distance import get_graph_distance
# from gene_trajectory.gene_distance_shared import cal_ot_mat
# from gene_trajectory.run_dm import run_dm
# from gene_trajectory.plot.gene_trajectory_plots import plot_gene_trajectory_3d, plot_gene_trajectory_umap
# from gene_trajectory.util.download_file import download_file_if_missing
# hvg_genes = select_top_genes(adata, layer='counts', n_variable_genes=500)
# neurodev_genes = ['ZIC2','ZIC1','OTX2','PENK','NR2F2','MAF','TSHZ1','ETV1',
#          'ZEB2','ERBB4','PAX6','FOXP2','FOXP1','TAC1','DRD2','DRD1',
#          'PPP1R1B','RASD2','BRS3','ZNF503','ISL1','MEIS2','EGR3',
#          'GBX2','GBX1','LHX8','LHX6','NKX2-1','DLX5','DLX1','FOXG1']
# neurodev_genes = [g for g in neurodev_genes if g in adata.var_names]
# all_genes = list(set(neurodev_genes).union(set(hvg_genes)))

# logging.info("Construct the cell-cell kNN graph and calculate cell-cell graph distances...")
# run_dm(adata)
# cell_graph_dist = get_graph_distance(adata, k=30)
# gene_expression_updated, graph_dist_updated = coarse_grain_adata(adata, graph_dist=cell_graph_dist, features=all_genes, n=500)
# gene_dist_mat = cal_ot_mat(gene_expr=gene_expression_updated, 
#                            ot_cost=graph_dist_updated,
#                            show_progress_bar=True)
# logging.info("Generate the gene embedding by employing Diffusion Map...")
# gene_embedding, _ = get_gene_embedding(gene_dist_mat, k = 10)
# logging.info(f'gene_dist_mat: {gene_dist_mat}')
# logging.info(f'gene_embedding: {gene_embedding}')
# pd.DataFrame(gene_dist_mat).to_csv(os.path.join(output_dir, f"{dataset}_gene_dist_mat.csv"))
# pd.DataFrame(gene_embedding).to_csv(os.path.join(output_dir, f"{dataset}_gene_embedding.csv"))

# # gene_dist_mat = pd.read_csv(os.path.join(output_dir, f"{dataset}_gene_dist_mat.csv"))
# # gene_embedding = pd.read_csv(os.path.join(output_dir, f"{dataset}_gene_embedding.csv"))
# gene_trajectory = extract_gene_trajectory(gene_embedding, gene_dist_mat, t_list = [4, 8, 7], gene_names=all_genes, k=10)
# logging.info("3D visualization...")
# plot_gene_trajectory_3d(gene_trajectory, label_genes=neurodev_genes)
# saveas(f"3D_{dataset}_GeneTrajectory.png")
# logging.info("GeneTrajectory inference complete!")

# %%
# # ========== PAGA ==========
# def run_PAGA(ad, descr = None):
#     ad.X = ad.layers['counts'].copy()
#     sc.pp.normalize_total(ad, target_sum=1e4)
#     sc.pp.log1p(ad)
#     sc.pp.scale(ad)
#     sc.pp.neighbors(ad, n_neighbors=5, n_pcs=30, use_rep='X_pca')
#     sc.tl.diffmap(ad, n_comps=30)
#     sc.pp.neighbors(ad, n_neighbors=10 ,use_rep='X_diffmap')
#     sc.tl.paga(ad, groups=ANNOTATION_KEY)
#     suffix = f"_{descr}" if descr is not None else ""
#     if dataset == 'salamander' and descr is None:
#         sc.set_figure_params(figsize=(8,6))
#         sc.pl.paga(ad, color=[ANNOTATION_KEY], fontsize=4, fontweight='bold', fontoutline=2,
#                    edge_width_scale=0.2, threshold=0.1, node_size_scale=0.5)
#     else:
#         sc.pl.paga(ad, color=[ANNOTATION_KEY])
#     saveas(f"PAGA_annot_nodes{suffix}_UMAP.png")
#     sc.tl.draw_graph(ad, init_pos='paga')
#     sc.pl.draw_graph(ad, color=ANNOTATION_KEY, legend_loc='on data', legend_fontsize = 'xx-small')
#     saveas(f"PAGA_annot_sc{suffix}_FA.png")
#     return ad

# adata = sc.read(os.path.join(data_dir, f"{FILE_NAME}_final.h5ad"))
# adata_hvg = adata[:, adata.var['highly_variable']==True].copy()
# adata_hvg = run_PAGA(adata_hvg)
# sc.pl.umap(adata_hvg, edges=True, color = ANNOTATION_KEY, legend_loc='on data', legend_fontsize= 'xx-small')
# saveas("PAGA_annot_sc_UMAP.png")

# root_cell = np.random.choice(adata_hvg.obs_names[(adata_hvg.obs[ANNOTATION_KEY] == PROG)])
# adata_hvg.uns['iroot'] = np.where(adata_hvg.obs_names == root_cell)[0][0]
# sc.tl.dpt(adata_hvg)
# sc.pl.umap(adata_hvg, color='dpt_pseudotime', title="DPT Pseudotime of HVG genes")
# saveas("dpt_pseudotime_hvg_UMAP.png")

# for descr, path in PATHS.items():
#     if descr == "Other":
#         continue
#     sc.pl.draw_graph(adata_hvg, color=ANNOTATION_KEY, legend_loc='on data', 
#                      legend_fontsize = 'xx-small', groups=path, title=f'{descr}')
#     saveas(f"PAGA_{descr}_FA.png")
#     subset = adata_hvg[adata_hvg.obs[ANNOTATION_KEY].isin(path)].copy()
#     _ = run_PAGA(subset, descr)

# adata.obs['dpt_pseudotime'] = adata_hvg.obs['dpt_pseudotime']
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)
# fig, axs = plt.subplots(ncols=4, gridspec_kw={'wspace': 0.05, 'left': 0.12})
# for ipath, (descr, path) in enumerate(PATHS.items()):
#     if descr == "Other":
#         continue
#     sc.pl.paga_path(adata=adata, nodes=path, keys=GENES,
#                     show_node_names=False, ax=axs[ipath],
#                     ytick_fontsize=12, left_margin=0.15, n_avg=50, 
#                     annotations=['dpt_pseudotime'], 
#                     show_yticks=True if ipath == 0 else False, 
#                     show_colorbar=False, color_map='Greys',
#                     groups_key=ANNOTATION_KEY,
#                     color_maps_annotations={'distance': 'viridis'},
#                     title='{}'.format(descr), return_data=True,
#                     use_raw=False, show=False)
# saveas("PAGA_path_gene_changes.png")

# adata_hvg.write(os.path.join(data_dir, f"{FILE_NAME}_hvg_PAGA.h5ad"), compression='gzip')
# logging.info("Algorithm PAGA complete!")

# %%
# output_dir = "./salamander_palantir_test"
# os.makedirs(output_dir, exist_ok=True)

# %%
# ========== CellRank: Combining MOSCOT with Pseudotime ==========
import pickle
import palantir
import scvelo as scv
import seaborn as sns
import cellrank as cr
from cellrank.estimators import GPCCA
from cellrank.kernels import PseudotimeKernel
from cellrank.kernels import CytoTRACEKernel
from cellrank.kernels import RealTimeKernel
from cellrank.kernels import PrecomputedKernel

# %%
# # ========== MOSCOT ==========
# from moscot.problems.time import TemporalProblem

# adata = sc.read(os.path.join(data_dir, f"{FILE_NAME}_final.h5ad"))
# sc.pp.normalize_total(adata)
# sc.pp.log1p(adata)
# sc.pp.neighbors(adata, use_rep='X_pca')

# logging.info("Computing RealTimeKernel...")
# adata.obs[TIME_KEY] = adata.obs[TIME_KEY].astype(float).astype("category")
# tp = TemporalProblem(adata)
# tp = tp.score_genes_for_marginals(gene_set_proliferation="human", gene_set_apoptosis="human")
# sc.pl.umap(adata, color=["proliferation", "apoptosis"])
# saveas("CellRank_birth_death_UMAP.png")
# tp = tp.prepare(time_key=TIME_KEY)
# if dataset == 'salamander':
#     epsilon = 1e-3
# elif dataset == 'chicken':
#     epsilon = 5e-3
# tp = tp.solve(epsilon=epsilon, tau_a=0.95, scale_cost="mean", device="cpu")
# tmk = RealTimeKernel.from_moscot(tp)
# tmk.compute_transition_matrix(self_transitions="all", conn_weight=0.5, threshold="auto")
# tmk.plot_random_walks(max_iter=500, start_ixs={TIME_KEY: START_TIME}, basis="umap", seed=0)
# saveas("CellRank_random_walks_stage_UMAP.png")

# tmk.write_to_adata(key="moscot_T", copy=False)
# adata.write(os.path.join(data_dir, f"{FILE_NAME}_moscot.h5ad"), compression='gzip')

# %%
# # ========== Palantir pseudotime ==========
# adata = sc.read(os.path.join(data_dir, f"{FILE_NAME}_moscot.h5ad"))
# logging.info("Computing PseudotimeKernel...")
# root_cell = np.random.choice(adata.obs_names[(adata.obs[ANNOTATION_KEY] == PROG)])
# adata.uns['iroot'] = np.where(adata.obs_names == root_cell)[0][0]
# dm_res = palantir.utils.run_diffusion_maps(adata, n_components=20)
# ms_data = palantir.utils.determine_multiscale_space(adata)
# palantir.plot.plot_diffusion_components(adata)
# saveas("palantir_diffusion_components.png", output_dir)
# imputed_X = palantir.utils.run_magic_imputation(adata, n_jobs=10)
# pr_res = palantir.core.run_palantir(adata, root_cell, num_waypoints=3000, n_jobs=10)
# palantir.plot.plot_palantir_results(adata)
# saveas("palantir_pseudotime_UMAP.png", output_dir)
# pk = PseudotimeKernel(adata, time_key="palantir_pseudotime")
# pk.compute_transition_matrix(threshold_scheme="soft")
# pk.write_to_adata(key="palantir_T", copy=False)

# %%
# # ========== CytoTRACE ==========
# adata = sc.read(os.path.join(data_dir, f"{FILE_NAME}_CellRank.h5ad"))
# adata.layers["spliced"] = adata.X
# adata.layers["unspliced"] = adata.X
# scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
# logging.info("Computing CytoTRACE Kernel...")
# ctk = CytoTRACEKernel(adata).compute_cytotrace()
# sc.pl.umap(adata, color="ct_score", color_map="gnuplot2")
# saveas("CytoTRACE_score_UMAP.png")
# sc.pl.umap(adata, color=["ct_pseudotime", "nFeature_RNA"], color_map="gnuplot2")
# saveas("CytoTRACE_pseudotime_UMAP.png")
# adata.write(os.path.join(data_dir, f"{FILE_NAME}_CellRank.h5ad"), compression='gzip')

# %%
# # ========== Combine kernels ==========
adata = sc.read(os.path.join(data_dir, f"{FILE_NAME}_CellRank.h5ad"))
# pk = PseudotimeKernel.from_adata(adata, key="palantir_T")
tmk = RealTimeKernel.from_adata(adata, key="moscot_T")
# combined_kernel = 0.5 * tmk + 0.5 * pk
# combined_kernel.write_to_adata(key="combined_kernel_T", copy=False)
# adata.write(os.path.join(data_dir, f"{FILE_NAME}_CellRank.h5ad"), compression='gzip')

# %%
# # ========== Computing initial and terminal states ==========
# adata = sc.read(os.path.join(data_dir, f"{FILE_NAME}_CellRank.h5ad"))
g = GPCCA.read(os.path.join(data_dir, f"{FILE_NAME}_gpcca_model.pkl"))

logging.info("Computing macrostates...")
if dataset == 'salamander':
    g = cr.estimators.GPCCA(combined_kernel)
elif dataset == 'chicken':
    g = cr.estimators.GPCCA(tmk)
# g.compute_schur(n_components=50, method='krylov', verbose=True)
# g.plot_spectrum(real_only=True, 
#                 save=os.path.join(output_dir, "CellRank_real_eigenvals.png"), dpi=300)
# g.fit(cluster_key=ANNOTATION_KEY, n_states=40)
# g.plot_macrostates(which="all", legend_loc="on data", size=80, legend_fontsize=4,
#                    save=os.path.join(output_dir, "CellRank_all_states_1.png"), dpi=300)
# g.plot_macrostates(which="all", legend_loc="right", size=80, legend_fontsize=4,
#                    save=os.path.join(output_dir, "CellRank_all_states_2.png"), dpi=300)
# g.plot_macrostate_composition(key=ANNOTATION_KEY, figsize=(10, 4), dpi=300,
#                               save=os.path.join(output_dir, "CellRank_macrostate_composition.png"))
# g.plot_coarse_T(annotate=False)
# saveas("CellRank_coarse_transitions.png", output_dir)

# logging.info("Computing fate probabilities and driver genes...")
# states_to_exclude = ['imm', 'prec', 'CGE', 'MGE', 'LGE', 'PROG', 'PREC']
# terminal_states = [x for x in adata.obs['macrostates_fwd'].unique().dropna().to_list() 
#                    if not any(x.startswith(prefix) for prefix in states_to_exclude)]
# g.set_terminal_states(states=terminal_states)
# g.plot_macrostates(which="terminal", s=80, legend_loc="on data", legend_fontoutline=2, 
#                    figsize=(8,6), legend_fontsize=4, legend_fontweight="normal")
# saveas("CellRank_terminal_states.png", output_dir)
g.set_initial_states(states=["CGE", "LGE"])
g.plot_macrostates(which="initial", legend_loc="on data", s=80)
saveas("CellRank_initial_states.png", output_dir)
g.compute_fate_probabilities()
g.plot_fate_probabilities(legend_loc="on data", legend_fontoutline=2, figsize=(8,6),
                          legend_fontsize=4, legend_fontweight="normal")
saveas("CellRank_fate_probabilities.png", output_dir)

g.write(os.path.join(data_dir, f"{FILE_NAME}_gpcca_model.pkl"))
adata.write(os.path.join(data_dir, f"{FILE_NAME}_CellRank.h5ad"), compression='gzip')

# %%
# g.write(os.path.join(output_dir, f"{FILE_NAME}_gpcca_model.pkl"))
# adata.write(os.path.join(output_dir, f"{FILE_NAME}_CellRank.h5ad"), compression='gzip')

# %%
# # ========== Visualizing and clustering gene expression trends ==========
# def find_cluster_list(terminal_state, PATHS):
#     if terminal_state.startswith("AMY"):
#         return None
#     for _, cell_types in PATHS.items():
#         for ct in cell_types:
#             if terminal_state.startswith(ct):
#                 return cell_types
#     return None

# model = cr.models.GAM(adata)
# gene_sets = {}
# for terminal_state in g.terminal_states.unique().dropna():
#     clusters = find_cluster_list(terminal_state, PATHS)
#     lineage_drivers = g.compute_lineage_drivers(lineages=terminal_state, n_jobs=16, 
#                                                 cluster_key=ANNOTATION_KEY, clusters=clusters)
#     driver_genes = lineage_drivers.head(30).index.to_list()
#     gene_sets[terminal_state] = driver_genes
#     cr.pl.heatmap(
#         adata,
#         model=model,
#         lineages=terminal_state,
#         cluster_key=ANNOTATION_KEY,
#         show_fate_probabilities=True,
#         data_key= "MAGIC_imputed_data",
#         genes=driver_genes,
#         time_key="palantir_pseudotime",
#         figsize=(12, 10),
#         show_all_genes=True,
#         weight_threshold=(1e-3, 1e-3),
#     )
#     saveas(f"CellRank_{terminal_state}_gene_trend.png", output_dir)
# pd.DataFrame(gene_sets).to_excel(os.path.join(output_dir, "CellRank_fates_top30_dexp.xlsx"), index=False)

# %%
# ========== Mapping progenitor and precursor states to terminal states ==========
# adata = sc.read(os.path.join(data_dir, f"{FILE_NAME}_CellRank.h5ad"))
# terminal_states = adata.obs['term_states_fwd'].unique().dropna().to_list()
# prec_clusters = [x for x in adata.obs[ANNOTATION_KEY].unique() if "prec" in str(x).lower()]
# sc.pl.umap(adata, color=ANNOTATION_KEY, groups=prec_clusters, legend_loc="right", 
#            legend_fontsize = 'x-small', frameon=True, size=10)
# saveas("CellRank_prec_clusters_UMAP.png", output_dir)
# cr.pl.aggregate_fate_probabilities(adata, mode="bar", lineages=terminal_states,
#                                    cluster_key=ANNOTATION_KEY, clusters=prec_clusters, ncols=1)
# fig = plt.gcf()
# for ax in fig.axes:
#     ax.tick_params(axis='x', labelsize=6)
#     ax.tick_params(axis='y', labelsize=6)
# saveas("CellRank_prec_fate_probs.png", output_dir)

# prog_clusters = adata[adata.obs['origin']=='progenitors'].obs['SCT_snn_res.2'].value_counts() > 30
# prog_clusters = prog_clusters[prog_clusters].index.to_list()
# prog_clusters = adata.obs['originalclusters_prog'].unique().tolist()
# prog_clusters.remove("NaN")
# sc.pl.umap(adata, color='originalclusters_prog', groups=prog_clusters, legend_fontsize='xx-small', frameon=True, size=10)
# saveas("CellRank_prog_clusters_UMAP.png", output_dir)
# cr.pl.aggregate_fate_probabilities(adata, mode="bar", lineages=terminal_states, ncols=1, 
#                                    cluster_key='originalclusters_prog', clusters=prog_clusters)
# fig = plt.gcf()
# for ax in fig.axes:
#     ax.tick_params(axis='x', labelsize=6)
#     ax.tick_params(axis='y', labelsize=6)
# saveas("CellRank_prog_fate_probs.png", output_dir)

# %%
# # ========== Generating a crosstable on two annotation columns ==========
# import scipy
# from scipy.optimize import linear_sum_assignment

# ct = pd.crosstab(adata.obs['macrostates_fwd'], adata.obs[ANNOTATION_KEY])
# pct = ct.div(ct.sum(axis=1), axis=0) * 100

# cost_matrix = -pct.values
# row_ind, col_ind = linear_sum_assignment(cost_matrix)
# pct_reordered = pct.iloc[row_ind, col_ind]
# pct_reordered.index = pct.index[row_ind]
# pct_reordered.columns = pct.columns[col_ind]

# plt.figure(figsize=(10, 8))
# sns.heatmap(pct_reordered, cmap="viridis")
# plt.xticks(fontsize=6)
# plt.yticks(fontsize=6)
# plt.xlabel(ANNOTATION_KEY)
# plt.ylabel('macrostates_fwd')
# plt.tight_layout()
# saveas("Cellrank_crosstable.png")
 
# %%
# # ========== scFates ==========
# os.environ['R_HOME'] = sys.exec_prefix+"/lib/R/"
# import seaborn
# import scFates as scf
# from anndata import AnnData
# seaborn.reset_orig()
# sc.set_figure_params()
# scf.set_figure_pubready()
# sc.settings.verbosity = 3
# sc.settings.logfile = sys.stdout

# adata = sc.read(os.path.join(data_dir, f"{FILE_NAME}_CellRank.h5ad"))

# sc.pl.draw_graph(adata, color=ANNOTATION_KEY, legend_loc='on data', legend_fontsize = 'xx-small')
# saveas("scFates_Cellrank_annot_sc_FA.png")
# scf.tl.cellrank_to_tree(adata, method="ppt", Nodes=400, time="palantir_pseudotime",
#                         device="gpu", seed=1, ppt_lambda=1, ppt_sigma=0.01, ppt_nsteps=200)
# sc.set_figure_params(figsize=(10,10))
# scf.pl.graph(adata, alpha_nodes=0.6, size_nodes=4, linewidth=1)
# saveas("scFates_Cellrank_tree_rep_UMAP.png")

# scf.tl.root(adata, 321)
# scf.tl.pseudotime(adata, n_jobs=16, n_map=200, seed=42)
# scf.pl.trajectory(adata, scale_path=0.6)
# saveas("scFates_Cellrank_trajectory_UMAP.png")
# sc.pl.umap(adata, color=[ANNOTATION_KEY, "seg"])
# saveas("scFates_Cellrank_annot_seg_UMAP.png")

# scf.tl.dendrogram(adata)
# sc.set_figure_params(figsize=(10,6))
# scf.pl.dendrogram(adata, color=ANNOTATION_KEY, legend_loc='on data', color_milestones=False, 
#                   legend_fontsize=4, legend_fontweight='normal')
# saveas("scFates_Cellrank_annot_dendrogram.svg", format='svg')
# scf.pl.dendrogram(adata, color='macrostates_fwd', legend_loc='on data', color_milestones=False, 
#                   legend_fontsize=4, legend_fontweight='normal')
# saveas("scFates_Cellrank_macrostates_dendrogram.svg", format='svg')
# scf.pl.dendrogram(adata, color="milestones", legend_loc="on data", color_milestones=False, 
#                   legend_fontsize=4, legend_fontweight='normal')
# saveas("scFates_Cellrank_milestones_dendrogram.png")

# scf.tl.test_association(adata, n_jobs=32)
# scf.pl.test_association(adata)
# saveas("genes_scFates_tree_associations.png")
# adata.write(os.path.join(data_dir, f"{FILE_NAME}_scFates.h5ad"), compression='gzip')
# adata = sc.read(os.path.join(data_dir, f"{FILE_NAME}_scFates.h5ad"))
# scf.tl.fit(adata, n_jobs=32, features = adata.var_names[adata.var["signi"] == True].tolist())
# adata.write(os.path.join(data_dir, f"{FILE_NAME}_scFates.h5ad"), compression='gzip')

# logging.info("scFates inference complete!")