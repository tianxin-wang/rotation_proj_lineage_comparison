Cross-species lineage signature detection


Problem Statement:
1. Given a species, how do we separate the trajectories of different developmental sub-processes?
2. Can we learn common developmental trajectories between species?
3. How do we quantify the similarities in developmental trajectories between species?


Dataset description:
* The Ggal_Kaess_devGABA_annotated dataset contains X cells and G genes and was obtained from https://www.science.org/doi/10.1126/science.adp5182 
* The salamander dataset contains X cells and G genes and is composed of developmental data from https://www.science.org/doi/10.1126/science.abp9186 and unpublished data from the lab




Guiding questions:
* Given a single cell RNA-seq data set, can we tell if the data is generated from a developmental trajectory? What makes a tree better than another tree? Assume we have a score for how “tree”-like a data set is!
* Can we find a group of genes to remove that would improve how likely a dataset is to be a tree?
* Given two trees/ branching trajectories, what is the best way to align them? Should we use optimal transport and pretend we do not know about the underlying tree? Should we use the knowledge regarding the underlying tree? Assume that given two trees, we have a score for how aligned they are?
* Can we remove genes that improve the score above?




Relevant Literature:
* Biology:
   * Useful single cell datasets and comparative analyses done to date: 
   * https://www.science.org/doi/10.1126/science.abp9186
   * https://www.science.org/doi/10.1126/science.adp5182 
   * https://www.nature.com/articles/s41586-022-04510-w
   * https://elifesciences.org/articles/71864 
   * https://www.nature.com/articles/s41586-020-2781-z
   * https://www.nature.com/articles/s41586-021-04237-0
   * Papers on evolution of telencephalic GABAergic neurons:
   * https://link.springer.com/chapter/10.1007/978-1-4419-0340-2_1
   * https://www.sciencedirect.com/science/article/pii/B978012369497310007X
   * https://www.sciencedirect.com/science/article/pii/S1084952117302410?via%3Dihub
   * https://www.nature.com/articles/s41583-019-0195-4
   * https://www.sciencedirect.com/science/article/pii/S1084952109000895?via%3Dihub
   * https://www.sciencedirect.com/science/article/pii/S0960982219305950
   * https://www.sciencedirect.com/science/article/pii/S0959438818302708


* Machine Learning/Comp Bio:
   * Revealing lineage-related signals in single-cell gene expression using random matrix theory
   * Revealing a coherent cell state landscape across single cell datasets with CONCORD | bioRxiv
   * Hierarchical Control of State Transitions in Dense Associative Memories
   * Toward universal cell embeddings: integrating single-cell RNA-seq datasets across species with SATURN | Nature Methods
   * Gene-level alignment of single-cell trajectories | Nature Methods




Relevant tools:
* Trajectory comparisons:
   * Genes2Genes
* Lineage reconstruction and simulations:
   * https://github.com/dynverse ; particularly https://www.nature.com/articles/s41467-021-24152-2