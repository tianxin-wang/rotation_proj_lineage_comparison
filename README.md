# Cross-Species Lineage Signature Detection

## Problem Statement

1. **Intra-species**: How do we disentangle trajectories of distinct developmental sub-processes within a species?
2. **Inter-species**: Can we learn common developmental trajectories across species?
3. **Comparative metrics**: How can we quantify the similarity between developmental trajectories of different species?

---

## Dataset Description

- **Chicken Dataset** (`Ggal_Kaess_devGABA_annotated`):   
  [Developmental origins and evolution of pallial cell types and structures in birds](https://www.science.org/doi/10.1126/science.adp5182)

- **Salamander Dataset** combinines:
  - Data from [Cell-type profiling in salamanders identifies innovations in vertebrate forebrain evolution](https://www.science.org/doi/10.1126/science.abp9186)
  - Unpublished developmental data from the Tosches lab
---

## Guiding Questions

- Given a single-cell RNA-seq dataset:
  - Can we assess whether it follows a developmental trajectory?
  - What criteria determine the "tree-likeness" of a trajectory?
  
- Can we:
  - Identify genes whose removal improves tree-likeness?
  - Align two branching trajectories effectively?
    - Use **optimal transport**, ignoring tree structure?
    - Use **tree-informed** methods?
  - Develop metrics to quantify trajectory alignment?
  - Remove genes to **enhance alignment scores**?

---

## Relevant Literature

### Biology
- Useful single cell datasets and comparative analyses done to date:
  - https://www.science.org/doi/10.1126/science.abp9186
  - https://www.science.org/doi/10.1126/science.adp5182
  - https://www.nature.com/articles/s41586-022-04510-w
  - https://elifesciences.org/articles/71864
  - https://www.nature.com/articles/s41586-020-2781-z
  - https://www.nature.com/articles/s41586-021-04237-0

- Evolution of telencephalic GABAergic neurons:
  - https://link.springer.com/chapter/10.1007/978-1-4419-0340-2_1
  - https://www.sciencedirect.com/science/article/pii/B978012369497310007X
  - https://www.sciencedirect.com/science/article/pii/S1084952117302410?via%3Dihub
  - https://www.nature.com/articles/s41583-019-0195-4
  - https://www.sciencedirect.com/science/article/pii/S1084952109000895?via%3Dihub
  - https://www.sciencedirect.com/science/article/pii/S0960982219305950
  - https://www.sciencedirect.com/science/article/pii/S0959438818302708

### Machine Learning & Computational Biology
- [Revealing lineage-related signals in single-cell gene expression using random matrix theory](https://www.pnas.org/doi/10.1073/pnas.1913931118)
- [*CONCORD*](https://qinzhu.github.io/Concord_documentation/)
- [Hierarchical Control of State Transitions in Dense Associative Memories](https://arxiv.org/html/2412.11336v1)
- [*SATURN*](https://www.nature.com/articles/s41592-024-02191-z#data-availability): Cross-species single-cell alignment
- [*Genes2Genes*](https://www.nature.com/articles/s41592-024-02378-4): Gene-level alignment of single-cell trajectories

---

## Relevant Tools

- **Trajectory Comparison**:
  - [Genes2Genes](https://github.com)

- **Lineage Reconstruction & Simulation**:
  - [Dynverse framework](https://github.com/dynverse)
    - Key paper: [Spearheading future omics analyses using dyngen, a multi-modal simulator of single cells](https://www.nature.com/articles/s41467-021-24152-2)
