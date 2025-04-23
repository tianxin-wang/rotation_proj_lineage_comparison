# Cross-Species Lineage Signature Detection

## Problem Statement

1. **Intra-species**: How do we disentangle trajectories of distinct developmental sub-processes within a species?
2. **Inter-species**: Can we learn common developmental trajectories across species?
3. **Comparative metrics**: How can we quantify the similarity between developmental trajectories of different species?

---

## Dataset Description

- **Chicken Dataset** (`Ggal_Kaess_devGABA_annotated`):  
  Contains **X cells** and **G genes**.  
  [Source](https://www.science.org/doi/10.1126/science.adp5182)

- **Salamander Dataset**:  
  Contains **X cells** and **G genes**, combining:
  - Data from [Science (2022)](https://www.science.org/doi/10.1126/science.abp9186)
  - Unpublished developmental data from the lab

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
- Core datasets & comparative analyses:
  - [Science (2022)](https://www.science.org/doi/10.1126/science.abp9186)
  - [Science (2023)](https://www.science.org/doi/10.1126/science.adp5182)
  - [Nature](https://www.nature.com/articles/s41586-022-04510-w)
  - [eLife](https://elifesciences.org/articles/71864)
  - [Nature (2020)](https://www.nature.com/articles
