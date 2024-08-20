# MALAT1-scRNAseq
A pipeline for pseudotime and differentiation potential analysis of human and mouse single cell RNA-sequencing datasets. 

## Overview
This repository contains code and analysis for single-cell RNA sequencing analysis of publicly available human and mouse hematopoeitic datasets to explore the role of TET2 and MALAT1 expression in normal hematopoesis as part of: 

Nana Adjoa Ben-Crentsil, Wazim Mohammed Ismail, Maria Balasis, Hannah Newman, Ariel Quintana, Moritz Binder, Traci Kruer, Surendra Neupane, Meghan Ferrall-Fairbanks, Jenna Fernandez, Terra Lasho, Christy Finke, Mohammed Ibrahim, Kathy Mcgraw, Michael Wysota, Amy Aldrich, Christopher Ryder, Christopher Letson, Joshua Traina, Amy McLemore, Nathalie Droin, Aditi Shastri, Seongseok Yun, Eric Solary, David Sallman, Amer Beg, Li Ma, Alexandre Gaspar-Maia, Mrinal Patnaik, and Eric Padron. "RNA shielding of P65 is required to potentiate oncogenic inflammation in TET2 mutated clonal hematopoiesis." (In Press)

## Requirements
- Python (v3.6)
- scanpy v1.4.4
- anndata v0.6.22.post1
- umap v0.3.7
- numpy v1.16.2
- scipy v1.3.1
- pandas v0.23.4
- scikit-learn v0.20.3
- statsmodels v0.10.1
- palantir v1.0.0

## Contents 
- [src/Dahlin-Tutorial+Palantir.ipynb](src/Dahlin-Tutorial+Palantir.ipynb): Implements the previous published Dahlin, Hamey et al (eBlood 2018) single-cell RNA sequencing analysis pipeline with a mouse hematopoietic dataset and then integrates the Palantir toolkid for pseudotime differentiation potential analysis in a Jupyter Notebook.
- [src/Dahlin-Tutorial+Palantir.py](src/Dahlin-Tutorial+Palantir.py): Implements the previous published Dahlin single-cell RNA sequencing analysis pipeline with a mouse hematopoietic dataset and then integrates the Palantir toolkid for pseudotime differentiation potential analysis in a Python Script.
- [src/Palantir-Tutorial-Adapted.ipynb](src/Palantir-Tutorial-Adapted.ipynb): Implements and adapts the Setty et al (Nat Biotech 2019) pipeline to explore single-cell expression of specific gene markers of hematopoiesis in a Jupyter Notebook. 
- [src/Palantir-Tutorial-Adapted.py](src/Palantir-Tutorial-Adapted.py): Implements and adapts the Setty et al (Nat Biotech 2019) pipeline to explore single-cell expression of specific gene markers of hematopoiesis in a Python Script. 
