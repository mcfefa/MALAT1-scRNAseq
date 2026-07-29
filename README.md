# MALAT1-scRNAseq
A pipeline for pseudotime and differentiation potential analysis of human and mouse single cell RNA-sequencing datasets. 

## Citation

If you use this code in your research, please cite this paper: 

NA Ben-Crentsil, WM Ismail, ME Balasis, H Newman, A Quintana, M Binder, T Kruer, S Neupane, MC Ferrall-Fairbanks, J Fernandez, TL Lasho, CM Finke, ML Ibrahim, KL McGraw, M Wysota, AL Aldrich, CB Ryder, CT Letson, J Traina, AF McLemore, N Droin, A Shastri, S Yun, E Solary, DA Sallman, AA Beg, L Ma, A Gaspar-Maia, MM Patnaik, and Eric Padron. (2024) "RNA shielding of P65 is required to potentiate oncogenic inflammation in TET2 mutated clonal hematopoiesis." _Cancer Discov._ 2024 Dec 2; 14(12):2509-2531. doi: [10.1158/2159-8920.CD-24-0093](https://aacrjournals.org/cancerdiscovery/article/14/12/2509/750138/RNA-Shielding-of-p65-Is-Required-to-Potentiate). PMID: [39189614](https://pubmed.ncbi.nlm.nih.gov/39189614/). PMCID: [PMC11611684](https://pmc.ncbi.nlm.nih.gov/articles/PMC11611684/).

## Overview
This repository contains code and analysis for single-cell RNA sequencing analysis of publicly available human and mouse hematopoeitic datasets to explore the role of TET2 and MALAT1 expression in normal hematopoesis as part of NA Ben-Crentsil et al's manuscript.

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
