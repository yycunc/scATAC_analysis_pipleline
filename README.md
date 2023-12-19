# scATAC analysis pipleline for zebrafish heart regeneration
## Overview
Cardiac regeneration requires coordinated participation of multiple cell types whereby their communications result in transient activation of pro-regenerative cell states. Although the molecular characteristics and lineage origins of these activated cell states and their contribution to cardiac regeneration have been studied, the extracellular signaling and the intrinsic genetic program underlying the activation of the transient functional cell states remain largely unexplored. In this study, we delineated the chromatin landscapes of the non-cardiomyocytes (nonCMs) of the regenerating heart at single cell level, and inferred the cis-regulatory architectures and trans-acting factors that control cell-type specific gene expression programs. Moreover, further motif analysis and cell-specific genetic manipulations suggest that the macrophage-derived inflammatory signal TNFα, acting via its downstream transcription factor complex AP-1, functions cooperatively with discrete transcription regulators to activate respective nonCM cell types critical for cardiac regeneration. Thus, our study defines the regulatory architectures and intercellular communication principles in zebrafish heart regeneration.

In this repository, we provide the code for our analysis pipeline of scATAC-seq data.

## Codes for analysis
1. scATACseq_analysis_pipeline_for_major_celltypes.R:
   The codes listed in this file were used for data pre-processing, clustering, peak calling and differential accessibility analysis for the scATAC-seq dataset

2. scATACseq_analysis_pipeline_for_subclusters.R:
   The codes listed in this file were used for data pre-processing, clustering, peak calling and differential accessibility analysis for the subclusters of fibroblasts, endothelial cells and macrophages

