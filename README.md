# Code repository for Sozza et al. "Combinatorial and Inducible CRISPRa/i Achieves Deterministic hiPSC Forward Programming and Iterative Refinement via Single-Cell Genomics"

This GitHub repository includes all the code aiming to reproduce the analysis presented in the publication:
1) scripts_fig5E contains the script of the scRNAseq data shown in Fig 5E. Primary data, input files, and processed data are deposited on Biostudies: Identification nb. XXXXXXX. Input files are also deposited on Zenodo Identification nb. XXXXXXX (Fig5_input_folder)
2) scripts_fig5F contains the script of the bulk RNAseq data shown in Fig 5F. Primary data, input files, and processed data are deposited on Biostudies: Identification nb. XXXXXXX. Input files are also deposited on Zenodo Identification nb. XXXXXXX (Fig5F_input_folder)
3) CIRI Analysis Pipeline is the pipeline complementing the CRISPR screening platform, aiming at identifying the genes and the combination of genes promoting forward programing. Figure 6 ans S6 show the graphs output of this analysis. Primary data, input files, and processed data are deposited on Biostudies: Identification nb. XXXXXXX XXXXXX. Input files are also deposited on Zenodo Identification nb. XXXXXXX (Fig6_dual_guides_INPUT_files and FigS6_single_guide_INPUT_files)
4) scripts_fig7 contains the script of the scRNAseq data shown in Fig 7 and S7. Primary data, input files, and processed data are deposited on Biostudies: Identification nb. XXXXXXX. Input files are also deposited on Zenodo Identification nb. XXXXXXX (Fig7_input_folder)

It is suggested to run scrips 1-2-4 in a Docker environment. The Docker images used for these scripts are:


hedgelab/rstudio-hedgelab:iPS2seq_new_cellranger7

hedgelab/rstudio-hedgelab:iPS2seq_new_cellranger9

These Images can run cellranger (either cellranger 7 or cellranger 9) on bash and the R codes on Rstudio server.

```docker run -d -itv /path/to/shared/folder:/scratch --privileged=true -p 8080:8787 --name container_name hedgelab/rstudio-hedgelab:iPS2seq_new_cellranger7```


```docker exec -idt container_name rstudio-server start```


then go on the browser http://localhost:8080/


#CIRI Analysis Pipeline

This repository contains an R-based pipeline for processing and analyzing single-cell RNA-seq data from CRISPRa/i screens. It handles perturbation deconvolution, quality control, dimensionality reduction (Monocle3), and trajectory analysis to assess the impact of guides on cell differentiation.

## Prerequisites

The pipeline runs in a Dockerized environment to ensure reproducibility (image available at docker.io/hedgelab/ciri_analysis).

Base Image: rocker/r-ver:4.4.2 (Ubuntu 24.04 LTS).

Key Libraries: Monocle3, Seurat, Tidyverse, Bioconductor.

### Build Docker Image

```docker build -t ciri_pipeline . ```

## Pipeline Overview
The analysis follows this sequence:

- Deconvolution: Assign CRISPR perturbations to cells based on UMI distributions.

- Filtering: Quality control and protein-coding gene selection.

- Loading: Normalization and initial dimensionality reduction.

- Validation: Check knockdown/activation efficiency.

- Cluster Analysis: Perturbation enrichment in target clusters.

- Trajectory Analysis: Pseudotime dynamics and statistical testing.

- Signatures: Biological scoring.

## Usage

### 1. Perturbation Deconvolution

Assigns guide identities based on "fixed" and "variable" guide capture. This is the first and mandatory step of the analysis.

The front end function were created with Baryon (see Fairflow-BioinformaticsFramework/baryon-lang). 

Command:

```./perturbation_assignment/CIRI_assign.sh  --target_directory "/home/AB11_screening_CRISPR/analysis_11/."  --script_directory "/home/CIRI_bioinfo/perturbation_assignment/."   --guide_strategy 1   --matrix_filename "filtered_feature_bc_matrix.h5"   --threshold_a -1   --threshold_i -1```

### --target_directory 

Required: Yes

Description: The working directory containing your input files (guides.csv and the .h5 matrix). Note: Output files will be saved here.

### --script_directory

Required: Yes

Description: The directory where the CIRI_unified.R script is located.

Example: /scripts/bioinfo

### --guide_strategy

Required: Yes

Description: The analysis logic mode, 1 or 2 guides experiment.

### --matrix_filename

Required: Yes

Description: The specific filename of the 10x Genomics H5 matrix inside the target directory.

Example: filtered_feature_bc_matrix.h5

### --threshold_a

Required: No

Description: Manual UMI threshold for CRISPRa. Set to -1 to auto-calculate.

Example: 23.5 (or -1)

### --threshold_i

Required: No

Description: Manual UMI threshold for CRISPRi. Set to -1 to auto-calculate.

## Methodology:

Fixed Guides: Fixed guides are assessed to establish the perturbation type (CRISPRa vs CRISPRi). The threshold is defined as the first derivative change of the Kernel Density Estimate curve of the UMI distribution, between the noise peak and the signal peak.

Variable Guides (Single): Variable guides are assigned if the top guide has ≥10 UMIs and the ratio between the top and second guide is ≥5.

Variable Guides (Dual): The sum of the top two guides must be ≥4 UMIs, and the ratio between this sum and the third guide must be ≥10. Cells are filtered out if the two variable guides do not target the same gene.

#### The following steps should be performed inside the docker after manually running it with: 

``` docker run -it  -v /home/CIRI_bioinfo/scripts:/scripts -v /home/AB11_screening_CRISPR/analysis_11/:/data docker.io/hedgelab/ciri_analysis bash ```

### 2. Annotation & Filtering

Performs quality control on the gene expression matrix, which should be customized based on the experiment.

Filters Applied:

- Percentage of ribosomal and mitochondrial genes is calculated.

- Cells with <250 total called genes are filtered out.

- Mitochondrial and ribosomal genes are removed, keeping only protein-coding genes.

- Genes with <3 total UMIs are removed.

Command:

``` Rscript /scripts/anno_filter.R /data/ filtered_feature_bc_matrix.h5 ```

### 3. Data Loading & Preprocessing

Initializes the Monocle3 object, performs normalization, and generates the initial UMAP.

Methodology:

Monocle3 clustering is performed. 

Command:

Args: <directory> <annotated_matrix_csv> <clustering_resolution>

``` Rscript /scripts/CIRI_load.R /data/ annotated_matrix.csv 0.5e-4 ```

### 4. Target Validation

Validates perturbation efficiency by comparing target gene expression in perturbed cells vs. non-targeting controls. Produces violin plots and tables.

Command:

``` Rscript /scripts/guide_genes_expr.R /data/ ```

### 5. Cluster Enrichment Analysis

Analyzes the distribution of perturbations across an user-defined target cluster.

Methodology:

- The cluster of interest is selected.

- The percentage of cells in the target cluster is calculated for each guide combination.

- Combinations are filtered for a minimum cell count (e.g., 10 or 20).

- Combinations present in both samples with sufficient cells are visualized on heatmaps and scatterplots.

Command:

Args: <directory> <cluster_ids> <control_name>

``` Rscript /scripts/CIRI_sec.R /data/ 5 "NTCa-NA" ```

### 6. Subclustering & Trajectory

Subsets the target lineage (e.g., muscle cluster), re-clusters at high resolution, and learns pseudotime trajectories.

Methodology:

- The cluster of interest is subclustered.

- Pseudotime is calculated using Monocle3, setting the root at the node with the highest gene expression of an user-specified gene.

- Pseudotime is produced for all samples together and for each sample separately.

Command:

Args: <directory> <cluster_ids> <root_gene> <group_name> <resolution>

``` Rscript /scripts/CIRI_sec_subclusters.R /data/ 5 SOX2 muscle 1e-3 ```

### 7. Pseudotime Statistics (KS Test)

Performs Kolmogorov-Smirnov tests to detect significant shifts in differentiation speed (pseudotime distribution) compared to controls. Can be run on each sample separately.

Methodology:

- generates ECDF (Cumulative Frequency) plots to show if perturbed cells are reaching differentiation milestones faster or slower than non-targeting controls.

- Kolmogorov-Smirnov (KS) test to quantify the shift in developmental pace, calculating a statistic to identify accelerated or stalled trajectories.

- produces volcano plots that map the magnitude of the developmental shift against its statistical significance (adjusted p-value).

- removes low-cell-count combinations (defaulting to a minimum of 8 cells) to ensure that only robust biological signals are reported.

Command:

Args: <dir> <cds_rdata> <pseudotime_csv> <analysis_name> <run_all> <control_grp> <min_cells>

``` Rscript /scripts/CIRI_pseudotime.R /data/ cds_sample_1_ordered.RData pseudotime_muscle_1.csv muscle_1 TRUE "NTCa-NA" 8 ```

### 8. Signature Scoring

Scores cells based on predefined gene sets (e.g., Cell Cycle, Sarcomere Core) defined in signatures.R.

``` Rscript /scripts/signatures.R ```


## Input/Output Structure
Input:

- filtered_feature_bc_matrix.h5 (10x Genomics output)

- guides.csv (Guide library definition)

- genes_of_interest.csv (custom gene list to visualize the expression)

Examples of input files and their structure can be found in the input_files folder of this repo. Filenames which are not provided in arguments are hardcoded and should not be changed. 

Output:

- processed_cds.RData: Monocle3 object.

- annotation_data.csv: Cell-guide assignments.

- pseudotime_*.csv: Pseudotime values per cell.

- *.pdf: UMAPs, Heatmaps, Violin plots, and Volcano plots.

- ...
