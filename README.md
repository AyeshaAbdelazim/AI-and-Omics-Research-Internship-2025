# AI-and-Omics-Research-Internship-2025

This repository contains the assignments and final project completed during my AI Omics internship. 

---


##  Assignments Overview

###  Assignment 1: R Project Setup & Data Cleaning
**Objective:** Establish a structured R environment and perform data cleaning on patient metadata.
* **Key Tasks:**
    * Fixed data types for variables: `age`, `gender`, `diagnosis`, and `smoker`.
    * Created a binary variable `smoker_binary` (1 = Yes, 0 = No).
    * Exported `patient_info_clean.csv` for downstream analysis.
* **Skills Covered:** Conditional statements for medical thresholds and automating data type conversions with `for` loops.

###  Assignment 2: DGE Classification
**Objective:** Automate the classification of genes based on expression significance.
* **Key Tasks:**
    * Developed a `classify_gene()` function using thresholds: $log_2FC > 1$ or $< -1$ and $padj < 0.05$.
    * Handled missing values by conservatively replacing `NA` $padj$ values with 1.
    * Processed multiple datasets using iterative loops and generated summary counts.

###  Assignment 4: Microarray Data Preprocessing
**Objective:** Preprocess a real-world GEO dataset (`E-GEOD-36980`) involving Alzheimer's Disease brain samples.
* **Workflow:**
    * Performed **RMA Normalization** and quality control.
    * Filtered low-intensity probes to refine the dataset.
* **Key Results:**
    * **Probes before filtering:** 32,321
    * **Transcripts after filtering:** 30,177
    * **Outlier Reduction:** Successfully reduced outliers from 5 samples to 3 post-normalization.

###  Assignment 5: DGE Analysis (Alzheimer's Disease)
**Objective:** Detailed differential expression analysis of the `E-GEOD-36980` dataset using Bioconductor.
* **Analysis Workflow:**
    * **Quality Control:** Used `arrayQualityMetrics` for assessment.
    * **Statistical Analysis:** Applied `limma` with empirical Bayes moderation for robust results.
    * **Visualization:** Generated Volcano plots (`ggplot2`) and Heatmaps (`pheatmap`) to identify significant gene clusters.

---
### My part of the final project: Co-expression & Network Analysis
**Objective:** Integrate ncRNA and mRNA expression data to identify regulatory hubs and disease-specific gene modules.

#### 1. ncRNA-mRNA Co-expression Profiling
* **Targeted Extraction:** Isolated expression data for 102 high-confidence genes identified from literature (Table S7).
* **Correlation Pipeline:** * Computed **Spearman correlations** between 706 ncRNAs and the mRNA matrix.
    * Applied **Benjamini-Hochberg (BH)** FDR correction to control for multiple testing.
    * **Results:** Identified significant ncRNA-mRNA regulatory pairs ($|r| > 0.5$, $FDR < 0.10$).

#### 2. Network Hub Identification
* **Hub Genes:** Focused on 20 critical hub genes (e.g., *CD69, CCR7, CXCL10*) from supplementary datasets.
* **Visualization:** * Generated **igraph** network plots to visualize the connectivity between ncRNAs and immune-related hub genes.
    * Produced annotated heatmaps distinguishing between Upregulated and Downregulated ncRNA clusters.

#### 3. WGCNA (Weighted Gene Co-expression Network Analysis)
* **Module Construction:** Performed hierarchical clustering and soft-thresholding (Power = 4) to group genes into functional modules.
* **Brain Metastasis Specificity:** * Conducted **Wilcoxon Rank Sum tests** on module eigengenes.
    * Identified modules significantly enriched in **Brain Metastasis (BM)** vs. Primary lung adenocarcinoma ($FDR < 0.05, |logFC| > 0.5$).
* **Immune Correlation:** Linked genomic modules to immune cell infiltration traits using adaptive correlation thresholds.

#### 4. Key Questions Answered
* **Hub Discovery:** Identified the top 20 network hubs based on intramodular connectivity ($kWithin$) and Module Membership ($kME > 0.7$).
* **Network Segregation:** Proved via **Chi-square testing** that Upregulated and Downregulated ncRNAs cluster into distinct functional modules.
* **Metastatic Drivers:** Isolated specific modules (e.g., MEblue, MEbrown) that act as primary drivers for the transition from primary tumor to brain metastasis.

---

