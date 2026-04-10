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

##  Tech Stack
* **Language:** R
* **Bioinformatics:** Bioconductor (`limma`, `affy`, `arrayQualityMetrics`)
* **Data Science:** `tidyverse`, `ggplot2`, `pheatmap`
