# AI-and-Omics-Research-Internship-2025
This repository contains the assignments and final project completed during my AI Omics internship.
Assignment 1 – R Project Setup & Data Cleaning
Objective
Set up a structured R project environment and perform basic data cleaning on a patient dataset.
What Was Done

Created an organized project folder structure (raw_data/, clean_data/, scripts/, results/, plots/)
Inspected and fixed data types for key variables (age, gender, diagnosis, smoker)
Created a binary smoking status variable (smoker_binary: 1 = Yes, 0 = No)
Exported the cleaned dataset as patient_info_clean.csv.

Practice Exercises
Supplementary exercises covering:

Conditional statements (if, if...else) for checking medical thresholds
Automating data type conversions using for loops
Binary encoding of categorical variables (e.g., Yes/No → 1/0)
Verifying dataset changes before and after transformation

Assignment 2 – Differential Gene Expression (DGE) Classification
Objective
Classify genes from DGE analysis results as Upregulated, Downregulated, or Not Significant based on log2FoldChange and adjusted p-value thresholds.
What Was Done

Wrote a classify_gene() function applying standard DGE thresholds (log2FC > 1 or < -1, padj < 0.05)
Processed two datasets (DEGs_Data_1.csv, DEGs_Data_2.csv) in a for-loop
Replaced missing padj values with 1 to handle NAs conservatively
Added a status column to each dataset and saved the results
Printed summary counts of upregulated, downregulated, and non-significant genes

Assignment 4 – Microarray Data Preprocessing
Objective
Full preprocessing workflow on a real microarray dataset from GEO.
Dataset
E-GEOD-36980: Post-mortem Alzheimer's disease brain samples (79 total: 47 Control, 32 Alzheimer's)
What Was Done

Quality control before and after normalization
RMA normalization
Low-intensity probe filtering
Group definition (Control vs Alzheimer's Disease)

Key Results

Probes before filtering: 32,321
Transcripts after filtering: 30,177
Outliers before normalization: 5 samples
Outliers after normalization: 3 samples (2 by distance, 1 by boxplot)
