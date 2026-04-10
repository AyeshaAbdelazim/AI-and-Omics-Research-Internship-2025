#How did I extract the mRNA_expr

library(dplyr)

# Define the 102 Paper Gene

# The 102 genes from the paper, they where mentioned in supplementary tables.doc.xls table S7
paper_genes <- c(
  "ABI3BP", "ADAMDEC1", "AEBP1", "ANGPTL2", "AOC3", "APOC1", "C1S", "C2CD4B",
  "C4BPA", "C7", "CCL14", "CCL18", "CCL19", "CCL5", "CCR7", "CD2", "CD247",
  "CD27", "CD3D", "CD48", "CD6", "CD69", "CD79A", "CD79B", "CD96", "CEBPA",
  "CES1", "CLDN11", "CLEC3B", "COL10A1", "COL1A2", "COL6A3", "CPA3", "CPXM2",
  "CRISPLD2", "CRYAB", "CTSK", "CXCL10", "CXCL13", "CXCL9", "DMBT1", "DPT",
  "DRAM1", "EEF1A2", "EGFL6", "FAM107A", "FCER1A", "FCN1", "FGL2", "FMO2",
  "FMO3", "FOSB", "FPR3", "GAS1", "GBP5", "GPM6B", "GZMA", "GZMK", "HLA-DPB1",
  "HLA-DQA1", "HLA-DQB2", "IL2RB", "INMT", "IRF8", "ITGBL1", "KIF5C", "LAX1",
  "LTB", "LUM", "LYZ", "MAP1B", "MFAP2", "MFAP4", "MMP11", "MXRA5", "MYH11",
  "OLFML1", "PDGFRA", "PDLIM3", "PGC", "PIGR", "PLP1", "PODN", "POSTN", "PRRX1",
  "RARRES2", "RGS1", "RNU4-1", "SCG2", "SCTR", "SFRP2", "SH2D1A", "SLAMF6",
  "SUSD2", "SVEP1", "TF", "TMEM119", "TPSAB1", "TPSB2", "UBD", "VCAM1", "ZIC2"
)

cat("Paper genes list loaded:", length(paper_genes), "genes\n\n")


# Check Prerequisites

if (!exists("master_norm_matrix")) {
  stop("ERROR: master_norm_matrix not found. Please run Parts 1-2 first.")
}

cat("✓ master_norm_matrix found\n")
cat("  Dimensions:", nrow(master_norm_matrix), "genes x", 
    ncol(master_norm_matrix), "samples\n\n")


# Extract Expression Data


# Find which paper genes are in our expression matrix
genes_available <- intersect(paper_genes, rownames(master_norm_matrix))
genes_missing <- setdiff(paper_genes, rownames(master_norm_matrix))

cat("Genes found in dataset:", length(genes_available), "out of 102\n")
cat("Genes missing:", length(genes_missing), "\n")

if (length(genes_missing) > 0) {
  cat("\nMissing genes:\n")
  print(genes_missing)
}

# Extract the expression matrix for available genes
mRNA_expr <- master_norm_matrix[genes_available, , drop = FALSE]

cat("\n✓ Extracted mRNA_expr:\n")
cat("  Dimensions:", nrow(mRNA_expr), "genes x", ncol(mRNA_expr), "samples\n")
cat("  Expression range:", round(range(mRNA_expr, na.rm=TRUE), 3), "\n\n")

# Quality Control

# Check for missing values
n_na <- sum(is.na(mRNA_expr))
cat("Missing values:", n_na, "\n")

if (n_na > 0) {
  cat("  Proportion:", round(n_na / (nrow(mRNA_expr) * ncol(mRNA_expr)) * 100, 2), "%\n")
}

# Check variance
gene_vars <- apply(mRNA_expr, 1, var, na.rm = TRUE)
zero_var <- sum(gene_vars == 0)

if (zero_var > 0) {
  cat("WARNING:", zero_var, "genes have zero variance\n")
  cat("  → Removing these genes...\n")
  mRNA_expr <- mRNA_expr[gene_vars > 0, , drop = FALSE]
  cat("  → Retained:", nrow(mRNA_expr), "genes\n")
}

cat("✓ QC complete\n\n")

####ncRNA-mRNA Co-expression starts from here####

#  Required packages 
library(Hmisc); library(dplyr); library(pheatmap); library(RColorBrewer); library(ggplot2); library(igraph)

#Inputs you must have 
# combined_expr_norm: expression matrix genes x samples (rownames = gene symbols)
# sig_results: data.frame with $hgnc_symbol (ncRNA list) and $logFC
# mRNA_expr: expression matrix for the 102 DEGs (rows = gene symbols, cols = samples)

# 1)  subsetting + numeric coercion

ncRNA_ids <- intersect(sig_results$hgnc_symbol, rownames(combined_expr_norm))
if(length(ncRNA_ids) == 0) stop("No ncRNA IDs found in combined_expr_norm. Check names.")

ncRNA_expr <- combined_expr_norm[ncRNA_ids, , drop = FALSE]
mRNA_expr <- as.matrix(mRNA_expr) # ensure matrix

# Ensure common samples and same order
common_samples <- intersect(colnames(ncRNA_expr), colnames(mRNA_expr))
if(length(common_samples) < 3) stop("Less than 3 common samples; cannot compute correlations reliably.")
ncRNA_expr <- as.matrix(ncRNA_expr[, common_samples, drop = FALSE])
mRNA_expr <- as.matrix(mRNA_expr[, common_samples, drop = FALSE])

# Coerce to numeric and remove genes with zero variance 
row_var <- function(x) { apply(x, 1, function(z) var(as.numeric(z), na.rm = TRUE)) }

v_nc <- row_var(ncRNA_expr)
v_m <- row_var(mRNA_expr)

ncRNA_expr <- ncRNA_expr[v_nc > 0, , drop = FALSE]
mRNA_expr <- mRNA_expr[v_m > 0, , drop = FALSE]

cat("ncRNAs:", nrow(ncRNA_expr), "mRNAs:", nrow(mRNA_expr), "Samples:", length(common_samples), "\n")


# 2) Compute pairwise Spearman correlations for ncRNA x mRNA only

cor_matrix <- cor(t(ncRNA_expr), t(mRNA_expr), method = "spearman", use = "pairwise.complete.obs")

# compute p-values per pair (vectorized)
n <- ncol(ncRNA_expr) # number of samples (same for both)
t_stat <- cor_matrix * sqrt((n - 2) / (1 - cor_matrix^2))
pval_matrix <- 2 * pt(-abs(t_stat), df = n - 2)

# Clean possible NaN/Inf from divisions
pval_matrix[!is.finite(pval_matrix)] <- NA
cor_matrix[!is.finite(cor_matrix)] <- NA


# 3) FDR correction (using column-wise BH correction)

fdr_matrix <- apply(pval_matrix, 2, function(x) p.adjust(x, method = "BH"))


# 4) Find significant pairs (General ncRNA-mRNA) 

cor_threshold <- 0.5
fdr_threshold <- 0.10

sig_idx <- which(abs(cor_matrix) > cor_threshold & fdr_matrix < fdr_threshold, arr.ind = TRUE)

cat("Significant correlations found:", nrow(sig_idx), "\n")

if(nrow(sig_idx) > 0) {
  cor_results <- data.frame(
    ncRNA = rownames(cor_matrix)[sig_idx[,1]],
    mRNA = colnames(cor_matrix)[sig_idx[,2]],
    Correlation = cor_matrix[sig_idx],
    Pvalue = pval_matrix[sig_idx],
    FDR = fdr_matrix[sig_idx],
    Regulation = ifelse(cor_matrix[sig_idx] > 0, "Positive", "Negative"),
    stringsAsFactors = FALSE
  ) %>% arrange(desc(abs(Correlation)))
  
  write.csv(cor_results, "ncRNA_mRNA_coexpression_results.csv", row.names = FALSE)
  print(head(cor_results, 20))
} else {
  cat("No significant correlations under thresholds |r| >", cor_threshold, "and FDR <", fdr_threshold, "\n")
}


# 5) Heatmap of top pairs (General ncRNA-mRNA)

if(exists("cor_results") && nrow(cor_results) > 0) {
  top_pairs <- head(cor_results, 50)
  unique_ncRNAs <- unique(top_pairs$ncRNA)
  unique_mRNAs <- unique(top_pairs$mRNA)
  heatmap_data <- matrix(NA, nrow = length(unique_ncRNAs), ncol = length(unique_mRNAs),
                         dimnames = list(unique_ncRNAs, unique_mRNAs))
  for(i in seq_len(nrow(top_pairs))) {
    heatmap_data[top_pairs$ncRNA[i], top_pairs$mRNA[i]] <- top_pairs$Correlation[i]
  }
  heatmap_data[is.na(heatmap_data)] <- 0
  if(nrow(heatmap_data) > 1 && ncol(heatmap_data) > 1) {
    pheatmap(heatmap_data,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             main = "Top 50 ncRNA-mRNA correlations",
             filename = "heatmap_top_correlations.pdf")
  }
}


# HUB GENE ANALYSIS SECTION 

up_ncRNAs <- sig_results$hgnc_symbol[sig_results$logFC > 0]
down_ncRNAs <- sig_results$hgnc_symbol[sig_results$logFC < 0]

# 1. Prepare hub gene matrices, hub genes were mentioned in supplementary tables.doc.xls Table S13  
hub_genes <- c("CD69", "CCR7", "CD27", "CD2", "CCL5", "CD247",
               "GZMA", "CD3D", "GZMK", "IL2RB", "CXCL9", "CCL19",
               "CXCL13", "CXCL10", "CD48", "VCAM1", "CD79B", "SLAMF6",
               "CD79A", "SH2D1A")

# Subset correlation and FDR matrices to hub genes
hub_cor_matrix <- cor_matrix[, hub_genes, drop = FALSE]
hub_fdr_matrix <- fdr_matrix[, hub_genes, drop = FALSE]

# 2. Align rows/columns
common_ncRNAs <- intersect(rownames(hub_cor_matrix), rownames(hub_fdr_matrix))
common_hubs <- intersect(colnames(hub_cor_matrix), colnames(hub_fdr_matrix))

hub_cor_matrix <- hub_cor_matrix[common_ncRNAs, common_hubs, drop = FALSE]
hub_fdr_matrix <- hub_fdr_matrix[common_ncRNAs, common_hubs, drop = FALSE]

# 3. Threshold correlations 
sig_idx <- which(abs(hub_cor_matrix) >= cor_threshold & hub_fdr_matrix <= fdr_threshold, arr.ind = TRUE)

# 4. Create results table
if(nrow(sig_idx) > 0) {
  hub_results <- data.frame(
    ncRNA = rownames(hub_cor_matrix)[sig_idx[, 1]],
    hub_gene = colnames(hub_cor_matrix)[sig_idx[, 2]],
    Correlation = hub_cor_matrix[sig_idx],
    FDR = hub_fdr_matrix[sig_idx],
    Regulation = ifelse(hub_cor_matrix[sig_idx] > 0, "Positive", "Negative")
  ) %>% arrange(desc(abs(Correlation)))
} else {
  stop("No significant ncRNA–hub gene correlations found with current thresholds.")
}

# 5. Rank ncRNAs by number of hub genes they regulate
ncRNA_rank <- hub_results %>%
  group_by(ncRNA) %>%
  summarise(n_hub_genes = n(),
            avg_abs_cor = mean(abs(Correlation))) %>%
  arrange(desc(n_hub_genes), desc(avg_abs_cor))

print(head(ncRNA_rank, 10))

# 6. Prepare heatmap matrix (ncRNA × hub gene)
heatmap_matrix <- matrix(0, nrow = length(unique(hub_results$ncRNA)),
                         ncol = length(unique(hub_results$hub_gene)),
                         dimnames = list(unique(hub_results$ncRNA), unique(hub_results$hub_gene)))

for(i in 1:nrow(hub_results)) {
  heatmap_matrix[hub_results$ncRNA[i], hub_results$hub_gene[i]] <- hub_results$Correlation[i]
}

# Remove rows/cols with all zeros
heatmap_matrix <- heatmap_matrix[rowSums(abs(heatmap_matrix)) > 0,
                                 colSums(abs(heatmap_matrix)) > 0, drop = FALSE]

# 7. Plot heatmap
if(nrow(heatmap_matrix) > 1 & ncol(heatmap_matrix) > 1) {
  pheatmap(heatmap_matrix,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           main = "ncRNA × Hub Gene Correlations (Relaxed)",
           cluster_rows = TRUE, cluster_cols = TRUE,
           filename = "ncRNA_HubGene_Heatmap_RELAXED.pdf")
}


# NETWORK AND FINAL ANNOTATED HEATMAP


# Network Plot
edges <- hub_results %>%
  dplyr::mutate(weight = abs(Correlation),
                color = ifelse(Correlation > 0, "red", "blue")) %>%
  dplyr::select(ncRNA, hub_gene, weight, color)

# Create graph
g <- graph_from_data_frame(edges, directed = FALSE)

# Set node types
V(g)$type <- ifelse(V(g)$name %in% hub_genes, "hub_gene", "ncRNA")
V(g)$color <- ifelse(V(g)$type == "hub_gene", "gold", "skyblue")
V(g)$size <- ifelse(V(g)$type == "hub_gene", 10, 6)
E(g)$color <- edges$color
E(g)$width <- edges$weight * 3

# Plot network
pdf("ncRNA_HubGene_Network_RELAXED.pdf", width = 10, height = 8)
plot(g, layout = layout_with_fr(g),
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     main = "ncRNA–Hub Gene Co-expression Network (Relaxed)")
dev.off()

plot(g, layout = layout_with_fr(g),
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     main = "ncRNA–Hub Gene Co-expression Network (Relaxed)")

# Check overlap with up/down regulated ncRNAs
up_overlap <- intersect(up_ncRNAs, rownames(hub_cor_matrix))
down_overlap <- intersect(down_ncRNAs, rownames(hub_cor_matrix))

cat("Upregulated ncRNAs in correlation matrix:", length(up_overlap), "\n")
cat("Downregulated ncRNAs in correlation matrix:", length(down_overlap), "\n")

# Combined annotated heatmap

# 1. Filter significant correlations - RELAXED THRESHOLDS
sig_idx <- which(abs(hub_cor_matrix) >= cor_threshold & hub_fdr_matrix <= fdr_threshold, arr.ind = TRUE)
if(nrow(sig_idx) == 0){
  stop("No significant correlations at current thresholds.")
}

sig_ncRNAs <- unique(rownames(hub_cor_matrix)[sig_idx[,1]])
sig_cor_matrix <- hub_cor_matrix[sig_ncRNAs, , drop = FALSE]

# 2. Annotate ncRNAs as Up or Down
ncRNA_status <- ifelse(sig_ncRNAs %in% up_overlap, "Upregulated", "Downregulated")
row_annotation <- data.frame(Status = ncRNA_status)
rownames(row_annotation) <- sig_ncRNAs

# 3. Define annotation colors
ann_colors <- list(Status = c(Upregulated = "red", Downregulated = "blue"))

# 4. Replace non-significant correlations with 0
sig_cor_matrix_filtered <- sig_cor_matrix
sig_cor_matrix_filtered[abs(sig_cor_matrix_filtered) < cor_threshold | hub_fdr_matrix[sig_ncRNAs, ] > fdr_threshold] <- 0

# 5. Plot heatmap
pheatmap(sig_cor_matrix_filtered,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_row = row_annotation,
         annotation_colors = ann_colors,
         main = "ncRNA-Hub Gene Correlations (Up/Downregulated) - Relaxed",
         fontsize_row = 8,
         fontsize_col = 10,
         filename = "Combined_ncRNA_HubGene_Heatmap_RELAXED.pdf")





#### WGCNA Network Analysis ####

library(WGCNA)
library(dplyr)
library(dynamicTreeCut)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# PREREQUISITE CHECK
if (!exists("matrix_A_immune")) {
  stop("Error: 'matrix_A_immune' not found. Please load immune trait data from Part 2.")
}
if (!exists("master_norm_matrix")) {
  stop("Error: 'master_norm_matrix' not found. Please run Parts 1-2 first.")
}

#  DATA PREPARATION AND QC


# Prepare expression data (samples as rows, genes as columns)
datExpr <- as.data.frame(t(master_norm_matrix))

# Check for good genes and samples
gsg <- goodSamplesGenes(datExpr, verbose = 3)

if (!gsg$allOK) {
  cat("Removing bad samples/genes...\n")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Check for outliers using hierarchical clustering
sampleTree <- hclust(dist(datExpr), method = "average")

# Plot sample tree 
pdf("WGCNA_Sample_Tree.pdf", width = 12, height = 6)
par(mar = c(0, 5, 2, 0))
plot(sampleTree, main = "Sample Clustering to Detect Outliers", 
     sub = "", xlab = "", cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.5)
abline(h = mean(sampleTree$height) + 2 * sd(sampleTree$height), col = "red")
dev.off()

# Identify outliers 
sample_h <- mean(sampleTree$height) + 2 * sd(sampleTree$height)
clusters <- cutreeStatic(sampleTree, cutHeight = sample_h, minSize = 10)
outlier_samples <- rownames(datExpr)[clusters == 0]

if (length(outlier_samples) > 0) {
  cat("Potential outliers detected:\n")
  print(outlier_samples)
  cat("Consider removing these samples if they represent technical issues.\n")
} else {
  cat("No outliers detected.\n")
}

cat("Using", nrow(datExpr), "samples and", ncol(datExpr), "genes for analysis.\n")


# SOFT-THRESHOLD POWER SELECTION


powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

#  Add diagnostic plots for soft threshold
pdf("WGCNA_Soft_Threshold_Selection.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))

# Scale-free topology fit index
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale Independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.85, col = "red")

# Mean connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9, col = "red")

dev.off()

# Select power (default to 6 if automatic selection fails)
softPower <- ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate)
cat("Selected soft power:", softPower, "\n") #soft power is 4


# NETWORK CONSTRUCTION (blockwiseModules)

net <- blockwiseModules(
  datExpr,
  power = softPower,
  maxBlockSize = 5000,
  networkType = "unsigned",
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "TOM_ncRNA_Network",
  verbose = 3
)

# Proper module color assignment
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
names(moduleColors) <- colnames(datExpr)

# Module eigengene naming
MEs <- net$MEs
unique_numeric <- sort(unique(net$colors))
unique_colors <- labels2colors(unique_numeric)
colnames(MEs) <- paste0("ME", unique_colors)

# Save module eigengenes
write.csv(MEs, "WGCNA_Module_Eigengenes.csv")

# Module summary
cat("\n=== Module Summary ===\n")
print(table(moduleColors))
cat("\nModule Eigengenes:\n")
print(colnames(MEs))


# DATA FILTERING (Common Samples Across All Datasets)


# Find common samples
common_samples <- intersect(rownames(datExpr), rownames(matrix_A_immune))

if (length(common_samples) == 0) {
  stop("FATAL ERROR: No common samples between expression and immune data.")
}

# Filter all datasets to common samples
datExpr_filtered <- datExpr[common_samples, , drop = FALSE]
matrix_A_immune_filtered <- matrix_A_immune[common_samples, , drop = FALSE]
MEs_filtered <- MEs[common_samples, , drop = FALSE]
nSamples <- nrow(datExpr_filtered)

cat("Using", nSamples, "common samples for all analyses.\n")

# SAVE WGCNA RESULTS

save(net, moduleLabels, moduleColors, MEs_filtered, 
     datExpr_filtered, matrix_A_immune_filtered, softPower,
     file = "WGCNA_Filtered_Results.RData")
cat("Saved WGCNA_Filtered_Results.RData\n")



# HUB ncRNA IDENTIFICATION (Custom Logic)


unique_modules <- setdiff(unique(moduleColors), "grey")

for (mod in unique_modules) {
  
  cat("\nProcessing module:", mod, "\n")
  
  genes_in_module <- names(moduleColors)[moduleColors == mod]
  ME_name <- paste0("ME", mod)
  
  if (!ME_name %in% colnames(MEs_filtered) || length(genes_in_module) == 0) {
    cat("Skipping module", mod, "- ME not found or no genes.\n")
    next
  }
  
  # Extract module eigengene
  ME_mod <- MEs_filtered[, ME_name]
  
  # Module membership (MM)
  MM <- cor(datExpr_filtered[, genes_in_module, drop = FALSE], 
            ME_mod, use = "pairwise.complete.obs")
  MM_abs <- abs(as.numeric(MM))
  names(MM_abs) <- genes_in_module
  
  # Calculate gene significance (GS) for each immune trait
  GS_matrix <- cor(datExpr_filtered[, genes_in_module, drop = FALSE], 
                   matrix_A_immune_filtered, 
                   use = "pairwise.complete.obs", method = "spearman")
  
  # Average immune correlation score
  immune_scores <- rowMeans(abs(GS_matrix), na.rm = TRUE)
  
  # Hub score = MM × ImmuneScore
  hub_results <- data.frame(
    ncRNA = genes_in_module,
    Module = mod,
    ModuleMembership = MM_abs,
    ImmuneScore = immune_scores,
    HubScore = MM_abs * immune_scores
  ) %>% arrange(desc(HubScore))
  
  # Save module-specific hub genes
  write.csv(hub_results, paste0("Hub_ncRNAs_", mod, ".csv"), row.names = FALSE)
  
  cat("Top 5 hub ncRNAs in", mod, "module:\n")
  print(head(hub_results, 5))
}


# MODULE-TRAIT ASSOCIATION ANALYSIS

# Calculate correlations
moduleTraitCor <- cor(MEs_filtered, matrix_A_immune_filtered, 
                      use = "pairwise.complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Create text matrix for heatmap
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", 
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

#  Proper heatmap visualization
pdf("Module_Trait_Heatmap.pdf", width = 14, height = 10)
par(mar = c(8, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(matrix_A_immune_filtered),
               yLabels = colnames(MEs_filtered),
               ySymbols = colnames(MEs_filtered),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = "Module-Trait Relationships")
dev.off()

# Save results
write.csv(moduleTraitCor, "Module_Trait_Correlations.csv")
write.csv(moduleTraitPvalue, "Module_Trait_Pvalues.csv")



# MODULE DENDROGRAM VISUALIZATION


pdf("WGCNA_Dendrogram_All_Blocks.pdf", width = 15, height = 8)
for (blockNum in 1:length(net$dendrograms)) {
  if (!is.null(net$dendrograms[[blockNum]])) {
    genesInBlock <- net$blockGenes[[blockNum]]
    colorsInBlock <- moduleColors[genesInBlock]
    
    plotDendroAndColors(net$dendrograms[[blockNum]],
                        colorsInBlock,
                        "Module colors",
                        dendroLabels = FALSE,
                        hang = 0.03,
                        addGuide = TRUE,
                        guideHang = 0.05,
                        main = paste("Gene Dendrogram - Block", blockNum))
  }
}
dev.off()

cat("\n=== WGCNA Analysis Complete! ===\n")
cat("Generated files:\n")
cat("  - WGCNA_Sample_Tree.pdf\n")
cat("  - WGCNA_Soft_Threshold_Selection.pdf\n")
cat("  - WGCNA_Module_Eigengenes.csv\n")
cat("  - Hub_ncRNAs_[module].csv (per module)\n")
cat("  - Module_Trait_Heatmap.pdf\n")
cat("  - Module_Trait_Correlations.csv\n")
cat("  - WGCNA_Dendrogram_All_Blocks.pdf\n")
cat("  - WGCNA_Filtered_Results.RData\n")



####Answering key questions####

library(WGCNA)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(igraph)
library(tidyr)
library(gridExtra)

# Load previous results
if (!exists("moduleColors")) {
  load("WGCNA_Filtered_Results.RData")
}

# Ensure sig_results is available (from Part 1)
if (!exists("sig_results")) {
  stop("sig_results not found. Please load your differential expression results from Part 1.")
}


# QUESTION 1: Which of 706 ncRNAs are network hubs?

# Step 1: Calculate intramodular connectivity (kWithin)
# This measures how connected each gene is within its own module
adjacency <- adjacency(datExpr_filtered, power = softPower, type = "unsigned")
kWithin <- intramodularConnectivity(adjacency, moduleColors)

# Step 2: Calculate module membership (kME) for all modules
MEs_all <- moduleEigengenes(datExpr_filtered, moduleColors)$eigengenes
kME <- signedKME(datExpr_filtered, MEs_all)

# Step 3: Identify significant ncRNAs from your 706 list
sig_ncRNA_names <- sig_results$hgnc_symbol
sig_ncRNAs_in_network <- intersect(sig_ncRNA_names, names(moduleColors))

cat("Significant ncRNAs found in network:", length(sig_ncRNAs_in_network), "out of 706\n")

# Step 4: Create comprehensive hub analysis
hub_analysis <- data.frame(
  ncRNA = sig_ncRNAs_in_network,
  Module = moduleColors[sig_ncRNAs_in_network],
  kWithin = kWithin[sig_ncRNAs_in_network, "kWithin"],
  kTotal = kWithin[sig_ncRNAs_in_network, "kOut"] + kWithin[sig_ncRNAs_in_network, "kWithin"],
  stringsAsFactors = FALSE
)

# Add module membership for each ncRNA's own module
for (i in 1:nrow(hub_analysis)) {
  mod <- hub_analysis$Module[i]
  kME_col <- paste0("kME", mod)
  if (kME_col %in% colnames(kME)) {
    hub_analysis$kME[i] <- kME[hub_analysis$ncRNA[i], kME_col]
  } else {
    hub_analysis$kME[i] <- NA
  }
}

# Add logFC and regulation status
hub_analysis <- hub_analysis %>%
  left_join(sig_results %>% select(hgnc_symbol, logFC, adj.P.Val), 
            by = c("ncRNA" = "hgnc_symbol")) %>%
  mutate(Regulation = ifelse(logFC > 0, "Upregulated", "Downregulated"))

# Calculate hub score (combination of connectivity and module membership)
hub_analysis$HubScore <- abs(hub_analysis$kME) * hub_analysis$kWithin

# Identify top hubs (kME > 0.7, kWithin > median)
median_kWithin <- median(hub_analysis$kWithin, na.rm = TRUE)
hub_analysis$IsHub <- (abs(hub_analysis$kME) > 0.7 & hub_analysis$kWithin > median_kWithin)

# Sort by hub score
hub_analysis <- hub_analysis %>% arrange(desc(HubScore))

hub_top20 <- hub_analysis %>% 
  filter(IsHub == TRUE) %>%
  distinct(ncRNA, .keep_all = TRUE) %>%
  arrange(desc(HubScore)) %>%
  head(20)

# Save results
write.csv(hub_analysis, "Q1_Hub_ncRNA_Analysis.csv", row.names = FALSE)

# Print top 20 hub ncRNAs
cat("\n=== TOP 20 NETWORK HUB ncRNAs ===\n")
print(hub_analysis %>% filter(IsHub == TRUE) %>% head(20) %>% 
        select(ncRNA, Module, kWithin, kME, HubScore, Regulation, logFC))

cat("\nTotal hub ncRNAs identified:", sum(hub_analysis$IsHub, na.rm = TRUE), "\n")

# Visualization: Hub ncRNAs distribution across modules
p1 <- ggplot(hub_analysis %>% filter(IsHub == TRUE), 
             aes(x = Module, fill = Regulation)) +
  geom_bar(position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = c("Upregulated" = "#E64B35FF", "Downregulated" = "#4DBBD5FF")) +
  theme_minimal(base_size = 12) +
  labs(title = "Q1: Hub ncRNAs by Module and Regulation",
       x = "Module", y = "Number of Hub ncRNAs")

ggsave("Q1_Hub_ncRNAs_by_Module.pdf", p1, width = 10, height = 6)
print(p1)

cat("\nHub Definition: kME > 0.7 (strong module membership) AND")
cat(" kWithin > median (highly connected within module)\n")


# QUESTION 2: Do Upregulated vs Downregulated ncRNAs form separate networks?

# Create regulation annotation
regulation_status <- sig_results %>%
  select(hgnc_symbol, logFC) %>%
  mutate(Regulation = ifelse(logFC > 0, "Upregulated", "Downregulated")) %>%
  filter(hgnc_symbol %in% sig_ncRNAs_in_network)

# Count ncRNAs per module by regulation
module_regulation <- data.frame(
  ncRNA = sig_ncRNAs_in_network,
  Module = moduleColors[sig_ncRNAs_in_network]
) %>%
  left_join(regulation_status, by = c("ncRNA" = "hgnc_symbol")) %>%
  group_by(Module, Regulation) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Regulation, values_from = Count, values_fill = 0) %>%
  mutate(
    Total = Upregulated + Downregulated,
    UpregPercent = round(Upregulated / Total * 100, 1),
    DownregPercent = round(Downregulated / Total * 100, 1),
    Dominance = ifelse(UpregPercent > 75, "Up-dominated",
                       ifelse(DownregPercent > 75, "Down-dominated", "Mixed"))
  )

write.csv(module_regulation, "Q2_Module_Regulation_Composition.csv", row.names = FALSE)

cat("\n=== MODULE COMPOSITION BY REGULATION ===\n")
print(module_regulation %>% arrange(desc(Total)))

# Create aligned data first
analysis_data <- data.frame(
  ncRNA = sig_ncRNAs_in_network,
  Module = moduleColors[sig_ncRNAs_in_network],
  stringsAsFactors = FALSE
) %>%
  left_join(regulation_status, by = c("ncRNA" = "hgnc_symbol")) %>%
  filter(!is.na(Regulation))

# Now create contingency table
contingency_table <- table(
  Module = factor(analysis_data$Module),
  Regulation = factor(analysis_data$Regulation)
)

chi_test <- chisq.test(contingency_table)

cat("\nChi-square test for module-regulation independence:\n")
cat("X-squared =", chi_test$statistic, ", p-value =", chi_test$p.value, "\n")
if (chi_test$p.value < 0.05) {
  cat("RESULT: Upregulated and downregulated ncRNAs show SIGNIFICANT clustering into separate modules.\n")
} else {
  cat("RESULT: No significant segregation detected.\n")
}

# Visualization: Stacked bar plot
p2 <- ggplot(module_regulation %>% filter(Total >= 5), 
             aes(x = reorder(Module, Total), y = Total)) +
  geom_bar(aes(fill = "Upregulated", y = Upregulated), stat = "identity") +
  geom_bar(aes(fill = "Downregulated", y = Downregulated), stat = "identity", 
           position = position_stack()) +
  scale_fill_manual(values = c("Upregulated" = "#E64B35FF", "Downregulated" = "#4DBBD5FF")) +
  coord_flip() +
  geom_text(aes(label = paste0(UpregPercent, "%"), y = Upregulated/2), 
            size = 3.5, color = "black", fontface = "bold") +
  geom_text(aes(label = paste0(DownregPercent, "%"), y = Upregulated + Downregulated/2), 
            size = 3.5, color = "black", fontface = "bold") +
  theme_minimal(base_size = 12) +
  labs(title = "Q2: Module Composition by Regulation Status",
       x = "Module", y = "Number of ncRNAs", fill = "Regulation")

ggsave("Q2_Module_Regulation_Composition.pdf", p2, width = 10, height = 8)
print(p2)


# QUESTION 3: Brain Metastasis-Specific Gene Modules

# Step 1: Define group labels (Primary vs BM)
if (!exists("group_all")) {
  stop("group_all not found. Please ensure it's defined from Part 2.")
}

# Filter to common samples
group_filtered <- group_all[match(rownames(datExpr_filtered), colnames(master_norm_matrix))]

# Step 2: Calculate module eigengene differences between groups
ME_by_group <- data.frame(MEs_filtered, Group = group_filtered)

me_comparison <- ME_by_group %>%
  pivot_longer(cols = starts_with("ME"), names_to = "Module", values_to = "Eigengene") %>%
  group_by(Module, Group) %>%
  summarise(
    Mean = mean(Eigengene),
    SD = sd(Eigengene),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Group, values_from = c(Mean, SD))

# Step 3: Wilcoxon test for each module
me_wilcox_results <- data.frame()

for (me_col in colnames(MEs_filtered)) {
  primary_vals <- ME_by_group %>% filter(Group == "Primary") %>% pull(me_col)
  bm_vals <- ME_by_group %>% filter(Group == "BM") %>% pull(me_col)
  
  test <- wilcox.test(primary_vals, bm_vals)
  
  me_wilcox_results <- rbind(me_wilcox_results, data.frame(
    Module = me_col,
    Primary_Mean = mean(primary_vals),
    BM_Mean = mean(bm_vals),
    LogFC = log2((mean(bm_vals) + 0.01) / (mean(primary_vals) + 0.01)),
    P_value = test$p.value
  ))
}

me_wilcox_results$FDR <- p.adjust(me_wilcox_results$P_value, method = "BH")
me_wilcox_results$BM_Specific <- (me_wilcox_results$FDR < 0.05 & abs(me_wilcox_results$LogFC) > 0.5)
me_wilcox_results <- me_wilcox_results %>% arrange(FDR)

write.csv(me_wilcox_results, "Q3_BM_Specific_Modules.csv", row.names = FALSE)

cat("\n=== BM-SPECIFIC MODULES (FDR < 0.05, |logFC| > 0.5) ===\n")
print(me_wilcox_results %>% filter(BM_Specific == TRUE))

# Visualization: Module eigengene differences
bm_specific_modules <- me_wilcox_results %>% filter(BM_Specific == TRUE) %>% pull(Module)

if (length(bm_specific_modules) > 0) {
  plot_data <- ME_by_group %>%
    select(Group, all_of(bm_specific_modules)) %>%
    pivot_longer(cols = -Group, names_to = "Module", values_to = "Eigengene")
  
  p3 <- ggplot(plot_data, aes(x = Group, y = Eigengene, fill = Group)) +
    geom_boxplot() +
    facet_wrap(~ Module, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = c("Primary" = "#4DBBD5FF", "BM" = "#E64B35FF")) +
    theme_minimal(base_size = 11) +
    labs(title = "Q3: Brain Metastasis-Specific Module Eigengenes",
         y = "Module Eigengene")
  
  ggsave("Q3_BM_Specific_Modules.pdf", p3, width = 12, height = 8)
  print(p3)
}


# QUESTION 4: Immune-correlated ncRNA Networks

if (!exists("ncRNA_vs_immune_rho")) {
  # Calculate correlations between ncRNAs and immune cells
  ncRNA_expr_for_immune <- combined_expr_norm[sig_ncRNAs_in_network, ]
  common_samples_immune <- intersect(colnames(ncRNA_expr_for_immune), 
                                     rownames(matrix_A_immune))
  
  if(length(common_samples_immune) > 10) {
    ncRNA_vs_immune_rho <- cor(t(ncRNA_expr_for_immune[, common_samples_immune]),
                               matrix_A_immune[common_samples_immune, ],
                               method = "spearman", use = "pairwise.complete.obs")
  } else {
    stop("Insufficient samples for immune correlation analysis")
  }
}
# Step 1: Diagnostic check of correlation distribution
cat("\n=== Correlation Distribution Check ===\n")
cat("Correlation matrix dimensions:", nrow(ncRNA_vs_immune_rho), "ncRNAs x", 
    ncol(ncRNA_vs_immune_rho), "immune cells\n")
cat("Correlation range:", round(min(ncRNA_vs_immune_rho, na.rm=TRUE), 3), "to", 
    round(max(ncRNA_vs_immune_rho, na.rm=TRUE), 3), "\n")
cat("Median |correlation|:", round(median(abs(ncRNA_vs_immune_rho), na.rm=TRUE), 3), "\n")
cat("90th percentile |correlation|:", 
    round(quantile(abs(ncRNA_vs_immune_rho), 0.90, na.rm=TRUE), 3), "\n")

# Step 2: Use adaptive threshold (90th percentile, minimum 0.3)
cor_90th <- quantile(abs(ncRNA_vs_immune_rho), 0.90, na.rm = TRUE)
cor_threshold <- max(0.3, cor_90th)  # At least 0.3, but higher if data supports it

cat("\nUsing adaptive threshold: |r| >", round(cor_threshold, 3), "\n")

# Step 3: Identify immune-correlated ncRNAs
# Use lower minimum cells requirement (2 instead of 3)
sig_mask <- abs(ncRNA_vs_immune_rho) > cor_threshold
immune_corr_ncRNAs <- rownames(ncRNA_vs_immune_rho)[rowSums(sig_mask, na.rm = TRUE) >= 2]

cat("Immune-correlated ncRNAs found:", length(immune_corr_ncRNAs), 
    "(|r| >", round(cor_threshold, 3), "with ≥2 immune cells)\n")

# Step 4: Fallback if still none found
if (length(immune_corr_ncRNAs) == 0) {
  cat("\nWARNING: No ncRNAs met threshold. Using top 50 by max correlation...\n")
  max_corr_per_ncRNA <- apply(abs(ncRNA_vs_immune_rho), 1, max, na.rm = TRUE)
  immune_corr_ncRNAs <- names(sort(max_corr_per_ncRNA, decreasing = TRUE)[1:50])
  cat("Selected top 50 ncRNAs (max |r| range:", 
      round(range(max_corr_per_ncRNA[immune_corr_ncRNAs]), 3), ")\n")
}

module_composition <- data.frame(
  ncRNA = immune_corr_ncRNAs,
  Module = moduleColors[immune_corr_ncRNAs]
) %>%
  group_by(Module) %>%
  summarise(ImmuneCorr_ncRNA = n()) %>%
  arrange(desc(ImmuneCorr_ncRNA))

write.csv(module_composition, "Q4_Module_Composition.csv", row.names = FALSE)

cat("\n=== MODULES WITH IMMUNE-CORRELATED ncRNAs ===\n")
print(head(module_composition, 10))


# SUMMARY REPORT

cat("Q1: Network Hub ncRNAs\n")
cat("  - Total hub ncRNAs identified:", sum(hub_analysis$IsHub, na.rm = TRUE), "\n")
cat("  - Top hub module:", hub_analysis %>% filter(IsHub == TRUE) %>% 
      count(Module) %>% arrange(desc(n)) %>% head(1) %>% pull(Module), "\n\n")

cat("Q2: Network Segregation\n")
cat("  - Chi-square p-value:", format(chi_test$p.value, scientific = TRUE), "\n")
cat("  - Up-dominated modules:", sum(module_regulation$Dominance == "Up-dominated"), "\n")
cat("  - Down-dominated modules:", sum(module_regulation$Dominance == "Down-dominated"), "\n\n")

cat("Q3: BM-Specific Modules\n")
cat("  - Significant BM-specific modules:", sum(me_wilcox_results$BM_Specific, na.rm = TRUE), "\n")

if (sum(me_wilcox_results$BM_Specific, na.rm = TRUE) > 0) {
  cat("  - Top BM-specific module:", me_wilcox_results %>% 
        filter(BM_Specific == TRUE) %>% 
        arrange(FDR) %>% 
        head(1) %>% 
        pull(Module), "\n\n")
}

cat("Q4: ncRNA-Protein Networks\n")
cat("  - Immune-correlated ncRNAs:", length(immune_corr_ncRNAs), "\n\n")


save.image(file = "Co-expression_WGCNA.RData")