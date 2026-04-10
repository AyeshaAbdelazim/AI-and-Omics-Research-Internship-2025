
classify_gene <- function(log2FC ,padj) {
  if(padj < 0.05 && log2FC > 1  ) {
    return("Upregulated")
  } else if(padj < 0.05 && log2FC < -1) {
      return("Downregulated")
  }else {
      return("Not_Significant")}
    }

input_dir <- "Raw_data_2"
output_dir <- "Results_2"
if(!dir.exists(input_dir)) {dir.create(input_dir)}
if(!dir.exists(output_dir)) {dir.create(output_dir)}

files_to_process <- c("DEGs_Data_1.csv" , "DEGs_Data_2.csv")

result_list <- list()

for (file_name in files_to_process) {
  cat("\nProcessing:", file_name, "\n")  
  input_file_path <- file.path(input_dir, file_name)
  
  data <- read.csv(input_file_path, header = TRUE)
  cat("File imported. Checking for missing values...\n")
  
  data$padj[is.na(data$padj)] <- 1
  data$status <- NA
  
  for (i in 1:nrow(data)) {
    data$status[i] <- classify_gene(data$logFC[i], data$padj[i])
  }
  
  output_file <- file.path(output_dir, paste0("Processed_", file_name))
  write.csv(data, output_file, row.names = FALSE)
  
  cat("Summary for", file_name, ":\n")
  up_count <- sum(data$status == "Upregulated")
  down_count <- sum(data$status == "Downregulated")
  not_sig_count <- sum(data$status == "Not_Significant")
  sig_count <- up_count + down_count
  
  cat("Total significant genes:", sig_count, "\n")
  cat("Upregulated genes:", up_count, "\n")
  cat("Downregulated genes:", down_count, "\n")
  cat("Not significant genes:", not_sig_count, "\n")
}

save.image("Assignment2.RData")

