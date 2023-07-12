#!/usr/bin/env Rscript
library(bigsnpr)
args <- commandArgs(trailingOnly = TRUE)
# Check if argument is provided
if (length(args) >= 3) {
    # Process the argumentS
    SAVE_DIR <- args[1]
    REF_FILE <- args[2]
    SAMPLE_FILE <- args[3]
    sample_file <- paste0(SAVE_DIR, "/pca/", SAMPLE_FILE, ".bed")
    ref_file <- paste0(SAVE_DIR, "/pca_reg/", REF_FILE, ".bed")

    # Load BED files
    ref.bed <- bed(ref_file)
    samples.bed <- bed(sample_file)

    #Run analysis
    sample_projected <- bed_projectPCA(ref.bed, sample.bed, k=20)

    #Save Sample PCs
    sample_pc <- paste0(SAVE_DIR, "/pca/samples_pcs.csv")
    write.csv(sample_projected$OADP_proj, file =sample_pc, row.names = FALSE)
    #Save Ref loadings
    save(obj.svd2, file = paste0(SAVE_DIR, "/pca_ref/obj_svd2.RData"))
    write.csv(obj.svd2$d, file = paste0(SAVE_DIR, "/pca_ref/d_values.csv"), row.names = FALSE)
    write.csv(obj.svd2$u, file = paste0(SAVE_DIR, "/pca_ref/u_values.csv"), row.names = FALSE)
    write.csv(obj.svd2$v, file = paste0(SAVE_DIR, "/pca_ref/v_loadings.csv"), row.names = FALSE)
    write.csv(obj.svd2$center, file = paste0(SAVE_DIR, "/pca_ref/center_values.csv"), row.names = FALSE)
    write.csv(obj.svd2$scale, file = paste0(SAVE_DIR, "/pca_ref/scale_values.csv"), row.names = FALSE)


} 
else {
  # Argument not provided, handle the case
  stop("All arguments not provided. Please provide save directory, name of reference bed file and samples bef file.")
}
