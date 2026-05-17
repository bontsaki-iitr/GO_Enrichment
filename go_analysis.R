if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
    !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Missing required R packages. Please ensure your Conda environment is set up correctly.")
}

library(clusterProfiler)
library(org.Hs.eg.db)

gene_file <- "hg38_GCGC_target_genes.tsv"

if (!file.exists(gene_file)) {
    stop(paste("Error: Input file", gene_file, "not found in the current directory!"))
}

genes <- read.table(gene_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)$V1
cat("Loaded", length(genes), "genes for analysis.\n")

cat("Running GO Enrichment Analysis...\n")
ego <- enrichGO(gene          = genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.01)

if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    warning("No significantly enriched GO terms found with the current thresholds.")
    write.csv(data.frame(Message = "No significant results"), "hg38_GCGC_GO_results.csv")
} else {
    write.csv(as.data.frame(ego), "hg38_GCGC_GO_results.csv")
    cat("Results saved to 'hg38_GCGC_GO_results.csv'.\n")

    pdf("hg38_GCGC_GO_dotplot.pdf", height = 8, width = 8)
    print(dotplot(ego, showCategory = 20, font.size = 6))
    dev.off()
    cat("Plot saved to 'hg38_GCGC_GO_dotplot.pdf'.\n")
}
