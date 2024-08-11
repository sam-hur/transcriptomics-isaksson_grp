# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# -- Dependencies
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("GenomicFeatures")
# BiocManager::install("apeglm")
# BiocManager::install("pheatmap")
# BiocManager::install("vsn")
# BiocManager::install("tximport")
# BiocManager::install("readr")

# install.packages("ggtext")
library(DESeq2)
library(GenomicFeatures)
library(stringr)
library(pheatmap)
library(ggplot2)
library(vsn)
library(purrr)
library(tibble)
library(reshape2)
library(tximport)
library(readr)
library(dplyr)
library(ggtext)
library(tidyr)
# --
# directory <- "D:\\Projects\\Bioinformatics\\Transcriptomics\\outputs\\Readcounts"
count_matrices <- "D:/Projects/Bioinformatics/Transcriptomics/outputs/counts_v2/stringtie-prepDE/_csv"
treatments <- c("ALAN", "noise", "soot", "control")

gene_count_matrix <- file.path(count_matrices, 'gene_count_matrix.csv')
transcript_count_matrix <- file.path(count_matrices, 'transcript_count_matrix.csv')

gene_matrix <- read.csv(gene_count_matrix, row.names=1)
transcript_matrix <- read.csv(transcript_count_matrix, row.names=1)

#========================================


sampleTable <- data.frame(
  sampleName=colnames(gene_matrix),
  condition=c(
    rep('ALAN', 5),
    rep('control', 5),
    rep('noise', 6),
    rep('soot', 6)
  )
)
colnames(sampleTable)
dds_data <- DESeqDataSetFromMatrix(countData = gene_matrix,
                                   colData = sampleTable,
                                   design = ~ condition)

#=============================================================================



# repeat DGE using each treatment as a reference base
for (treatment in treatments){
  print(treatment)
  
  
  # By default, R will choose a reference level for factors based on alphabetical order
  # Explicitly set reference level to the current treatment:
  dds_data$condition <- relevel(dds_data$condition, ref = treatment) 
  
  
  dim(dds_data)
  dds_data <- dds_data[rowSums(counts(dds_data, normalized=F)) >= 10, ] # filter out where count < 10
  dim(dds_data)
  
  dds <- DESeq(dds_data)
  
  # dds <- dds[rowSums(counts(dds, normalized=F)) >= 10, ] # filter out where count < 10
  dim(dds)
  
  read_counts <- counts(dds, normalized = FALSE)
  head(sort(read_counts, decreasing = T), 25)
  resNames <- resultsNames(dds)
  treatment.conditions <- c(
    resNames[2],
    resNames[3],
    resNames[4]
  )
  
  get_contrast <- function (cond) {
    tmp <- unlist(strsplit(cond, "_"))
    return(c(tmp[2], tmp[4]))
  }
  cond1 <- get_contrast(treatment.conditions[1])
  cond2 <- get_contrast(treatment.conditions[2])
  cond3 <- get_contrast(treatment.conditions[3])
  
  c_conds <- list(cond1, cond2, cond3)
  
  
  
  comparisons.pairwise <- list(
    # Comparisons vs. control
    results(dds, contrast = c("condition", cond1[1], cond1[2]), alpha = 0.05),
    results(dds, contrast = c("condition", cond2[1], cond2[2]), alpha = 0.05),
    results(dds, contrast = c("condition", cond3[1], cond3[2]), alpha = 0.05)
  )
  
  #__________ lists for diff comparisons
  MAplots <- list()
  MAplots.significant <- list()
  comparisons.significant <- list()
  
  resLFC <- list()
  #_____________________________________
  
  
  draw_MAplot <- function(res, name, padj=0.05){
    title <- paste("MA Plot of DESeq results for", name)
    sig_label <- paste("Significant padj (", padj, ")")
    plot.ma <- plotMA(res, cex=0.75, main = title)
    abline(h=c(-1,1), col="dodgerblue", lwd=2)
    legend("topright", legend = c(sig_label, "Not Significant"), 
           col = c("blue", "grey"), pch = 16, cex = 0.5)
  }

      
  for(i in 1:length(treatment.conditions)){
    res <- comparisons.pairwise[[i]]
    
    cat("Summary for", c_conds[[i]][1], "vs.", c_conds[[i]][2],":\n")
    print(capture.output(summary(res)), sep="\n")
    
    name <- paste(c_conds[[i]][1], "vs.", c_conds[[i]][2], "(ref)")
    
    # plot results
    draw_MAplot(res, name)
    MAplots[[name]] <- recordPlot()
    MAplots
    # subset to only include significant genes
    sig <- subset(res, padj < 0.001 & abs(log2FoldChange) > 1) # keep only the v. significant results
    comparisons.significant[[name]] <- sig
    head(sort(sig), 10)
    summary(sig)
    
    draw_MAplot(sig, name)
    MAplots.significant[[name]] <- recordPlot()
  
    
    # there are a lot of NA values in the padj of the fold shrinkage. Not sure why.
    # but they are most probably insignificant, so I'll just omit them
    resLFC[[treatment.conditions[1]]] <- na.omit(lfcShrink(dds, coef=treatment.conditions[1], type="apeglm"))
    resLFC[[treatment.conditions[2]]] <- na.omit(lfcShrink(dds, coef=treatment.conditions[2], type="apeglm"))
    resLFC[[treatment.conditions[3]]] <- na.omit(lfcShrink(dds, coef=treatment.conditions[3], type="apeglm"))

    
    for (j in 1:length(resLFC)){
      res <- resLFC[[j]]
      res <- res[res$padj < 0.001, ]
      res_name <- names(resLFC)[j]
      summary(res)
      head(res$pvalue, 10)
      dim(res)
      upreg <- rownames(res[order(-res$log2FoldChange),])[1:10]
      downreg <- rownames(res[order(res$log2FoldChange),])[1:10]
      top_genes <- c(upreg, downreg)
      top_genes
      
      ntd <- normTransform(dds)
      
      colnames(ntd)
      
      matrix <- assay(ntd)[top_genes,]
      matrix

      predefined_order <- c(treatment, setdiff(
        c("control", "soot", "noise", "ALAN"), 
        treatment)
      )
      predefined_order
      
      
      annotation <- as.data.frame(colData(dds)[, "condition", drop = FALSE])
      # add a type condition, which is rep as we only have 1 type
      annotation["type"] <- factor(rep("paired-end", 22))
      annotation <- annotation[order(annotation$condition), ]
      annotation$condition <- factor(annotation$condition, levels = predefined_order)
      annotation
      
      
      dim(annotation)
      
      
      title <- paste("Pheatmap of", res_name, "(ref)")
      pheatmap(
        matrix, 
        cluster_rows=T, 
        show_rownames=T,
        cluster_cols=F,
        annotation_col=annotation,
        main = title
      )
     
      filename <- paste0(res_name,"(ref)-top_genes.csv")
      filename
      
      write.csv(matrix, filename)
      
      rownames(matrix)
      matrix_filtered <- matrix[grep("\\|", rownames(matrix)), ]
      rownames(matrix_filtered) <- sapply(strsplit(rownames(matrix_filtered), "\\|"), `[`, 2)
      matrix_filtered
      
      title_filtered <- paste("Pheatmap of", res_name, "(ref) - filtered")
      pheatmap(
        matrix_filtered, 
        cluster_rows=T, 
        show_rownames=T,
        cluster_cols=F,
        annotation_col=annotation,
        main = title_filtered,
        filename = paste0("pheatmap-", res_name, ".png")
      )
      
      
      filename <- paste0(res_name,"(ref)-top_genes.csv")
      filename
      
      write.csv(matrix, filename)
      
      cnts <- as.data.frame(read_counts[rownames(read_counts) %in% top_genes, ])
      cnts
      
      df <- data.frame(Gene=rownames(cnts), Frequency=rowSums(cnts), row.names = NULL)
      df
      
      cnts_per_treatment <- as.data.frame(cnts[rownames(cnts) %in% top_genes, ])
      
      head(cnts_per_treatment, 10)
      df_treatment <- as.data.frame(cnts[])
      df_treatment
      
      
      ALAN <- rowSums(cnts_per_treatment[, grepl("ALAN", colnames(cnts_per_treatment))])
      control <- rowSums(cnts_per_treatment[, grepl("control", colnames(cnts_per_treatment))])
      noise <- rowSums(cnts_per_treatment[, grepl("noise", colnames(cnts_per_treatment))])
      soot <- rowSums(cnts_per_treatment[, grepl("soot", colnames(cnts_per_treatment))])
      
      cnts_all <- data.frame(
        Gene = rownames(cnts_per_treatment),
        ALAN = ALAN,
        control = control,
        noise = noise,
        soot = soot,
        row.names = NULL
      )
      cnts_all
      
      
      
      #  Removing novel transcripts without gene name
      cnts_all_filtered <- cnts_all[grep("\\|", cnts_all$Gene), ]
      cnts_all_filtered$Gene <- sapply(strsplit(cnts_all_filtered$Gene, "\\|"), `[`, 2)
      cnts_all_filtered
      
      
      cnts_all$total_counts <- rowSums(cnts_all[, treatments])
      cnts_all$Gene <- reorder(cnts_all$Gene, cnts_all$total_counts)
      
      cnts_all_long <- cnts_all %>%
        select(-total_counts) %>%
        pivot_longer(cols = c(
          predefined_order[1],
          predefined_order[2],
          predefined_order[3],
          predefined_order[4]
        ),
        names_to = "Condition",
        values_to = "Count")
      
      cnts_all_filtered$total_counts <- rowSums(cnts_all_filtered[, treatments])
      cnts_all_filtered$Gene <- reorder(cnts_all_filtered$Gene, cnts_all_filtered$total_counts)
      cnts_all_filtered_long <- cnts_all_filtered %>%
        select(-total_counts) %>%
        pivot_longer(cols = c(
          predefined_order[1],
          predefined_order[2],
          predefined_order[3],
          predefined_order[4]
        ),
        names_to = "Condition",
        values_to = "Count")
      cnts_all_filtered_long
      
      legend_labels <- c(
        ALAN = ifelse("ALAN" == predefined_order[1], "**ALAN**", "ALAN"),
        control = ifelse("control" == predefined_order[1], "**control**", "control"),
        soot = ifelse("soot" == predefined_order[1], "**soot**", "soot"),
        noise = ifelse("noise" == predefined_order[1], "**noise**", "noise")
      )
      
      
      plt <- ggplot(cnts_all_long, aes(x = Gene, y=Count, fill=Condition)) +
        geom_bar(stat = "identity", position = "stack") +
        labs(title = paste0("Sum of Read Counts for Each Gene by Condition\n",res_name),
             x = "Gene",
             y = "Sum of Read Counts",
             fill = "Condition") +
        theme_minimal() +
        scale_fill_discrete(labels = legend_labels) +  # Only modify labels
        theme(legend.text = element_markdown(),
              axis.text.x = element_text(angle = 90, hjust = 1),
              plot.title = element_text(hjust = 0.5, size = 12, margin = margin(t = 20, b = 20))) + coord_flip()
      plt
      ggsave(paste0(res_name,"-top_genes_hist.png"), plt, height = 5, width = 10, bg="white")
      
      
      plt <- ggplot(cnts_all_filtered_long, aes(x=Gene, y=Count, fill=Condition)) +
        geom_bar(stat = "identity", position = "stack") +
        labs(title = paste0("Sum of Read Counts for Each Gene by Condition\n",res_name, " (filtered)"),
             x = "Gene",
             y = "Sum of Read Counts",
             fill = "Condition") +
        theme_minimal() +
        scale_fill_discrete(labels = legend_labels) +  # Use manual labels with markdown
        theme(legend.text = element_markdown(),
              axis.text.y = element_text(angle = 0, hjust = 1),  # Adjust for horizontal bars
              plot.title = element_text(hjust = 0.5, size = 12, margin = margin(t = 20, b = 20))) +
        coord_flip()
      plt
      ggsave(paste0(res_name,"-top_genes_hist-filtered.png"), plt, height = 5, width = 10, bg="white")
    }
  }
}
