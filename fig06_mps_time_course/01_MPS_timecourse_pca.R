#Project: Live-cell Transcriptomics Manuscript
#Experiment: Microphysiological Systems Transcriptional Time-course
#Author: Mohamad Najia
#Objective: perform QC analysis on VLP cell media samples from Gag-EIF4E HUVEC MPS experiments with Griffith lab 

library(DESeq2)
library(jsonlite)
library(dplyr)
library(data.table)
library(tximport)
library(rhdf5)
library(ggplot2)
library(scales)
library(rjson)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(circlize)
library(irlba)
library(matrixStats)
library(M3C)
library(clusterProfiler)
library(enrichplot)


####################################################
#Function Declarations
####################################################

scaleRows <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
}

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

subset_txi <- function(txi, samples) {
  # check samples exist
  missing <- setdiff(samples, colnames(txi$counts))
  if (length(missing) > 0) {
    stop("These samples are not in txi: ", paste(missing, collapse = ", "))
  }
  
  # subset core matrices
  txi_subset <- txi
  for (nm in c("abundance", "counts", "length")) {
    if (!is.null(txi[[nm]])) {
      txi_subset[[nm]] <- txi[[nm]][, samples, drop = FALSE]
    }
  }
  
  # carry over countsFromAbundance if present
  if (!is.null(txi$countsFromAbundance)) {
    txi_subset$countsFromAbundance <- txi$countsFromAbundance
  }
  
  return(txi_subset)
}


####################################################
#Set-up Environment
####################################################

# Declare colormaps
pal_rna <- colorRampPalette(c("#352A86","#343DAE","#0262E0","#1389D2",
                              "#2DB7A3","#A5BE6A","#F8BA43","#F6DA23","#F8FA0D"))(100)
pal_atac <- colorRampPalette(c('#3361A5', '#248AF3', '#14B3FF', 
                               '#88CEEF', '#C1D5DC', '#EAD397', 
                               '#FDB31A','#E42A2A', '#A31D1D'))(100)

# Initialize variables
project_dir <- "/Volumes/mnajia/fig06_mps_time_course/"
kallisto_dir <- paste0(project_dir, "kallisto_output/")
output_dir <- paste0(project_dir, "analysis/")

# Import sample metadata
fn <- paste0(project_dir, "sample_metadata.txt")
df.samples <- fread(fn, header = TRUE, data.table = FALSE)
rownames(df.samples) <- df.samples$sample_name
df.samples$condition <- paste0(df.samples$timepoint, "_", df.samples$environment, "_", df.samples$HUVEC_condition, "_", df.samples$stroma)

# Import kallisto matrices
fn <- paste0(kallisto_dir, "HUVEC_MPS_timecourse_kallisto_txi.rds")
txi <- readRDS(fn)



####################################################
#Plot signal-to-background of cell media measurements
####################################################

#quantify and plot number of expressed genes
df.samples$expressed_genes <- colSums(txi$counts > 1e-1)
df.samples$grouping <- paste0(df.samples$HUVEC_condition, "_", df.samples$stroma)

df.plot <- df.samples 

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(df.plot, varname="expressed_genes",
                    groupnames=c("timepoint", "HUVEC_condition", "grouping"))
df2$timepoint <- as.numeric(str_split_fixed(df2$timepoint, pattern = "Day0", n=2)[,2])
df2$grouping <- factor(df2$grouping, levels = c("HUVEC-Gag-EIF4E_NHLF","HUVEC-Gag-EIF4E_ESC","Control_NHLF"))

pdf(paste0(output_dir, "QC_signal_to_background_expressed_genes.pdf"), width = 5.5, height = 3.5, useDingbats = FALSE)

ggplot(df2, aes(x=timepoint, y=expressed_genes, group=grouping, color=grouping)) +
  geom_line() +
  geom_pointrange(aes(ymin=expressed_genes-sd, ymax=expressed_genes+sd)) +
  geom_errorbar(aes(ymin=expressed_genes-sd, ymax=expressed_genes+sd), width = 0.25) +
  scale_color_manual(
    values = c(
      #"HUVEC-Gag-EIF4E_ESC_IL1b" = "#1B9E77",
      "HUVEC-Gag-EIF4E_NHLF" = "#D95F02",
      "HUVEC-Gag-EIF4E_ESC" = "#7570B3",
      "Control_NHLF" = "grey")
  ) + 
  scale_x_continuous(limits = c(0.5, 7.5), breaks = c(1,2,3,4,5,6,7)) +
  scale_y_continuous(limits = c(0, 15000), breaks = c(0,3e3,6e3,9e3,12e3,15e3)) +
  xlab("Day") +
  ylab("Genes detected") +
  ggtitle("Signal-to-background") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12),
        axis.line = element_line(colour = 'black', linewidth = 0.25),
        panel.border = element_rect(fill=NA, colour = "black", linewidth = 0.25),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

dev.off()



####################################################
#Principal Components Analysis on Stroma Conditions
####################################################

# Focus on Gag-EIF4E expressing HUVECs for all stromal conditions
df.samples.subset <- df.samples[df.samples$HUVEC_condition == "HUVEC-Gag-EIF4E",]
txi <- subset_txi(txi, df.samples.subset$sample_name)

# Set up DESeq 
sampleTable <- data.frame(condition = df.samples.subset$condition)
rownames(sampleTable) <- df.samples.subset$sample_name

dds <- DESeqDataSetFromTximport(txi, sampleTable, design = ~ condition)
dds <- DESeq(dds)

# Perform PCA using variance-stabilizing transformation
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds.keep <- dds[keep,]
nrow(dds.keep)

vsd <- vst(dds.keep, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)

pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE) 
pcaData$grouping <- df.samples.subset$grouping
pcaData$timepoint <- df.samples.subset$timepoint
pcaData$timepoint <- factor(
  pcaData$timepoint,
  levels = sort(unique(pcaData$timepoint))
)
pcaData$sampleID <- df.samples.subset$sample_ID
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$grouping <- factor(pcaData$grouping, levels = c("HUVEC-Gag-EIF4E_NHLF","HUVEC-Gag-EIF4E_ESC"))

gg <- ggplot(pcaData, aes(x=PC1, y=PC2, fill=timepoint, color=timepoint)) + 
  geom_path(aes(group = sampleID), color = "grey70", linewidth = 0.25, alpha = 0.5) + 
  geom_point() + 
  scale_x_continuous(limits = c(-30, 20)) +
  scale_y_continuous(limits = c(-15, 15)) +
  xlab(paste0("PC1: (", percentVar[1], "% Variance Explained)")) +
  ylab(paste0("PC2: (", percentVar[2], "% Variance Explained)")) +
  ggtitle("PCA on VLP-RNAs from HUVEC-Gag-EIF4E cells co-cultured with various stroma") + 
  facet_wrap(~grouping) + 
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        legend.title = element_blank())

pdf(file = paste(output_dir, "QC_pca_vst_Pump_HUVEC-Gag-EIF4E_NHLF_ESC_stroma_timecourses.pdf", sep = "/"), width = 6, height = 3, useDingbats = FALSE)  
print(gg)
dev.off()

#investigate gene loadings for PC1 (time) and PC2 (stroma condition)
mat <- assay(vsd)
pca <- prcomp(t(mat), scale. = FALSE)

pc1_loadings <- pca$rotation[, 1]
pc2_loadings <- pca$rotation[, 2]

head(sort(pc1_loadings, decreasing = TRUE), 20)
head(sort(pc1_loadings, decreasing = FALSE), 20)

head(sort(pc2_loadings, decreasing = TRUE), 20)
head(sort(pc2_loadings, decreasing = FALSE), 20)

#plot top genes by PC loadings
df.pc1 <- data.frame(
  gene = names(pc1_loadings),
  loading = as.numeric(pc1_loadings)
)

df.pc1.plot <- df.pc1 %>%
  arrange(loading) %>%
  slice(c(1:15, (n() - 14):n())) %>%
  mutate(
    gene = reorder(gene, loading),
    direction = ifelse(loading > 0, "Positive", "Negative")
  )

df.pc2 <- data.frame(
  gene = names(pc2_loadings),
  loading = as.numeric(pc2_loadings)
)

df.pc2.plot <- df.pc2 %>%
  arrange(loading) %>%
  slice(c(1:15, (n() - 14):n())) %>%
  mutate(
    gene = reorder(gene, loading),
    direction = ifelse(loading > 0, "Positive", "Negative")
  )

gg1 <- ggplot(df.pc1.plot, aes(x = gene, y = loading, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Negative" = "steelblue", "Positive" = "firebrick")) +
  #facet_wrap(~PC, scales = "free") + 
  xlab("Gene") + 
  ylab("PC Loading") +
  ggtitle("PC1") + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

gg2 <- ggplot(df.pc2.plot, aes(x = gene, y = loading, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Negative" = "steelblue", "Positive" = "firebrick")) +
  #facet_wrap(~PC, scales = "free") + 
  xlab("") + 
  ylab("PC Loading") +
  ggtitle("PC2") + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf(file = paste(output_dir, "QC_pca_vst_Pump_HUVEC-Gag-EIF4E_NHLF_ESC_stroma_PCA_PC1_PC2_gene_loadings.pdf", sep = "/"), width = 6, height = 4, useDingbats = FALSE)  
gg1+gg2
dev.off()

#save PC1 and PC2 gene loadings
ranked_genes <- sort(pc1_loadings, decreasing = TRUE)

ranked_genes_df <- data.frame(
  gene = names(ranked_genes),
  loading = as.numeric(ranked_genes)
)

write.table(
  ranked_genes_df,
  file = paste0(output_dir, "QC_pca_vst_Pump_HUVEC-Gag-EIF4E_NHLF_ESC_PC1_gene_loadings.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


ranked_genes <- sort(pc2_loadings, decreasing = TRUE)

ranked_genes_df <- data.frame(
  gene = names(ranked_genes),
  loading = as.numeric(ranked_genes)
)

write.table(
  ranked_genes_df,
  file = paste0(output_dir, "QC_pca_vst_Pump_HUVEC-Gag-EIF4E_NHLF_ESC_PC2_gene_loadings.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

#perform pathway enrichment on genes ranked by PC1 (time) and PC2 (condition) loadings
library(msigdbr)
library(fgsea)

#hallmark
hallmark <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)

set.seed(123)
fgsea_res <- fgsea(
  pathways = hallmark_list,
  stats    = sort(pc1_loadings, decreasing = TRUE),
  minSize  = 15,
  maxSize  = 500,
  nperm    = 10000
)

pdf(file = paste(output_dir, "QC_pca_vst_Pump_HUVEC-Gag-EIF4E_NHLF_ESC_PC1_loadings_Hallmark_GSEA.pdf", sep = "/"), width = 6, height = 4, useDingbats = FALSE)  

fgsea_res %>%
  filter(padj < 0.05) %>%
  mutate(pathway = reorder(pathway, NES)) %>%
  ggplot(aes(NES, pathway, color = NES)) +
  geom_point(size = 3) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  ylab("Hallmark gene set (Padj < 0.05)") + 
  ggtitle("GSEA using PC1 gene loadings") + 
  theme_bw()

dev.off()

tmp <- fgsea_res %>% filter(padj < 0.05) %>% mutate(pathway = reorder(pathway, NES)) %>% as.data.frame()
tmp$leadingEdge <- vapply(tmp$leadingEdge, paste, collapse = "/", character(1L))

write.table(
  tmp,
  file = paste0(output_dir, "QC_pca_vst_Pump_HUVEC-Gag-EIF4E_NHLF_ESC_PC1_loadings_Hallmark_GSEA_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)



fgsea_res <- fgsea(
  pathways = hallmark_list,
  stats    = sort(pc2_loadings, decreasing = TRUE),
  minSize  = 15,
  maxSize  = 500,
  nperm    = 10000
)

pdf(file = paste(output_dir, "QC_pca_vst_Pump_HUVEC-Gag-EIF4E_NHLF_ESC_PC2_loadings_Hallmark_GSEA.pdf", sep = "/"), width = 6, height = 4, useDingbats = FALSE)  

fgsea_res %>%
  filter(padj < 0.05) %>%
  mutate(pathway = reorder(pathway, NES)) %>%
  ggplot(aes(NES, pathway, color = NES)) +
  geom_point(size = 3) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  ylab("Hallmark gene set (Padj < 0.05)") + 
  ggtitle("GSEA using PC2 gene loadings") + 
  theme_bw()

dev.off()

tmp <- fgsea_res %>% filter(padj < 0.05) %>% mutate(pathway = reorder(pathway, NES)) %>% as.data.frame()
tmp$leadingEdge <- vapply(tmp$leadingEdge, paste, collapse = "/", character(1L))

write.table(
  tmp,
  file = paste0(output_dir, "QC_pca_vst_Pump_HUVEC-Gag-EIF4E_NHLF_ESC_PC2_loadings_Hallmark_GSEA_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


#reactome
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

ranked_genes <- sort(pc1_loadings, decreasing = TRUE)

gene_df <- bitr(
  names(ranked_genes),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
ranked_genes_entrez <- ranked_genes[gene_df$SYMBOL]
names(ranked_genes_entrez) <- gene_df$ENTREZID

reactome_gsea <- gsePathway(
  geneList     = ranked_genes_entrez,
  organism     = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize    = 15,
  maxGSSize    = 500,
  verbose      = FALSE
)

pdf(file = paste(output_dir, "QC_pca_vst_Pump_HUVEC-Gag-EIF4E_NHLF_ESC_PC1_loadings_Reactome_GSEA.pdf", sep = "/"), width = 10, height = 5, useDingbats = FALSE)  

reactome_gsea@result %>%
  filter(p.adjust < 0.05) %>%
  arrange(NES) %>%
  slice(c(1:10, (n() - 19):n())) %>%   
  mutate(Description = reorder(Description, NES)) %>%
  ggplot(aes(NES, Description, color = NES)) +
  geom_point(size = 3) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  ylab("Reactome gene set (Padj < 0.05)") + 
  ggtitle("GSEA using PC1 gene loadings") + 
  theme_bw()

dev.off()



ranked_genes <- sort(pc2_loadings, decreasing = TRUE)

gene_df <- bitr(
  names(ranked_genes),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
ranked_genes_entrez <- ranked_genes[gene_df$SYMBOL]
names(ranked_genes_entrez) <- gene_df$ENTREZID

reactome_gsea <- gsePathway(
  geneList     = ranked_genes_entrez,
  organism     = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize    = 15,
  maxGSSize    = 500,
  verbose      = FALSE
)

pdf(file = paste(output_dir, "QC_pca_vst_Pump_HUVEC-Gag-EIF4E_NHLF_ESC_PC2_loadings_Reactome_GSEA.pdf", sep = "/"), width = 10, height = 5, useDingbats = FALSE)  

reactome_gsea@result %>%
  filter(p.adjust < 0.05) %>%
  arrange(NES) %>%
  slice(c(1:10, (n() - 19):n())) %>%   
  mutate(Description = reorder(Description, NES)) %>%
  ggplot(aes(NES, Description, color = NES)) +
  geom_point(size = 3) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  ylab("Reactome gene set (Padj < 0.05)") + 
  ggtitle("GSEA using PC2 gene loadings") + 
  theme_bw()

dev.off()



