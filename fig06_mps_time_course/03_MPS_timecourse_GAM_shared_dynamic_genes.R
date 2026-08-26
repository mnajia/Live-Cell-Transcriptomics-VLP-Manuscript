#Project: Live-cell Transcriptomics Manuscript
#Experiment: Microphysiological Systems Transcriptional Time-course
#Author: Mohamad Najia
#Objective: identify *shared* temporal gene dynamics across stroma conditions by 
#           applying joint Generalized Additive Models (GAM) to VLP cell media 
#           samples produced from Gag-EIF4E expressing HUVECs

library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(scales)
library(cluster)
library(stringr)
library(tximport)
library(DESeq2)
library(mgcv)
library(ComplexHeatmap)
library(circlize)



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
#Fit a Joint GAM to Identify Temporally Significant Genes
####################################################

# Focus on HUVEC Gag-EIF4E samples
df.samples <- df.samples %>% filter(HUVEC_condition == "HUVEC-Gag-EIF4E")
txi <- subset_txi(txi, df.samples$sample_name)

# Create DESeq2 object to normalize kallisto counts
sampleTable <- data.frame(sample = df.samples$sample_name,
                          subject = factor(df.samples$sample_ID),
                          condition = factor(df.samples$stroma),
                          time = as.numeric(str_split_fixed(df.samples$timepoint, pattern = "Day0", n=2)[,2])
)
rownames(sampleTable) <- df.samples$sample_name

dds <- DESeqDataSetFromTximport(
  txi = txi,
  colData = sampleTable,
  design = ~ subject + time  
)

# Run vst to generate variance-stabilized expression counts 
#Note: GAMs assume homoskedastic, approximately Gaussian errors... RNA-seq counts violate this, but VST-normalized values satisfy the assumptions reasonably well
vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)

# Prepare data for GAM fitting
expr_df <- vst_mat %>%
  as.data.frame() %>%
  mutate(gene = rownames(.)) %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "expr"
  ) %>%
  left_join(sampleTable, by = "sample")

expr_df <- expr_df %>%
  mutate(
    subject = factor(subject),
    condition = factor(condition),
    time = as.numeric(time)
  )

# Fit a joint GAM with random subject effects for each gene
fit_gene <- function(g) {
  
  df_g <- expr_df[expr_df$gene == g, ]
  
  failed_df <- data.frame(gene = g, 
                          p_condition = as.numeric(NA),
                          p_time = as.numeric(NA),
                          p_time_condition = as.numeric(NA)
  )
  
  # Safety checks
  if (nrow(df_g) < 5 || 
      sd(df_g$expr) == 0 || 
      length(unique(df_g$time)) < 2) {
    
    return(failed_df)
  }
  
  # Fit model with try()
  fit <- try(
    gam(expr ~ condition + s(time, k=5) + s(time, condition, bs="fs", k=5) + s(subject, bs="re"),
        data=df_g, method="REML"),
    silent = TRUE
  )
  
  # If fit failed → return NA
  if (inherits(fit, "try-error")) {
    return(failed_df)
  }
  
  # Extract Parametric coefficients
  p_table <- summary(fit)$p.table
  
  # Extract smooth table safely
  s_table <- summary(fit)$s.table
  
  # Compile results 
  success_df <- data.frame(gene = g, 
                           p_condition = p_table[grepl("condition",rownames(p_table)),"Pr(>|t|)"],
                           p_time = s_table["s(time)", "p-value"],
                           p_time_condition = s_table["s(time,condition)", "p-value"]
  )
  
  #Interpertation of model parameters:
  #condition = baseline differences in gene expression between conditions independent of time
  #s(time) = shared temporal program between conditions over time
  #s(time,condition) = condition-specific temporal deviation
  
  return(success_df)
}

genes <- unique(expr_df$gene)
results <- do.call(rbind, lapply(genes, fit_gene))

# Correct for multiple hypothesis testing
results <- results %>% filter(complete.cases(.))
results$padj_condition <- p.adjust(results$p_condition, method = "BH")
results$padj_time <- p.adjust(results$p_time, method = "BH")
results$padj_time_condition <- p.adjust(results$p_time_condition, method = "BH")

# Classify genes by model parameters
results <- results %>%
  mutate(
    classification = case_when(
      padj_time >= 0.05 & padj_time_condition >= 0.05 ~ "Not dynamic",
      padj_time < 0.05 & padj_time_condition >= 0.05  ~ "Shared dynamics",
      TRUE                                            ~ "Not shared dynamics",
    )
  )

table(results$classification)
#Total significant genes: 1798
#Shared dynamics: 1472
#Not shared dynamics: 326

filter(results, padj_condition < 0.05 | padj_time < 0.05 | padj_time_condition < 0.05)
filter(results, classification == "Not shared dynamics")
filter(results, padj_condition < 0.05)



####################################################
#Identify Shared Temporal Programs
####################################################

# Fit temporal trajectories for all significant genes for each stroma condition

#isolate significant genes
df_sig_genes <- results %>% filter(classification == "Shared dynamics")
rownames(df_sig_genes) <- df_sig_genes$gene
sig_genes <- df_sig_genes$gene

#create a time grid for predictions
time_grid <- seq(min(expr_df$time), max(expr_df$time), length.out = 100)

#initialize an empty prediction matrix for each condition
conds <- c("NHLF", "ESC")

smooth_mat <- list(
  NHLF  = matrix(NA, nrow = length(sig_genes), ncol = length(time_grid)),
  ESC = matrix(NA, nrow = length(sig_genes), ncol = length(time_grid))
)

for (cc in conds) {
  rownames(smooth_mat[[cc]]) <- sig_genes
  colnames(smooth_mat[[cc]]) <- time_grid
}

#get the predicted trajectory per gene and for each stroma condition
for (i in seq_along(sig_genes)) {
  g <- sig_genes[i]
  df_g <- expr_df %>% filter(gene == g)
  
  fit <- try(
    gam(expr ~ condition + s(time, k=5) + s(time, condition, bs="fs", k=5) + s(subject, bs="re"),
        data = df_g, method = "REML"),
    silent = TRUE)
  
  if (inherits(fit, "try-error")) next
  
  for (cc in conds) {
    pred_df <- data.frame(
      time = time_grid,
      condition = cc,
      subject = df_g$subject[1]  # dummy subject
    )
    
    smooth_mat[[cc]][g, ] <- predict(fit, newdata = pred_df, exclude = "s(subject)", type = "response")
    
  }
  
}

#standardize the rows to compare across genes (genes x time grid, row-wise z-scored)
scale_minmax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

smooth_mat_nhlf_scaled <- t(apply(smooth_mat$NHLF, 1, scale_minmax)) #t(scale(t(smooth_mat$NHLF)))
smooth_mat_esc_scaled <- t(apply(smooth_mat$ESC, 1, scale_minmax)) #t(scale(t(smooth_mat$ESC)))


# Plot a heatmap of shared dynamic genes between NHLF and ESC conditions

#perform clustering and assess robustness of clusters
pca_res <- prcomp(smooth_mat_nhlf_scaled, scale. = TRUE)
scores <- pca_res$x[,1:3]
row_hclust <- hclust(dist(scores), method = "ward.D2")

library(fpc)
k <- 6
stab <- clusterboot(
  scores,
  B = 100,
  clustermethod = hclustCBI,
  k = k,
  method = "ward.D2"
)
stab$bootmean

gene_clusters <- cutree(row_hclust, k = k)
gene_clusters <- factor(gene_clusters)

#prepare Complexheatmap of trajectories of significant genes
row_ha <- rowAnnotation(
  Cluster = gene_clusters,
  col = list(Cluster = structure(
    RColorBrewer::brewer.pal(k, "Set1"),
    names = levels(gene_clusters)
  )),
  show_annotation_name = TRUE
)

highlight <- c("AXL", "AKT1", "ATF5", "BRD4", "CALM1", "CBX6", "CCDC85B", 
               "CCND1", "CHD4", "CLDN11", "COL1A1", "CXCL14", "CDH6", "DDX5", 
               "E2F7", "EDN1", "ENG", "EPHB4", "FGD5", "FGF23", "FN1", "FZD4", 
               "ID1", "ITGA5", "JAG2", "KDR", "NOTCH1", "KLF2", "PDGFB", "PGF", 
               "ROCK2", "SMAD5", "SOX17", "SOX11", "STAT5B", "TEAD1", "TCF4", 
               "ZBTB16", "VWF", "ZNF175")

ha = rowAnnotation(foo = anno_mark(at = match(highlight, rownames(smooth_mat_nhlf_scaled)), 
                                   labels = highlight,
                                   labels_gp = gpar(fontsize = 6) ) )

hm_nhlf <- Heatmap(
  smooth_mat_nhlf_scaled,
  cluster_rows = row_hclust,       # use precomputed clustering
  cluster_columns = FALSE,         # keep time order
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_split = k, 
  row_title_gp = gpar(fontsize = 0),
  show_row_dend = FALSE,
  col = pal_atac, 
  top_annotation = NULL,
  left_annotation = row_ha,
  #right_annotation = ha,
  column_title = "NHLF",
  column_gap = unit(0, "mm"),
  #row_gap = unit(0, "mm"),
  heatmap_legend_param = list(title = "Scaled expression", direction = "horizontal"),
  border = TRUE
)

hm_esc <- Heatmap(
  smooth_mat_esc_scaled,
  cluster_rows = row_hclust,       # use precomputed clustering so rows match with NHLF heatmap
  cluster_columns = FALSE,         # keep time order
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_split = k, 
  row_title_gp = gpar(fontsize = 0),
  show_row_dend = FALSE,
  col = pal_atac, 
  top_annotation = NULL,
  #left_annotation = row_ha,
  right_annotation = ha,
  column_title = "ESC",
  column_gap = unit(0, "mm"),
  #row_gap = unit(0, "mm"),
  heatmap_legend_param = list(title = "Scaled expression", direction = "horizontal"),
  border = TRUE
)

pdf(paste0(output_dir, "Joint_GAM_HUVEC-Gag-EIF4E_NHLF_ESC_Pump_shared_dynamic_genes_temporal_heatmap.pdf"), width = 4, height = 6, useDingbats = FALSE)
draw(hm_nhlf + hm_esc, heatmap_legend_side = "bottom", merge_legend = TRUE) 
dev.off()

#export shared dynamic genes with the trajectory cluster identity
df_sig_genes$trajectory_cluster <- gene_clusters

write.table(df_sig_genes, 
            file = paste0(output_dir, "Joint_GAM_HUVEC-Gag-EIF4E_NHLF_ESC_Pump_shared_dynamic_genes_temporal_results.tsv"), 
            sep = "\t", 
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)


# Compute the average trajectory and standard dev across genes in each temporal cluster
time_grid_char <- as.character(time_grid)

cluster_stats <- lapply(levels(gene_clusters), function(cl) {
  rows <- which(gene_clusters == cl)
  mean_traj <- colMeans(smooth_mat_nhlf_scaled[rows, , drop = FALSE])
  se_traj <- apply(smooth_mat_nhlf_scaled[rows, , drop = FALSE], 2, sd) 
  data.frame(
    Cluster = cl,
    Time = time_grid,
    Mean = mean_traj,
    SEM = se_traj
  )
}) %>% bind_rows()

cluster_stats <- cluster_stats %>%
  mutate(ymin = Mean - SEM,
         ymax = Mean + SEM)

pdf(paste0(output_dir, "Joint_GAM_HUVEC-Gag-EIF4E_NHLF_ESC_Pump_shared_dynamic_genes_avg_trajectory_per_cluster.pdf"), width = 1.5, height = 6, useDingbats = FALSE)

ggplot(cluster_stats, aes(x = Time, y = Mean)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = Cluster), alpha = 0.2) +
  geom_line(aes(color = Cluster), size = 1.2) +
  facet_wrap(~ Cluster, ncol = 1) + 
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(title = "Average trajectory",
       x = "Day", y = "Scaled expression") +
  theme(
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none")

dev.off()


# Perform pathway enrichment within each cluster
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

genes_by_module <- split(df_sig_genes$gene, df_sig_genes$trajectory_cluster)

gene_map <- bitr(
  df_sig_genes$gene,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

module_df2 <- merge(df_sig_genes, gene_map,
                    by.x = "gene", by.y = "SYMBOL")

genes_by_module <- split(module_df2$ENTREZID, module_df2$trajectory_cluster)

#GO BP
clusterComp <- compareCluster(
  genes_by_module,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP"
)

dotplot(clusterComp, showCategory = 5, label_format = 100) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10))

#reactome
clusterComp <- compareCluster(
  genes_by_module,
  fun = "enrichPathway",
  organism="human", 
  pvalueCutoff=0.05
)

dotplot(clusterComp, showCategory = 10, label_format = 100) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10))

#KEGG
clusterComp <- compareCluster(
  geneClusters = genes_by_module, 
  fun = "enrichKEGG",
  organism = "hsa", 
  pvalueCutoff = 0.05
)

dotplot(clusterComp, showCategory = 10, label_format = 100) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10))

#hallmark
library(msigdbr)

h_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol) 

genes_by_module <- split(df_sig_genes$gene, df_sig_genes$trajectory_cluster)

clusterComp <- compareCluster(
  geneClusters = genes_by_module, 
  fun = "enricher",
  TERM2GENE = h_t2g,
  pvalueCutoff = 0.1
)

dotplot(clusterComp, showCategory = 10, label_format = 100) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10))








