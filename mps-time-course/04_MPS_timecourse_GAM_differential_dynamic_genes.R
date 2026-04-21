#Project: Live-cell Transcriptomics Manuscript
#Experiment: Microphysiological Systems Transcriptional Time-course
#Author: Mohamad Najia
#Objective: identify differentially dynamic genes between stroma conditions
#           by applying joint Generalized Additive Models (GAMs) to VLP
#           cell media samples produced from Gag-EIF4E expressing HUVECs

library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggrepel)
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
project_dir <- "/Volumes/mnajia/2026-01-14_CSR/mps-time-course/"
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
#NHLF and ESC Stroma Condition Temporal DEG Analysis
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
    condition = factor(condition, levels = c("NHLF", "ESC"), ordered = TRUE),
    time = as.numeric(time)
  )
contrasts(expr_df$condition) <- "contr.treatment"

# Fit a joint GAM with random subject effects for each gene to determine if stroma conditions vary differently over time
#Use s(time, by=condition) GAM model to explicitly test if the difference between the two condition curves is non-zero 
fit_gene <- function(g) {
  
  df_g <- expr_df[expr_df$gene == g, ]
  
  failed_df <- data.frame(gene = g, 
                          edf_ref = as.numeric(NA),
                          edf_diff = as.numeric(NA),
                          p_ref = as.numeric(NA),
                          p_diff = as.numeric(NA)
  )
  
  # Safety checks
  if (nrow(df_g) < 5 || 
      sd(df_g$expr) == 0 || 
      length(unique(df_g$time)) < 2) {
    
    return(failed_df)
  }
  
  fit <- try(
    gam(expr ~ condition + s(time, k=5) + s(time, by=condition, k=5) + s(subject, bs="re"),
        data=df_g, method="REML"),
    silent = TRUE
  )
  
  # If fit failed → return NA
  if (inherits(fit, "try-error")) {
    return(failed_df)
  }
  
  # Extract p-value and EDF statistics
  tidy_stats <- broom::tidy(fit)
  p_ref <- tidy_stats$p.value[tidy_stats$term == "s(time)"]
  p_diff <- tidy_stats$p.value[grepl("s\\(time\\):condition", tidy_stats$term)]
  edf_ref <- tidy_stats$edf[tidy_stats$term == "s(time)"]
  edf_diff <- tidy_stats$edf[grepl("s\\(time\\):condition", tidy_stats$term)]
  
  # Compile results 
  success_df <- data.frame(gene = g, 
                           edf_ref = edf_ref,
                           edf_diff = edf_diff,
                           p_ref = p_ref,
                           p_diff = p_diff
  )
  
  #Interpertation of model parameters:
  #p_ref = s(time) = shared temporal component between conditions over time
  #p_diff = s(time,condition) = difference between the ESC trajectory and the NHLF trajectory
  
  return(success_df)
}

genes <- unique(expr_df$gene)
results <- do.call(rbind, lapply(genes, fit_gene))
rownames(results) <- results$gene

# Correct for multiple hypothesis testing
results <- results %>% filter(complete.cases(.))
results$padj_ref <- p.adjust(results$p_ref, method = "BH")
results$padj_diff <- p.adjust(results$p_diff, method = "BH")

# Classify genes by model parameters into temporal categories
results <- results %>%
  mutate(
    classification = case_when(
      padj_ref >= 0.05 & padj_diff >= 0.05 ~ "Not dynamic",
      padj_ref < 0.05 & padj_diff >= 0.05  ~ "Shared dynamics",
      padj_ref >= 0.05 & padj_diff < 0.05  ~ "Condition-specific dynamics",
      padj_ref < 0.05 & padj_diff < 0.05   ~ "Differential dynamics",
      TRUE                                 ~ "Ambiguous"
    )
  )

table(results$classification)
#Total significant genes: 1653
#Shared dynamics: 1421
#Condition-specific dynamics: 184
#Differential dynamics: 48

filter(results, classification == "Shared dynamics")
filter(results, classification == "Condition-specific dynamics")
filter(results, classification == "Differential dynamics")

# Export per-gene model parameter statistics from GAM
#write.table(results, 
#            file = paste0(output_dir, "Joint_GAM_HUVEC-Gag-EIF4E_NHLF_ESC_Pump_diff_model_results.tsv"), 
#            sep = "\t", 
#            quote = FALSE,
#            row.names = FALSE,
#            col.names = TRUE)


# Isolate significant genes
df_sig_genes <- results %>% filter(classification != "Not dynamic")
rownames(df_sig_genes) <- df_sig_genes$gene
sig_genes <- df_sig_genes$gene

# Create a time grid for predictions
time_grid <- seq(min(expr_df$time), max(expr_df$time), length.out = 100)

# Initialize an empty prediction matrix for each condition
conds <- c("NHLF", "ESC")

smooth_mat <- list(
  NHLF  = matrix(NA, nrow = length(sig_genes), ncol = length(time_grid)),
  ESC = matrix(NA, nrow = length(sig_genes), ncol = length(time_grid))
)

for (cc in conds) {
  rownames(smooth_mat[[cc]]) <- sig_genes
  colnames(smooth_mat[[cc]]) <- time_grid
}

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

# Get the predicted temporal trajectories for all significant genes in each stroma condition
for (i in seq_along(sig_genes)) {
  g <- sig_genes[i]
  df_g <- expr_df %>% filter(gene == g)
  
  fit <- gam(expr ~ condition + s(time, k=5) + s(time, condition, bs="fs", k=5) + s(subject, bs="re"),
             data = df_g, method = "REML")
  
  for (cc in conds) {
    pred_df <- data.frame(
      time = time_grid,
      condition = cc, #factor(cc, levels = c("NHLF", "ESC"), ordered = TRUE),
      subject = df_g$subject[1]  # dummy subject
    )
    #contrasts(pred_df$condition) <- "contr.treatment"
    
    smooth_mat[[cc]][g, ] <- predict(fit, newdata = pred_df, exclude = "s(subject)", type = "response")
    
  }
  
}

# Standardize the rows to compare across genes (genes x time grid, row-wise z-scored)
scale_minmax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

smooth_mat_nhlf_scaled <- t(apply(smooth_mat$NHLF, 1, scale_minmax)) #t(scale(t(smooth_mat$NHLF)))
smooth_mat_esc_scaled <- t(apply(smooth_mat$ESC, 1, scale_minmax)) #t(scale(t(smooth_mat$ESC)))

# Plot the regulatory landscape of condition-specific/differentially dynamic genes
df_sig_cond_genes <- results %>% filter(classification %in% c("Condition-specific dynamics", "Differential dynamics") ) 
sig_cond_genes <- df_sig_cond_genes %>% pull(gene)
pred_nhlf <- smooth_mat$NHLF[sig_cond_genes,]
pred_esc <- smooth_mat$ESC[sig_cond_genes,]
diff_vector <- pred_nhlf - pred_esc

df_sig_cond_genes$peak_time <- as.numeric(NA)
df_sig_cond_genes$peak_val <- as.numeric(NA)

for (i in 1:nrow(df_sig_cond_genes)) {
  peak_idx <- which.max(abs(diff_vector[i,]))
  df_sig_cond_genes[i,"peak_time"] <- as.numeric(names(peak_idx))
  df_sig_cond_genes[i,"peak_val"]  <- diff_vector[i,peak_idx]
}

filter(df_sig_cond_genes, peak_time >3 & peak_time <6)

label <- c("SUZ12", "ROCK1", "RNF181", "DBNDD2", "UBXN6", "KCNJ2",
           "DUSP3", "ATOH8", "FOXO3B", "VTI1B", "BACE2", "COX8A",
           "MMP2", "RNF219", "ITGA10", "BAHD1", "PLVAP",
           "S100A4", "IGF2", "ZEB1", "MIF", "LGALS1")

df_sig_cond_genes$highlight <- df_sig_cond_genes$gene %in% label
df_sig_cond_genes$label <- ""
df_sig_cond_genes[df_sig_cond_genes$highlight, "label"] <- df_sig_cond_genes[df_sig_cond_genes$highlight, "gene"]

df_sig_cond_genes$response <- ""
df_sig_cond_genes[df_sig_cond_genes$peak_time <= 2, "response"] <- "Early"
df_sig_cond_genes[df_sig_cond_genes$peak_time >= 6, "response"] <- "Late"
df_sig_cond_genes[df_sig_cond_genes$peak_time > 2 & df_sig_cond_genes$peak_time < 6, "response"] <- "Intermediate"

gg <- ggplot(df_sig_cond_genes, aes(x=peak_time, y=peak_val, color=response)) + 
  geom_point() +
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  scale_x_continuous(breaks = c(1:7)) + 
  scale_y_continuous(limits = c(-5, 5)) + 
  xlab("Time of peak difference (Day)") + 
  ylab("Peak expression difference (NHLF - ESC)") + 
  ggtitle("Temporal response of differentially dynamic genes") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_text_repel(data = df_sig_cond_genes,
                  aes(label = label), 
                  segment.color = "black",
                  colour = "black",
                  size = 4,
                  box.padding = 1,
                  max.overlaps = 200,
                  show.legend = FALSE)

pdf(paste0(output_dir, "Joint_GAM_HUVEC-Gag-EIF4E_NHLF_ESC_Pump_differential_dynamic_genes_temporal_response.pdf"), width = 6, height = 4.5, useDingbats = FALSE)
print(gg)
dev.off()


# Identify temporal modules of differential dynamic genes
smooth_mat_nhlf_scaled_cond <- smooth_mat_nhlf_scaled[sig_cond_genes,]
smooth_mat_esc_scaled_cond <- smooth_mat_esc_scaled[sig_cond_genes,]

z_nhlf <- t(scale(t(smooth_mat$NHLF[sig_cond_genes,])))
z_esc <- t(scale(t(smooth_mat$ESC[sig_cond_genes,])))
diff_curves <- z_nhlf - z_esc

corr_mat_concat <- cor(t(diff_curves), method = "pearson")

dist_mat <- as.dist(1 - corr_mat_concat)
row_hclust <- hclust(dist_mat, method = "ward.D2")

sil_width <- sapply(2:10, function(k) {
  cl <- cutree(row_hclust, k = k)
  mean(silhouette(cl, dist_mat)[, 3])
})

plot(2:10, sil_width, type = "b")

k<-4 

#cluster stability under resampling
library(fpc)
stab <- clusterboot(
  corr_mat_concat,
  B = 100,
  clustermethod = hclustCBI,
  k = k,
  method = "ward.D2"
)
stab$bootmean

modules <- cutree(row_hclust, k = k)

module_df <- data.frame(
  gene = names(modules),
  module = factor(modules)
)

row_ha <- rowAnnotation(
  Modules = module_df$module,
  col = list(Cluster = structure(
    RColorBrewer::brewer.pal(k, "Set3"),
    names = levels(module_df$module)
  )),
  show_annotation_name = TRUE
)

highlight <- c("ATOH8", "CALD1", "COL7A1", "COX8A", "DBNDD2", "EN1", "EPHB3", "FOXO3B",
               "GPX7", "IFI27L2", "IFITM3", "IGF2", "ITGA10", "LGALS1", "LGR4", "MIF",
               "MMP2", "PFN1", "ROCK1", "RNF219", "S100A11", "S100A4", "S100A6", "SH3BGRL3",
               "SLC17A8", "SUZ12", "TMSB4X", "VIM", "ZEB1", "VTI1B", "ARPC3", "COX6B1",
               "ATP5E", "ATP5L", "ATP5G1", "PLVAP")

ha = rowAnnotation(foo = anno_mark(at = match(highlight, rownames(corr_mat_concat)), 
                                   labels = highlight,
                                   side = "left",
                                   labels_gp = gpar(fontsize = 6, just = "left")) )

hm <- Heatmap(
  corr_mat_concat,
  col = pal_rna, 
  cluster_rows = row_hclust,
  cluster_columns = row_hclust,
  row_title_gp = gpar(fontsize = 0),
  column_title_gp = gpar(fontsize = 0),
  row_split = k, 
  column_split = k,
  row_gap = unit(1, "mm"),
  column_gap = unit(1, "mm"),
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  left_annotation = ha, #row_ha,
  right_annotation = row_ha,
  heatmap_legend_param = list(title = "Pearson", direction = "horizontal"),
  border = TRUE
)

pdf(paste0(output_dir, "Joint_GAM_HUVEC-Gag-EIF4E_NHLF_ESC_Pump_differential_dynamic_genes_pearson_cor_heatmap.pdf"), width = 5.25, height = 5, useDingbats = FALSE)
draw(hm, heatmap_legend_side = "bottom", merge_legend = TRUE) 
dev.off()

# Compute the average trajectory of genes per module 
library(purrr)

z_nhlf <- smooth_mat_nhlf_scaled_cond 
z_esc <- smooth_mat_esc_scaled_cond 

traj_summary <- module_df %>%
  group_by(module) %>%
  summarise(
    esc_mat  = list(z_esc[gene, , drop = FALSE]),
    nhlf_mat = list(z_nhlf[gene, , drop = FALSE]),
    n_genes  = n(),
    .groups = "drop"
  ) %>%
  mutate(
    esc_mean  = map(esc_mat, colMeans),
    nhlf_mean = map(nhlf_mat, colMeans),
    esc_sem   = map2(esc_mat, n_genes, ~ apply(.x, 2, sd)/ sqrt(.y) ),
    nhlf_sem  = map2(nhlf_mat, n_genes, ~ apply(.x, 2, sd)/ sqrt(.y) )
  )

plot_df <- traj_summary %>%
  dplyr::select(module, esc_mean, nhlf_mean, esc_sem, nhlf_sem) %>%
  tidyr::pivot_longer(
    cols = -module,
    names_to = c("condition", ".value"),
    names_pattern = "(esc|nhlf)_(mean|sem)"
  ) %>%
  dplyr::mutate(
    time = list(time_grid)
  ) %>%
  tidyr::unnest(c(mean, sem, time))

module_sizes <- table(module_df$module) %>% as.data.frame()
colnames(module_sizes) <- c("module", "n")
plot_df$condition <- factor(plot_df$condition, levels = c("nhlf", "esc"))

gg <- ggplot(plot_df, aes(x = time, y = mean, color = condition, fill = condition)) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.25, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(
    ~ module,
    labeller = labeller(
      module = function(m)
        paste0("Module ", m, " (n=", module_sizes$n[module_sizes$module == m], ")")
    ), ncol = 1
  ) + 
  scale_color_manual(values = c(nhlf = "#E41A1C", esc = "#377EB8"), labels = c("NHLF", "ESC")) + 
  scale_fill_manual(values = c(nhlf = "#E41A1C", esc = "#377EB8"), labels = c("NHLF", "ESC")) +
  xlab("Day") + 
  ylab("Scaled Expression") + 
  labs(color = "Stroma Condition", fill = "Stroma Condition") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

pdf(paste0(output_dir, "Joint_GAM_HUVEC-Gag-EIF4E_NHLF_ESC_Pump_differential_dynamic_genes_avg_trajectory_per_module.pdf"), width = 3.5, height = 6, useDingbats = FALSE)
print(gg)
dev.off()


# Plot gene trajectory heatmaps for each differential time module 
module_df <- module_df[
  match(rownames(z_nhlf), module_df$gene),  #smooth_mat_nhlf_scaled_cond
]

module_colors <- structure(
  RColorBrewer::brewer.pal(length(levels(module_df$module)), "Set3"),
  names = levels(module_df$module)
)

row_ha <- rowAnnotation(
  Module = module_df$module,
  col = list(Module = module_colors),
  show_annotation_name = FALSE
)

ht <- draw(hm, heatmap_legend_side = "bottom", merge_legend = TRUE) 
row_order_vec <- unlist(row_order(ht))

ht_nhlf <- Heatmap(
  z_nhlf,
  col = pal_atac,
  cluster_rows = FALSE, 
  cluster_columns = FALSE,
  row_order = row_order_vec,
  row_split = module_df$module,
  show_row_names = FALSE,
  show_column_names = FALSE,
  left_annotation = row_ha,
  column_title = "NHLF",
  heatmap_legend_param = list(title = "Scaled expression", direction = "horizontal"),
  border = TRUE
)

ht_esc <- Heatmap(
  z_esc, 
  col = pal_atac,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = row_order_vec,
  row_split = module_df$module,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_title = "ESC",
  heatmap_legend_param = list(title = "Scaled expression", direction = "horizontal"),
  border = TRUE
)

pdf(paste0(output_dir, "Joint_GAM_HUVEC-Gag-EIF4E_NHLF_ESC_Pump_differential_dynamic_genes_temporal_heatmap.pdf"), width = 4, height = 5, useDingbats = FALSE)
draw(ht_nhlf + ht_esc, heatmap_legend_side = "bottom", annotation_legend_side = "right")
dev.off()


# Export GAM results for significant genes
tmp <- df_sig_cond_genes
tmp$module <- module_df[rownames(tmp), "module"]

write.table(tmp, 
            file = paste0(output_dir, "Joint_GAM_HUVEC-Gag-EIF4E_NHLF_ESC_Pump_diff_model_results.tsv"), 
            sep = "\t", 
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)






# Identify enriched pathways in each temporal module
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

genes_by_module <- split(module_df$gene, module_df$module)

gene_map <- bitr(
  module_df$gene,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

module_df2 <- merge(module_df, gene_map,
                    by.x = "gene", by.y = "SYMBOL")

genes_by_module <- split(module_df2$ENTREZID, module_df2$module)

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

genes_by_module <- split(module_df$gene, module_df$module)

clusterComp <- compareCluster(
  geneClusters = genes_by_module, 
  fun = "enricher",
  TERM2GENE = h_t2g,
  pvalueCutoff = 0.1
)

dotplot(clusterComp, showCategory = 10, label_format = 100) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 10))









