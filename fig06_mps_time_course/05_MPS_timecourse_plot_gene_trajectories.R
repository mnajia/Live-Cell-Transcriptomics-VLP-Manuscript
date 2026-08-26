#Project: Live-cell Transcriptomics Manuscript
#Experiment: Microphysiological Systems Transcriptional Time-course
#Author: Mohamad Najia
#Objective: plot temporal trajectory of genes identified from GAMs

library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(scales)
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
df.samples <- df.samples %>% filter(HUVEC_condition == "HUVEC-Gag-EIF4E")
txi <- subset_txi(txi, df.samples$sample_name)



####################################################
#Temporal DEG Analysis 
####################################################

# Create DESeq2 object for vst transformation
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

# Run vst
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

# Plot variance-stabilized expression trajectories over time
time_grid <- seq(1,7, length.out = 100)

plot_gene_trajectory <- function(expr_df, g, time_grid, ymin, ymax) {
  
  # Get the desired gene
  df_g <- expr_df[expr_df$gene == g, ]
  
  # Fit Joint GAM
  fit <- gam(expr ~ condition + s(time, k=5) + s(time, condition, bs="fs", k=5) + s(subject, bs="re"), data=df_g, method="REML")
  
  # Predict the temporal trajectory of the desired gene from the GAM
  grid.nhlf <- data.frame(
    time = time_grid,
    condition = "NHLF",
    subject = df_g$subject[1]  # dummy subject
  )
  
  pred.nhlf <- predict(fit, newdata = grid.nhlf, exclude = "s(subject)", type = "response", se.fit = TRUE)
  grid.nhlf$fit <- pred.nhlf$fit
  grid.nhlf$se  <- pred.nhlf$se.fit
  
  grid.esc <- data.frame(
    time = time_grid,
    condition = "ESC",
    subject = df_g$subject[1]  # dummy subject
  )
  
  pred.esc <- predict(fit, newdata = grid.esc, exclude = "s(subject)", type = "response", se.fit = TRUE)
  grid.esc$fit <- pred.esc$fit
  grid.esc$se  <- pred.esc$se.fit
  
  grid_all <- rbind(grid.nhlf, grid.esc)
  grid_all$condition <- factor(grid_all$condition, levels = c("NHLF", "ESC"))
  
  # Plot trajectory for both conditions
  ggplot() +
    geom_ribbon(data = grid_all, aes(x = time, y = fit, ymin = fit - (1.96 * se), ymax = fit + (1.96 * se), fill = condition), 
                alpha = 0.15, color = NA) +
    geom_line(data = grid_all, aes(x = time, y = fit, color = condition), size = 1.25) +
    geom_point(data = df_g, aes(x = time, y = expr, color = condition), alpha = 0.6, size = 2) +
    geom_line(data = df_g, aes(x = time, y = expr, color = condition, group = subject), 
              alpha = 0.3, size = 0.5, linetype = "dashed") +
    scale_color_manual(values = c("NHLF" = "#E41A1C", "ESC" = "#377EB8")) +
    scale_fill_manual(values = c("NHLF" = "#E41A1C", "ESC" = "#377EB8")) +
    scale_y_continuous(limits = c(ymin, ymax)) + 
    scale_x_continuous(breaks = 1:7) +
    labs(
      title = g,
      x = "Day",
      y = "VST-normalized expression",
      color = "Condition",
      fill = "Condition"
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(color = "black"),
      axis.text  = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5)
    )
  
}


# Genes with differential dynamics
gene <- "TMSB4X"
pdf(paste0(output_dir, "HUVEC-Gag-EIF4E_NHLF_ESC_Pump_temporal_trajectory_", gene, ".pdf"), width = 4, height = 3.2, useDingbats = FALSE)
plot_gene_trajectory(expr_df, gene, time_grid, 12, 16)
dev.off()

gene <- "VIM"
pdf(paste0(output_dir, "HUVEC-Gag-EIF4E_NHLF_ESC_Pump_temporal_trajectory_", gene, ".pdf"), width = 4, height = 3.2, useDingbats = FALSE)
plot_gene_trajectory(expr_df, gene, time_grid, 12, 16)
dev.off()

gene <- "IGF2"
pdf(paste0(output_dir, "HUVEC-Gag-EIF4E_NHLF_ESC_Pump_temporal_trajectory_", gene, ".pdf"), width = 4, height = 3.2, useDingbats = FALSE)
plot_gene_trajectory(expr_df, gene, time_grid, 6, 12)
dev.off()

gene <- "PLVAP"
pdf(paste0(output_dir, "HUVEC-Gag-EIF4E_NHLF_ESC_Pump_temporal_trajectory_", gene, ".pdf"), width = 4, height = 3.2, useDingbats = FALSE)
plot_gene_trajectory(expr_df, gene, time_grid, 6, 9)
dev.off()

gene <- "ZEB1"
pdf(paste0(output_dir, "HUVEC-Gag-EIF4E_NHLF_ESC_Pump_temporal_trajectory_", gene, ".pdf"), width = 4, height = 3.2, useDingbats = FALSE)
plot_gene_trajectory(expr_df, gene, time_grid, 6, 16)
dev.off()

gene <- "MIF"
pdf(paste0(output_dir, "HUVEC-Gag-EIF4E_NHLF_ESC_Pump_temporal_trajectory_", gene, ".pdf"), width = 4, height = 3.2, useDingbats = FALSE)
plot_gene_trajectory(expr_df, gene, time_grid, 8, 15)
dev.off()

gene <- "CALD1"
pdf(paste0(output_dir, "HUVEC-Gag-EIF4E_NHLF_ESC_Pump_temporal_trajectory_", gene, ".pdf"), width = 4, height = 3.2, useDingbats = FALSE)
plot_gene_trajectory(expr_df, gene, time_grid, 5, 16)
dev.off()

gene <- "ITGA10"
pdf(paste0(output_dir, "HUVEC-Gag-EIF4E_NHLF_ESC_Pump_temporal_trajectory_", gene, ".pdf"), width = 4, height = 3.2, useDingbats = FALSE)
plot_gene_trajectory(expr_df, gene, time_grid, 5, 10)
dev.off()

gene <- "ATP5E"
pdf(paste0(output_dir, "HUVEC-Gag-EIF4E_NHLF_ESC_Pump_temporal_trajectory_", gene, ".pdf"), width = 4, height = 3.2, useDingbats = FALSE)
plot_gene_trajectory(expr_df, gene, time_grid, 8, 13)
dev.off()

gene <- "COX6B1"
pdf(paste0(output_dir, "HUVEC-Gag-EIF4E_NHLF_ESC_Pump_temporal_trajectory_", gene, ".pdf"), width = 4, height = 3.2, useDingbats = FALSE)
plot_gene_trajectory(expr_df, gene, time_grid, 4, 13)
dev.off()


# Genes with shared dynamics 
gene <- "CCDC85B"
pdf(paste0(output_dir, "HUVEC-Gag-EIF4E_NHLF_ESC_Pump_temporal_trajectory_", gene, ".pdf"), width = 4, height = 3.2, useDingbats = FALSE)
plot_gene_trajectory(expr_df, gene, time_grid, 8, 15)
dev.off()

gene <- "EDN1"
pdf(paste0(output_dir, "HUVEC-Gag-EIF4E_NHLF_ESC_Pump_temporal_trajectory_", gene, ".pdf"), width = 4, height = 3.2, useDingbats = FALSE)
plot_gene_trajectory(expr_df, gene, time_grid, 4, 14)
dev.off()

gene <- "CCND1"
pdf(paste0(output_dir, "HUVEC-Gag-EIF4E_NHLF_ESC_Pump_temporal_trajectory_", gene, ".pdf"), width = 4, height = 3.2, useDingbats = FALSE)
plot_gene_trajectory(expr_df, gene, time_grid, 10, 18)
dev.off()

gene <- "TCF4"
pdf(paste0(output_dir, "HUVEC-Gag-EIF4E_NHLF_ESC_Pump_temporal_trajectory_", gene, ".pdf"), width = 4, height = 3.2, useDingbats = FALSE)
plot_gene_trajectory(expr_df, gene, time_grid, 6, 14)
dev.off()

gene <- "FGF23"
pdf(paste0(output_dir, "HUVEC-Gag-EIF4E_NHLF_ESC_Pump_temporal_trajectory_", gene, ".pdf"), width = 4, height = 3.2, useDingbats = FALSE)
plot_gene_trajectory(expr_df, gene, time_grid, 6, 10)
dev.off()

gene <- "FZD4"
pdf(paste0(output_dir, "HUVEC-Gag-EIF4E_NHLF_ESC_Pump_temporal_trajectory_", gene, ".pdf"), width = 4, height = 3.2, useDingbats = FALSE)
plot_gene_trajectory(expr_df, gene, time_grid, 6, 14)
dev.off()







