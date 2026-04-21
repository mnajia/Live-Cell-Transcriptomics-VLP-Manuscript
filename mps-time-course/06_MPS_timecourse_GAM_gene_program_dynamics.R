#Project: Live-cell Transcriptomics Manuscript
#Experiment: Microphysiological Systems Transcriptional Time-course
#Author: Mohamad Najia
#Objective: investigate gene programs over time for each stroma condition
#           using single-sample GSEA on VLP cell media samples produced 
#           from Gag-EIF4E expressing HUVECs

library(matrixStats)
library(circlize)
library(data.table)
library(tximport)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggExtra)
library(scales)
library(lme4)
library(lmerTest)



####################################################
#Function Declarations
####################################################

scaleRows <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
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

#Source: https://rpubs.com/pranali018/SSGSEA
#Analysis inspired by: https://www.biorxiv.org/content/10.1101/2025.11.13.687739v1.full#sec-10
ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  
  # Ranks for genes
  R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
  
  # Calculate enrichment score (es) for each sample (column)
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)
    
    # Calc es for each gene set
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos
      
      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
      
      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
      
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      
      # Normalize by gene number
      if (scale) step_cdf_diff = step_cdf_diff / num_genes
      
      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })
  
  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
  
  # Normalize by absolute diff between max and min
  if (norm) es = es / diff(range(es))
  
  # Prepare output
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}

### Helper Functions ###
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])) )
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
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
df.samples <- df.samples %>% filter(HUVEC_condition == "HUVEC-Gag-EIF4E")
txi <- subset_txi(txi, df.samples$sample_name)



####################################################
# Perform ssGSEA on gene sets over time
####################################################

# log2 transform the kallisto TPM matrix
#df.tpm <- txi$abundance
#df.tpm <- log2(df.tpm + 1)
#data <- as.matrix(df.tpm)

# Generate vst-transformed counts
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

vsd <- vst(dds, blind = TRUE)
data <- as.matrix(assay(vsd))


# Import Hallmark and vascular gene sets motivated by condition-specific temporal genes
library(msigdbr)
library(mgcv)

h_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

gs_name <- list.files(path = paste0(project_dir, "gene_sets/"), pattern = "\\.gmt$", full.names = FALSE)

res <- lapply(gs_name, function(x) {
  fn <- paste0(project_dir, "gene_sets/", x)
  gene_set <- clusterProfiler::read.gmt(fn)
  colnames(gene_set) <- c("gs_name", "gene_symbol")
  return(gene_set)
})

gs_t2g <- rbind(h_t2g, do.call(rbind, res))

# Perform ssGSEA on gene sets over each timepoint
df_gs <- data.frame(gs = unique(gs_t2g$gs_name),
                    p_time = as.numeric(NA),
                    p_time_condition = as.numeric(NA))
rownames(df_gs) <- df_gs$gs

time_grid <- seq(1,7, length.out = 100)

for (i in unique(gs_t2g$gs_name)) {
  
  #get hallmark gene set and format for ssGSEA
  gene_set <- filter(gs_t2g, gs_name == i)
  gene_set$gs_name <- NULL
  colnames(gene_set) <- c("gene_set")
  gene_sets = as.list(as.data.frame(gene_set))
  
  #perform ssGSEA on current gene set
  res <- ssgsea(data, gene_sets, scale = TRUE, norm = FALSE)
  
  #format ssGSEA result
  df.plot <- t(res) %>% as.data.frame()
  colnames(df.plot) <- "score"
  tmp <- str_split_fixed(rownames(df.plot), pattern = "_", n=5)
  df.plot$timepoint <- str_split_fixed( tmp[,4], pattern = "Day0", n=2)[,2] %>% as.numeric()
  df.plot$sampleID <- factor(tmp[,1])
  df.plot$condition <- factor(tmp[,3], levels = c("NHLF", "ESC"))
  
  #fit GAM to assess if gene set varies over time and between conditions
  fit_gsea <- gam(score ~ condition + s(timepoint, k=5) + s(timepoint, condition, bs="fs", k=5) + s(sampleID, bs="re"), data = df.plot)
  s_table <- summary(fit_gsea)$s.table
  
  df_gs[i,"p_time"] <- s_table["s(timepoint)","p-value"]
  df_gs[i,"p_time_condition"] <- s_table["s(timepoint,condition)","p-value"]
  
  #test with a linear mixed model if the conditions are different
  fit_lmm <- lmer(score ~ condition * timepoint + (1 | sampleID), data = df.plot)
  summary(fit_lmm)
  
  #predict the temporal trajectory of the gene set from the GAM
  grid.nhlf <- data.frame(
    timepoint = time_grid,
    condition = "NHLF",
    sampleID = df.plot$sampleID[1]  # dummy subject
  )
  
  pred.nhlf <- predict(fit_gsea, newdata = grid.nhlf, exclude = "s(sampleID)", type = "response", se.fit = TRUE)
  grid.nhlf$fit <- pred.nhlf$fit
  grid.nhlf$se  <- pred.nhlf$se.fit
  
  grid.esc <- data.frame(
    timepoint = time_grid,
    condition = "ESC",
    sampleID = df.plot$sampleID[1]  # dummy subject
  )
  
  pred.esc <- predict(fit_gsea, newdata = grid.esc, exclude = "s(sampleID)", type = "response", se.fit = TRUE)
  grid.esc$fit <- pred.esc$fit
  grid.esc$se  <- pred.esc$se.fit
  
  grid_all <- rbind(grid.nhlf, grid.esc)
  grid_all$condition <- factor(grid_all$condition, levels = c("NHLF", "ESC"))
  
  gg <- ggplot(grid_all, aes(x = timepoint, y = fit, color = condition, fill = condition)) +
    geom_ribbon(aes(ymin = fit - (1.96 * se), ymax = fit + (1.96 * se)), alpha = 0.2, color = NA) + #95% confidence bands
    geom_line(size = 1.2) +
    scale_color_manual(values = c("NHLF" = "#E41A1C", "ESC" = "#377EB8")) +
    scale_fill_manual(values = c("NHLF" = "#E41A1C", "ESC" = "#377EB8")) +
    xlab("Day") + 
    ylab("ssGSEA score") + 
    ggtitle(paste0("HUVEC Gag-EIF4E\n",i)) +
    theme_bw() +
    theme(#legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(color = "black"),
      axis.text  = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5))
  
  pdf(paste0(output_dir, "ssGSEA_HUVEC-Gag-EIF4E_NHLF_ESC_Pump_GAM_fit_", i, ".pdf"), width = 4, height = 3.5, useDingbats = FALSE)
  print(gg)
  dev.off()
  
}

df_gs$padj_time <- p.adjust(df_gs$p_time, method = "BH")
df_gs$padj_time_condition <- p.adjust(df_gs$p_time_condition, method = "BH")

filter(df_gs, padj_time < 0.05)
filter(df_gs, padj_time_condition < 0.1)

# Export results
write.table(df_gs, 
            file = paste0(output_dir, "ssGSEA_HUVEC-Gag-EIF4E_NHLF_ESC_Pump_results.tsv"), 
            sep = "\t", 
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)



