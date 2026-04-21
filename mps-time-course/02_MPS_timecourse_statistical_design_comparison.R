#Project: Live-cell Transcriptomics Manuscript
#Experiment: Microphysiological Systems Transcriptional Time-course
#Author: Mohamad Najia
#Objective: compare independent versus repeated-measures statistical models for identifying temporally dynamic genes

library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(scales)
library(stringr)
library(tximport)
library(DESeq2)
library(mgcv)
library(circlize)



####################################################
#Function Declarations
####################################################

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
#Repeated-measures v independent measures comparsion
####################################################

# Focus on HUVEC Gag-EIF4E samples
df.samples <- df.samples %>% filter(HUVEC_condition == "HUVEC-Gag-EIF4E")
txi <- subset_txi(txi, df.samples$sample_name)

# Create DESeq2 object for vst
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
#Note: GAMs assume homoskedastic, approximately Gaussian errors...RNA-seq counts violate this, but VST-normalized values satisfy the assumptions reasonably well
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

# Fit a repeated-measures or independent-measures GAM for each gene
genes <- unique(expr_df$gene)

fit_gene <- function(g) {
  
  df_g <- expr_df[expr_df$gene == g, ]
  
  # Safety checks
  if (nrow(df_g) < 5 || 
      sd(df_g$expr) == 0 || 
      length(unique(df_g$time)) < 2) {
    
    rep_resVar <- as.numeric(NA)
    ind_resVar <- as.numeric(NA)
    p_time_rep <- as.numeric(NA)
    p_time_condition_rep <- as.numeric(NA)
    p_time_ind <- as.numeric(NA)
    p_time_condition_ind <- as.numeric(NA)
    
  } else {
    
    # Fit repeated-measures model accounting for subject-level random effects variation
    fitRep <- try(
      #gam(expr ~ s(time, k=5) + s(subject, bs="re"), data=df_g, method="REML"),
      gam(expr ~ condition + s(time, k=5) + s(time, condition, bs="fs", k=5) + s(subject, bs="re"),
          data=df_g, method="REML"),
      silent = TRUE
    )
    
    # If fit failed → return NA
    if (inherits(fitRep, "try-error")) {
      rep_resVar <- as.numeric(NA)
      p_time_rep <- as.numeric(NA)
      p_time_condition_rep <- as.numeric(NA)
    } else {
      rep_resVar <- fitRep$sig2
      
      # Extract smooth table and p-values
      s_table <- summary(fitRep)$s.table
      p_time_rep <- s_table["s(time)", "p-value"]
      p_time_condition_rep <- s_table["s(time,condition)", "p-value"]
    }
    
    # Fit independent-measures model
    fitInd <- try(
      #gam(expr ~ s(time, k=5), data=df_g, method="REML"),
      gam(expr ~ condition + s(time, k=5) + s(time, condition, bs="fs", k=5),
          data=df_g, method="REML"),
      silent = TRUE
    )
    
    # If fit failed → return NA
    if (inherits(fitInd, "try-error")) {
      ind_resVar <- as.numeric(NA)
      p_time_ind <- as.numeric(NA)
      p_time_condition_ind <- as.numeric(NA)
    } else {
      ind_resVar <- fitInd$sig2
      
      # Extract smooth table and p-values
      s_table <- summary(fitInd)$s.table
      p_time_ind <- s_table["s(time)", "p-value"]
      p_time_condition_ind <- s_table["s(time,condition)", "p-value"]
    }
    
  }
  
  # Create data frame to output residual variance and p-values from both models
  df <- data.frame(gene = g, 
                   residual_var_repeated = rep_resVar,
                   residual_var_independent = ind_resVar,
                   p_time_rep = p_time_rep,
                   p_time_condition_rep = p_time_condition_rep,
                   p_time_ind = p_time_ind,
                   p_time_condition_ind = p_time_condition_ind)
  
  return(df)
}

results <- do.call(rbind, lapply(genes, fit_gene))
rownames(results) <- results$gene

# Correct for multiple hypothesis testing
results <- results %>% filter(complete.cases(.))
results$padj_time_rep <- p.adjust(results$p_time_rep, method = "BH")
results$padj_time_condition_rep <- p.adjust(results$p_time_condition_rep, method = "BH")
results$padj_time_ind <- p.adjust(results$p_time_ind, method = "BH")
results$padj_time_condition_ind <- p.adjust(results$p_time_condition_ind, method = "BH")

alpha = 0.05
table(results$padj_time_rep < alpha & results$padj_time_condition_rep >= alpha)
table(results$padj_time_ind < alpha & results$padj_time_condition_ind >= alpha)

#significant shared dynamic genes following multiple hypothesis correction of p-values:
#repeated-measures: 1472
#independent-measures: 1363

setdiff(results[results$padj_time_rep < alpha & results$padj_time_condition_rep >= alpha,"gene"], 
        results[results$padj_time_ind < alpha & results$padj_time_condition_ind >= alpha,"gene"])


# Plot residual variance from each model
results$significant <- results$padj_time_rep < alpha & results$padj_time_condition_rep >= alpha

df2 <- results %>%
  mutate(on_top = significant == TRUE) %>%
  arrange(on_top)

pdf(file = paste0(output_dir, "GAM_Pump_HUVEC-Gag-EIF4E_NHLF_residual_var_comparison.pdf"), width = 5, height = 4, useDingbats = FALSE)

ggplot(df2, aes(x=residual_var_independent, y=residual_var_repeated, color=significant, alpha=significant)) + 
  geom_point() + 
  scale_color_manual(values=c("#999999", "#E69F00")) + 
  scale_alpha_manual(values=c(0.2,1)) + 
  scale_size_manual(values=c(0.75,1)) + 
  geom_abline(slope=1, intercept=0, linetype="dashed", color = "black") + 
  scale_x_continuous(limits = c(0, 3)) + 
  scale_y_continuous(limits = c(0, 3)) + 
  xlab("Residual Variance: Independent Measures Model") + 
  ylab("Residual Variance: Repeated-Measures Model") + 
  ggtitle("Repeated-measures experiments reduce residual variance") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

dev.off()

# Export results
write.table(results, 
            file = paste0(output_dir, "GAM_Pump_HUVEC-Gag-EIF4E_NHLF_residual_var_comparison.tsv"), 
            sep = "\t", 
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)




