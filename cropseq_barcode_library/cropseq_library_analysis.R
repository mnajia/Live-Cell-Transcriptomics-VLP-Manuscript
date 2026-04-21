#Project: Live-cell Transcriptomics Manuscript
#Experiment: VLP export of CROPseq barcode libraries
#Author: Mohamad Najia
#Objective: Analyze barcode count distributions from the plasmid library, 293T genomic DNA, and VLP RNA

library(systemPipeR)
library(Biostrings)
library(parallel)
library(LSD)
library(data.table)
library(ggplot2)
library(ggExtra)
library(edgeR)
library(ineq)



####################################################
#Function Declarations
####################################################

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}



####################################################
#Set-up Environment
####################################################

# Declare color pallets
ScatterPlotColorPanel = rev(c("#F8FA0D", "#F6DA23", "#F8BA43","#A5BE6A","#2DB7A3","#1389D2","#0262E0","#343DAE","#352A86"))

# Initialize variables
project_dir <- "/Volumes/mnajia/cropseq_barcode_library/"
output_dir <- paste0(project_dir, "analysis/")
counts_dir <- paste0(project_dir, "barcode_counts/")
barcode_library_fn <- paste0(project_dir, "barcode_library.tsv")
samplesheet <- paste0(project_dir, "samples.txt")

# Import the samplesheet
df.samples <- fread(samplesheet, data.table = FALSE, header = FALSE)
colnames(df.samples) <- c("sample_name", "fastq")

# Import the barcode library
df.lib <- fread(barcode_library_fn, data.table = FALSE)
df.counts <- df.lib[, c("sgRNAID", "Protospacer")]

# Import and compile barcode counts for each sample
for (i in df.samples$sample_name) {
  fn <- paste0(counts_dir, i, "_counts.tsv")
  tmp <- fread(fn, data.table = FALSE)
  sample_name <- colnames(tmp)[3]
  df.counts[,sample_name] <- tmp[,c(3)]
}

# Save the compiled sample counts matrix to disk
write.table(df.counts, 
            file = paste0(output_dir, "CROPseq-24xMS2_barcode_library_counts_matrix.tsv"), 
            quote = FALSE, 
            sep = " \t", 
            row.names = FALSE, 
            col.names = TRUE)



####################################################
#CROPseq Barcode Distribution Analysis
####################################################

# Perform counts-per-million normalization 
cpm_mat <- cpm(df.counts[, df.samples$sample_name])
df.cpms <- cbind(df.counts[,c(1:2)], as.data.frame(cpm_mat))

# Save the CPM-normalized sample counts matrix to disk
write.table(df.cpms, 
            file = paste0(output_dir, "CROPseq-24xMS2_barcode_library_CPM_matrix.tsv"), 
            quote = FALSE, 
            sep = " \t", 
            row.names = FALSE, 
            col.names = TRUE)

# Generate Lorenz curves for each sample
for (i in df.samples$sample_name) {
  
  neworder <- df.cpms[order(df.cpms[[i]]),]
  
  lcolc <- Lc(df.cpms[[i]])
  lcdf <- data.frame(L = rev(1-lcolc$L), p = lcolc$p, Uprob = c(1:length(lcolc$L)/length(lcolc$L)))
  
  gini <- ineq(df.cpms[[i]], type = "Gini")
  auc <- (gini+1)*0.5
  
  gg <- ggplot(lcdf, aes(x = Uprob, y = L)) + 
    geom_line(colour = hcl(h=15, l=65, c=100)) + 
    geom_abline(slope=1, intercept=0, linetype="dotted") + 
    scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1)) + scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
    theme_bw() + 
    ggtitle(paste0(i, "\nLorenz curve of barcode representation")) + 
    xlab("Descending gRNA Abundance Rank") + 
    ylab("Cumulative Fraction of Reads") + 
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black")) +
    geom_text(aes(x=Inf, 
                  y=-Inf,
                  hjust=1, 
                  vjust=0, 
                  label=paste0("Gini AUC = ", round(auc, digits = 4))
                  )
              )
  
  pdf(paste0(output_dir, i, "_lorenz_curve.pdf"), width = 4, height = 4)
  print(gg)
  dev.off()
}


# Plot CPM concordance of barcodes between plasmid and genomic DNA libraries 
df.cpms$density <- get_density(df.cpms$`CROPseq-24xMS2_plasmid_library`, df.cpms$`CROPseq-24xMS2_293T_gDNA_library`, n = 100)

lib_cor <- cor(df.cpms$`CROPseq-24xMS2_plasmid_library`, df.cpms$`CROPseq-24xMS2_293T_gDNA_library`, method = "pearson")
print(lib_cor)
#r = 0.8543432

pdf(paste0(output_dir, "CROPseq-24xMS2_plasmid_v_gDNA_cpm_concordance.pdf"), width = 4, height = 4)

gg <- ggplot(df.cpms) + 
  geom_point(aes(log10(`CROPseq-24xMS2_plasmid_library`), log10(`CROPseq-24xMS2_293T_gDNA_library`), color = density)) + 
  scale_color_gradientn(colours = ScatterPlotColorPanel) +
  scale_x_continuous(limits = c(0, 3.25)) + 
  scale_y_continuous(limits = c(0, 3.25)) + 
  geom_abline(slope=1, intercept=0, linetype="dotted") + 
  ggtitle("CROPseq-24xMS2 barcode library") + 
  xlab("Plasmid Library log10(CPM)") + 
  ylab("293T Genomic DNA log10(CPM)") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none") + 
  geom_text(aes(x=Inf, y=-Inf, hjust=1, vjust=0, label=paste0("r = ", round(lib_cor, digits = 3))))
ggMarginal(gg, type = "histogram", fill = "grey")

dev.off()


# Replicate reproducibility
cor(df.cpms$`CROPseq-24xMS2_293T_gDNA_library`, df.cpms$`CROPseq-24xMS2_293T_Gag_VLPs_rep1`)
cor(df.cpms$`CROPseq-24xMS2_293T_gDNA_library`, df.cpms$`CROPseq-24xMS2_293T_Gag_VLPs_rep2`)
cor(df.cpms$`CROPseq-24xMS2_293T_gDNA_library`, df.cpms$`CROPseq-24xMS2_293T_Gag_VLPs_rep3`)
cor(df.cpms$`CROPseq-24xMS2_293T_Gag_VLPs_rep1`, df.cpms$`CROPseq-24xMS2_293T_Gag_VLPs_rep2`)
cor(df.cpms$`CROPseq-24xMS2_293T_Gag_VLPs_rep1`, df.cpms$`CROPseq-24xMS2_293T_Gag_VLPs_rep3`)

cor(df.cpms$`CROPseq-24xMS2_293T_gDNA_library`, df.cpms$`CROPseq-24xMS2_293T_Gag-MCP_VLPs_rep1`)
cor(df.cpms$`CROPseq-24xMS2_293T_gDNA_library`, df.cpms$`CROPseq-24xMS2_293T_Gag-MCP_VLPs_rep2`)
cor(df.cpms$`CROPseq-24xMS2_293T_gDNA_library`, df.cpms$`CROPseq-24xMS2_293T_Gag-MCP_VLPs_rep3`)
cor(df.cpms$`CROPseq-24xMS2_293T_Gag-MCP_VLPs_rep1`, df.cpms$`CROPseq-24xMS2_293T_Gag-MCP_VLPs_rep2`)
cor(df.cpms$`CROPseq-24xMS2_293T_Gag-MCP_VLPs_rep1`, df.cpms$`CROPseq-24xMS2_293T_Gag-MCP_VLPs_rep3`)


# Plot replicate concordance
df.cpms$density <- get_density(df.cpms$`CROPseq-24xMS2_293T_Gag_VLPs_rep1`, df.cpms$`CROPseq-24xMS2_293T_Gag_VLPs_rep2`, n = 100)

lib_cor <- cor(df.cpms$`CROPseq-24xMS2_293T_Gag_VLPs_rep1`, df.cpms$`CROPseq-24xMS2_293T_Gag_VLPs_rep2`, method = "pearson")
print(lib_cor)
#r = 0.1854755

pdf(paste0(output_dir, "CROPseq-24xMS2_Gag_VLPs_rep1_v_rep2_cpm_concordance.pdf"), width = 4, height = 4)

gg <- ggplot(df.cpms) + 
  geom_point(aes(log10(`CROPseq-24xMS2_293T_Gag_VLPs_rep1`), log10(`CROPseq-24xMS2_293T_Gag_VLPs_rep2`), color = density)) + 
  scale_color_gradientn(colours = ScatterPlotColorPanel) +
  scale_x_continuous(limits = c(0, 3.25)) + 
  scale_y_continuous(limits = c(0, 3.25)) + 
  geom_abline(slope=1, intercept=0, linetype="dotted") + 
  xlab("Barcode Abundance in Rep 1, log10(CPM)") + 
  ylab("Barcode Abundance in Rep 2, log10(CPM)") + 
  ggtitle("MLV Gag\nCROPseq-24xMS2 barcodes") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none") + 
  geom_text(aes(x=Inf, y=-Inf, hjust=1, vjust=0, label=paste0("r = ", round(lib_cor, digits = 3))))
ggMarginal(gg, type = "histogram", fill = "grey")

dev.off()



df.cpms$density <- get_density(df.cpms$`CROPseq-24xMS2_293T_Gag-MCP_VLPs_rep1`, df.cpms$`CROPseq-24xMS2_293T_Gag-MCP_VLPs_rep2`, n = 100)

lib_cor <- cor(df.cpms$`CROPseq-24xMS2_293T_Gag-MCP_VLPs_rep1`, df.cpms$`CROPseq-24xMS2_293T_Gag-MCP_VLPs_rep2`, method = "pearson")
print(lib_cor)
#r = 0.7686647

pdf(paste0(output_dir, "CROPseq-24xMS2_Gag-MCP_VLPs_rep1_v_rep2_cpm_concordance.pdf"), width = 4, height = 4)

gg <- ggplot(df.cpms) + 
  geom_point(aes(log10(`CROPseq-24xMS2_293T_Gag-MCP_VLPs_rep1`), log10(`CROPseq-24xMS2_293T_Gag-MCP_VLPs_rep2`), color = density)) + 
  scale_color_gradientn(colours = ScatterPlotColorPanel) +
  scale_x_continuous(limits = c(0, 3.25)) + 
  scale_y_continuous(limits = c(0, 3.25)) + 
  geom_abline(slope=1, intercept=0, linetype="dotted") + 
  xlab("Barcode Abundance in Rep 1, log10(CPM)") + 
  ylab("Barcode Abundance in Rep 2, log10(CPM)") + 
  ggtitle("MLV Gag-MCP\nCROPseq-24xMS2 barcodes") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none") + 
  geom_text(aes(x=Inf, y=-Inf, hjust=1, vjust=0, label=paste0("r = ", round(lib_cor, digits = 3))))
ggMarginal(gg, type = "histogram", fill = "grey")

dev.off()

# Average CPMs across replicates
df.cpms$CROPseq_24xMS2_Gag_VLP_avg <- rowMeans(df.cpms[,c(5:7)])
df.cpms$CROPseq_24xMS2_Gag_MCP_VLP_avg <- rowMeans(df.cpms[,c(8:10)])

# Plot CPM concordance of barcodes between gDNA and Gag VLP libraries 
df.cpms$density <- get_density(df.cpms$`CROPseq-24xMS2_293T_gDNA_library`, df.cpms$CROPseq_24xMS2_Gag_VLP_avg, n = 100)

lib_cor <- cor(df.cpms$`CROPseq-24xMS2_293T_gDNA_library`, df.cpms$CROPseq_24xMS2_Gag_VLP_avg, method = "pearson")
print(lib_cor)
#r = 0.4773324

pdf(paste0(output_dir, "CROPseq-24xMS2_293T_gDNA_v_Gag_VLPs_cpm_concordance.pdf"), width = 4, height = 4)

gg <- ggplot(df.cpms) + 
  geom_point(aes(log10(`CROPseq-24xMS2_293T_gDNA_library`), log10(`CROPseq_24xMS2_Gag_VLP_avg`), color = density)) + 
  scale_color_gradientn(colours = ScatterPlotColorPanel) +
  scale_x_continuous(limits = c(0, 3.25)) + 
  scale_y_continuous(limits = c(0, 3.25)) + 
  geom_abline(slope=1, intercept=0, linetype="dotted") + 
  xlab("Barcode Abundance in Genomic DNA, log10(CPM)") + 
  ylab("Barcode Abundance in VLPs, log10(CPM)") + 
  ggtitle("MLV Gag\nCROPseq-24xMS2 barcodes") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none") + 
  geom_text(aes(x=Inf, y=-Inf, hjust=1, vjust=0, label=paste0("r = ", round(lib_cor, digits = 3))))
ggMarginal(gg, type = "histogram", fill = "grey")

dev.off()


# Plot CPM concordance of barcodes between gDNA and Gag-MCP VLP libraries 
df.cpms$density <- get_density(df.cpms$`CROPseq-24xMS2_293T_gDNA_library`, df.cpms$CROPseq_24xMS2_Gag_MCP_VLP_avg, n = 100)

lib_cor <- cor(df.cpms$`CROPseq-24xMS2_293T_gDNA_library`, df.cpms$CROPseq_24xMS2_Gag_MCP_VLP_avg, method = "pearson")
print(lib_cor)
#r = 0.8219897

pdf(paste0(output_dir, "CROPseq-24xMS2_293T_gDNA_v_Gag-MCP_VLPs_cpm_concordance.pdf"), width = 4, height = 4)

gg <- ggplot(df.cpms) + 
  geom_point(aes(log10(`CROPseq-24xMS2_293T_gDNA_library`), log10(`CROPseq_24xMS2_Gag_MCP_VLP_avg`), color = density)) + 
  scale_color_gradientn(colours = ScatterPlotColorPanel) +
  scale_x_continuous(limits = c(0, 3.25)) + 
  scale_y_continuous(limits = c(0, 3.25)) + 
  geom_abline(slope=1, intercept=0, linetype="dotted") + 
  xlab("Barcode Abundance in Genomic DNA, log10(CPM)") + 
  ylab("Barcode Abundance in VLPs, log10(CPM)") + 
  ggtitle("MLV Gag-MCP\nCROPseq-24xMS2 barcodes") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none") + 
  geom_text(aes(x=Inf, y=-Inf, hjust=1, vjust=0, label=paste0("r = ", round(lib_cor, digits = 3))))
ggMarginal(gg, type = "histogram", fill = "grey")

dev.off()





