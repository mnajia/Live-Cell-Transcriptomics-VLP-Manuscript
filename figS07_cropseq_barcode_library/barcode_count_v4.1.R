#Experiment: CROPseq Barcode Library
#Objective: Quantify barcode counts from NGS reads
#Author: Mohamad Najia
#barcode_count_v4.1.R

library(ShortRead)
library(Biostrings)
library(stringr)
library(data.table)


# Parse input arguments 
args <- commandArgs(TRUE)
sample_name <- args[1] 
fq_file <- args[2]
barcode_lib_fn <- args[3]
output_dir <- args[4]

# Import the barcode library
df.lib <- fread(barcode_lib_fn, data.table = FALSE)
df.lib$Protospacer <- toupper(df.lib$Protospacer)
barcode_list <- as.list(df.lib$Protospacer)

# Import the sample's fastq file
fq <- readFastq(fq_file)
reads <- sread(fq)
reads <- as.character(reads)
total_reads <- length(reads)
print(paste("Total reads: ", total_reads, sep = ""))

# Find the "CACCG" TSS motif at position 20 of the read 
nums <- str_locate(reads, "CACCG")
reads <- reads[!is.na(nums[,1])]
nums <- nums[ !is.na(nums[,1]) ,]
inds <- nums[,1] == 20
print(paste0("Reads with detectable CACCG TSS motif: ", sum(inds), " (", sum(inds)/total_reads*100, "%)"))
nums <- nums[inds,]
reads <- reads[inds]

# Extract the expected barcode sequence following the TSS motif 
reads <- substr(reads, nums[,2], nums[,2]+20)

# Map reads to barcodes with exact matches only 
print("Mapping reads to barcode library (exact matches only)")
sample_counts <- df.lib[,c("sgRNAID","Protospacer")]
truths <- rep(FALSE, length(reads))
counts_exact <- lapply(barcode_list, function(x) {
  hits <- str_detect(reads, x)
  truths <<- truths | hits
  return( sum(hits) )
})

reads_map_exact <- sum(truths)
print(paste0("Reads with exact barcode matches: ", reads_map_exact, " (", reads_map_exact/total_reads*100, "%)"))

# Export barcode counts to disk 
sample_counts[[sample_name]] <- unlist(counts_exact)
write.table(sample_counts,
            file = paste0(output_dir, "/", sample_name, "_counts.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)



