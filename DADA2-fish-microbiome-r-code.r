library(dada2)
library(ggplot2)
library(phyloseq)
library(dplyr)
library(reshape2)
library(tibble)


###################### Quality-filter reads and create Amplicon Sequence Variant tables 

###adjust path name as appropriate

getwd() 

path <- "C:/Users/eweli/PC-Documents/1-UOG-research/Fish-microbiome"



list.files(path)

#you have to see your fasta.gz files -- if you get character() it means dada2 cannot see your sequencing data 

# Samplename is everything before the first underscore

fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq.gz", full.names = TRUE))


# Extract sample names, assuming filenames have format: SAMPLENAME_*.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Perform filtering and trimming
# Assign the filenames for the filtered fastq.gz files.
# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "filtered") 

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter the forward and reverse reads
# WINDOWS USERS: set multithread=FALSE

#notes on trimming - see folder in Geneious - if you want to trimm both fwd and rev PCR1 primers you only need to trim left side on both fwd and rev 
#because the reverse read starts with the reverse primer, 
#the fwd primer is 19 bp long plus 4 bp of index equals 23 bp 
#the rev primer is 20 bp long plus 4 bp of index equals 24 bp 

#note I tried truncQ=2 and many ASVs were assinged to "Bacteria" only -no lower tax level after manual blastn check their were Porites mitochondrial seq
# re-do everthing wiht truncQ-10



out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, truncQ=10, rm.phix=TRUE, compress=TRUE, multithread=FALSE, 
            trimLeft = 24, trimRight = 24, minLen = 50)
            
write.table(out,file="filtered-information.txt")   

#selecting option for trimming 

#  compress = TRUE, compress files to .gz 
#  truncQ = 10 (Optional). Default 2. Truncate reads at the first instance of a quality score less than or equal to truncQ.
#  truncLen = 0,(Optional). Default 0 (no truncation). Truncate reads after truncLen bases. Reads shorter than this are discarded.
# trimLeft = 0,(Optional). Default 0. The number of nucleotides to remove from the start of each read. If both truncLen and trimLeft are provided, filtered reads will have length truncLen-trimLeft.
# trimRight = 0,(Optional). Default 0. The number of nucleotides to remove from the end of each read. If both truncLen and trimRight are provided, truncation will be performed after trimRight is enforced.
# maxLen = Inf,
# minLen = 50,(Optional). Default 20. Remove reads with length less than minLen. minLen is enforced after trimming and truncation.
# maxN = 0,
# minQ = 0, (Optional). Default 0. After truncation, reads contain a quality score less than minQ will be discarded.
# maxEE = Inf,(Optional). Default Inf (no EE filtering). After truncation, reads with higher than maxEE "expected errors" will be discarded. Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))
# rm.phix = TRUE,



# Learn the Error Rates, it TAKES TIME! do first forward and then reverse
# Forward reads

errF <- learnErrors(filtFs, multithread=FALSE)

# Reverse reads

errR <- learnErrors(filtRs, multithread=FALSE)


# Dereplicate the filtered fastq files

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)


# Name the derep-class objects by the sample names

names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Infer the sequence variants in each sample

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merge the denoised forward and reverse reads and 
# save so you can merge runs

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)


seqtab <- makeSequenceTable(mergers)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE)

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("raw-reads", "filtered-reads", "denoised-reads", "merged-reads", "tabled-reads", "nonchim-reads")
rownames(track) <- sample.names
head(track)

write.table(track, "read_stats.txt",sep="\t",col.names=NA)

saveRDS(seqtab.nochim, "seqtabnochim.rds") 

---------------------------------------------------------------------


AVSseqlengthtab <- table(nchar(getSequences(seqtab.nochim)))
write.table(AVSseqlengthtab,file="AVSseqlengthtab.txt")




# Assign taxonomy


taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=FALSE)


#updated taxonomy database was downloaded from here 
https://benjjneb.github.io/dada2/training.html

https://zenodo.org/records/14169026

# FIX the NAs in the taxa table


taxon <- as.data.frame(taxa,stringsAsFactors=FALSE)


taxon$Phylum[is.na(taxon$Phylum)] <- taxon$Kingdom[is.na(taxon$Phylum)]
taxon$Class[is.na(taxon$Class)] <- taxon$Phylum[is.na(taxon$Class)]
taxon$Order[is.na(taxon$Order)] <- taxon$Class[is.na(taxon$Order)]
taxon$Family[is.na(taxon$Family)] <- taxon$Order[is.na(taxon$Family)]
taxon$Genus[is.na(taxon$Genus)] <- taxon$Family[is.na(taxon$Genus)]


write.table(taxon,"silva_taxa_table.txt",sep="\t",col.names=NA)
write.table(seqtab.nochim, "silva_ASVs_table.txt",sep="\t",col.names=NA)

seqtab.nochim.tranposed <- t(seqtab.nochim) 

write.table(seqtab.nochim.transposed, "silva_ASVs_table_transposed.txt",sep="\t",col.names=NA)



##################### REMOVING MITOCHONDIRAL AND CHLOROPLAST READS

# convert the above files to excel and examine - removed NC - whihc should have  zero reads after chimera removal 

#you can replace the V4 sequences with ASV1- to ASV_no , but you have to do it taxon file and ASs frequency 

# Create phyloseq object from otu and taxonomy tables from dada2, along with the sample metadata

otu <- read.table("silva_ASVs_table.txt",sep="\t",header=TRUE, row.names=1)

taxon <- read.table("silva_taxa_table.txt",sep="\t",header=TRUE,row.names=1)

#make metadata file - the names of the samples must match 

otu2 <- read.table("silva_ASVs_table-2.txt",sep="\t",header=FALSE)
otu2T <- t(otu2) 
write.table(otu2T, "ASVs_table_transposed.txt",sep="\t",col.names=NA)





----------------------------------------------------------------------------------

otu <- read.table("ASVs_table_transposed.txt",sep="\t",header=TRUE, row.names=1)
otu <- t(otu) 
taxon <- read.table("silva_taxa_table.txt",sep="\t",header=TRUE,row.names=1)

samples<-read.table("metadata.txt",sep="\t",header=T,row.names=1)


OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))

ps

#output
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 44869 taxa and 66 samples ]
sample_data() Sample Data:       [ 66 samples by 2 sample variables ]
tax_table()   Taxonomy Table:    [ 44869 taxa by 6 taxonomic ranks ]


# you have to removed ASVs that are designated mitochondria or chloroplast or eukaryotes and NA for kingdom - usuall and artifact 

ps <- subset_taxa(ps, Family !="Mitochondria")
ps <- subset_taxa(ps, Order !="Chloroplast")
ps <- subset_taxa(ps, Kingdom !="Eukaryota")
ps <- subset_taxa(ps, Kingdom !="NA")

ps 
# filtered taxa with phyloseq, now export cleaned otu and taxa tables from phyloseq
otu = as(otu_table(ps), "matrix")
taxon = as(tax_table(ps), "matrix")

write.table(otu,"nochloronomito_otu_table.txt",sep="\t",col.names=NA)
write.table(taxon,"nochloronomito_taxa_table.txt",sep="\t",col.names=NA)


# you can remove very low abundance ASVs -otherwise you will have a lot of taxa to plot 

ntaxa(ps) 

#keep ASVs with minium of 5 reads on average


psmin5 <-filter_taxa(ps, function(x) mean(x) >5, TRUE)

ntaxa(psmin5)

psmin5


length(get_taxa_unique(psmin5, "Phylum"))
length(get_taxa_unique(psmin5, "Class"))
length(get_taxa_unique(psmin5, "Order"))
length(get_taxa_unique(psmin5, "Family"))
length(get_taxa_unique(psmin5, "Genus")) 

 
psmin5

#output
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 2821 taxa and 66 samples ]
sample_data() Sample Data:       [ 66 samples by 2 sample variables ]
tax_table()   Taxonomy Table:    [ 2821 taxa by 6 taxonomic ranks ]

otu = as(otu_table(psmin5), "matrix")
#you can exports this filtered data as you final data for plotting 

taxon = as(tax_table(psmin5), "matrix")
write.table(otu,"psmin5_ASV_table.txt",sep="\t",col.names=NA)
write.table(taxon,"psmin5_taxa_table.txt",sep="\t",col.names=NA)

