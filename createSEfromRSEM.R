#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(affy))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("dplyr"))

option_list <- list(
  make_option(c("-a", "--genesresultsfile"), type="character", help="genesresultsfile"),
  make_option(c("-b", "--isoformsresultsfile"), type="character", help="isoformfile"),
  make_option(c("-c", "--demographicfile"), type="character", default = NA, help="demographicsfile"),
  make_option(c("-d", "--inputfile"), type="character", help="inputfile"),
  make_option(c("-e", "--fastqfile"), type="character", help="fastqfile"),
  make_option(c("-f", "--picardmarkeddupfile"), type="character", help="picardmarkeddupfile"),
  make_option(c("-g", "--rseqcbamstatfile"), type="character", help="rseqcbamstatfile"),
  make_option(c("-z", "--multiqcgeneralstatfile"), type="character", help="genstatfile"),
  make_option(c("-i", "--multiqcinferexpfile"), type="character", help="multiqcinferexpfile"),
  make_option(c("-x", "--rseqcjunctionannfile"), type="character", help="rseqcjunctionannfile"),
  make_option(c("-k", "--rseqcreaddistfile"), type="character", help="rseqcreaddistfile"),
  make_option(c("-n", "--gtffile"), type="character", help="GTFfile"),
  make_option(c("-m", "--inferfile"), type="character", help="inferfile"),
  make_option(c("-o", "--prefix"), type="character", help="prefix"),
  make_option(c("-v", "--starfile"), type="character", help="starfiles"),
  make_option(c("-w", "--tinfile"), type="character", help="tinfiles")
)

opt <- parse_args(OptionParser(option_list=option_list))
genefile <- opt$genesresultsfile
trfile <- opt$isoformsresultsfile

removalOfChars <- function(file){
  colnames(file) <- gsub("[' \"\\.]","_",colnames(file))
}

if (!is.na(opt$demographicfile)){
  pfile1 <- opt$demographicfile
}
pfile2 <- opt$inputfile
pfile3 <- scan(opt$fastqfile,what="")
pfile4 <- scan(opt$picardmarkeddupfile,what="")
pfile5 <- scan(opt$rseqcbamstatfile,what="")
pfile6 <- scan(opt$multiqcgeneralstatfile,what="")
pfile7 <- scan(opt$multiqcinferexpfile,what="")
pfile8 <- scan(opt$rseqcjunctionannfile,what="")
pfile9 <- scan(opt$rseqcreaddistfile,what="")
pfile10<- opt$gtffile
pfile11<- scan(opt$inferfile,what="")
pfile12<- scan(opt$starfile,what="")
pfile13<- scan(opt$tinfile,what="")
##Reading in Phenotype Data Files


#Input file
file2 =  read.table(pfile2, header=TRUE, stringsAsFactors=FALSE, quote="")
inputfile <- file2 %>% group_by_(.dots=c("INDIVIDUAL_ID","SAMPLE_ID","LIBRARY_ID")) %>% 
  summarise(RG_ID = paste(RG_ID, collapse = ","),
            PLATFORM_UNIT = paste(PLATFORM_UNIT, collapse = ","),
            PLATFORM = paste(PLATFORM, collapse = ","),
            PLATFORM_MODEL = paste(PLATFORM_MODEL, collapse = ","),
            RUN_DATE = paste(RUN_DATE, collapse = ","),
            CENTER = paste(CENTER, collapse = ","), 
            R1 = paste(R1, collapse = ","), 
            R2 = paste(R2, collapse = ","))
phenofile = as.data.frame(inputfile)
rownames(phenofile) = inputfile$SAMPLE_ID
colnames(phenofile) = removalOfChars(phenofile)

#Demographics
if (!is.na(opt$demographicfile)){
  demo = read.table(pfile1, header=TRUE, stringsAsFactors=FALSE, quote="",row.names = 1, sep = "\t")
  indices <- match(row.names(phenofile),row.names(demo))
  file1 <- demo[indices,]
  file1 <- file1[(rownames(phenofile)),]
  colnames(file1) <- paste("Annot",colnames(file1),sep="_")
  colnames(file1) <-removalOfChars(file1)
  phenofile = cbind(phenofile,(file1))
}

#FASTQ
file3 = read.table(pfile3, header=TRUE,sep = "\t", stringsAsFactors=FALSE, quote="")

#Modify file3$Sample so it can be grouped by lib.ID via dplyr
file3$Sample = substr(file3$Sample,1,8)
for(x in 1:nrow(file3)){
  if(x%%2 != 0){
    file3$Sample[x] <- paste0(file3$Sample[x],"R1") 
  }else{
    file3$Sample[x] <- paste0(file3$Sample[x],"R2")
  }
}

#New "file3"
file3.new <- file3 %>% group_by_(.dots="Sample") %>% 
  summarise(adapter_content = paste(adapter_content, collapse = ","),
            Sequences.flagged.as.poor.quality = sum(Sequences.flagged.as.poor.quality),
            sequence_duplication_levels = paste(sequence_duplication_levels, collapse = ","),
            avg_sequence_length = mean(avg_sequence_length),
            Encoding = paste(Encoding, collapse = ","),
            kmer_content = paste(kmer_content, collapse = ","),
            per_base_sequence_quality = paste(per_base_sequence_quality, collapse = ","),
            sequence_length_distribution = paste(sequence_length_distribution, collapse = ","),
            Sequence.length = paste(Sequence.length, collapse = ","),
            File.type = paste(File.type, collapse = ","),
            basic_statistics = paste(basic_statistics, collapse = ","),
            per_sequence_gc_content = paste(per_sequence_gc_content, collapse = ","),
            Total.Sequences = sum(Total.Sequences),
            per_base_n_content = paste(per_base_n_content, collapse = ","),
            per_base_sequence_content = paste(per_base_sequence_content, collapse = ","),
            overrepresented_sequences = paste(overrepresented_sequences, collapse = ","),
            X.GC = mean(X.GC),
            total_deduplicated_percentage = mean(total_deduplicated_percentage),
            per_tile_sequence_quality = paste(per_tile_sequence_quality, collapse = ","),
            per_sequence_quality_scores = paste(per_sequence_quality_scores, collapse = ","))

a = 0
colnames(file3.new) = removalOfChars(file3.new)
matrix <- c()
while(a <nrow(file3.new)){
  a = a+2
  x <- file3.new[a-1,]
  colnames(x) <- paste("FastQC_R1",colnames(x),sep="_")
  y <- file3.new[a,]
  colnames(y) <- paste("FastQC_R2",colnames(y),sep="_")
  filex <- cbind(x,y)
  matrix <- rbind(matrix,filex)
}

phenofile <- cbind(phenofile,matrix)

#Picard Marked Duplicates
file4 = read.table(pfile4,header=TRUE, stringsAsFactors=FALSE, quote="",row.names = 1)
colnames(file4) <- paste("Picard_MarkDuplicates",colnames(file4),sep="_")
rownames(file4) <- gsub("_.*rgid$", "",rownames(file4))

file4 <- file4[(rownames(phenofile)),]
colnames(file4) = removalOfChars(file4)
phenofile = cbind(phenofile,file4)
#RSeqC BamStat
file5 = read.table(pfile5,header=TRUE, stringsAsFactors=FALSE, quote="", row.names = 1) 
colnames(file5) <- paste("RSeQC_BamStat",colnames(file5),sep="_")
file5 <- file5[(rownames(phenofile)),]
colnames(file5) = removalOfChars(file5)
phenofile = cbind(phenofile,file5)

#MultiQC, General Stats
file6 = read.table(pfile6, header=TRUE, stringsAsFactors=FALSE, quote="",row.names = 1)
colnames(file6) <- paste("MultiQC_General_Stats",colnames(file6),sep="_")
file6$F <- "F"
file6 <- file6[(rownames(phenofile)),]
file6$F <- NULL
colnames(file6) = removalOfChars(file6)
phenofile = cbind(phenofile,file6)

#MultiQC, InferredExperiment
file7 = read.table(pfile7, header=TRUE, stringsAsFactors=FALSE, quote="", row.names = 1)
rownames(file7) <-  gsub(".inferred_experiment","",row.names(file7))
colnames(file7) <- paste("RSeQC_Infer_Experiment", colnames(file7), sep="_")
file7 <- file7[(rownames(phenofile)),]
colnames(file7) = removalOfChars(file7)
phenofile = cbind(phenofile,file7)
#RSeqC Junction Annotation
file8 = read.table(pfile8, header=TRUE, stringsAsFactors=FALSE, quote="", row.names = 1)
rownames(file8) <-  gsub(".junction_annotation","",row.names(file8))
colnames(file8) <- paste("RSeQC_Junction_Annotation", colnames(file8), sep="_")
file8 <- file8[(rownames(phenofile)),]
colnames(file8) = removalOfChars(file8)
phenofile = cbind(phenofile,file8)

#RSeqC Read Distribution
file9  = read.table(pfile9, header=TRUE, stringsAsFactors=FALSE, quote="", row.names = 1)
rownames(file9) <-  gsub(".read_distribution","",row.names(file9))
colnames(file9) <- paste("RSeQC_Read_Distribution", colnames(file9), sep="_")
file9 <- file9[(rownames(phenofile)),]
colnames(file9) = removalOfChars(file9)
phenofile = cbind(phenofile,file9)

#Reading in PCA data
file11 = read.table(pfile11, header=TRUE, stringsAsFactors=FALSE, quote="", row.names=2)
file11 = file11[2:3]
colnames(file11) <- c("SNPRelate_PCA_PC1","SNPRelate_PCA_PC2")
file11 <- file11[(rownames(phenofile)),]
colnames(file11) = removalOfChars(file11)
phenofile = cbind(phenofile,file11)

file12 = read.table(pfile12, row.names = 1, header = TRUE)
file12 = file12[(rownames(phenofile)),]
colnames(file12) = removalOfChars(file12)
colnames(file12) = paste("STAR_2pass",colnames(file12),sep="_")
phenofile = cbind(phenofile,file12)


forFile13 <- c()
for (i in pfile13){
  file13 = read.table(i, row.names = 1, header = TRUE)
  forFile13 <- rbind(forFile13,file13)
}
rownames(forFile13) <- gsub(".bam","",rownames(forFile13))

forFile13$C = substr(x = rownames(forFile13),start = 1, stop =34)

file13.new <- forFile13 %>% group_by_(.dots="C") %>% 
  summarise(TIN_mean = mean(TIN.mean.),
            TIN_median = mean(TIN.median.),
            TIN_stdev = mean(TIN.stdev.))

class(file13.new) <- "data.frame"
rownames(file13.new) <- file13.new$C
file13.new$C <- NULL


file13.new = file13.new[(rownames(phenofile)),]
colnames(file13.new) = removalOfChars(file13.new)
colnames(file13.new) = gsub("_$","",colnames(file13.new))
colnames(file13.new) = paste("RSeQC",colnames(file13.new),sep="_")
phenofile = cbind(phenofile,file13.new)

print("File loading complete")
##Reading in assay data

#gene data
files = scan(genefile, what="")
posterior_mean_count = c()
title.matrix = c()
fpkm = c()
exp_count = c()
tpm = c()
posterior_stdev_count = c()
pme_tpm = c()
pme_fpkm = c()
tpm_lowerbound = c()
tpm_upperbound = c()
tpm_coeff_quartile_variation = c()
fpkm_ci_lowerbound = c()
fpkm_ci_upperbound = c()
fpkm_coefficient_quartile_variation = c()

for(i in files) {
  temp = read.table(i, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
  posterior_mean_count = cbind(posterior_mean_count, temp$posterior_mean_count)
  fpkm = cbind(fpkm, temp$FPKM)
  exp_count = cbind(exp_count, temp$expected_count)
  tpm = cbind(tpm, temp$TPM)
  posterior_stdev_count = cbind(posterior_stdev_count, temp$posterior_standard_deviation_of_count)
  pme_tpm = cbind(pme_tpm, temp$pme_TPM)
  pme_fpkm = cbind(pme_fpkm, temp$pme_FPKM)
  tpm_lowerbound = cbind(tpm_lowerbound, temp$TPM_ci_lower_bound)
  tpm_upperbound = cbind(tpm_upperbound, temp$TPM_ci_upper_bound)
  tpm_coeff_quartile_variation = cbind(tpm_coeff_quartile_variation, temp$TPM_coefficient_of_quartile_variation)
  fpkm_ci_lowerbound = cbind(fpkm_ci_lowerbound, temp$FPKM_ci_lower_bound)
  fpkm_ci_upperbound = cbind(fpkm_ci_upperbound, temp$FPKM_ci_upper_bound)
  fpkm_coefficient_quartile_variation = cbind(fpkm_coefficient_quartile_variation, temp$FPKM_coefficient_of_quartile_variation)
  title = basename(i)
  title <- gsub(".genes.results","",title)
  title.matrix = cbind(title.matrix, title)
}

name.matrix <- temp$gene_id

naming <- function(name,title,matrix,phenofile){
  colnames(matrix) <- as.character(title)
  rownames(matrix) <- as.character(name)
  matrix <- matrix[,(rownames(phenofile))]
  return(matrix)
}

#Naming the rows/matrices
posterior_mean_count <- naming(name.matrix,title.matrix,posterior_mean_count,phenofile)
fpkm <- naming(name.matrix,title.matrix,fpkm,phenofile)
exp_count <- naming(name.matrix,title.matrix,exp_count,phenofile)
tpm <- naming(name.matrix,title.matrix,tpm,phenofile)
posterior_stdev_count <- naming(name.matrix,title.matrix,posterior_stdev_count,phenofile)
pme_tpm <- naming(name.matrix,title.matrix,pme_tpm,phenofile)
pme_fpkm <- naming(name.matrix,title.matrix,pme_fpkm,phenofile)
tpm_lowerbound <- naming(name.matrix,title.matrix,tpm_lowerbound,phenofile)
tpm_upperbound <- naming(name.matrix,title.matrix,tpm_upperbound,phenofile)
tpm_coeff_quartile_variation <- naming(name.matrix,title.matrix,tpm_coeff_quartile_variation,phenofile)
fpkm_ci_lowerbound <- naming(name.matrix,title.matrix,fpkm_ci_lowerbound,phenofile)
fpkm_ci_upperbound <- naming(name.matrix,title.matrix,fpkm_ci_upperbound,phenofile)
fpkm_coefficient_quartile_variation <- naming(name.matrix,title.matrix,fpkm_coefficient_quartile_variation,phenofile)


#trancsript data
trfiles = scan(trfile, what="")
titletr.matrix = c()

posterior_mean_count_tr = c()
fpkm_tr = c()
exp_count_tr = c()
tpm_tr = c()
posterior_stdev_count_tr = c()
pme_tpm_tr = c()
pme_fpkm_tr = c()
tpm_lowerbound_tr = c()
tpm_upperbound_tr = c()
tpm_coeff_quartile_variation_tr = c()
fpkm_ci_lowerbound_tr = c()
fpkm_ci_upperbound_tr = c()
fpkm_coefficient_quartile_variation_tr = c()
isopct_tr = c()
isopct_pmetpm_tr = c()

#Creating assaydata for transcripts
for(i in trfiles) {
  temp = read.table(i, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
  posterior_mean_count_tr = cbind(posterior_mean_count_tr, temp$posterior_mean_count)
  fpkm_tr =cbind(fpkm_tr,temp$FPKM)
  exp_count_tr = cbind(exp_count_tr, temp$expected_count)
  tpm_tr = cbind(tpm_tr, temp$TPM)
  isopct_tr = cbind(isopct_tr, temp$IsoPct)
  posterior_stdev_count_tr = cbind(posterior_stdev_count_tr, temp$posterior_standard_deviation_of_count)
  pme_tpm_tr = cbind(pme_tpm_tr, temp$pme_TPM)
  pme_fpkm_tr = cbind(pme_fpkm_tr, temp$pme_FPKM)
  isopct_pmetpm_tr = cbind(isopct_pmetpm_tr, temp$IsoPct_from_pme_TPM)
  tpm_lowerbound_tr = cbind(tpm_lowerbound_tr, temp$TPM_ci_lower_bound)
  tpm_upperbound_tr = cbind(tpm_upperbound_tr, temp$TPM_ci_upper_bound)
  tpm_coeff_quartile_variation_tr = cbind(tpm_coeff_quartile_variation_tr, temp$TPM_coefficient_of_quartile_variation)
  fpkm_ci_lowerbound_tr = cbind(fpkm_ci_lowerbound_tr, temp$FPKM_ci_lower_bound)
  fpkm_ci_upperbound_tr = cbind(fpkm_ci_upperbound_tr, temp$FPKM_ci_upper_bound)
  fpkm_coefficient_quartile_variation_tr = cbind(fpkm_coefficient_quartile_variation_tr,temp$FPKM_coefficient_of_quartile_variation)
  title = basename(i)
  title <- gsub(".isoforms.results","",title)
  titletr.matrix = cbind(titletr.matrix, title)
}

nametr.matrix <- temp$transcript_id

#Naming the rows/columns
posterior_mean_count_tr <- naming(nametr.matrix,titletr.matrix,posterior_mean_count_tr,phenofile)
fpkm_tr <- naming(nametr.matrix,titletr.matrix,fpkm_tr,phenofile)
exp_count_tr <- naming(nametr.matrix,titletr.matrix,exp_count_tr,phenofile)
tpm_tr <- naming(nametr.matrix,titletr.matrix,tpm_tr,phenofile)
posterior_stdev_count_tr <- naming(nametr.matrix,titletr.matrix,posterior_stdev_count_tr,phenofile)
pme_tpm_tr <- naming(nametr.matrix,titletr.matrix,pme_tpm_tr,phenofile)
pme_fpkm_tr <- naming(nametr.matrix,titletr.matrix,pme_fpkm_tr,phenofile)
tpm_lowerbound_tr <- naming(nametr.matrix,titletr.matrix,tpm_lowerbound_tr,phenofile)
tpm_upperbound_tr <- naming(nametr.matrix,titletr.matrix,tpm_upperbound_tr,phenofile)
tpm_coeff_quartile_variation_tr <- naming(nametr.matrix,titletr.matrix,tpm_coeff_quartile_variation_tr,phenofile)
fpkm_ci_lowerbound_tr <- naming(nametr.matrix,titletr.matrix,fpkm_ci_lowerbound_tr,phenofile)
fpkm_ci_upperbound_tr <- naming(nametr.matrix,titletr.matrix,fpkm_ci_upperbound_tr,phenofile)
fpkm_coefficient_quartile_variation_tr <- naming(nametr.matrix,titletr.matrix,fpkm_coefficient_quartile_variation_tr,phenofile)
isopct_tr <- naming(nametr.matrix,titletr.matrix,isopct_tr,phenofile)
isopct_pmetpm_tr <- naming(nametr.matrix,titletr.matrix,isopct_pmetpm_tr,phenofile)
print("Finished with making counts matrix")
##Adding in rowRanges data
txdb = makeTxDbFromGFF(pfile10)
genedata = genes(txdb)
trdata = transcripts(txdb)

##Adding in rowData data from Biomart
mart = useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl", GRCh=37)
gene.attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "description","band", "external_gene_name", "external_gene_source", "transcript_count","entrezgene")
gene.info = getBM(gene.attributes, mart = mart)

transcript.attributes = c("ensembl_transcript_id", "ensembl_gene_id", "hgnc_symbol", "gene_biotype", "transcript_biotype", "description","band", "external_gene_name", "external_gene_source", "external_transcript_name", "external_transcript_source_name", "entrezgene")
transcript.info = getBM(transcript.attributes, mart = mart)

g = gene.info[gene.info[,1]%in%row.names(posterior_mean_count),]
df <- aggregate(g[2:9],g[1],unique)
mcols(genedata) <- c(mcols(genedata),df[2:9])

t = transcript.info[transcript.info[,1]%in%row.names(posterior_mean_count_tr),]
dft <- aggregate(t[2:12],t[1],unique)
mcols(trdata) <- c(mcols(trdata),dft[2:12])


##Creation of SummarizedExperiment object
se <- SummarizedExperiment(assays = list(posterior_mean_count=posterior_mean_count,fpkm=fpkm,exp_count=exp_count,tpm=tpm,posterior_stdev_count=posterior_stdev_count,pme_tpm=pme_tpm,pme_fpkm=pme_fpkm,tpm_lowerbound=tpm_lowerbound,tpm_upperbound=tpm_upperbound,tpm_coeff_quartile_variation=tpm_coeff_quartile_variation,fpkm_ci_lowerbound=fpkm_ci_lowerbound,fpkm_ci_upperbound=fpkm_ci_upperbound,fpkm_coefficient_quartile_variation=fpkm_coefficient_quartile_variation), colData=phenofile,rowRanges = genedata)

set <- SummarizedExperiment(assays = list(posterior_mean_count_tr=posterior_mean_count_tr,fpkm_tr=fpkm_tr,exp_count_tr=exp_count_tr,tpm_tr=tpm_tr,posterior_stdev_count_tr=posterior_stdev_count_tr,pme_tpm_tr=pme_tpm_tr,pme_fpkm_tr=pme_fpkm_tr,tpm_lowerbound_tr=tpm_lowerbound_tr,tpm_upperbound_tr=tpm_upperbound_tr,tpm_coeff_quartile_variation_tr=tpm_coeff_quartile_variation_tr,fpkm_ci_lowerbound_tr=fpkm_ci_lowerbound_tr,fpkm_ci_upperbound_tr=fpkm_ci_upperbound_tr,fpkm_coefficient_quartile_variation_tr=fpkm_coefficient_quartile_variation_tr,isopct_tr=isopct_tr,isopct_pmetpm_tr=isopct_pmetpm_tr), colData=phenofile,rowRanges = trdata)

saveRDS(se,paste(opt$prefix,"_Gene_Expression.rds",sep=""))
saveRDS(set,paste(opt$prefix,"_Isoform_Expression.rds",sep=""))

