#!/usr/bin/env Rscript

## how to run 
# Rscript updog_Rscript.R vcf_file ploidy cpus

### NOTE: this will only work if names in vcf are in particular order
## SpeciesPloidy_Pop_ID: CD2_AT_6 

# packages needed: data.table, vcfR, updog (https://github.com/dcgerard/updog)

## accept argument from the command line. 
args <- commandArgs(trailingOnly=TRUE)

## check correct
if (length(args) != 3) {
  stop("Arguments for vcf, ploidy, and number of cpus must be provided", call.=FALSE)
}

## set all necessary file/variable names for this run using command args
vcf_file <- args[1]
ploidy <- as.numeric(args[2])
cpus <- as.numeric(args[3])

cat('vcf file: ',vcf_file,'\n')
cat('ploidy: ',ploidy,'\n')
cat('cpus: ',cpus,'\n','\n')

## load packages 
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(vcfR)))
suppressWarnings(suppressPackageStartupMessages(library(updog)))

## read in vcf
cat('Reading in vcf.... ','\n')
vcf <- read.vcfR(vcf_file, verbose = FALSE)
cat('vcf successfully read in.... ','\n','\n')

#get positions
chrom <- getCHROM(vcf)
pos <- getPOS(vcf)
pos_ID <- paste(chrom,pos,sep = ':')

## get AD: allelic depth
cat('calculating allelic depth....','\n')
ad <- extract.gt(vcf, element = 'AD')

#make Pop_ID
indv <- colnames(ad)
Sp <- rep(NA,times=length(indv))
Ploidy <- rep(NA,times=length(indv))
Pop <- rep(NA,times=length(indv))
ID <- rep(NA,times=length(indv))
All <- rep(NA,times=length(indv))
for (i in 1:length(indv)){
  SpP <- unlist(strsplit(as.character(indv[i]),"_"))[1]
  Sp[i] <- gsub('\\d','',SpP,perl=TRUE)
  Ploidy[i] <-  gsub('(\\D)','',SpP,perl=TRUE)
  Pop[i] <- unlist(strsplit(as.character(indv[i]),"_"))[2]
  ID[i] <- unlist(strsplit(as.character(indv[i]),"_"))[3]
  All[i] <- as.character(indv[i])
}
Pop_ID <- data.frame(Sp=Sp,Ploidy=Ploidy,Pop=Pop,ID=ID,All=All,
                     SpPloidy=paste0(Sp,Ploidy))

## split ad by ploidy
ploidy_index <- which(Pop_ID$Ploidy == ploidy)
ad_ploidy <- ad[,ploidy_index]

cat('Number of individuals: ',dim(ad_ploidy)[2],'....','\n')
cat('Number of loci: ',dim(ad_ploidy)[1],'....','\n','\n')

## get total and reference ad
cat('gathering total and reference allelic depth....','\n')
tot_ad <- apply(ad_ploidy, c(1,2), function(df) sum(as.numeric(unlist(strsplit(as.character(df),',')))))
ref_ad <- apply(ad_ploidy, c(1,2), function(df) as.numeric(unlist(strsplit(as.character(df),','))[1]))

## run updog 
cat('running updog....','\n')
updog_out <- multidog(refmat=ref_ad,
                       sizemat=tot_ad,
                       ploidy=ploidy,
                       model = "hw",
                       nc = cpus)

#write it out
outfile <- paste0('updog',ploidy,'_out.RDS')
cat('\n','saving output as: ',outfile,'\n')
saveRDS(updog_out,outfile)


