#Set parameters
args <- commandArgs(trailingOnly = TRUE)
bed.fn <- args[1]
bim.fn <- args[2]
fam.fn <- args[3]
gds.fn <- args[4]
rdata.out.fn <- args[5]
log.fn <- args[6]
slide.max.bp <- as.numeric(args[7])
r2.ld.threshold <- as.numeric(args[8])

#Load librariess
library(SNPRelate)
library(GWASTools)
library(GENESIS)

#Open the log file
sink(log.fn)

#Convert PLINK files to GDS files
snpgdsBED2GDS(bed.fn = bed.fn, 
              bim.fn = bim.fn, 
              fam.fn = fam.fn, 
              out.gdsfn = gds.fn)

#Set seed for LD pruning
set.seed(100)

#Perform LD pruning
gds <- snpgdsOpen(gds.fn)
snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=slide.max.bp, 
                          ld.threshold=sqrt(r2.ld.threshold), verbose=FALSE)
pruned <- unlist(snpset, use.names=FALSE)
print(paste0("Number of SNPs after LD pruning: ", length(pruned)))

#Get KING relatedness estimates
ibd <- snpgdsIBDKING(gds, snp.id=pruned)
colnames(ibd$kinship) <- ibd$sample.id
rownames(ibd$kinship) <- ibd$sample.id

#Close the GDS file so it can be read again for PC-AiR
snpgdsClose(gds)

#Run PC-AiR
geno <- GdsGenotypeReader(filename=gds.fn)
geno.data <- GenotypeData(geno)
pca <- pcair(geno.data, 
                       kinobj=ibd$kinship,
                       divobj=ibd$kinship,
                       snp.include=pruned)

#Write PC-AiR output
save(pca, file=rdata.out.fn)

#Close the GDS file and log file
close(geno.data)
sink()
