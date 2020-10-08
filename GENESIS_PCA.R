#Set parameters
args <- commandArgs(trailingOnly = TRUE)
bed.fn <- args[1]
bim.fn <- args[2]
fam.fn <- args[3]
gds.fn <- args[4]
pcvec.out.fn <- args[5]
eigenval.out.fn <- args[6]
slide.max.bp <- as.numeric(args[7])
r2.ld.threshold <- as.numeric(args[8])

#Load librariess
library(SNPRelate)
library(GWASTools)
library(GENESIS)

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

#Get KING relatedness estimates
ibd <- snpgdsIBDKING(gds, snp.id=pruned)
colnames(ibd$kinship) <- ibd$sample.id
rownames(ibd$kinship) <- ibd$sample.id

#Close the GDS file so it can be read again for PC-AiR
snpgdsClose(gds)

#Run PC-AiR
geno <- GdsGenotypeReader(filename=gds.fn)
geno.data <- GenotypeData(geno)
pcair.results <- pcair(geno.data, 
                       kinobj=ibd$kinship,
                       divobj=ibd$kinship,
                       snp.include=pruned)

#Write PC-AiR output
pc.frame <- as.data.frame(pcair.results$vector)
colnames(pc.frame) <- paste0("PC", 1:ncol(pc.frame))
pc.frame$SampleID <- rownames(pc.frame)
pc.frame <- pc.frame[,c(ncol(pc.frame), 1:(ncol(pc.frame)-1))]
write.table(pc.frame, pcvec.out.fn,  sep="\t", quote=F, row.names=F, col.names=T)
write.table(pcair.results$values, eigenval.out.fn,  sep="\t", quote=F, row.names=F, col.names=F)

#Close the GDS file
close(geno.data)
