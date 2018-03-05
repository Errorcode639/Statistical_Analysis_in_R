# R Linkage Analysis Test Code
# Ada Madejska
# March 5th, 2018
# The purpose of this code is to do analysis on genotype and phenotype data

library(GenABEL)
library(MASS)
library(MatrixEQTL)
data <- load.gwaa.data(phenofile="phenotype.dat", genofile="genotype.raw")

# Look at the distribution of the continuous trait
ct <- phdata(data)$ct

hist(ct, col="slateblue", main="Distribution of ct")
# Check wheter the ct depends on sex col of phenotypic dataset
boxplot(ct~phdata(data)$sex, col=c("blue", "red"), xlab="sex", names=c("F", "M"), main="ct by sex")

# Look at the quality of markers
qc <- check.marker(data, call=0.95, perid.call=0.95, maf=1e-08, p.lev=1e-08)

# Create a new dataset based on the qc results
data.qc <-data[qc$idok, qc$snpok]

# Test using  Cochran-Armitage Trend Test
an <- qtscore(ct~1, data=data.qc, trait="gaussian")

# See if there's a strong correlation btwn the continous trait and any chromosome
plot(an, col=c("olivedrab", "slateblue"), cex=.5, main="Manhattan plot")

# Make sure each p-value is actually significant
estlambda(an[, "P1df"], plot=T)

plot(an, col=c("olivedrab", "slateblue"), cex=.5, main="Manhattan plot", df="Pc1df")

# Use Bonferroni correction to correct for false positives
pval.threshold <- 0.05
bonferroni <- -log10(pval.threshold / nids(data.qc))

plot(an, col=c("olivedrab", "slateblue"), cex=.5, main="Manhattan plot", df="Pc1df")
abline(h=bonferroni, lty=3, col="red")

useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
pvOutputThershold = 1e-2
errorCovariance = numeric()

# Initialize new dataframe called snps
snps = SlicedData$new()
snps$fileDelimiter = "\t"
snps$fileOmitCharacters = "NA"
snps$fileSkipRows = 1
snps$fileSkipColumns = 1
snps$fileSliceSize = 2000
snps$LoadFile("SNP.txt")

# Initialize new dataframe called gene
gene = SlicedData$new()
gene$fileDelimiter = "\t"      
gene$fileOmitCharacters = "NA" 
gene$fileSkipRows = 1          
gene$fileSkipColumns = 1       
gene$fileSliceSize = 2000      
gene$LoadFile("GE.txt")

# Initialize new dataframe called cvrt
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      
cvrt$fileOmitCharacters = "NA" 
cvrt$fileSkipRows = 1          
cvrt$fileSkipColumns = 1       
cvrt$LoadFile("Covariates.txt")

# Start qQTL Analysis - get a list of SNPs that exceed the thresholds 
me = Matrix_eQTL_engine(snps=snps, gene=gene, cvrt=cvrt, output_file_name="outfile.txt",
  pvOutputThreshold = pvOutputThershold, useModel = useModel, errorCovariance = errorCovariance, 
  verbose = TRUE, pvalue.hist = FALSE,min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE)

cat('Analysis done in:', me$time.in.sec, 'seconds', '\n')
cat('Detected eQTLs:', '\n')
show(me$all$eqtls)

# Make temporary files to store cis vs trans information
output_file_name_cis = tempfile()
output_file_name_tra = tempfile()

pvOutputThershold_cis = 1e-2
pvOutputThershold_tra = 1e-2
cisDist = 1e6

snpspos = read.table("snpsloc.txt", header=TRUE, stringsAsFactors=FALSE)
genepos = read.table("geneloc.txt", header=TRUE, stringsAsFactors=FALSE)

me = Matrix_eQTL_main(snps=snps, gene=gene, cvrt=cvrt, output_file_name=output_file_name_tra,
                      pvOutputThreshold=pvOutputThershold_tra, useModel=useModel, errorCovariance=errorCovariance,
                      verbose=TRUE, output_file_name.cis=output_file_name_cis, pvOutputThreshold.cis=pvOutputThershold_cis,
                      snpspos=snpspos, genepos=genepos, cisDist=cisDist, pvalue.hist="qqplot",
                      min.pv.by.genesnp=FALSE, noFDRsaveMemory=FALSE)