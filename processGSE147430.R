
install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("minfi")

options(stringsAsFactors=FALSE)

library(GEOquery)
library(minfi)

# specify directory to which the data will be downloaded and results will be saved
outpath = ""

gse.id = "GSE147430"
gset <- getGEO(
  gse.id, GSEMatrix=TRUE, getGPL=FALSE
)
gset1 <- gset[[1]]
gset2 <- gset[[2]]

covariates1 = pData(gset1)
covariates2 = pData(gset2)
covariates = rbind(covariates1,covariates2)
covariates$wellID = unlist(regmatches(covariates$title, gregexpr("\\[.*?\\]", covariates$title)))
covariates$wellID = gsub("\\[", "", covariates$wellID)
covariates$wellID = gsub("\\]", "", covariates$wellID)
covariates$wellID = paste0("X",covariates$wellID)
row.names(covariates) = covariates$wellID
covariates$smoking = 1.0*(covariates$"smoking_status:ch1" == "Smoker")

download.file(
  url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147430/suppl/GSE147430_Matrix_AverageBeta.txt.gz", 
  destfile = paste0(outpath,"/",gse.id,"_beta_values.txt.gz")
)
gunzip(
  filename = paste0(outpath,"/",gse.id,"_beta_values.txt.gz"),
  destname = paste0(outpath,"/",gse.id,"_beta_values.txt")
)
betas = read.table(
  paste0(outpath,"/",gse.id,"_beta_values.txt"),
  sep="\t", header=TRUE
)
row.names(betas) = betas$ID_REF

betas = betas[,colnames(betas) %in% row.names(covariates)]
covariates = covariates[row.names(covariates) %in% colnames(betas),]
covariates = covariates[match(colnames(betas),row.names(covariates)),]
betas = as.matrix(betas)

# analysis
testForEffectOfSmoking = function(Y,meta){
	return( summary(lm(Y~smoking,data=meta))$coefficients[2,] )
}
isolate_tests = apply(betas,MARGIN=1,FUN=testForEffectOfSmoking,meta=covariates)
isolate_tests = t(isolate_tests)
isolate_tests = data.frame(isolate_tests)
isolate_tests$p.adj = p.adjust(isolate_tests$Pr,method="fdr")

save(isolate_tests,file=paste0(outpath,"GSE147430_testing.RDATA"))
