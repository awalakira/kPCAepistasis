
#####################################################################################################
############    code by Walakira Andrew. adwalakira@gmail.com. For Walakira et al 2020   ############
#####################################################################################################

#####Imputing missing genotypes using kNN
##### Used the bnstruct package in R


### Imputing missing data for each gene file

library(bnstruct)

setwd("your path here")

#Specify where the gene files are
Path2GeneFiles = "your path here" # original file source

files <- list.files(path = Path2GeneFiles, pattern="\\.csv$")

get_pheno <- read.csv("your path here /ABO.csv") # get pheno continuous
PHENOTYPE = get_pheno$PHENOTYPE

for(file in files)
{
  dat = read.csv(file.path(Path2GeneFiles, file), header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
  print("******************************************************************************************")
  print(head(dat))
  print(sum(is.na(dat)))
  print("====================================================================")
  datx = as.matrix(dat[,2:ncol(dat)])
  print(head(datx))
  print(sum(is.na(datx)))
  
  
  dat_imputed = knn.impute(datx, k = 10, cat.var = 1:ncol(datx), to.impute = 1:nrow(datx), using = 1:nrow(datx))
  
  dat_imputedOut = cbind(PHENOTYPE, dat_imputed)
  
  print("====================================================================")
  print(sum(is.na(dat_imputedOut)))
  print(head(dat_imputedOut))
  print("******************************************************************************************")
  
  
  write.table(dat_imputedOut, file, row.names=F, col.names=TRUE, quote=FALSE, sep=",")
  
  
}

