#####################################################################################################
############    code by Walakira Andrew. adwalakira@gmail.com. For Walakira et al 2020   ############
#####################################################################################################


####################### make data frame of 1st kpca per gene with gene name as file name.



########################################################################################
######### making the data frame

setwd("your path here")

#Specify where the gene files are
Path2GeneFiles = "your path here" # original file source

files <- list.files(path = Path2GeneFiles, pattern="\\.csv$")


datF <- NULL
names <-c()   # get file names as gene names
for(file in files)
{
  data_in <- read.table(file.path(Path2GeneFiles, file), header=F, check.names=FALSE, stringsAsFactors=FALSE)
  print(data_in[1:5,1:4])
  data = data_in[,1]
  dim(data)
  
  names[[file]] = file
  datF = cbind(datF, data)
  
}

colnames(datF) = names
pheno = read.csv("/massstorage/URT/GEN/BIO3/Andrew/analysis/analysis_final/geneFiles_imputed/ABO.csv")
PHENOTYPE = pheno$PHENOTYPE
newdat = cbind(PHENOTYPE, datF)
newdat = data.frame(newdat)

write.csv(newdat, "kpca1_genes_final.csv")


# ==================================================================
######################################################################################
