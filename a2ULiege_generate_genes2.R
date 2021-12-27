
#####################################################################################################
############    code by Walakira Andrew. adwalakira@gmail.com. For Walakira et al 2020   ############
#####################################################################################################

### GENERATING GENE FILES


library(data.table)
library(dplyr)
library(stringi)
library(tidyverse)

setwd("your path here")

print(" saving files in folder: genes ***************************")

dir.create("genes_tabdelim")   #store tab delim files
Path2ResultsG <- "./genes_tabdelim"
outdirG=Path2ResultsG


#read in the files
data = read.csv("your path here/dataFull_ALL.csv", header=TRUE)
print(data[1:8,1:30])
data = data.frame(data)

datx = data[1:100,1:50]
write.csv(datx, "datx.csv")

# remove the part "_T*"
names(data) = gsub(pattern = "_.*", replacement="", x=names(data))
print("################## =========================================================================== data")
print(data[1:8,1:30])


## read in the snp-gene mapping file   snp_geneMapCorrect.xlsx
snp_map = read.delim("your path here/snp_gene_physicalmapping.txt", header=TRUE)

print("################## ========================================================================== genes")
print(head(snp_map))
print(dim(snp_map))


############ data
dataCols = data[,18:ncol(data)]
print(dataCols[1:8,1:8])
print("dimension of snps only ####################################### ----------------------------")
print(dim(dataCols))
print(class(dataCols))

cnames = colnames(dataCols)
cnames = as.vector(cnames)
print(length(cnames))

print(sum(is.na(dataCols)))


print("reduced mapping file ================================================ ********************************")
#### *****************************    Alternative code: take note - NAs cause downstream error
#mapped_subx = setDT(snp_map, key="rsID")[.(cnames)]
##mapped_subDF = setDF(mapped_subx)
#mapped_sub = mapped_subDF %>% select(rsID,ensg)

#gene_list = list(mapped_sub$ensg)  #1
#genes_list = as.vector(unlist(gene_list))  #1 gene names appearing alongside mapped snp
#genes = unique(genes_list)   # gene names appearing only once
#print(length(genes))
#print(genes[5])
#print(length(unique(mapped_sub$rsID)))

#snps_add_list = list(mapped_sub$rsID)   # list of snps in reduced map file
#snps_add = as.vector(unlist(snps_add_list))
#print(snps_add[5])


######### end of alternative code


#mapped_sub = filter(snp_map, rsID %in% cnames)  # can work
mapped_sub = snp_map[as.vector(snp_map$rsID) %in%cnames,]  # working


print("################## ================================================================ UNIQUE genes and SNPs")
print(class(mapped_sub))
print(head(mapped_sub))
print(dim(mapped_sub))


genes_list = as.vector(mapped_sub$ensg)  #1
genes = unique(genes_list)   # gene names appearing only once
print(length(genes))
print(genes[5])
print(length(unique(mapped_sub$rsID)))
print(length(unique(cnames)))

snps_add = as.vector(mapped_sub$rsID)   # list of snps in reduced map file
print(snps_add[5])


# ****************************************************************

print("sample subset ##################################")


### names as a vector
keepcovs <-c("rs9729550", "rs1815606", "rs7515488")
print(keepcovs)
dat.sub = dataCols[, names(dataCols)%in%keepcovs]
write.csv(dat.sub, "covsfSamp.csv")



for(i in 1:length(genes)){
  genex = genes[i]
  gene = as.character(genex)
  #snps = c()
  snps <- NULL
  for(j in 1:length(genes_list)){
    gene_chkx = as.vector(genes_list)
    gene_chkxv = gene_chkx[j]
    gene_chk = as.character(gene_chkxv)
    #(gene == gene_chk)==TRUE
    if(gene == gene_chk){
      #snps_add = as.vector(mapped_sub$rsID)
      snpx = snps_add[j]
      snpxc = as.character(snpx)
      snps=c(snps,snpxc)
      # names(dataCols)%in%snpx
      #if(grepl(snpx, cnames)==TRUE){
      #  snps=c(snps,snpxc)
      #}else{
      #  print("not in cnames")
      #  }
    }
  }
  # extract sub data set for specific gene 
  lensnp = length(snps)
  print(paste0(i, "_ ", gene, "_ ", lensnp))
  print(gene)
  print(snps)

  if(lensnp>1){
    gene_filex = dataCols[, names(dataCols)%in%snps]
    gene_filex = data.frame(gene_filex)
  
    print("##################### ==================== ************ gene file with extracted snps only") 
    print(class(gene_filex))
    print(head(gene_filex))
    print(dim(gene_filex))

    PHENOTYPE = data$PHENOTYPE
    print(length(PHENOTYPE))
    gene_file = cbind(PHENOTYPE, gene_filex)
  
    #write.csv  
    write.csv(gene_file, paste0(gene, ".csv"))
  }
}







#########################################################################




