
#####################################################################################################
############    code by Walakira Andrew. adwalakira@gmail.com. For Walakira et al 2020   ############
#####################################################################################################

#### Generation of kernel principal components. Includes calculation of MI (information gain)


library(stringi)
library(infotheo)
library(expm)
library(data.table)
library(RSpectra)
library(foreach)
library(parallel)
library(doParallel)
library(dplyr)
library(infotheo)
library(snow)
library(SMUT)
library(microbenchmark)
library(filesstrings)

################################################################
################ Load data
################################################################ 

setwd("your path here")
dir.create("Results_Transfo_LDNet")
dir.create("Results_MI")   #store MI

#Specify where the gene files are
Path2GeneFiles = "your path here/geneFiles_imputed" # original file source

Path2Results <- "./Results_Transfo_LDNet"
outdir=Path2Results

Path2ResultsMI <- "./Results_MI"
outdirMI=Path2ResultsMI


files <- list.files(path = Path2GeneFiles, pattern="\\.csv$")

datpheno <- read.csv("your path here/ABO.csv") # get continupous phenotype
#datpEF = discretize(datpheno$PHENOTYPE, disc="equalfreq", nbins=2)
#print("equal freq")
#print(table(datpEF))


datp = discretize(datpheno$PHENOTYPE, disc="equalwidth", nbins=2)
print("Discretising by equal width")
print(table(datp))



ncores = detectCores()
print("########################################## ****************** number of cores")
print(ncores)



for(file in files)
{
  print(file)
  dat = read.csv(file.path(Path2GeneFiles, file), header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
  colms = ncol(dat)
  if(colms>3) {
    print("#### ====================================================================")
    print(head(dat))
    print(dim(dat))
    dataxy = cbind(datp, dat[,-1])  # dat[,-1] changed for dat

    print(head(dataxy))
    print("#### ====================================================================")

    nc = ncol(dataxy)
    #######################################################################
    ############### calculating Information gain
    #######################################################################
      	
    #Normalize by dividing by the entropy
    H <- entropy(datp, method="emp") #calculate entropy
        
    inter2=c()
    SNP1=c()
    SNP2=c()
    for(i in 2 : nc){
      for(j in 2 : nc){
        ii <- interinformation(dataxy[,c(1,i,j)], method = "emp")/H
        SNP1=c(SNP1,colnames(dataxy[i]))
        SNP2=c(SNP2,colnames(dataxy[j]))
        inter2=c(inter2,-ii)
      }
    }
    infoGain2=data.frame(SNP1,SNP2,inter2)
    #length(infoGain2$inter2)
          
        
    ###########################################################################################
    ## Making a square matrix from the dataframe of interinformation
    ###########################################################################################
    # get names for row and columns
    nameVals <- sort(unique(unlist(infoGain2[1:2])))
    # construct 0 matrix of correct dimensions with row and column names
    myMat <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
    # fill in the matrix with matrix indexing on row and column names
    myMat[as.matrix(infoGain2[c("SNP1", "SNP2")])] <- infoGain2[["inter2"]] 
    LD<-myMat
    diag(LD) = 0     #make diagonal zero i.e no info btn the same snp
    print(LD)



    ##################################################################
    ###### Make the Diffusion Laplacian matrix
    ################################################################## 
    #LD<-myMat 
    D<-diag(rowSums(LD))
    print("D")
    print(D)
    Laplacian<-as.matrix(D-LD)
    print("Laplacian ===========================================")
    print(Laplacian)   
        
    ########################################
    ## Make the Diffusion kernel for the various betas and get the linear combination using STATIS-UMKL -- no cosinus scaling
    ########################################
    pcs = 10 
    Beta=seq(0,10,0.1)
        
    mydata <- list()
    for(i in 1:length(Beta))
    {
      beta <- Beta[[i]]
      KG1<-expm(-beta*Laplacian)
      KG1<-as.matrix(KG1)
      mydata[[i]] <- KG1
      names(mydata)[[i]] <- paste0("KG_beta_", i) #Each kernel should have a unique name as required in mixKernel package
    }
        
    X <- mydata
    ind_beta <- 1/length(X)
    #Method2: Calculate every matrix with the equal weight and sum it up
    meta.kernel2 <- lapply(as.list(1:length(X)), function(x){
    X[[x]]*ind_beta })
    #Get a linear combination of all the matrix to get the final kernel
    beta_avg <- Reduce("+",meta.kernel2) #This is our new diffusion parameter
        
        
    ##************************************************************************
    ## Construct the diffusion kernels with the average diffusion parameter
    ##************************************************************************
    N<-dim(dataxy)[1] # # nb of individuals               
    KG<-as.matrix(beta_avg)
    datax<-as.matrix(dataxy[,-1])
        
    trdata <- t(datax)

    ### *** ORIGINAL    K0<-datax%*%KG%*%trdata  # Calculates the similarity between individuals 
    ###################### Parallelising1
    print("big mat print dataxy   Parallelising1 ************************************************")
    cl = makeSOCKcluster(ncores)
    registerDoParallel(cl)
    K0<-datax%*%KG%*%trdata
    stopCluster(cl)
    print(dim(K0))

    ### centering the matrix
    cl = makeSOCKcluster(ncores)
    registerDoParallel(cl)
    K = scale(K0, center=TRUE, scale=FALSE)
    stopCluster(cl)

    print("Centered matrix")
    print(dim(K))
    print(K[1:5, 1:5])


    # eigen value analysis
    #dim(K) ----> n by n
    cl = makeSOCKcluster(ncores)
    registerDoParallel(cl)
    eig <- eigs_sym(K/N,k=pcs, which="LM", sigma=NULL, lower=TRUE, retvec=TRUE)
    stopCluster(cl)
    
    print("eigen done")    

    eigVector = eig$vectors
    cl = makeSOCKcluster(ncores)
    registerDoParallel(cl)
    Yx <- K%*%eigVector
    stopCluster(cl)
         
    #save
    outputfile=paste0(outdir,"/",file)
    write.table(Yx,outputfile,row.names=F, col.names=FALSE,quote=FALSE,sep="\t")

    outputfile_MI=paste0(outdirMI,"/",file)
    write.table(infoGain2,outputfile_MI,row.names=F, col.names=FALSE,quote=FALSE,sep="\t")
    
  }
}






